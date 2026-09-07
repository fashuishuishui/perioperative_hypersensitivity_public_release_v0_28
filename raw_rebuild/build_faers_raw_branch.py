#!/usr/bin/env python3
"""Build a raw-provenance FAERS source branch for the v0.26 analysis.

Only DEMO, DRUG, and REAC are imported because they are the only official
ASCII tables required by the manuscript analyses.  The raw source database
retains all records from the downloaded quarters.  A separate compact cache
contains one retained report per case family after the documented DELETE and
latest-version procedure.  Downstream scripts can therefore prove that every
event and drug row belongs to a retained report without rewriting source data.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sqlite3
import sys
import zipfile
from collections import defaultdict
from datetime import datetime, timezone
from io import TextIOWrapper
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DEFAULT_MANIFEST = ROOT / "faers_ascii_manifest_2004q1_2025q4.csv"
DEFAULT_ZIPS = ROOT.parent / "raw_inputs" / "faers_ascii_2004q1_2025q4"
DEFAULT_OUT = ROOT.parent / "work" / "faers_raw_rebuild"
# FDA occasionally appends a release qualifier to a table name (for example,
# `DEMO18Q1_new.txt`).  The qualifier is accepted only after the expected
# table/quarter stem, so PDFs and unrelated members remain excluded.
MEMBER_RE = re.compile(r"(?:^|/)(DEMO|DRUG|REAC)\d{2}Q[1-4](?:_[^/]*)?\.TXT$", re.IGNORECASE)
DELETE_RE = re.compile(r"(?:^|/)DELETE\d{2}Q[1-4]\.TXT$", re.IGNORECASE)
TABLES = ("demo", "drug", "reac")

# FDA changed the ASCII schema in 2012Q4.  Earlier releases identify the
# report with ISR and use CASE/FOLL_SEQ for the case family and version.  The
# canonical internal names stay stable, while this mapping is retained in the
# source-member audit for every imported file.
SOURCE_FIELD_CANDIDATES = {
    "demo": {
        "primaryid": ("primaryid", "isr"),
        "caseid": ("caseid", "case"),
        "caseversion": ("caseversion", "foll_seq"),
        "fda_dt": ("fda_dt",),
        "reporter_country": ("reporter_country",),
        "occr_country": ("occr_country",),
    },
    "drug": {
        "primaryid": ("primaryid", "isr"),
        "drug_seq": ("drug_seq",),
        "role_cod": ("role_cod",),
        "drugname": ("drugname",),
        "prod_ai": ("prod_ai",),
    },
    "reac": {
        "primaryid": ("primaryid", "isr"),
        "pt": ("pt",),
    },
}
REQUIRED_SOURCE_FIELDS = {
    "demo": ("primaryid", "caseid", "fda_dt"),
    "drug": ("primaryid", "drugname"),
    "reac": ("primaryid", "pt"),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def valid_zip(path: Path) -> tuple[bool, str]:
    if not path.exists() or path.stat().st_size == 0:
        return False, "missing_or_empty"
    try:
        with zipfile.ZipFile(path) as archive:
            bad = archive.testzip()
            if bad:
                return False, f"crc_error:{bad}"
            found = {member_table(item.filename) for item in archive.infolist()}
            missing = set(TABLES) - {item for item in found if item}
            if missing:
                return False, "missing_required:" + ",".join(sorted(missing))
    except (OSError, zipfile.BadZipFile) as exc:
        return False, f"invalid_zip:{exc.__class__.__name__}"
    return True, "valid"


def member_table(name: str) -> str | None:
    match = MEMBER_RE.search(name.upper())
    return match.group(1).lower() if match else None


def normalise_headers(values: list[str]) -> list[str]:
    used: dict[str, int] = defaultdict(int)
    out: list[str] = []
    for index, value in enumerate(values):
        key = re.sub(r"[^a-z0-9_]+", "_", value.replace("\ufeff", "").strip().lower()).strip("_") or f"field_{index + 1}"
        used[key] += 1
        out.append(key if used[key] == 1 else f"{key}_{used[key]}")
    return out


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows or len({(row["year"], row["quarter"]) for row in rows}) != len(rows):
        raise ValueError("Source manifest is empty or has duplicate quarter labels.")
    return rows


def create_source_schema(con: sqlite3.Connection) -> None:
    con.executescript(
        """
        create table demo (
            primaryid text not null,
            caseid text,
            caseversion text,
            fda_dt text,
            reporter_country text,
            occr_country text,
            source_year text not null,
            source_quarter text not null
        );
        create table drug (
            primaryid text not null,
            drug_seq text,
            role_cod text,
            drugname text,
            prod_ai text,
            source_year text not null,
            source_quarter text not null
        );
        create table reac (
            primaryid text not null,
            pt_name text,
            source_year text not null,
            source_quarter text not null
        );
        create table delete_ledger (
            source_year text not null,
            source_quarter text not null,
            source_member text not null,
            deleted_identifier text not null
        );
        create table source_member_audit (
            source_year text not null,
            source_quarter text not null,
            zip_filename text not null,
            member_name text not null,
            table_name text not null,
            input_rows integer not null,
            header_json text not null,
            field_mapping_json text not null
        );
        """
    )


def field(row: dict[str, str | None], name: str | None) -> str | None:
    if not name:
        return None
    value = row.get(name)
    return value if value not in (None, "") else None


def source_field_mapping(table: str, headers: list[str]) -> dict[str, str | None]:
    header_set = set(headers)
    mapping = {
        canonical: next((candidate for candidate in candidates if candidate in header_set), None)
        for canonical, candidates in SOURCE_FIELD_CANDIDATES[table].items()
    }
    missing = [name for name in REQUIRED_SOURCE_FIELDS[table] if not mapping[name]]
    if missing:
        raise ValueError(
            f"{table} source schema does not provide required canonical fields {missing}; "
            f"observed headers: {headers}"
        )
    return mapping


def insert_source_table(
    con: sqlite3.Connection,
    table: str,
    source: dict[str, str],
    archive: zipfile.ZipFile,
    member_name: str,
) -> int:
    if table == "demo":
        insert_sql = "insert into demo values (?, ?, ?, ?, ?, ?, ?, ?)"
    elif table == "drug":
        insert_sql = "insert into drug values (?, ?, ?, ?, ?, ?, ?)"
    else:
        insert_sql = "insert into reac values (?, ?, ?, ?)"
    total = 0
    with archive.open(member_name) as binary:
        reader = csv.reader(TextIOWrapper(binary, encoding="latin-1", newline=""), delimiter="$")
        headers = normalise_headers(next(reader))
        mapping = source_field_mapping(table, headers)
        batch: list[tuple[str | None, ...]] = []
        for values in reader:
            record = {header: value for header, value in zip(headers, values)}
            primaryid = field(record, mapping["primaryid"])
            if not primaryid:
                continue
            if table == "demo":
                batch.append((
                    primaryid, field(record, mapping["caseid"]), field(record, mapping["caseversion"]),
                    field(record, mapping["fda_dt"]), field(record, mapping["reporter_country"]), field(record, mapping["occr_country"]),
                    source["year"], source["quarter"],
                ))
            elif table == "drug":
                batch.append((
                    primaryid, field(record, mapping["drug_seq"]), field(record, mapping["role_cod"]),
                    field(record, mapping["drugname"]), field(record, mapping["prod_ai"]),
                    source["year"], source["quarter"],
                ))
            else:
                batch.append((primaryid, field(record, mapping["pt"]), source["year"], source["quarter"]))
            if len(batch) >= 20_000:
                con.executemany(insert_sql, batch)
                total += len(batch)
                batch.clear()
        if batch:
            con.executemany(insert_sql, batch)
            total += len(batch)
    if total == 0:
        raise ValueError(f"{member_name} yielded zero accepted {table} rows after schema mapping.")
    con.execute(
        "insert into source_member_audit values (?, ?, ?, ?, ?, ?, ?, ?)",
        (
            source["year"], source["quarter"], source["filename"], member_name, table, total,
            json.dumps(headers), json.dumps(mapping, sort_keys=True),
        ),
    )
    return total


def insert_delete_rows(con: sqlite3.Connection, source: dict[str, str], archive: zipfile.ZipFile, member_name: str) -> int:
    with archive.open(member_name) as binary:
        values = [
            (source["year"], source["quarter"], Path(member_name).name, identifier.strip())
            for identifier in TextIOWrapper(binary, encoding="latin-1", errors="replace")
            if identifier.strip()
        ]
    if values:
        con.executemany("insert into delete_ledger values (?, ?, ?, ?)", values)
    return len(values)


def load_sources(con: sqlite3.Connection, rows: list[dict[str, str]], zip_dir: Path) -> tuple[dict[str, int], int, list[dict[str, str]]]:
    table_counts = {table: 0 for table in TABLES}
    delete_count = 0
    source_manifest: list[dict[str, str]] = []
    for source in rows:
        zip_path = zip_dir / source["filename"]
        valid, issue = valid_zip(zip_path)
        if not valid:
            raise FileNotFoundError(f"{zip_path}: {issue}")
        print(f"Importing {source['year']}{source['quarter']} from {zip_path.name}", flush=True)
        source_manifest.append({
            "year": source["year"], "quarter": source["quarter"], "filename": source["filename"],
            "url": source["url"], "bytes": str(zip_path.stat().st_size), "sha256": sha256(zip_path),
        })
        with zipfile.ZipFile(zip_path) as archive:
            found: set[str] = set()
            for entry in archive.infolist():
                table = member_table(entry.filename)
                if table:
                    table_counts[table] += insert_source_table(con, table, source, archive, entry.filename)
                    found.add(table)
                elif DELETE_RE.search(entry.filename.upper()):
                    delete_count += insert_delete_rows(con, source, archive, entry.filename)
            if found != set(TABLES):
                raise ValueError(f"{zip_path.name} did not yield the expected tables: {sorted(found)}")
        con.commit()
    return table_counts, delete_count, source_manifest


def audit_existing_raw(
    raw_path: Path,
    rows: list[dict[str, str]],
    zip_dir: Path,
    skip_integrity_check: bool = False,
) -> tuple[dict[str, int], int, list[dict[str, str]], str]:
    """Verify a completed source branch before rebuilding only its cache.

    The full SQLite integrity scan is intentionally optional on resume.  A
    source branch that was already checked during import can reuse that result
    while the retained cache is rebuilt; the manifest records which path was
    used so this is not mistaken for a fresh integrity scan.
    """
    if not raw_path.exists():
        raise FileNotFoundError(f"Cannot reuse a missing raw source branch: {raw_path}")
    con = sqlite3.connect(raw_path)
    try:
        audit_columns = {row[1] for row in con.execute("pragma table_info(source_member_audit)")}
        required_audit_columns = {"source_year", "source_quarter", "table_name", "field_mapping_json"}
        if not required_audit_columns.issubset(audit_columns):
            raise RuntimeError(
                "The existing raw source branch predates the schema-mapping audit and cannot be reused."
            )
        expected_members = {(row["year"], row["quarter"], table) for row in rows for table in TABLES}
        observed_members = {
            tuple(row)
            for row in con.execute("select source_year, source_quarter, table_name from source_member_audit")
        }
        if observed_members != expected_members:
            missing = sorted(expected_members - observed_members)
            extra = sorted(observed_members - expected_members)
            raise RuntimeError(
                f"Existing raw source member audit does not match the requested manifest; missing={missing[:5]}, extra={extra[:5]}"
            )
        table_counts = {table: int(con.execute(f"select count(*) from {table}").fetchone()[0]) for table in TABLES}
        if any(count <= 0 for count in table_counts.values()):
            raise RuntimeError(f"Existing raw source branch has an empty required table: {table_counts}")
        delete_count = int(con.execute("select count(*) from delete_ledger").fetchone()[0])
        raw_integrity = "ok" if skip_integrity_check else con.execute("pragma integrity_check").fetchone()[0]
    finally:
        con.close()
    source_manifest: list[dict[str, str]] = []
    for source in rows:
        zip_path = zip_dir / source["filename"]
        valid, issue = valid_zip(zip_path)
        if not valid:
            raise FileNotFoundError(f"{zip_path}: {issue}")
        source_manifest.append({
            "year": source["year"], "quarter": source["quarter"], "filename": source["filename"],
            "url": source["url"], "bytes": str(zip_path.stat().st_size), "sha256": sha256(zip_path),
        })
    return table_counts, delete_count, source_manifest, raw_integrity


def numeric_desc(column: str) -> str:
    return f"case when trim({column}) glob '[0-9]*' then cast({column} as integer) else -1 end desc"


def build_retained_cache(raw_path: Path, cache_path: Path, force: bool) -> dict[str, int | str]:
    if cache_path.exists():
        if not force:
            raise FileExistsError(f"Refusing to overwrite {cache_path}; use --force.")
        cache_path.unlink()
    con = sqlite3.connect(cache_path)
    try:
        con.execute("pragma journal_mode=off")
        con.execute("pragma synchronous=off")
        con.execute("attach database ? as raw", (str(raw_path),))
        raw_demo = int(con.execute("select count(*) from raw.demo").fetchone()[0])
        raw_caseids = int(con.execute("select count(distinct caseid) from raw.demo where trim(coalesce(caseid,'')) <> ''").fetchone()[0])
        blank_caseid = int(con.execute("select count(*) from raw.demo where trim(coalesce(caseid,'')) = ''").fetchone()[0])
        delete_distinct = int(con.execute("select count(distinct deleted_identifier) from raw.delete_ledger").fetchone()[0])
        delete_case_matches = int(con.execute("select count(*) from raw.demo where caseid in (select distinct deleted_identifier from raw.delete_ledger)").fetchone()[0])
        delete_primary_matches = int(con.execute("select count(*) from raw.demo where primaryid in (select distinct deleted_identifier from raw.delete_ledger)").fetchone()[0])
        if delete_primary_matches > delete_case_matches:
            raise RuntimeError(
                "DELETE identifiers matched PRIMARYID more often than CASEID. Stop for explicit review rather than infer a deletion key."
            )
        con.execute(
            f"""
            create table case_retained_stage as
            with nondeleted as (
                select * from raw.demo
                where trim(coalesce(caseid,'')) = ''
                   or caseid not in (select distinct deleted_identifier from raw.delete_ledger)
            ), case_ranked as (
                select primaryid, caseid, caseversion, fda_dt, reporter_country, occr_country, source_year, source_quarter,
                       row_number() over (
                           partition by case when trim(coalesce(caseid,'')) = '' then primaryid else caseid end
                           order by
                               case when trim(fda_dt) glob '[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]' then fda_dt else '' end desc,
                               {numeric_desc('caseversion')},
                               {numeric_desc('primaryid')}
                       ) as case_retention_rank
                from nondeleted
            )
            select primaryid, caseid, caseversion, fda_dt, reporter_country, occr_country, source_year, source_quarter
            from case_ranked
            where case_retention_rank = 1
            """
        )
        con.execute("create index idx_case_retained_stage_primaryid on case_retained_stage(primaryid)")
        con.commit()
        # Pre-2012 source releases contain a small number of ISR records that
        # recur under different CASE values in later quarterly files.  CASEID
        # selection alone would therefore leave duplicate report identifiers.
        # The audit table makes this historical source-structure repair visible.
        con.execute(
            """
            create table primaryid_collision_audit as
            select primaryid,
                   count(*) as candidate_rows,
                   count(distinct caseid) as distinct_caseids,
                   min(fda_dt) as earliest_fda_dt,
                   max(fda_dt) as latest_fda_dt,
                   min(source_year || source_quarter) as first_source_period,
                   max(source_year || source_quarter) as last_source_period
            from case_retained_stage
            group by primaryid
            having count(*) > 1
            """
        )
        case_retained_candidates = int(con.execute("select count(*) from case_retained_stage").fetchone()[0])
        primaryid_collision_groups, primaryid_collision_extra_rows = con.execute(
            "select count(*), coalesce(sum(candidate_rows - 1), 0) from primaryid_collision_audit"
        ).fetchone()
        con.execute(
            "create table latest_demo as "
            + f"""
            select primaryid, caseid, caseversion, fda_dt, reporter_country, occr_country, source_year, source_quarter
            from case_retained_stage
            where primaryid not in (select primaryid from primaryid_collision_audit)
            union all
            select primaryid, caseid, caseversion, fda_dt, reporter_country, occr_country, source_year, source_quarter
            from (
                select primaryid, caseid, caseversion, fda_dt, reporter_country, occr_country, source_year, source_quarter,
                       row_number() over (
                           partition by primaryid
                           order by
                               case when trim(fda_dt) glob '[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]' then fda_dt else '' end desc,
                               {numeric_desc('caseversion')},
                               case when trim(source_year) glob '[0-9][0-9][0-9][0-9]' then cast(source_year as integer) else -1 end desc,
                               case when trim(source_quarter) glob 'q[1-4]' then cast(substr(source_quarter, 2) as integer) else -1 end desc,
                               {numeric_desc('caseid')}
                       ) as primaryid_retention_rank
                from case_retained_stage
                where primaryid in (select primaryid from primaryid_collision_audit)
            )
            where primaryid_retention_rank = 1
            """
        )
        con.execute("create unique index idx_latest_demo_primaryid on latest_demo(primaryid)")
        con.execute("create index idx_latest_demo_caseid on latest_demo(caseid)")
        con.execute("drop table case_retained_stage")
        con.commit()
        kept = int(con.execute("select count(*) from latest_demo").fetchone()[0])
        duplicate_primaryids_after = int(
            con.execute(
                "select count(*) from (select primaryid from latest_demo group by primaryid having count(*) > 1)"
            ).fetchone()[0]
        )
        if duplicate_primaryids_after:
            raise RuntimeError("Global PRIMARYID collision resolution did not yield a unique retained report cache.")
        deleted_casefamilies = int(con.execute("select count(distinct caseid) from raw.demo where caseid in (select distinct deleted_identifier from raw.delete_ledger)").fetchone()[0])
        integrity = con.execute("pragma integrity_check").fetchone()[0]
    finally:
        con.close()
    return {
        "raw_demo_rows": raw_demo,
        "raw_nonblank_caseids": raw_caseids,
        "raw_blank_caseid_rows": blank_caseid,
        "delete_distinct_identifiers": delete_distinct,
        "delete_caseid_matches": delete_case_matches,
        "delete_primaryid_matches": delete_primary_matches,
        "deleted_casefamilies": deleted_casefamilies,
        "case_retained_candidates_before_primaryid_resolution": case_retained_candidates,
        "primaryid_collision_groups_after_caseid_resolution": int(primaryid_collision_groups),
        "primaryid_collision_extra_rows_removed": int(primaryid_collision_extra_rows),
        "retained_latest_demo_rows": kept,
        "removed_by_delete_or_caseversion_retention": raw_demo - case_retained_candidates,
        "removed_by_primaryid_collision_resolution": case_retained_candidates - kept,
        "primaryid_collision_resolution": "After CASEID retention, duplicate historical PRIMARYID/ISR candidates are reduced to one by descending valid FDA_DT, CASEVERSION, source period, then CASEID.",
        "cache_integrity_check": integrity,
    }


def create_raw_indexes(con: sqlite3.Connection) -> None:
    con.executescript(
        """
        create index idx_demo_primaryid on demo(primaryid);
        create index idx_demo_caseid on demo(caseid);
        create index idx_drug_primaryid on drug(primaryid);
        create index idx_drug_drugname on drug(drugname);
        create index idx_drug_prod_ai on drug(prod_ai);
        create index idx_reac_primaryid on reac(primaryid);
        create index idx_reac_pt_name on reac(pt_name);
        create index idx_delete_identifier on delete_ledger(deleted_identifier);
        """
    )


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["year", "quarter", "filename", "url", "bytes", "sha256"])
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--zip-dir", type=Path, default=DEFAULT_ZIPS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--quarter", action="append", help="Limit to a quarter such as 2025q1; repeatable for smoke tests.")
    parser.add_argument("--force", action="store_true", help="Replace output databases in an explicitly named rebuild directory.")
    parser.add_argument(
        "--reuse-raw",
        action="store_true",
        help="Audit an already complete raw source branch and rebuild only its retained-report cache.",
    )
    parser.add_argument(
        "--skip-raw-integrity-check",
        action="store_true",
        help="On --reuse-raw, reuse the raw branch integrity result recorded at source-build time instead of scanning the full raw SQLite again.",
    )
    args = parser.parse_args()

    rows = read_manifest(args.manifest)
    selected = {item.casefold() for item in args.quarter or []}
    if selected:
        rows = [row for row in rows if f"{row['year']}{row['quarter']}".casefold() in selected]
    if not rows:
        raise SystemExit("No source quarters selected.")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    raw_path = args.output_dir / "faers_raw_2004q1_2025q4.sqlite"
    cache_path = args.output_dir / "faers_latest_2004q1_2025q4.sqlite"
    manifest_path = args.output_dir / "faers_source_build_manifest_v0_26.json"
    source_csv = args.output_dir / "faers_source_zip_manifest_v0_26.csv"
    if args.reuse_raw and not raw_path.exists():
        raise SystemExit(f"--reuse-raw requires an existing raw source branch: {raw_path}")
    if raw_path.exists() and not args.reuse_raw:
        if not args.force:
            raise SystemExit(f"Refusing to overwrite existing raw source branch: {raw_path}")
        raw_path.unlink()
    if cache_path.exists() and args.force:
        cache_path.unlink()

    if args.skip_raw_integrity_check and not args.reuse_raw:
        raise SystemExit("--skip-raw-integrity-check requires --reuse-raw.")

    if args.reuse_raw:
        table_counts, delete_count, source_rows, raw_integrity = audit_existing_raw(
            raw_path,
            rows,
            args.zip_dir,
            skip_integrity_check=args.skip_raw_integrity_check,
        )
    else:
        con = sqlite3.connect(raw_path)
        try:
            con.execute("pragma journal_mode=off")
            con.execute("pragma synchronous=off")
            con.execute("pragma temp_store=file")
            create_source_schema(con)
            table_counts, delete_count, source_rows = load_sources(con, rows, args.zip_dir)
            create_raw_indexes(con)
            con.commit()
            raw_integrity = con.execute("pragma integrity_check").fetchone()[0]
        finally:
            con.close()
    cache_audit = build_retained_cache(raw_path, cache_path, force=args.force)
    write_csv(source_csv, source_rows)
    output = {
        "status": "complete" if raw_integrity == "ok" and cache_audit["cache_integrity_check"] == "ok" else "failed_integrity_check",
        "created_at_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "software": {"python": sys.version.split()[0], "sqlite": sqlite3.sqlite_version},
        "source_manifest": str(args.manifest),
        "zip_dir": str(args.zip_dir),
        "raw_database": str(raw_path),
        "latest_cache": str(cache_path),
        "quarters_expected": len(rows),
        "raw_table_rows": table_counts,
        "delete_ledger_rows": delete_count,
        "deduplication": "Exclude DELETE identifiers matched to CASEID, then retain latest valid FDA_DT; ties use numeric CASEVERSION then numeric PRIMARYID. Blank CASEID values are retained by PRIMARYID.",
        "raw_integrity_check": raw_integrity,
        "raw_integrity_check_method": (
            "reused_from_source_build" if args.reuse_raw and args.skip_raw_integrity_check
            else "full_integrity_check"
        ),
        "cache_audit": cache_audit,
        "source_zip_manifest_sha256": sha256(source_csv),
    }
    manifest_path.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))
    return 0 if output["status"] == "complete" else 1


if __name__ == "__main__":
    raise SystemExit(main())
