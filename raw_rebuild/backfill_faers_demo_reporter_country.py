#!/usr/bin/env python3
"""Backfill `REPORTER_COUNTRY` from the retained official DEMO source ZIPs.

This is a deterministic schema repair for the first v0.26 raw build, which
already imported every DEMO row but omitted a field required by the prespecified
country fallback audit.  The final source builder imports this column directly;
the backfill leaves a separate audit trail for this execution.
"""

from __future__ import annotations

import argparse
import csv
import json
import sqlite3
import zipfile
from datetime import datetime, timezone
from io import TextIOWrapper
from pathlib import Path

from build_faers_raw_branch import MEMBER_RE, normalise_headers


ROOT = Path(__file__).resolve().parent
DEFAULT_MANIFEST = ROOT / "faers_ascii_manifest_2004q1_2025q4.csv"
DEFAULT_ZIPS = ROOT.parent / "raw_inputs" / "faers_ascii_2004q1_2025q4"
DEFAULT_RAW = ROOT.parent / "work" / "faers_raw_2004q1_2025q4" / "faers_raw_2004q1_2025q4.sqlite"
DEFAULT_OUT = ROOT.parent / "run_logs" / "faers_reporter_country_backfill_v0_26.json"


def demo_member(names: list[str]) -> str:
    choices = [name for name in names if MEMBER_RE.search(name.upper()) and "DEMO" in Path(name).name.upper()]
    if len(choices) != 1:
        raise ValueError(f"Expected one DEMO member, found: {choices}")
    return choices[0]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--zip-dir", type=Path, default=DEFAULT_ZIPS)
    parser.add_argument("--raw-db", type=Path, default=DEFAULT_RAW)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    with args.manifest.open(encoding="utf-8-sig", newline="") as handle:
        sources = list(csv.DictReader(handle))
    con = sqlite3.connect(args.raw_db)
    try:
        columns = {row[1] for row in con.execute("pragma table_info(demo)")}
        if "reporter_country" not in columns:
            con.execute("alter table demo add column reporter_country text")
        con.execute("drop table if exists reporter_country_stage")
        con.execute("create table reporter_country_stage (primaryid text not null, source_year text not null, source_quarter text not null, reporter_country text)")
        total = 0
        per_quarter: list[dict[str, object]] = []
        for source in sources:
            path = args.zip_dir / source["filename"]
            with zipfile.ZipFile(path) as archive:
                member = demo_member(archive.namelist())
                with archive.open(member) as binary:
                    reader = csv.reader(TextIOWrapper(binary, encoding="latin-1", newline=""), delimiter="$")
                    headers = normalise_headers(next(reader))
                    primary_index = headers.index("primaryid")
                    country_index = headers.index("reporter_country")
                    batch: list[tuple[str, str, str, str | None]] = []
                    count = 0
                    for row in reader:
                        primaryid = row[primary_index].strip() if primary_index < len(row) else ""
                        if not primaryid:
                            continue
                        country = row[country_index].strip() if country_index < len(row) else ""
                        batch.append((primaryid, source["year"], source["quarter"], country or None))
                        count += 1
                        if len(batch) >= 20_000:
                            con.executemany("insert into reporter_country_stage values (?, ?, ?, ?)", batch)
                            batch.clear()
                    if batch:
                        con.executemany("insert into reporter_country_stage values (?, ?, ?, ?)", batch)
            con.commit()
            total += count
            per_quarter.append({"year": source["year"], "quarter": source["quarter"], "demo_rows_staged": count})
        con.execute("create index idx_reporter_country_stage_key on reporter_country_stage(primaryid, source_year, source_quarter)")
        con.execute(
            """
            update demo
            set reporter_country = (
                select s.reporter_country from reporter_country_stage s
                where s.primaryid = demo.primaryid
                  and s.source_year = demo.source_year
                  and s.source_quarter = demo.source_quarter
            )
            where exists (
                select 1 from reporter_country_stage s
                where s.primaryid = demo.primaryid
                  and s.source_year = demo.source_year
                  and s.source_quarter = demo.source_quarter
            )
            """
        )
        filled = int(con.execute("select count(*) from demo where trim(coalesce(reporter_country,'')) <> ''").fetchone()[0])
        demo_total = int(con.execute("select count(*) from demo").fetchone()[0])
        missing_stage = int(
            con.execute(
                """select count(*) from demo d where not exists (
                    select 1 from reporter_country_stage s
                    where s.primaryid=d.primaryid and s.source_year=d.source_year and s.source_quarter=d.source_quarter
                )"""
            ).fetchone()[0]
        )
        con.execute("drop table reporter_country_stage")
        con.commit()
        integrity = con.execute("pragma integrity_check").fetchone()[0]
    finally:
        con.close()
    result = {
        "status": "complete" if integrity == "ok" and missing_stage == 0 else "failed",
        "completed_at_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "raw_database": str(args.raw_db), "quarters": len(sources),
        "staged_demo_rows": total, "raw_demo_rows": demo_total,
        "rows_with_nonblank_reporter_country": filled,
        "demo_rows_without_stage_match": missing_stage,
        "integrity_check": integrity,
        "per_quarter": per_quarter,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0 if result["status"] == "complete" else 1


if __name__ == "__main__":
    raise SystemExit(main())
