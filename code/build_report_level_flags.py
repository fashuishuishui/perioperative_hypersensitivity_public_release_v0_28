from __future__ import annotations

import argparse
import csv
import itertools
import os
import re
import sqlite3
import time
import unicodedata
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
WORK = REPO_ROOT / "work"
CONFIG = REPO_ROOT / "config"
OUT = Path(os.environ.get("PV_AGGREGATE_OUTPUT_DIR", REPO_ROOT / "aggregate_outputs"))
WORK.mkdir(parents=True, exist_ok=True)
OUT.mkdir(parents=True, exist_ok=True)

CANADA_DIR = Path(os.environ.get("PV_CANADA_DIR", REPO_ROOT / "external_data" / "canada_20250731"))
JADER_DIR = Path(os.environ.get("PV_JADER_DIR", REPO_ROOT / "external_data" / "jader_202511"))
DB_PATH = Path(os.environ.get("PV_REPORT_FLAGS_DB", WORK / "report_level_flags_canada_jader.sqlite"))

EVENT_GROUPS = ["anaphylaxis_core", "periop_hypersensitivity_broad", "kounis_takotsubo_probe"]
PRIMARY_CULPRIT_CLASSES = ["NMBA", "antibiotic_anchor", "chlorhexidine", "dye_anchor"]
BACKGROUND_CLASSES = [
    "sedative_hypnotic",
    "opioid",
    "volatile_anesthetic",
    "local_anesthetic",
    "vasoactive_rescue",
]
DRUG_CLASSES = PRIMARY_CULPRIT_CLASSES + BACKGROUND_CLASSES


def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def norm(value: str | None) -> str:
    return unicodedata.normalize("NFKC", value or "").casefold().strip()


def read_config(name: str) -> list[dict[str, str]]:
    with (CONFIG / name).open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def read_dollar_file(path: Path):
    with path.open("r", encoding="utf-8", errors="replace", newline="") as f:
        yield from csv.reader(f, delimiter="$", quotechar='"')


def resolve_source_file(directory: Path, canonical_name: str) -> Path:
    """Resolve a verified local source file without silently changing inputs."""
    candidates = [
        directory / canonical_name,
        directory / f"{canonical_name}.baiduyun.p.downloading",
    ]
    for candidate in candidates:
        if candidate.is_file() and candidate.stat().st_size > 0:
            return candidate
    raise FileNotFoundError(
        f"Required source file {canonical_name} was not found under {directory}; "
        "a local download-suffix variant is accepted only when explicitly audited."
    )


def safe_get(row: list[str], index: int) -> str:
    return row[index] if index < len(row) else ""


def canonical_canada_report_id(value: str) -> str:
    value = (value or "").strip()
    stripped = value.lstrip("0")
    return stripped or value


def alias_matches(raw_text: str, alias: str, database: str) -> bool:
    text = norm(raw_text)
    needle = norm(alias)
    if not needle:
        return False
    if database in {"FAERS", "Canada"}:
        pattern = rf"(?<![a-z0-9]){re.escape(needle)}(?![a-z0-9])"
        return re.search(pattern, text) is not None
    return needle in text


def event_term_maps(event_rows: list[dict[str, str]], database: str) -> dict[str, set[str]]:
    column = "faers_canada_pt" if database == "Canada" else "jader_pt_ja"
    by_term: dict[str, set[str]] = {}
    for row in event_rows:
        term = row[column]
        if not term:
            continue
        by_term.setdefault(term, set()).add(row["event_group"])
    return by_term


def aliases_for(alias_rows: list[dict[str, str]], database: str) -> list[dict[str, str]]:
    out = []
    for row in alias_rows:
        scope = row["database_scope"]
        if scope == database or scope == "All" or (database in {"FAERS", "Canada"} and scope == "FAERS_Canada"):
            out.append(row)
    return out


def open_output_db(rebuild: bool) -> sqlite3.Connection:
    if rebuild and DB_PATH.exists():
        DB_PATH.unlink()
    con = sqlite3.connect(DB_PATH)
    con.execute("pragma journal_mode=off")
    con.execute("pragma synchronous=off")
    con.execute("pragma temp_store=memory")
    con.executescript(
        """
        create table if not exists reports (
            database text not null,
            report_id text not null,
            report_round text,
            sex text,
            age text,
            age_unit text,
            year_quarter text,
            serious text,
            report_type text,
            reporter_qualification text,
            e2b text,
            primary key (database, report_id)
        );

        create table if not exists event_hits (
            database text not null,
            report_id text not null,
            event_group text not null,
            primary key (database, report_id, event_group)
        );

        create table if not exists class_hits (
            database text not null,
            report_id text not null,
            drug_class text not null,
            generic_or_group text not null,
            any_role integer not null default 0,
            suspect_role integer not null default 0,
            primary key (database, report_id, drug_class, generic_or_group)
        );
        """
    )
    return con


def insert_many(con: sqlite3.Connection, sql: str, rows, chunk_size: int = 10000) -> int:
    count = 0
    batch = []
    for row in rows:
        batch.append(row)
        if len(batch) >= chunk_size:
            con.executemany(sql, batch)
            count += len(batch)
            batch.clear()
    if batch:
        con.executemany(sql, batch)
        count += len(batch)
    con.commit()
    return count


def build_canada(con: sqlite3.Connection, event_rows: list[dict[str, str]], alias_rows: list[dict[str, str]]) -> None:
    database = "Canada"
    reports_file = resolve_source_file(CANADA_DIR, "reports.txt")
    reactions_file = resolve_source_file(CANADA_DIR, "reactions.txt")
    report_drug_file = resolve_source_file(CANADA_DIR, "report_drug.txt")
    log("Canada: inserting report metadata")

    def report_iter():
        for row in read_dollar_file(reports_file):
            report_id = canonical_canada_report_id(safe_get(row, 0))
            if not report_id:
                continue
            yield (
                database,
                report_id,
                safe_get(row, 2),
                safe_get(row, 10),
                safe_get(row, 12),
                safe_get(row, 14),
                "",
                safe_get(row, 26),
                safe_get(row, 7),
                safe_get(row, 34),
                "",
            )

    n_reports = insert_many(
        con,
        """
        insert or ignore into reports
        (database, report_id, report_round, sex, age, age_unit, year_quarter, serious, report_type, reporter_qualification, e2b)
        values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        report_iter(),
    )
    log(f"Canada: reports inserted/seen = {n_reports}")

    term_to_groups = event_term_maps(event_rows, database)
    log("Canada: inserting event hits")

    def event_iter():
        for row in read_dollar_file(reactions_file):
            report_id = canonical_canada_report_id(safe_get(row, 1))
            pt = safe_get(row, 5)
            for group in term_to_groups.get(pt, ()):
                yield (database, report_id, group)

    n_events = insert_many(
        con,
        "insert or ignore into event_hits (database, report_id, event_group) values (?, ?, ?)",
        event_iter(),
    )
    log(f"Canada: event-hit rows inserted/seen = {n_events}")

    aliases = aliases_for(alias_rows, database)
    log("Canada: inserting drug-class hits")

    def class_iter():
        for row in read_dollar_file(report_drug_file):
            report_id = canonical_canada_report_id(safe_get(row, 1))
            drug_name = safe_get(row, 3)
            role = safe_get(row, 4)
            is_suspect = 1 if "suspect" in norm(role) else 0
            for cfg in aliases:
                if alias_matches(drug_name, cfg["alias"], database):
                    yield (
                        database,
                        report_id,
                        cfg["drug_class"],
                        cfg["generic_or_group"],
                        1,
                        is_suspect,
                    )

    n_classes = insert_many(
        con,
        """
        insert into class_hits (database, report_id, drug_class, generic_or_group, any_role, suspect_role)
        values (?, ?, ?, ?, ?, ?)
        on conflict(database, report_id, drug_class, generic_or_group) do update set
            any_role = max(any_role, excluded.any_role),
            suspect_role = max(suspect_role, excluded.suspect_role)
        """,
        class_iter(),
    )
    log(f"Canada: class-hit rows inserted/seen = {n_classes}")


def latest_jader_demo() -> dict[str, dict[str, str]]:
    latest: dict[str, dict[str, str]] = {}
    with resolve_source_file(JADER_DIR, "demo202511.csv").open(
        "r", encoding="cp932", errors="replace", newline=""
    ) as f:
        reader = csv.DictReader(f)
        for row in reader:
            report_id = row.get("識別番号", "")
            round_value = row.get("報告回数", "")
            if not report_id:
                continue
            current = latest.get(report_id)
            if current is None or int(round_value or "0") > int(current.get("報告回数") or "0"):
                latest[report_id] = row
    return latest


def build_jader(con: sqlite3.Connection, event_rows: list[dict[str, str]], alias_rows: list[dict[str, str]]) -> None:
    database = "JADER"
    reactions_file = resolve_source_file(JADER_DIR, "reac202511.csv")
    drug_file = resolve_source_file(JADER_DIR, "drug202511.csv")
    log("JADER: reading latest demo rows")
    latest = latest_jader_demo()
    latest_keys = {(rid, row.get("報告回数", "")) for rid, row in latest.items()}
    log(f"JADER: latest report units = {len(latest):,}")

    def report_iter():
        for report_id, row in latest.items():
            yield (
                database,
                report_id,
                row.get("報告回数", ""),
                row.get("性別", ""),
                row.get("年齢", ""),
                "",
                row.get("報告年度・四半期", ""),
                "",
                row.get("報告の種類", ""),
                row.get("報告者の資格", ""),
                row.get("E2B", ""),
            )

    n_reports = insert_many(
        con,
        """
        insert or ignore into reports
        (database, report_id, report_round, sex, age, age_unit, year_quarter, serious, report_type, reporter_qualification, e2b)
        values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        report_iter(),
    )
    log(f"JADER: reports inserted/seen = {n_reports}")

    term_to_groups = event_term_maps(event_rows, database)
    log("JADER: inserting event hits for latest rounds")

    def event_iter():
        with reactions_file.open(
            "r", encoding="cp932", errors="replace", newline=""
        ) as f:
            reader = csv.DictReader(f)
            for row in reader:
                report_id = row.get("識別番号", "")
                report_round = row.get("報告回数", "")
                if (report_id, report_round) not in latest_keys:
                    continue
                pt = row.get("有害事象", "")
                for group in term_to_groups.get(pt, ()):
                    yield (database, report_id, group)

    n_events = insert_many(
        con,
        "insert or ignore into event_hits (database, report_id, event_group) values (?, ?, ?)",
        event_iter(),
    )
    log(f"JADER: event-hit rows inserted/seen = {n_events}")

    aliases = aliases_for(alias_rows, database)
    log("JADER: inserting drug-class hits for latest rounds")

    def class_iter():
        with drug_file.open(
            "r", encoding="cp932", errors="replace", newline=""
        ) as f:
            reader = csv.DictReader(f)
            for row in reader:
                report_id = row.get("識別番号", "")
                report_round = row.get("報告回数", "")
                if (report_id, report_round) not in latest_keys:
                    continue
                generic = row.get("医薬品（一般名）", "")
                brand = row.get("医薬品（販売名）", "")
                role = row.get("医薬品の関与", "")
                text = f"{generic} {brand}"
                is_suspect = 1 if "被疑薬" in role else 0
                for cfg in aliases:
                    if alias_matches(text, cfg["alias"], database):
                        yield (
                            database,
                            report_id,
                            cfg["drug_class"],
                            cfg["generic_or_group"],
                            1,
                            is_suspect,
                        )

    n_classes = insert_many(
        con,
        """
        insert into class_hits (database, report_id, drug_class, generic_or_group, any_role, suspect_role)
        values (?, ?, ?, ?, ?, ?)
        on conflict(database, report_id, drug_class, generic_or_group) do update set
            any_role = max(any_role, excluded.any_role),
            suspect_role = max(suspect_role, excluded.suspect_role)
        """,
        class_iter(),
    )
    log(f"JADER: class-hit rows inserted/seen = {n_classes}")


def q(name: str) -> str:
    return '"' + name.replace('"', '""') + '"'


def build_aggregates(con: sqlite3.Connection) -> None:
    log("Building aggregate and wide flag tables")
    con.executescript(
        """
        drop table if exists event_agg;
        drop table if exists report_class_flags;
        drop table if exists class_agg;
        drop table if exists report_level_flags;

        create table report_class_flags as
        select database, report_id, drug_class,
               max(any_role) as any_role,
               max(suspect_role) as suspect_role
        from class_hits
        group by database, report_id, drug_class;

        create index if not exists idx_reports_db_id on reports(database, report_id);
        create index if not exists idx_event_hits_db_id on event_hits(database, report_id);
        create index if not exists idx_class_hits_db_id on class_hits(database, report_id);
        create index if not exists idx_report_class_flags_db_id on report_class_flags(database, report_id);
        """
    )

    event_exprs = [
        f"max(case when event_group = '{group}' then 1 else 0 end) as event_{group}"
        for group in EVENT_GROUPS
    ]
    con.execute(
        f"""
        create table event_agg as
        select database, report_id, {", ".join(event_exprs)}
        from event_hits
        group by database, report_id
        """
    )
    con.execute("create index if not exists idx_event_agg_db_id on event_agg(database, report_id)")

    class_exprs: list[str] = []
    for cls in DRUG_CLASSES:
        class_exprs.append(f"max(case when drug_class = '{cls}' and any_role = 1 then 1 else 0 end) as any_{cls}")
        class_exprs.append(
            f"max(case when drug_class = '{cls}' and suspect_role = 1 then 1 else 0 end) as suspect_{cls}"
        )
    con.execute(
        f"""
        create table class_agg as
        select database, report_id, {", ".join(class_exprs)}
        from report_class_flags
        group by database, report_id
        """
    )
    con.execute("create index if not exists idx_class_agg_db_id on class_agg(database, report_id)")

    select_cols = [
        "r.database",
        "r.report_id",
        "r.report_round",
        "r.sex",
        "r.age",
        "r.age_unit",
        "r.year_quarter",
        "r.serious",
        "r.report_type",
        "r.reporter_qualification",
        "r.e2b",
    ]
    for group in EVENT_GROUPS:
        select_cols.append(f"coalesce(e.event_{group}, 0) as event_{group}")
    for cls in DRUG_CLASSES:
        select_cols.append(f"coalesce(c.any_{cls}, 0) as any_{cls}")
        select_cols.append(f"coalesce(c.suspect_{cls}, 0) as suspect_{cls}")

    any_count = " + ".join([f"coalesce(c.any_{cls}, 0)" for cls in PRIMARY_CULPRIT_CLASSES])
    suspect_count = " + ".join([f"coalesce(c.suspect_{cls}, 0)" for cls in PRIMARY_CULPRIT_CLASSES])
    anesthesia_terms = ["NMBA", "chlorhexidine", "dye_anchor", *BACKGROUND_CLASSES]
    anesthesia_background = "max(" + ", ".join([f"coalesce(c.any_{cls}, 0)" for cls in anesthesia_terms]) + ")"
    broad_background = f"max({anesthesia_background}, coalesce(c.any_antibiotic_anchor, 0))"
    select_cols.append(f"({any_count}) as any_primary_culprit_class_count")
    select_cols.append(f"({suspect_count}) as suspect_primary_culprit_class_count")
    select_cols.append(f"{anesthesia_background} as anesthesia_background_any")
    select_cols.append(f"{broad_background} as broad_periop_background_any")

    for a, b in itertools.combinations(PRIMARY_CULPRIT_CLASSES, 2):
        select_cols.append(f"(coalesce(c.any_{a}, 0) * coalesce(c.any_{b}, 0)) as any_pair_{a}__{b}")
        select_cols.append(
            f"(coalesce(c.suspect_{a}, 0) * coalesce(c.suspect_{b}, 0)) as suspect_pair_{a}__{b}"
        )

    con.execute(
        f"""
        create table report_level_flags as
        select {", ".join(select_cols)}
        from reports r
        left join event_agg e on e.database = r.database and e.report_id = r.report_id
        left join class_agg c on c.database = r.database and c.report_id = r.report_id
        """
    )
    con.execute("create index if not exists idx_report_level_flags_db_id on report_level_flags(database, report_id)")
    con.execute("create index if not exists idx_report_level_flags_db_event on report_level_flags(database, event_anaphylaxis_core)")
    con.commit()


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def query_dicts(con: sqlite3.Connection, sql: str, params: tuple = ()) -> list[dict[str, object]]:
    cur = con.execute(sql, params)
    names = [d[0] for d in cur.description]
    return [dict(zip(names, row)) for row in cur.fetchall()]


def export_summaries(con: sqlite3.Connection) -> None:
    log("Exporting summary CSV files")
    summary_rows: list[dict[str, object]] = []
    for database in ["Canada", "JADER"]:
        base = query_dicts(
            con,
            """
            select database, count(*) as report_n,
                   sum(event_anaphylaxis_core) as anaphylaxis_core_n,
                   sum(event_periop_hypersensitivity_broad) as broad_hypersensitivity_n,
                   sum(event_kounis_takotsubo_probe) as kounis_takotsubo_n,
                   sum(anesthesia_background_any) as anesthesia_background_n,
                   sum(broad_periop_background_any) as broad_periop_background_n
            from report_level_flags
            where database = ?
            group by database
            """,
            (database,),
        )
        summary_rows.extend(base)
    write_rows(OUT / "report_level_database_summary_canada_jader.csv", summary_rows)

    attribution_rows: list[dict[str, object]] = []
    for database in ["Canada", "JADER"]:
        for event_group in EVENT_GROUPS:
            event_col = f"event_{event_group}"
            event_n = con.execute(
                f"select sum({q(event_col)}) from report_level_flags where database = ?", (database,)
            ).fetchone()[0]
            for cls in PRIMARY_CULPRIT_CLASSES:
                row = con.execute(
                    f"""
                    select
                        sum(case when {q(event_col)} = 1 and {q('any_' + cls)} = 1 then 1 else 0 end),
                        sum(case when {q(event_col)} = 1 and {q('suspect_' + cls)} = 1 then 1 else 0 end)
                    from report_level_flags
                    where database = ?
                    """,
                    (database,),
                ).fetchone()
                attribution_rows.append(
                    {
                        "database": database,
                        "event_group": event_group,
                        "event_n": event_n or 0,
                        "drug_class": cls,
                        "event_with_class_any_n": row[0] or 0,
                        "event_with_class_suspect_n": row[1] or 0,
                    }
                )
    write_rows(OUT / "attribution_counts_canada_jader.csv", attribution_rows)

    coexposure_rows: list[dict[str, object]] = []
    for database in ["Canada", "JADER"]:
        coexposure_rows.extend(
            query_dicts(
                con,
                """
                select database, 'anaphylaxis_core' as event_group,
                       suspect_primary_culprit_class_count as culprit_class_count,
                       count(*) as report_n
                from report_level_flags
                where database = ? and event_anaphylaxis_core = 1
                group by database, suspect_primary_culprit_class_count
                order by database, suspect_primary_culprit_class_count
                """,
                (database,),
            )
        )
        for a, b in itertools.combinations(PRIMARY_CULPRIT_CLASSES, 2):
            row = con.execute(
                f"""
                select
                    sum(case when event_anaphylaxis_core = 1 and {q('any_pair_' + a + '__' + b)} = 1 then 1 else 0 end),
                    sum(case when event_anaphylaxis_core = 1 and {q('suspect_pair_' + a + '__' + b)} = 1 then 1 else 0 end)
                from report_level_flags
                where database = ?
                """,
                (database,),
            ).fetchone()
            coexposure_rows.append(
                {
                    "database": database,
                    "event_group": "anaphylaxis_core",
                    "pair": f"{a}__{b}",
                    "any_pair_n": row[0] or 0,
                    "suspect_pair_n": row[1] or 0,
                }
            )
    write_rows(OUT / "coexposure_counts_canada_jader.csv", coexposure_rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rebuild", action="store_true", help="Rebuild the SQLite cache from scratch.")
    parser.add_argument(
        "--databases",
        default="Canada,JADER",
        help="Comma-separated databases to build. Supported: Canada,JADER.",
    )
    args = parser.parse_args()

    selected = {x.strip() for x in args.databases.split(",") if x.strip()}
    event_rows = read_config("event_pt_groups.csv")
    alias_rows = read_config("drug_class_aliases.csv")
    con = open_output_db(args.rebuild)

    if "Canada" in selected:
        con.execute("delete from reports where database = 'Canada'")
        con.execute("delete from event_hits where database = 'Canada'")
        con.execute("delete from class_hits where database = 'Canada'")
        con.commit()
        build_canada(con, event_rows, alias_rows)
    if "JADER" in selected:
        con.execute("delete from reports where database = 'JADER'")
        con.execute("delete from event_hits where database = 'JADER'")
        con.execute("delete from class_hits where database = 'JADER'")
        con.commit()
        build_jader(con, event_rows, alias_rows)

    build_aggregates(con)
    export_summaries(con)
    con.close()
    log(f"Done: {DB_PATH}")


if __name__ == "__main__":
    main()
