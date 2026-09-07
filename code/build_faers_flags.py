from __future__ import annotations

import csv
import itertools
import math
import os
import sqlite3
import time
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
WORK = REPO_ROOT / "work"
CONFIG = REPO_ROOT / "config"
OUT = Path(os.environ.get("PV_AGGREGATE_OUTPUT_DIR", REPO_ROOT / "aggregate_outputs"))
WORK.mkdir(parents=True, exist_ok=True)
OUT.mkdir(parents=True, exist_ok=True)

RAW_FAERS = Path(os.environ.get("PV_FAERS_RAW_DB", REPO_ROOT / "external_data" / "faers.sqlite"))
LATEST_CACHE = Path(os.environ.get("PV_FAERS_CACHE_DB", REPO_ROOT / "external_data" / "faers_cache.sqlite"))
FAERS_OUT = Path(os.environ.get("PV_FAERS_FLAGS_DB", WORK / "faers_periop_hypersensitivity_prelim.sqlite"))

EVENT_GROUPS = ["anaphylaxis_core", "periop_hypersensitivity_broad", "kounis_takotsubo_probe"]
PRIMARY_CLASSES = ["NMBA", "antibiotic_anchor", "chlorhexidine", "dye_anchor"]
BACKGROUND_CLASSES = [
    "sedative_hypnotic",
    "opioid",
    "volatile_anesthetic",
    "local_anesthetic",
    "vasoactive_rescue",
]
DRUG_CLASSES = PRIMARY_CLASSES + BACKGROUND_CLASSES
ROLE_SCOPES = ["suspect", "any"]
BACKGROUNDS = {
    "whole_database": "1 = 1",
    "anesthesia_background": "anesthesia_background_any = 1",
    "broad_periop_background": "broad_periop_background_any = 1",
}


def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def read_config(name: str) -> list[dict[str, str]]:
    with (CONFIG / name).open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def q(name: str) -> str:
    return '"' + name.replace('"', '""') + '"'


def upper_prefix(alias: str) -> str:
    return alias.strip().upper() + "*"


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
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


def metric_2x2(a: int, b: int, c: int, d: int) -> dict[str, object]:
    n = a + b + c + d
    out: dict[str, object] = {"a": a, "b": b, "c": c, "d": d, "n": n}
    if min(a, b, c, d) > 0:
        ror = (a * d) / (b * c)
        se = math.sqrt((1 / a) + (1 / b) + (1 / c) + (1 / d))
    else:
        aa, bb, cc, dd = a + 0.5, b + 0.5, c + 0.5, d + 0.5
        ror = (aa * dd) / (bb * cc)
        se = math.sqrt((1 / aa) + (1 / bb) + (1 / cc) + (1 / dd))
        out["zero_cell_correction"] = "haldane_0.5"
    out["ror"] = ror
    out["ror025"] = math.exp(math.log(ror) - 1.96 * se)
    out["ror975"] = math.exp(math.log(ror) + 1.96 * se)

    exposed_total = a + b
    unexposed_total = c + d
    event_total = a + c
    prr = None
    if exposed_total and unexposed_total and c:
        prr = (a / exposed_total) / (c / unexposed_total)
    out["prr"] = prr
    denom = exposed_total * unexposed_total * event_total * (b + d)
    out["chi_square"] = ((a * d - b * c) ** 2 * n / denom) if denom else None

    if n and exposed_total and event_total:
        aa = a + 0.5
        expected = (exposed_total * event_total) / n
        ic = math.log2(aa / (expected + 0.5))
        variance_ln = max((1 / aa) - (1 / exposed_total) - (1 / event_total) + (1 / (n + 0.5)), 0)
        out["ic_approx"] = ic
        out["ic025_approx"] = ic - 1.96 * math.sqrt(variance_ln) / math.log(2)
    out["positive_composite_prelim"] = int(
        a >= 3
        and float(out["ror025"]) > 1
        and prr is not None
        and prr >= 2
        and out["chi_square"] is not None
        and float(out["chi_square"]) >= 4
        and out.get("ic025_approx") is not None
        and float(out["ic025_approx"]) > 0
    )
    return out


def setup_db(rebuild: bool) -> sqlite3.Connection:
    if rebuild and FAERS_OUT.exists():
        FAERS_OUT.unlink()
    con = sqlite3.connect(FAERS_OUT)
    con.execute("pragma journal_mode=off")
    con.execute("pragma synchronous=off")
    con.execute("pragma temp_store=memory")
    con.execute("attach database ? as src", (str(RAW_FAERS),))
    con.execute("attach database ? as ld", (str(LATEST_CACHE),))
    con.executescript(
        """
        create table if not exists event_hits (
            primaryid text not null,
            event_group text not null,
            primary key (primaryid, event_group)
        );

        create table if not exists class_hits (
            primaryid text not null,
            drug_class text not null,
            generic_or_group text not null,
            any_role integer not null default 0,
            suspect_role integer not null default 0,
            primary key (primaryid, drug_class, generic_or_group)
        );
        """
    )
    return con


def build_events(con: sqlite3.Connection, event_rows: list[dict[str, str]]) -> None:
    con.execute("delete from event_hits")
    con.commit()
    seen_terms: set[tuple[str, str]] = set()
    for row in event_rows:
        group = row["event_group"]
        term = row["faers_canada_pt"]
        if not term or (group, term) in seen_terms:
            continue
        seen_terms.add((group, term))
        log(f"FAERS event: {group} / {term}")
        con.execute(
            """
            insert or ignore into event_hits (primaryid, event_group)
            select r.primaryid, ?
            from src.reac r indexed by idx_reac_pt_name
            join ld.latest_demo x on x.primaryid = r.primaryid
                                and x.source_year = r.source_year
                                and x.source_quarter = r.source_quarter
            where r.pt_name = ?
            """,
            (group, term),
        )
        con.commit()
    con.execute("create index if not exists idx_faers_event_hits_group on event_hits(event_group)")
    con.commit()


def faers_alias_rows(alias_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    out = []
    for row in alias_rows:
        if row["database_scope"] in {"FAERS", "FAERS_Canada", "All"}:
            out.append(row)
    return out


def build_classes(con: sqlite3.Connection, alias_rows: list[dict[str, str]]) -> None:
    con.execute("delete from class_hits")
    con.commit()
    rows = faers_alias_rows(alias_rows)
    for row in rows:
        alias = row["alias"]
        if not alias:
            continue
        pattern = upper_prefix(alias)
        for column, index_name in [("drugname", "idx_drug_drugname"), ("prod_ai", "idx_drug_prod_ai")]:
            log(f"FAERS drug: {row['drug_class']} / {row['generic_or_group']} / {alias} / {column}")
            con.execute(
                f"""
                insert into class_hits (primaryid, drug_class, generic_or_group, any_role, suspect_role)
                select d.primaryid, ?, ?, 1,
                       case when upper(coalesce(d.role_cod, '')) in ('PS', 'SS') then 1 else 0 end
                from src.drug d indexed by {index_name}
                join ld.latest_demo x on x.primaryid = d.primaryid
                                    and x.source_year = d.source_year
                                    and x.source_quarter = d.source_quarter
                where d.{column} glob ?
                on conflict(primaryid, drug_class, generic_or_group) do update set
                    any_role = max(any_role, excluded.any_role),
                    suspect_role = max(suspect_role, excluded.suspect_role)
                """,
                (row["drug_class"], row["generic_or_group"], pattern),
            )
            con.commit()
    con.execute("create index if not exists idx_faers_class_hits_class on class_hits(drug_class)")
    con.execute("create index if not exists idx_faers_class_hits_pid on class_hits(primaryid)")
    con.commit()


def build_aggregates(con: sqlite3.Connection) -> None:
    log("FAERS: building aggregate tables")
    con.executescript(
        """
        drop table if exists event_agg;
        drop table if exists report_class_flags;
        drop table if exists class_agg;
        drop table if exists context_flags;

        create table report_class_flags as
        select primaryid, drug_class, max(any_role) as any_role, max(suspect_role) as suspect_role
        from class_hits
        group by primaryid, drug_class;
        create index idx_faers_report_class_flags_pid on report_class_flags(primaryid);
        """
    )
    event_expr = ", ".join(
        [f"max(case when event_group = '{group}' then 1 else 0 end) as event_{group}" for group in EVENT_GROUPS]
    )
    con.execute(
        f"""
        create table event_agg as
        select primaryid, {event_expr}
        from event_hits
        group by primaryid
        """
    )
    con.execute("create index idx_faers_event_agg_pid on event_agg(primaryid)")

    class_exprs = []
    for cls in DRUG_CLASSES:
        class_exprs.append(f"max(case when drug_class = '{cls}' and any_role = 1 then 1 else 0 end) as any_{cls}")
        class_exprs.append(
            f"max(case when drug_class = '{cls}' and suspect_role = 1 then 1 else 0 end) as suspect_{cls}"
        )
    con.execute(
        f"""
        create table class_agg as
        select primaryid, {", ".join(class_exprs)}
        from report_class_flags
        group by primaryid
        """
    )
    con.execute("create index idx_faers_class_agg_pid on class_agg(primaryid)")

    any_count = " + ".join([f"coalesce(any_{cls}, 0)" for cls in PRIMARY_CLASSES])
    suspect_count = " + ".join([f"coalesce(suspect_{cls}, 0)" for cls in PRIMARY_CLASSES])
    anesthesia_terms = ["NMBA", "chlorhexidine", "dye_anchor", *BACKGROUND_CLASSES]
    anesthesia_background = "max(" + ", ".join([f"coalesce(any_{cls}, 0)" for cls in anesthesia_terms]) + ")"
    broad_background = f"max({anesthesia_background}, coalesce(any_antibiotic_anchor, 0))"
    con.execute(
        f"""
        create table context_flags as
        select *,
               ({any_count}) as any_primary_culprit_class_count,
               ({suspect_count}) as suspect_primary_culprit_class_count,
               {anesthesia_background} as anesthesia_background_any,
               {broad_background} as broad_periop_background_any
        from class_agg
        """
    )
    con.execute("create index idx_faers_context_flags_pid on context_flags(primaryid)")
    con.commit()


def scalar(con: sqlite3.Connection, sql: str, params: tuple = ()) -> int:
    return int(con.execute(sql, params).fetchone()[0] or 0)


def event_total(con: sqlite3.Connection, event_group: str, background: str) -> int:
    event_col = f"event_{event_group}"
    if background == "whole_database":
        return scalar(con, f"select count(*) from event_agg where {q(event_col)} = 1")
    context_col = "anesthesia_background_any" if background == "anesthesia_background" else "broad_periop_background_any"
    return scalar(
        con,
        f"""
        select count(*)
        from event_agg e
        join context_flags c on c.primaryid = e.primaryid
        where e.{q(event_col)} = 1 and c.{q(context_col)} = 1
        """,
    )


def exposure_total(con: sqlite3.Connection, role_scope: str, drug_class: str, background: str) -> int:
    exposure_col = f"{role_scope}_{drug_class}"
    if background == "whole_database":
        return scalar(con, f"select count(*) from context_flags where {q(exposure_col)} = 1")
    context_col = "anesthesia_background_any" if background == "anesthesia_background" else "broad_periop_background_any"
    return scalar(
        con,
        f"select count(*) from context_flags where {q(exposure_col)} = 1 and {q(context_col)} = 1",
    )


def background_total(con: sqlite3.Connection, background: str) -> int:
    if background == "whole_database":
        return scalar(con, "select count(*) from ld.latest_demo")
    context_col = "anesthesia_background_any" if background == "anesthesia_background" else "broad_periop_background_any"
    return scalar(con, f"select count(*) from context_flags where {q(context_col)} = 1")


def exposed_event_total(con: sqlite3.Connection, event_group: str, role_scope: str, drug_class: str, background: str) -> int:
    event_col = f"event_{event_group}"
    exposure_col = f"{role_scope}_{drug_class}"
    where = [f"e.{q(event_col)} = 1", f"c.{q(exposure_col)} = 1"]
    if background != "whole_database":
        context_col = "anesthesia_background_any" if background == "anesthesia_background" else "broad_periop_background_any"
        where.append(f"c.{q(context_col)} = 1")
    return scalar(
        con,
        f"""
        select count(*)
        from event_agg e
        join context_flags c on c.primaryid = e.primaryid
        where {" and ".join(where)}
        """,
    )


def export_summaries(con: sqlite3.Connection) -> None:
    log("FAERS: exporting summaries")
    database = "FAERS"
    db_summary = [
        {
            "database": database,
            "report_n": background_total(con, "whole_database"),
            "anaphylaxis_core_n": event_total(con, "anaphylaxis_core", "whole_database"),
            "broad_hypersensitivity_n": event_total(con, "periop_hypersensitivity_broad", "whole_database"),
            "kounis_takotsubo_n": event_total(con, "kounis_takotsubo_probe", "whole_database"),
            "anesthesia_background_n": background_total(con, "anesthesia_background"),
            "broad_periop_background_n": background_total(con, "broad_periop_background"),
        }
    ]
    write_csv(OUT / "faers_report_level_database_summary.csv", db_summary)

    attribution_rows = []
    for event_group in EVENT_GROUPS:
        event_n = event_total(con, event_group, "whole_database")
        for cls in PRIMARY_CLASSES:
            attribution_rows.append(
                {
                    "database": database,
                    "event_group": event_group,
                    "event_n": event_n,
                    "drug_class": cls,
                    "event_with_class_any_n": exposed_event_total(con, event_group, "any", cls, "whole_database"),
                    "event_with_class_suspect_n": exposed_event_total(con, event_group, "suspect", cls, "whole_database"),
                }
            )
    write_csv(OUT / "faers_attribution_counts.csv", attribution_rows)

    coexposure_rows = []
    for count_value in range(len(PRIMARY_CLASSES) + 1):
        coexposure_rows.append(
            {
                "database": database,
                "event_group": "anaphylaxis_core",
                "culprit_class_count": count_value,
                "report_n": scalar(
                    con,
                    """
                    select count(*)
                    from event_agg e
                    left join context_flags c on c.primaryid = e.primaryid
                    where e.event_anaphylaxis_core = 1
                      and coalesce(c.suspect_primary_culprit_class_count, 0) = ?
                    """,
                    (count_value,),
                ),
            }
        )
    for a, b in itertools.combinations(PRIMARY_CLASSES, 2):
        coexposure_rows.append(
            {
                "database": database,
                "event_group": "anaphylaxis_core",
                "pair": f"{a}__{b}",
                "any_pair_n": scalar(
                    con,
                    f"""
                    select count(*)
                    from event_agg e join context_flags c on c.primaryid = e.primaryid
                    where e.event_anaphylaxis_core = 1
                      and c.{q('any_' + a)} = 1 and c.{q('any_' + b)} = 1
                    """,
                ),
                "suspect_pair_n": scalar(
                    con,
                    f"""
                    select count(*)
                    from event_agg e join context_flags c on c.primaryid = e.primaryid
                    where e.event_anaphylaxis_core = 1
                      and c.{q('suspect_' + a)} = 1 and c.{q('suspect_' + b)} = 1
                    """,
                ),
            }
        )
    write_csv(OUT / "faers_coexposure_counts.csv", coexposure_rows)


def run_da(con: sqlite3.Connection) -> None:
    log("FAERS: running preliminary DA")
    rows = []
    for event_group in EVENT_GROUPS:
        for role_scope in ROLE_SCOPES:
            for drug_class in PRIMARY_CLASSES:
                for background in BACKGROUNDS:
                    n = background_total(con, background)
                    e = event_total(con, event_group, background)
                    x = exposure_total(con, role_scope, drug_class, background)
                    a = exposed_event_total(con, event_group, role_scope, drug_class, background)
                    b = x - a
                    c = e - a
                    d = n - a - b - c
                    row: dict[str, object] = {
                        "database": "FAERS",
                        "event_group": event_group,
                        "role_scope": role_scope,
                        "drug_class": drug_class,
                        "background": background,
                    }
                    row.update(metric_2x2(a, b, c, d))
                    rows.append(row)
    write_csv(OUT / "faers_preliminary_da.csv", rows)


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--rebuild", action="store_true")
    args = parser.parse_args()

    event_rows = read_config("event_pt_groups.csv")
    alias_rows = read_config("drug_class_aliases.csv")
    con = setup_db(args.rebuild)
    build_events(con, event_rows)
    build_classes(con, alias_rows)
    build_aggregates(con)
    export_summaries(con)
    run_da(con)
    con.close()
    log(f"Done: {FAERS_OUT}")


if __name__ == "__main__":
    main()
