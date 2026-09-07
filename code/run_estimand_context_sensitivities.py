from __future__ import annotations

"""Recompute the v0.14 estimand and context-background sensitivity.

The primary agent contrast is deliberately role-symmetric: within reports that
carry a culprit class at suspect role, the selected agent's suspect-role share
among core-anaphylaxis reports is compared with its share among non-core
reports. The latter excludes the outcome from the reporting-volume baseline.

The context sensitivity uses only non-culprit perioperative medication markers
(sedatives, opioids, volatile anesthetics, local anesthetics, vasoactive
rescue medicines). It therefore does not make eligibility depend on the target
culprit class. It is a sensitivity analysis, not a purified comparator.
"""

import argparse
import csv
import math
import os
import sqlite3
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
WORK = REPO_ROOT / "work"
OUT = Path(os.environ.get("PV_ANALYSIS_OUTPUT_DIR", REPO_ROOT / "aggregate_outputs"))
FAERS_DB = Path(os.environ.get("PV_FAERS_FLAGS_DB", WORK / "faers_periop_hypersensitivity_prelim.sqlite"))
CJ_DB = Path(os.environ.get("PV_CANADA_JADER_FLAGS_DB", WORK / "report_level_flags_canada_jader.sqlite"))

PRIMARY_AGENTS = (
    ("NMBA", "rocuronium"),
    ("NMBA", "succinylcholine"),
    ("antibiotic_anchor", "cefazolin"),
    ("antibiotic_anchor", "teicoplanin"),
)
PRIMARY_CLASSES = ("NMBA", "antibiotic_anchor", "chlorhexidine", "dye_anchor")
NONCULPRIT_CONTEXT_CLASSES = (
    "sedative_hypnotic",
    "opioid",
    "volatile_anesthetic",
    "local_anesthetic",
    "vasoactive_rescue",
)


def scalar(con: sqlite3.Connection, sql: str, params: tuple[object, ...] = ()) -> int:
    return int(con.execute(sql, params).fetchone()[0] or 0)


def wilson_interval(successes: int, total: int, z: float = 1.959963984540054) -> tuple[float | None, float | None]:
    if total == 0:
        return None, None
    p = successes / total
    denom = 1 + z * z / total
    centre = (p + z * z / (2 * total)) / denom
    half = z * math.sqrt((p * (1 - p) + z * z / (4 * total)) / total) / denom
    return centre - half, centre + half


def newcombe_wilson_difference_interval(
    p_core: float, core_low: float, core_high: float, p_noncore: float, noncore_low: float, noncore_high: float
) -> tuple[float, float]:
    """Newcombe method 10 (MOVER-Wilson) interval for p_core - p_noncore."""
    difference = p_core - p_noncore
    lower = difference - math.sqrt((p_core - core_low) ** 2 + (noncore_high - p_noncore) ** 2)
    upper = difference + math.sqrt((core_high - p_core) ** 2 + (p_noncore - noncore_low) ** 2)
    return lower, upper


def ror_with_ci(a: int, b: int, c: int, d: int) -> dict[str, float | str | None]:
    correction = "none"
    aa, bb, cc, dd = a, b, c, d
    if min(a, b, c, d) == 0:
        aa, bb, cc, dd = a + 0.5, b + 0.5, c + 0.5, d + 0.5
        correction = "Haldane-Anscombe 0.5"
    ror = (aa * dd) / (bb * cc)
    se = math.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    return {
        "ror": ror,
        "ror_ci_low": math.exp(math.log(ror) - 1.96 * se),
        "ror_ci_high": math.exp(math.log(ror) + 1.96 * se),
        "zero_cell_correction": correction,
    }


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def add_primary_agent_row(
    con: sqlite3.Connection,
    database: str,
    id_col: str,
    db_where: str,
    db_params: tuple[object, ...],
    drug_class: str,
    agent: str,
) -> dict[str, object]:
    raw_db_where = db_where.replace("h.", "")
    class_suspect_n = scalar(
        con,
        f"""select count(distinct {id_col}) from class_hits
            where {raw_db_where} drug_class=? and suspect_role=1""",
        db_params + (drug_class,),
    )
    core_class_suspect_n = scalar(
        con,
        f"""select count(distinct h.{id_col}) from class_hits h
            join event_agg e on {('e.primaryid=h.primaryid' if id_col == 'primaryid' else 'e.database=h.database and e.report_id=h.report_id')}
            where {db_where} h.drug_class=? and h.suspect_role=1 and e.event_anaphylaxis_core=1""",
        db_params + (drug_class,),
    )
    agent_core_n = scalar(
        con,
        f"""select count(distinct h.{id_col}) from class_hits h
            join event_agg e on {('e.primaryid=h.primaryid' if id_col == 'primaryid' else 'e.database=h.database and e.report_id=h.report_id')}
            where {db_where} h.drug_class=? and h.generic_or_group=?
              and h.suspect_role=1 and e.event_anaphylaxis_core=1""",
        db_params + (drug_class, agent),
    )
    agent_total_n = scalar(
        con,
        f"""select count(distinct {id_col}) from class_hits
            where {raw_db_where} drug_class=? and generic_or_group=? and suspect_role=1""",
        db_params + (drug_class, agent),
    )
    noncore_class_suspect_n = class_suspect_n - core_class_suspect_n
    agent_noncore_n = agent_total_n - agent_core_n

    # Rows are core versus non-core class-suspect reports. Agent presence is a
    # binary report-level attribute, so co-reported same-class agents remain
    # in the agent-present cell rather than being double-counted across cells.
    a = agent_core_n
    b = core_class_suspect_n - agent_core_n
    c = agent_noncore_n
    d = noncore_class_suspect_n - agent_noncore_n
    p_core = a / core_class_suspect_n if core_class_suspect_n else None
    p_noncore = c / noncore_class_suspect_n if noncore_class_suspect_n else None
    core_low, core_high = wilson_interval(a, core_class_suspect_n)
    noncore_low, noncore_high = wilson_interval(c, noncore_class_suspect_n)
    sparse_agent_contrast = min(a, c) < 3
    # The percentage-point contrast remains estimable with a zero agent-event
    # cell. A continuity-corrected odds ratio can reverse direction in that
    # setting, so it is deliberately withheld rather than treated as evidence.
    pp_low, pp_high = newcombe_wilson_difference_interval(
        p_core, core_low, core_high, p_noncore, noncore_low, noncore_high
    ) if p_core is not None and core_low is not None and core_high is not None and p_noncore is not None and noncore_low is not None and noncore_high is not None else (None, None)
    ror = (
        ror_with_ci(a, c, b, d)
        if not sparse_agent_contrast
        else {
            "ror": None,
            "ror_ci_low": None,
            "ror_ci_high": None,
            "zero_cell_correction": "not reported; sparse agent-event or agent-noncore cell",
        }
    )
    return {
        "database": database,
        "drug_class": drug_class,
        "agent": agent,
        "estimand": "suspect-role agent share: core anaphylaxis minus non-core same-class reports",
        "class_suspect_reports_n": class_suspect_n,
        "core_class_suspect_reports_n": core_class_suspect_n,
        "noncore_class_suspect_reports_n": noncore_class_suspect_n,
        "agent_suspect_all_n": agent_total_n,
        "agent_suspect_core_n": a,
        "agent_suspect_noncore_n": c,
        "all_same_class_agent_report_share_proxy": agent_total_n / class_suspect_n if class_suspect_n else None,
        "core_agent_share": p_core,
        "core_agent_share_ci_low": core_low,
        "core_agent_share_ci_high": core_high,
        "noncore_agent_share_leave_event_out": p_noncore,
        "noncore_agent_share_ci_low": noncore_low,
        "noncore_agent_share_ci_high": noncore_high,
        "percentage_point_difference": (p_core - p_noncore) * 100 if p_core is not None and p_noncore is not None else None,
        "core_minus_all_report_share_proxy_pp": (p_core - agent_total_n / class_suspect_n) * 100 if p_core is not None and class_suspect_n else None,
        "newcombe_score_ci_low_pp": pp_low * 100 if pp_low is not None else None,
        "newcombe_score_ci_high_pp": pp_high * 100 if pp_high is not None else None,
        "a_agent_core": a,
        "b_other_agent_status_core": b,
        "c_agent_noncore": c,
        "d_other_agent_status_noncore": d,
        "sparse_agent_contrast": int(sparse_agent_contrast),
        **ror,
    }


def add_report_volume_proxy_rows(
    con: sqlite3.Connection,
    database: str,
    id_col: str,
    db_where: str,
    db_params: tuple[object, ...],
) -> list[dict[str, object]]:
    """Describe same-class reporting shares without treating them as use denominators."""
    rows: list[dict[str, object]] = []
    event_join = (
        "e.primaryid=h.primaryid"
        if id_col == "primaryid"
        else "e.database=h.database and e.report_id=h.report_id"
    )
    raw_db_where = db_where.replace("h.", "")
    for drug_class in ("NMBA", "antibiotic_anchor"):
        class_total = scalar(
            con,
            f"select count(distinct {id_col}) from class_hits where {raw_db_where} drug_class=? and suspect_role=1",
            db_params + (drug_class,),
        )
        class_core = scalar(
            con,
            f"""select count(distinct h.{id_col}) from class_hits h
                join event_agg e on {event_join}
                where {db_where} h.drug_class=? and h.suspect_role=1 and e.event_anaphylaxis_core=1""",
            db_params + (drug_class,),
        )
        agents = con.execute(
            f"""select h.generic_or_group, count(distinct h.{id_col})
                from class_hits h
                where {db_where} h.drug_class=? and h.suspect_role=1
                group by h.generic_or_group order by h.generic_or_group""",
            db_params + (drug_class,),
        ).fetchall()
        for agent, agent_total in agents:
            agent_core = scalar(
                con,
                f"""select count(distinct h.{id_col}) from class_hits h
                    join event_agg e on {event_join}
                    where {db_where} h.drug_class=? and h.generic_or_group=?
                      and h.suspect_role=1 and e.event_anaphylaxis_core=1""",
                db_params + (drug_class, agent),
            )
            agent_noncore = int(agent_total) - agent_core
            class_noncore = class_total - class_core
            all_share = int(agent_total) / class_total if class_total else None
            core_share = agent_core / class_core if class_core else None
            noncore_share = agent_noncore / class_noncore if class_noncore else None
            rows.append(
                {
                    "database": database,
                    "drug_class": drug_class,
                    "agent": agent,
                    "class_suspect_reports_n": class_total,
                    "core_class_suspect_reports_n": class_core,
                    "agent_suspect_all_n": int(agent_total),
                    "agent_suspect_core_n": agent_core,
                    "agent_suspect_noncore_n": agent_noncore,
                    "all_same_class_agent_report_share_proxy": all_share,
                    "core_agent_share": core_share,
                    "noncore_agent_share_leave_event_out": noncore_share,
                    "core_minus_all_report_share_proxy_pp": (core_share - all_share) * 100 if core_share is not None and all_share is not None else None,
                    "core_minus_noncore_agent_share_pp": (core_share - noncore_share) * 100 if core_share is not None and noncore_share is not None else None,
                    "sparse_agent_event_or_noncore": int(min(agent_core, agent_noncore) < 5),
                    "interpretation": "Same-class suspect-report share is a reporting-volume proxy, not a drug-utilization or exposure denominator; core reports are a subset of all reports.",
                }
            )
    return rows


def add_context_row(
    con: sqlite3.Connection,
    database: str,
    id_col: str,
    db_where: str,
    db_params: tuple[object, ...],
    drug_class: str,
) -> dict[str, object]:
    join = "join v14_context f on " + ("f.report_id=h.primaryid" if id_col == "primaryid" else "f.report_id=h.report_id")
    event_join = "join event_agg e on " + (
        "e.primaryid=h.primaryid" if id_col == "primaryid" else "e.database=h.database and e.report_id=h.report_id"
    )
    context_n = scalar(con, "select count(*) from v14_context")
    context_event_n = scalar(
        con,
        """select count(*) from v14_context f join event_agg e on """
        + ("e.primaryid=f.report_id" if id_col == "primaryid" else "e.report_id=f.report_id and e.database=?")
        + " where e.event_anaphylaxis_core=1",
        (() if id_col == "primaryid" else db_params),
    )
    a = scalar(
        con,
        f"""select count(distinct h.{id_col}) from class_hits h {join} {event_join}
            where {db_where} h.drug_class=? and h.suspect_role=1
              and e.event_anaphylaxis_core=1""",
        db_params + (drug_class,),
    )
    exposed_n = scalar(
        con,
        f"""select count(distinct h.{id_col}) from class_hits h {join}
            where {db_where} h.drug_class=? and h.suspect_role=1""",
        db_params + (drug_class,),
    )
    b = exposed_n - a
    c = context_event_n - a
    d = context_n - a - b - c
    return {
        "database": database,
        "drug_class": drug_class,
        "background": "nonculprit_medication_marker_context",
        "context_definition": "at least one sedative, opioid, volatile anesthetic, local anesthetic, or vasoactive-rescue marker; no culprit-class marker used for eligibility",
        "context_reports_n": context_n,
        "context_core_anaphylaxis_n": context_event_n,
        "a_target_suspect_core": a,
        "b_target_suspect_noncore": b,
        "c_non_target_core": c,
        "d_non_target_noncore": d,
        **ror_with_ci(a, b, c, d),
    }


def create_context_temp(
    con: sqlite3.Connection, id_col: str, database: str, flags_table: str
) -> None:
    predicate = " OR ".join(f"any_{name}=1" for name in NONCULPRIT_CONTEXT_CLASSES)
    con.execute("drop table if exists temp.v14_context")
    if id_col == "primaryid":
        con.execute(
            f"create temp table v14_context as select primaryid as report_id from {flags_table} where {predicate}"
        )
    else:
        con.execute(
            f"create temp table v14_context as select report_id from {flags_table} where database=? and ({predicate})",
            (database,),
        )
    con.execute("create unique index temp.idx_v14_context_report_id on v14_context(report_id)")


def add_context_agent_row(
    con: sqlite3.Connection,
    database: str,
    id_col: str,
    db_where: str,
    db_params: tuple[object, ...],
    drug_class: str,
    agent: str,
) -> dict[str, object]:
    """Repeat the role-consistent agent contrast inside target-independent context."""
    join = "join v14_context x on " + ("x.report_id=h.primaryid" if id_col == "primaryid" else "x.report_id=h.report_id")
    event_join = "join event_agg e on " + (
        "e.primaryid=h.primaryid" if id_col == "primaryid" else "e.database=h.database and e.report_id=h.report_id"
    )
    class_total = scalar(
        con,
        f"select count(distinct h.{id_col}) from class_hits h {join} where {db_where} h.drug_class=? and h.suspect_role=1",
        db_params + (drug_class,),
    )
    class_core = scalar(
        con,
        f"select count(distinct h.{id_col}) from class_hits h {join} {event_join} where {db_where} h.drug_class=? and h.suspect_role=1 and e.event_anaphylaxis_core=1",
        db_params + (drug_class,),
    )
    agent_total = scalar(
        con,
        f"select count(distinct h.{id_col}) from class_hits h {join} where {db_where} h.drug_class=? and h.generic_or_group=? and h.suspect_role=1",
        db_params + (drug_class, agent),
    )
    agent_core = scalar(
        con,
        f"select count(distinct h.{id_col}) from class_hits h {join} {event_join} where {db_where} h.drug_class=? and h.generic_or_group=? and h.suspect_role=1 and e.event_anaphylaxis_core=1",
        db_params + (drug_class, agent),
    )
    class_noncore, agent_noncore = class_total - class_core, agent_total - agent_core
    p_core = agent_core / class_core if class_core else None
    p_noncore = agent_noncore / class_noncore if class_noncore else None
    core_low, core_high = wilson_interval(agent_core, class_core)
    noncore_low, noncore_high = wilson_interval(agent_noncore, class_noncore)
    pp_low, pp_high = newcombe_wilson_difference_interval(
        p_core, core_low, core_high, p_noncore, noncore_low, noncore_high
    ) if p_core is not None and core_low is not None and core_high is not None and p_noncore is not None and noncore_low is not None and noncore_high is not None else (None, None)
    a, b, c, d = agent_core, class_core - agent_core, agent_noncore, class_noncore - agent_noncore
    sparse = min(a, c) < 3
    ror = ror_with_ci(a, c, b, d) if not sparse else {
        "ror": None, "ror_ci_low": None, "ror_ci_high": None,
        "zero_cell_correction": "not reported; sparse agent-event or agent-noncore cell",
    }
    return {
        "database": database,
        "drug_class": drug_class,
        "agent": agent,
        "context_definition": "target-independent non-culprit medication-marker context",
        "context_class_suspect_reports_n": class_total,
        "context_core_class_suspect_reports_n": class_core,
        "context_noncore_class_suspect_reports_n": class_noncore,
        "agent_suspect_core_n": agent_core,
        "agent_suspect_noncore_n": agent_noncore,
        "core_agent_share": p_core,
        "noncore_agent_share_leave_event_out": p_noncore,
        "percentage_point_difference": (p_core - p_noncore) * 100 if p_core is not None and p_noncore is not None else None,
        "newcombe_score_ci_low_pp": pp_low * 100 if pp_low is not None else None,
        "newcombe_score_ci_high_pp": pp_high * 100 if pp_high is not None else None,
        "a_agent_core": a, "b_other_agent_status_core": b,
        "c_agent_noncore": c, "d_other_agent_status_noncore": d,
        "sparse_agent_contrast": int(sparse),
        **ror,
    }


def run() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--databases",
        default="FAERS,Canada,JADER",
        help="Comma-separated subset of FAERS, Canada, and JADER for independent source-branch rebuilding.",
    )
    args = parser.parse_args()
    selected = {item.strip() for item in args.databases.split(",") if item.strip()}
    unknown = selected - {"FAERS", "Canada", "JADER"}
    if unknown:
        raise SystemExit(f"Unknown database selection: {sorted(unknown)}")
    primary_rows: list[dict[str, object]] = []
    context_rows: list[dict[str, object]] = []
    context_agent_rows: list[dict[str, object]] = []
    report_volume_proxy_rows: list[dict[str, object]] = []
    sources = (
        ("FAERS", FAERS_DB, "primaryid", "1=1 and", (), "context_flags"),
        ("Canada", CJ_DB, "report_id", "h.database=? and", ("Canada",), "report_level_flags"),
        ("JADER", CJ_DB, "report_id", "h.database=? and", ("JADER",), "report_level_flags"),
    )
    for database, path, id_col, db_where, db_params, flags_table in sources:
        if database not in selected:
            continue
        if not path.exists():
            raise FileNotFoundError(f"Selected {database} source flags database is missing: {path}")
        con = sqlite3.connect(path)
        create_context_temp(con, id_col, database, flags_table)
        for drug_class, agent in PRIMARY_AGENTS:
            primary_rows.append(
                add_primary_agent_row(con, database, id_col, db_where, db_params, drug_class, agent)
            )
            context_agent_rows.append(
                add_context_agent_row(con, database, id_col, db_where, db_params, drug_class, agent)
            )
        for drug_class in PRIMARY_CLASSES:
            context_rows.append(
                add_context_row(con, database, id_col, db_where, db_params, drug_class)
            )
        report_volume_proxy_rows.extend(
            add_report_volume_proxy_rows(con, database, id_col, db_where, db_params)
        )
        con.close()
    write_csv(OUT / "agent_primary_suspect_suspect_v0_26.csv", primary_rows)
    write_csv(OUT / "nonculprit_context_background_v0_26.csv", context_rows)
    write_csv(OUT / "agent_within_nonculprit_context_v0_26.csv", context_agent_rows)
    write_csv(OUT / "agent_same_class_report_volume_proxy_v0_26.csv", report_volume_proxy_rows)


if __name__ == "__main__":
    run()
