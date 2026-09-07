from __future__ import annotations

"""Audit JADER reporting structure without treating it as clinical heterogeneity.

The audit identifies vaccine-like reports from the local raw DRUG file, measures
their contribution to the three-PT core-outcome density, and reruns the
chlorhexidine reporting contrast after predefined structural exclusions.  The
same report unit used by the main JADER branch is retained throughout.
"""

import csv
import math
import os
import sqlite3
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
WORK = REPO_ROOT / "work"
JADER_DIR = Path(os.environ.get("PV_JADER_DIR", REPO_ROOT / "external_data" / "jader_202511"))
CJ_DB = Path(os.environ.get("PV_CANADA_JADER_FLAGS_DB", WORK / "report_level_flags_canada_jader.sqlite"))
OUT = Path(os.environ.get("PV_ANALYSIS_OUTPUT_DIR", REPO_ROOT / "aggregate_outputs"))

VACCINE_PATTERNS = (
    "covid", "sars", "\u30b3\u30ed\u30ca", "\u30ef\u30af\u30c1\u30f3", "\u30b3\u30df\u30ca\u30c6\u30a3",
    "\u30b9\u30d1\u30a4\u30af\u30d0\u30c3\u30af\u30b9", "\u30e2\u30c7\u30eb\u30ca", "\u30d5\u30a1\u30a4\u30b6\u30fc",
    "\u30d0\u30ad\u30b9\u30bc\u30d6\u30ea\u30a2", "\u30cc\u30d0\u30ad\u30bd\u30d3\u30c3\u30c9",
)
EARLY_2021_TOKENS = ("2021\u30fb\u7b2c\u4e00", "2021\u30fb\u7b2c\u4e8c")
CONTEXT_PREDICATE = " OR ".join(
    f"any_{name}=1"
    for name in ("sedative_hypnotic", "opioid", "volatile_anesthetic", "local_anesthetic", "vasoactive_rescue")
)


def resolve_source_file(canonical_name: str) -> Path:
    for candidate in (JADER_DIR / canonical_name, JADER_DIR / f"{canonical_name}.baiduyun.p.downloading"):
        if candidate.is_file() and candidate.stat().st_size > 0:
            return candidate
    raise FileNotFoundError(f"Missing JADER source file: {canonical_name} under {JADER_DIR}")


def write_csv(name: str, rows: list[dict[str, object]]) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with (OUT / name).open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def latest_jader_rounds() -> set[tuple[str, str]]:
    latest: dict[str, str] = {}
    with resolve_source_file("demo202511.csv").open("r", encoding="cp932", errors="strict", newline="") as handle:
        reader = csv.reader(handle)
        next(reader, None)
        for row in reader:
            report_id = row[0].strip() if len(row) > 0 else ""
            report_round = row[1].strip() if len(row) > 1 else ""
            if not report_id:
                continue
            previous = latest.get(report_id)
            if previous is None or int(report_round or "0") > int(previous or "0"):
                latest[report_id] = report_round
    return set(latest.items())


def vaccine_like_ids() -> set[str]:
    latest_pairs = latest_jader_rounds()
    identified: set[str] = set()
    with resolve_source_file("drug202511.csv").open("r", encoding="cp932", errors="strict", newline="") as handle:
        reader = csv.reader(handle)
        next(reader, None)
        for row in reader:
            report_id = row[0].strip() if len(row) > 0 else ""
            report_round = row[1].strip() if len(row) > 1 else ""
            if (report_id, report_round) not in latest_pairs:
                continue
            generic = row[4] if len(row) > 4 else ""
            brand = row[5] if len(row) > 5 else ""
            name = f"{generic} {brand}".casefold()
            if any(pattern.casefold() in name for pattern in VACCINE_PATTERNS):
                identified.add(report_id)
    return identified


def scalar(con: sqlite3.Connection, sql: str) -> int:
    return int(con.execute(sql).fetchone()[0] or 0)


def contrast_metrics(a: int, b: int, c: int, d: int) -> dict[str, object]:
    sparse = min(a, c) < 5 or min(a, b, c, d) == 0
    if sparse:
        return {
            "ror": None,
            "ror_ci_low": None,
            "ror_ci_high": None,
            "prr": None,
            "chi_square": None,
            "evans_plus_ror_positive": 0,
            "sparse_or_zero_support": 1,
            "display_note": "ROR withheld because an agent-event/nonevent cell is <5 or a 2x2 cell is zero.",
        }
    ror = a * d / (b * c)
    se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    prr = (a / (a + b)) / (c / (c + d))
    n = a + b + c + d
    chi_square = (a * d - b * c) ** 2 * n / ((a + b) * (c + d) * (a + c) * (b + d))
    low = math.exp(math.log(ror) - 1.96 * se)
    return {
        "ror": ror,
        "ror_ci_low": low,
        "ror_ci_high": math.exp(math.log(ror) + 1.96 * se),
        "prr": prr,
        "chi_square": chi_square,
        "evans_plus_ror_positive": int(a >= 3 and low > 1 and prr >= 2 and chi_square >= 4),
        "sparse_or_zero_support": 0,
        "display_note": "Descriptive reporting contrast; not a risk estimate or truth test.",
    }


def audit_structure(con: sqlite3.Connection, vaccine_ids: set[str]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    con.execute("drop table if exists temp.vaccine_like")
    con.execute("create temp table vaccine_like (report_id text primary key)")
    con.executemany("insert into vaccine_like(report_id) values (?)", [(value,) for value in vaccine_ids])
    con.execute("drop table if exists temp.early_2021")
    con.execute("create temp table early_2021 (report_id text primary key)")
    con.execute(
        """insert into early_2021(report_id)
           select report_id from report_level_flags
           where database='JADER' and (year_quarter like '2021%\u7b2c\u4e00%' or year_quarter like '2021%\u7b2c\u4e8c%')"""
    )
    scenarios = {
        "baseline": "1=1",
        "exclude_vaccine_like": "report_id not in (select report_id from vaccine_like)",
        "exclude_2021_q1_q2": "report_id not in (select report_id from early_2021)",
        "exclude_vaccine_like_and_2021_q1_q2": "report_id not in (select report_id from vaccine_like) and report_id not in (select report_id from early_2021)",
        "vaccine_like_only": "report_id in (select report_id from vaccine_like)",
    }
    rows: list[dict[str, object]] = []
    for scenario, condition in scenarios.items():
        n = scalar(con, f"select count(*) from report_level_flags where database='JADER' and {condition}")
        core = scalar(con, f"select count(*) from report_level_flags where database='JADER' and {condition} and event_anaphylaxis_core=1")
        for scope, scope_condition in (("whole_database", "1=1"), ("nonculprit_medication_marker_context", CONTEXT_PREDICATE)):
            where = f"database='JADER' and ({condition}) and ({scope_condition})"
            context_n = scalar(con, f"select count(*) from report_level_flags where {where}")
            core_n = scalar(con, f"select count(*) from report_level_flags where {where} and event_anaphylaxis_core=1")
            target_n = scalar(con, f"select count(*) from report_level_flags where {where} and suspect_chlorhexidine=1")
            a = scalar(con, f"select count(*) from report_level_flags where {where} and event_anaphylaxis_core=1 and suspect_chlorhexidine=1")
            b, c, d = target_n - a, core_n - a, context_n - target_n - core_n + a
            row = {
                "database": "JADER",
                "scenario": scenario,
                "scope": scope,
                "scenario_report_n": n,
                "scenario_core_anaphylaxis_n": core,
                "scenario_core_anaphylaxis_percent": 100 * core / n if n else None,
                "scope_report_n": context_n,
                "scope_core_anaphylaxis_n": core_n,
                "a_chlorhexidine_suspect_core": a,
                "b_chlorhexidine_suspect_noncore": b,
                "c_nonchlorhexidine_core": c,
                "d_nonchlorhexidine_noncore": d,
                "vaccine_like_report_n_total": len(vaccine_ids),
                "early_2021_report_n_total": scalar(con, "select count(*) from early_2021"),
                "interpretation": "Structural sensitivity; removal changes the reporting reference set and is not a causal adjustment.",
            }
            row.update(contrast_metrics(a, b, c, d))
            rows.append(row)
    distribution_rows: list[dict[str, object]] = []
    for report_round, year_quarter, report_n, core_n in con.execute(
        """select report_round, year_quarter, count(*), sum(event_anaphylaxis_core)
           from report_level_flags where database='JADER'
           group by report_round, year_quarter
           order by year_quarter, report_round"""
    ):
        distribution_rows.append(
            {
                "database": "JADER",
                "report_round": report_round,
                "year_quarter": year_quarter,
                "report_n": report_n,
                "core_anaphylaxis_n": core_n,
                "core_anaphylaxis_percent": 100 * core_n / report_n if report_n else None,
                "purpose": "Source-structure density audit; not a time-trend or incidence analysis.",
            }
        )
    return rows, distribution_rows


def main() -> None:
    if not CJ_DB.exists():
        raise FileNotFoundError(f"JADER/Canada flags database is missing: {CJ_DB}")
    vaccine_ids = vaccine_like_ids()
    con = sqlite3.connect(CJ_DB)
    try:
        rows, distribution_rows = audit_structure(con, vaccine_ids)
    finally:
        con.close()
    write_csv("jader_structure_and_vaccine_audit_v0_26.csv", rows)
    write_csv("jader_round_year_core_density_v0_26.csv", distribution_rows)
    print(f"JADER vaccine-like reports identified: {len(vaccine_ids):,}")


if __name__ == "__main__":
    main()
