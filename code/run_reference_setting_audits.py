from __future__ import annotations

"""Audit medication-burden selection, teicoplanin geography, and source counts.

Medication-burden strata are diagnostic sensitivity analyses, not confounding
adjustment. Upstream agent-contrast inputs and this script's outputs use the
same configured aggregate-output directory unless explicitly overridden.
"""

import csv
import math
import os
import sqlite3
from collections import defaultdict
from pathlib import Path
from statistics import NormalDist


REPO_ROOT = Path(__file__).resolve().parents[1]
WORK = REPO_ROOT / "work"
FLAG_DB = Path(os.environ.get("PV_FAERS_FLAGS_DB", WORK / "faers_periop_hypersensitivity_prelim.sqlite"))
RAW_DB = Path(os.environ.get("PV_FAERS_RAW_DB", REPO_ROOT / "external_data" / "faers.sqlite"))
LATEST_CACHE = Path(os.environ.get("PV_FAERS_CACHE_DB", REPO_ROOT / "external_data" / "faers_cache.sqlite"))
OUT = Path(os.environ.get("PV_ANALYSIS_OUTPUT_DIR", REPO_ROOT / "aggregate_outputs"))
ESTIMAND_OUT = Path(os.environ.get("PV_ESTIMAND_OUTPUT_DIR", OUT))
PRIMARY_CLASSES = ("NMBA", "antibiotic_anchor", "chlorhexidine", "dye_anchor")


def percentile(values: list[int], q: float) -> float | None:
    if not values:
        return None
    values = sorted(values)
    position = (len(values) - 1) * q
    low, high = int(position), min(int(position) + 1, len(values) - 1)
    return values[low] + (values[high] - values[low]) * (position - low)


def burden_stratum(value: int) -> str:
    if value < 2:
        return "<2"
    if value <= 4:
        return "2-4"
    if value <= 8:
        return "5-8"
    return ">=9"


def ror_with_ci(a: int, b: int, c: int, d: int) -> dict[str, float | str]:
    aa, bb, cc, dd = a, b, c, d
    correction = "none"
    if min(a, b, c, d) == 0:
        aa, bb, cc, dd = a + 0.5, b + 0.5, c + 0.5, d + 0.5
        correction = "Haldane-Anscombe 0.5"
    estimate = (aa * dd) / (bb * cc)
    se = math.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    return {
        "ror": estimate,
        "ror_ci_low": math.exp(math.log(estimate) - 1.96 * se),
        "ror_ci_high": math.exp(math.log(estimate) + 1.96 * se),
        "zero_cell_correction": correction,
    }


def wilson_interval(successes: int, total: int, z: float) -> tuple[float | None, float | None]:
    if total == 0:
        return None, None
    proportion = successes / total
    denominator = 1 + z * z / total
    centre = (proportion + z * z / (2 * total)) / denominator
    half = z * math.sqrt((proportion * (1 - proportion) + z * z / (4 * total)) / total) / denominator
    return centre - half, centre + half


def newcombe_difference_interval(a: int, n1: int, c: int, n0: int, z: float) -> tuple[float | None, float | None]:
    if n1 == 0 or n0 == 0:
        return None, None
    p1, p0 = a / n1, c / n0
    l1, u1 = wilson_interval(a, n1, z)
    l0, u0 = wilson_interval(c, n0, z)
    if None in (l1, u1, l0, u0):
        return None, None
    difference = p1 - p0
    lower = difference - math.sqrt((p1 - l1) ** 2 + (u0 - p0) ** 2)
    upper = difference + math.sqrt((u1 - p1) ** 2 + (p0 - l0) ** 2)
    return lower * 100, upper * 100


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


def create_temp_tables(con: sqlite3.Connection) -> None:
    predicate = " OR ".join(
        f"c.any_{name}=1"
        for name in ("sedative_hypnotic", "opioid", "volatile_anesthetic", "local_anesthetic", "vasoactive_rescue")
    )
    con.executescript(
        """
        drop table if exists temp.v16_context;
        drop table if exists temp.v16_context_burden;
        drop table if exists temp.v16_target_burden;
        drop table if exists temp.v16_teicoplanin;
        """
    )
    con.execute(
        f"""
        create temp table v16_context as
        select c.primaryid, coalesce(e.event_anaphylaxis_core, 0) as core
        from context_flags c
        left join event_agg e on e.primaryid=c.primaryid
        where {predicate}
        """
    )
    con.execute("create unique index temp.idx_v16_context_primaryid on v16_context(primaryid)")
    con.execute(
        """
        create temp table v16_context_burden as
        select x.primaryid, x.core, count(distinct d.drug_seq) as drug_sequence_n
        from v16_context x
        join raw.drug d indexed by idx_drug_primaryid on d.primaryid=x.primaryid
        join latest.latest_demo l on l.primaryid=d.primaryid
                                and l.source_year=d.source_year
                                and l.source_quarter=d.source_quarter
        group by x.primaryid, x.core
        """
    )
    con.execute("create unique index temp.idx_v16_context_burden_primaryid on v16_context_burden(primaryid)")
    placeholders = ",".join("?" for _ in PRIMARY_CLASSES)
    con.execute(
        f"""
        create temp table v16_target_burden as
        select distinct h.primaryid, h.drug_class, coalesce(e.event_anaphylaxis_core, 0) as core,
               case when x.primaryid is null then 0 else 1 end as in_context,
               coalesce(b.drug_sequence_n, raw_n.drug_sequence_n, 0) as drug_sequence_n
        from class_hits h
        left join event_agg e on e.primaryid=h.primaryid
        left join v16_context x on x.primaryid=h.primaryid
        left join v16_context_burden b on b.primaryid=h.primaryid
        left join (
            select h2.primaryid, count(distinct d2.drug_seq) as drug_sequence_n
            from class_hits h2
            join raw.drug d2 indexed by idx_drug_primaryid on d2.primaryid=h2.primaryid
            join latest.latest_demo l2 on l2.primaryid=d2.primaryid
                                     and l2.source_year=d2.source_year
                                     and l2.source_quarter=d2.source_quarter
            where h2.suspect_role=1 and h2.drug_class in ({placeholders})
            group by h2.primaryid
        ) raw_n on raw_n.primaryid=h.primaryid
        where h.suspect_role=1 and h.drug_class in ({placeholders})
        """,
        PRIMARY_CLASSES + PRIMARY_CLASSES,
    )
    con.execute("create index temp.idx_v16_target_burden_class on v16_target_burden(drug_class, core, in_context)")
    con.execute(
        """
        create temp table v16_teicoplanin as
        select distinct h.primaryid, coalesce(e.event_anaphylaxis_core, 0) as core,
               case when x.primaryid is null then 0 else 1 end as in_context
        from class_hits h
        left join event_agg e on e.primaryid=h.primaryid
        left join v16_context x on x.primaryid=h.primaryid
        where h.suspect_role=1 and h.drug_class='antibiotic_anchor'
          and h.generic_or_group='teicoplanin'
        """
    )


def medication_burden_audit(con: sqlite3.Connection) -> None:
    records = con.execute(
        "select drug_class, core, in_context, drug_sequence_n from v16_target_burden"
    ).fetchall()
    grouped: dict[tuple[str, int, int], list[int]] = defaultdict(list)
    for drug_class, core, in_context, burden in records:
        grouped[(drug_class, core, in_context)].append(int(burden))
    rows: list[dict[str, object]] = []
    for drug_class in PRIMARY_CLASSES:
        for core in (0, 1):
            for in_context in (0, 1):
                values = grouped[(drug_class, core, in_context)]
                rows.append(
                    {
                        "database": "FAERS",
                        "drug_class": drug_class,
                        "outcome_group": "core anaphylaxis" if core else "non-core reports",
                        "context_status": "inside" if in_context else "outside",
                        "report_n": len(values),
                        "drug_sequence_median": percentile(values, 0.5),
                        "drug_sequence_q1": percentile(values, 0.25),
                        "drug_sequence_q3": percentile(values, 0.75),
                        "drug_sequence_mean": sum(values) / len(values) if values else None,
                        "definition": "distinct raw FAERS DRUG_SEQ records per suspect-class report; descriptive selection audit",
                    }
                )
    write_csv("faers_context_medication_burden_by_outcome_v0_16.csv", rows)


def burden_stratified_class_rors(con: sqlite3.Connection) -> None:
    context_rows = con.execute(
        "select primaryid, core, drug_sequence_n from v16_context_burden"
    ).fetchall()
    target_rows = con.execute(
        "select primaryid, drug_class from v16_target_burden where in_context=1"
    ).fetchall()
    target_by_class: dict[str, set[str]] = defaultdict(set)
    for primaryid, drug_class in target_rows:
        target_by_class[str(drug_class)].add(str(primaryid))
    by_stratum: dict[str, list[tuple[str, int]]] = defaultdict(list)
    for primaryid, core, burden in context_rows:
        by_stratum[burden_stratum(int(burden))].append((str(primaryid), int(core)))
    rows: list[dict[str, object]] = []
    for drug_class in PRIMARY_CLASSES:
        target_ids = target_by_class[drug_class]
        for stratum in ("<2", "2-4", "5-8", ">=9"):
            reports = by_stratum[stratum]
            a = sum(1 for primaryid, core in reports if core and primaryid in target_ids)
            b = sum(1 for primaryid, core in reports if not core and primaryid in target_ids)
            c = sum(1 for primaryid, core in reports if core and primaryid not in target_ids)
            d = sum(1 for primaryid, core in reports if not core and primaryid not in target_ids)
            sparse = min(a, b, c, d) < 5
            estimate = ror_with_ci(a, b, c, d) if not sparse else {
                "ror": None,
                "ror_ci_low": None,
                "ror_ci_high": None,
                "zero_cell_correction": "not reported; sparse or zero-support burden stratum",
            }
            rows.append(
                {
                    "database": "FAERS",
                    "drug_class": drug_class,
                    "drug_sequence_stratum": stratum,
                    "context_definition": "target-independent non-culprit medication-marker context",
                    "a_target_suspect_core": a,
                    "b_target_suspect_noncore": b,
                    "c_non_target_core": c,
                    "d_non_target_noncore": d,
                    "sparse": int(sparse),
                    "interpretation": "descriptive when sparse; burden-stratified sensitivity, not adjustment",
                    **estimate,
                }
            )
    write_csv("faers_context_class_ror_by_drug_burden_v0_16.csv", rows)


def teicoplanin_country_audit(con: sqlite3.Connection) -> None:
    rows = con.execute(
        """
        select 'all suspect reports' as scope, t.core,
               case
                 when trim(coalesce(d.occr_country,'')) <> '' then upper(trim(d.occr_country))
                 when trim(coalesce(d.reporter_country,'')) <> '' then upper(trim(d.reporter_country))
                 else 'MISSING'
               end as country,
               count(*) as report_n
        from v16_teicoplanin t
        left join latest.latest_demo d indexed by idx_latest_demo_primaryid on d.primaryid=t.primaryid
        group by scope, t.core, country
        union all
        select 'inside context' as scope, t.core,
               case
                 when trim(coalesce(d.occr_country,'')) <> '' then upper(trim(d.occr_country))
                 when trim(coalesce(d.reporter_country,'')) <> '' then upper(trim(d.reporter_country))
                 else 'MISSING'
               end as country,
               count(*) as report_n
        from v16_teicoplanin t
        left join latest.latest_demo d indexed by idx_latest_demo_primaryid on d.primaryid=t.primaryid
        where t.in_context=1
        group by scope, t.core, country
        union all
        select 'outside context' as scope, t.core,
               case
                 when trim(coalesce(d.occr_country,'')) <> '' then upper(trim(d.occr_country))
                 when trim(coalesce(d.reporter_country,'')) <> '' then upper(trim(d.reporter_country))
                 else 'MISSING'
               end as country,
               count(*) as report_n
        from v16_teicoplanin t
        left join latest.latest_demo d indexed by idx_latest_demo_primaryid on d.primaryid=t.primaryid
        where t.in_context=0
        group by scope, t.core, country
        order by scope, t.core, report_n desc, country
        """
    ).fetchall()
    aliases = {"GB": "United Kingdom", "UNITED KINGDOM": "United Kingdom", "JAPAN": "Japan", "JP": "Japan"}
    output = [
        {
            "database": "FAERS",
            "agent": "teicoplanin",
            "scope": scope,
            "outcome_group": "core anaphylaxis" if core else "non-core reports",
            "country_field": "OCCR_COUNTRY when present, otherwise REPORTER_COUNTRY",
            "country": country,
            "country_normalized_for_display": aliases.get(country, country),
            "report_n": report_n,
            "purpose": "geographic structure audit; not an exposure denominator",
        }
        for scope, core, country, report_n in rows
    ]
    write_csv("faers_teicoplanin_country_audit_v0_16.csv", output)


def agent_scale_and_multiplicity_audits() -> None:
    """Expose all 12 context changes and a labelled post hoc multiplicity sensitivity."""
    unrestricted_path = ESTIMAND_OUT / "agent_primary_suspect_suspect_v0_26.csv"
    context_path = ESTIMAND_OUT / "agent_within_nonculprit_context_v0_26.csv"
    with unrestricted_path.open(encoding="utf-8-sig", newline="") as handle:
        unrestricted = {
            (row["database"], row["drug_class"], row["agent"]): row
            for row in csv.DictReader(handle)
        }
    with context_path.open(encoding="utf-8-sig", newline="") as handle:
        context = list(csv.DictReader(handle))

    comparison_rows: list[dict[str, object]] = []
    multiplicity_rows: list[dict[str, object]] = []
    family_size = len(context)
    adjusted_z = NormalDist().inv_cdf(1 - (0.05 / family_size) / 2)
    for row in context:
        key = (row["database"], row["drug_class"], row["agent"])
        unrestricted_row = unrestricted[key]
        context_pp = float(row["percentage_point_difference"])
        unrestricted_pp = float(unrestricted_row["percentage_point_difference"])
        context_ror = float(row["ror"]) if row["ror"] else None
        unrestricted_ror = float(unrestricted_row["ror"]) if unrestricted_row["ror"] else None
        comparison_rows.append(
            {
                "database": row["database"],
                "drug_class": row["drug_class"],
                "agent": row["agent"],
                "unrestricted_percentage_point_difference": unrestricted_pp,
                "context_percentage_point_difference": context_pp,
                "percentage_point_change": context_pp - unrestricted_pp,
                "unrestricted_ror": unrestricted_ror,
                "context_ror": context_ror,
                "comparison_note": "Both measures are descriptive reporting contrasts. Different movements can occur when context restriction changes baseline agent shares.",
            }
        )
        a, b, c, d = (int(row[name]) for name in (
            "a_agent_core", "b_other_agent_status_core", "c_agent_noncore", "d_other_agent_status_noncore"
        ))
        adjusted_low, adjusted_high = newcombe_difference_interval(a, a + b, c, c + d, adjusted_z)
        multiplicity_rows.append(
            {
                "database": row["database"],
                "drug_class": row["drug_class"],
                "agent": row["agent"],
                "family_size": family_size,
                "familywise_alpha": 0.05,
                "method": "post hoc Bonferroni familywise Newcombe method 10 sensitivity",
                "context_percentage_point_difference": context_pp,
                "bonferroni_newcombe_ci_low_pp": adjusted_low,
                "bonferroni_newcombe_ci_high_pp": adjusted_high,
                "includes_zero": int(adjusted_low <= 0 <= adjusted_high) if adjusted_low is not None and adjusted_high is not None else None,
                "sparse_agent_contrast": row["sparse_agent_contrast"],
            }
        )
    write_csv("agent_unrestricted_vs_context_scale_audit_v0_16.csv", comparison_rows)
    write_csv("agent_context_bonferroni_newcombe_sensitivity_v0_16.csv", multiplicity_rows)


def observed_source_counts(con: sqlite3.Connection) -> None:
    raw = sqlite3.connect(f"file:{RAW_DB}?mode=ro", uri=True)
    try:
        counts = {
            "demo_rows": raw.execute("select count(*) from demo").fetchone()[0],
            "unique_primaryid": raw.execute("select count(distinct primaryid) from demo").fetchone()[0],
            "unique_caseid": raw.execute("select count(distinct caseid) from demo").fetchone()[0],
            "caseids_with_multiple_demo_rows": raw.execute(
                "select count(*) from (select caseid from demo group by caseid having count(*) > 1)"
            ).fetchone()[0],
            "coverage": raw.execute("select min(source_year || source_quarter), max(source_year || source_quarter) from demo").fetchone(),
        }
    finally:
        raw.close()
    write_csv(
        "faers_observed_source_branch_counts_v0_16.csv",
        [
            {
                "source_branch": "official quarterly ASCII raw-source rebuild with per-ZIP SHA-256 manifest and DELETE ledger",
                "coverage_start": counts["coverage"][0],
                "coverage_end": counts["coverage"][1],
                "demo_rows": counts["demo_rows"],
                "unique_primaryid": counts["unique_primaryid"],
                "unique_caseid": counts["unique_caseid"],
                "caseids_with_multiple_demo_rows": counts["caseids_with_multiple_demo_rows"],
                "deleted_case_application_status": "DELETE identifiers were audited against CASEID and PRIMARYID before retained-report construction",
            }
        ],
    )


def main() -> None:
    if not FLAG_DB.exists() or not RAW_DB.exists():
        raise FileNotFoundError("Required local FAERS analysis database is missing")
    con = sqlite3.connect(FLAG_DB)
    try:
        con.execute("attach database ? as raw", (str(RAW_DB),))
        con.execute("attach database ? as latest", (str(LATEST_CACHE),))
        create_temp_tables(con)
        medication_burden_audit(con)
        burden_stratified_class_rors(con)
        teicoplanin_country_audit(con)
    finally:
        con.close()
    observed_source_counts(sqlite3.connect(FLAG_DB))
    agent_scale_and_multiplicity_audits()


if __name__ == "__main__":
    main()
