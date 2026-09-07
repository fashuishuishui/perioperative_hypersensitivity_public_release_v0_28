"""Validate and lock the schema-fixed v0.26 aggregate artifacts.

This script does not rebuild report-level flags or modify source data.  It
checks that the independently generated aggregate outputs reconcile to their
reported 2 x 2 tables and to each other before they are used in manuscript
tables or figures.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import sqlite3
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = Path(os.environ.get("PV_ANALYSIS_OUTPUT_DIR", ROOT / "aggregate_outputs" / "schema_fixed"))
FAERS_FLAGS = Path(os.environ.get("PV_FAERS_FLAGS_DB", ROOT / "work" / "faers_flags_schema_fixed.sqlite"))
CJ_FLAGS = Path(
    os.environ.get(
        "PV_CANADA_JADER_FLAGS_DB",
        ROOT / "work" / "report_level_flags_canada_jader_schema_fixed.sqlite",
    )
)
FAERS_LATEST = Path(
    os.environ.get(
        "PV_FAERS_CACHE_DB",
        ROOT / "work" / "faers_raw_2004q1_2025q4_schema_fixed" / "faers_latest_2004q1_2025q4.sqlite",
    )
)

FILES = {
    "agent_primary": "agent_primary_suspect_suspect_v0_26.csv",
    "agent_context": "agent_within_nonculprit_context_v0_26.csv",
    "agent_volume_proxy": "agent_same_class_report_volume_proxy_v0_26.csv",
    "class_context": "nonculprit_context_background_v0_26.csv",
    "broad_pt": "faers_context_broad_pt_sensitivity_v0_20.csv",
    "country": "faers_agent_occr_country_stratified_v0_16.csv",
    "role_profile": "faers_context_suspect_role_profile_v0_16.csv",
    "role_contrast": "faers_agent_ps_only_role_contrast_v0_16.csv",
    "burden": "faers_context_class_ror_by_drug_burden_v0_16.csv",
    "jader_structure": "jader_structure_and_vaccine_audit_v0_26.csv",
    "smq": "official_anaphylactic_reaction_smq_v28_pilot_v0_26.csv",
    "breslow_day": "breslow_day_variant_validation_v0_26.csv",
    "database_summary": "report_level_database_summary_canada_jader.csv",
}


def read_csv(name: str) -> list[dict[str, str]]:
    path = OUT / FILES[name]
    if not path.is_file():
        raise FileNotFoundError(f"Required v0.26 aggregate output is missing: {path}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def as_int(row: dict[str, str], key: str) -> int:
    return int(row[key])


def as_float(row: dict[str, str], key: str) -> float | None:
    value = row.get(key, "")
    return float(value) if value not in (None, "") else None


def close(actual: float | None, expected: float, label: str, tolerance: float = 1e-9) -> None:
    if actual is None or not math.isclose(actual, expected, rel_tol=tolerance, abs_tol=tolerance):
        raise ValueError(f"{label}: observed={actual!r}; expected={expected!r}")


def ror(a: int, b: int, c: int, d: int) -> float | None:
    if min(a, b, c, d) == 0:
        return None
    return a * d / (b * c)


def validate_agent_rows(rows: list[dict[str, str]], label: str, context: bool) -> list[dict[str, object]]:
    locked: list[dict[str, object]] = []
    for row in rows:
        prefix = f"{label}: {row['database']} {row['drug_class']} {row['agent']}"
        if context:
            class_total_key = "context_class_suspect_reports_n"
            core_total_key = "context_core_class_suspect_reports_n"
            noncore_total_key = "context_noncore_class_suspect_reports_n"
        else:
            class_total_key = "class_suspect_reports_n"
            core_total_key = "core_class_suspect_reports_n"
            noncore_total_key = "noncore_class_suspect_reports_n"
        a, b, c, d = (as_int(row, key) for key in (
            "a_agent_core",
            "b_other_agent_status_core",
            "c_agent_noncore",
            "d_other_agent_status_noncore",
        ))
        core_n, noncore_n = a + b, c + d
        if core_n != as_int(row, core_total_key) or noncore_n != as_int(row, noncore_total_key):
            raise ValueError(f"{prefix}: 2 x 2 counts do not sum to reported core/non-core totals")
        if core_n + noncore_n != as_int(row, class_total_key):
            raise ValueError(f"{prefix}: class total does not equal core plus non-core reports")
        close(as_float(row, "core_agent_share"), a / core_n, f"{prefix} core share")
        close(as_float(row, "noncore_agent_share_leave_event_out"), c / noncore_n, f"{prefix} non-core share")
        close(as_float(row, "percentage_point_difference"), 100 * (a / core_n - c / noncore_n), f"{prefix} percentage-point difference")
        expected_ror = ror(a, b, c, d)
        sparse = as_int(row, "sparse_agent_contrast")
        if sparse:
            if as_float(row, "ror") is not None:
                raise ValueError(f"{prefix}: sparse agent contrast has a displayed ROR")
        elif expected_ror is not None:
            close(as_float(row, "ror"), expected_ror, f"{prefix} ROR")
        locked.append({
            "database": row["database"],
            "drug_class": row["drug_class"],
            "agent": row["agent"],
            "context": "target-independent medication context" if context else "role-symmetric baseline",
            "core_class_n": core_n,
            "noncore_class_n": noncore_n,
            "a": a,
            "b": b,
            "c": c,
            "d": d,
            "core_agent_share": a / core_n,
            "noncore_agent_share": c / noncore_n,
            "difference_pp": 100 * (a / core_n - c / noncore_n),
            "ror": expected_ror if not sparse else None,
            "sparse_display": int(bool(sparse)),
        })
    return locked


def validate_volume_proxy(
    rows: list[dict[str, str]], primary: list[dict[str, str]]
) -> None:
    primary_index = {
        (row["database"], row["drug_class"], row["agent"]): row
        for row in primary
    }
    for row in rows:
        key = (row["database"], row["drug_class"], row["agent"])
        primary_row = primary_index.get(key)
        if primary_row is None:
            continue
        prefix = f"report-volume proxy: {' '.join(key)}"
        class_n = as_int(row, "class_suspect_reports_n")
        agent_all = as_int(row, "agent_suspect_all_n")
        agent_core = as_int(row, "agent_suspect_core_n")
        agent_noncore = as_int(row, "agent_suspect_noncore_n")
        if agent_all != agent_core + agent_noncore:
            raise ValueError(f"{prefix}: all-agent count does not sum to core plus non-core")
        if class_n != as_int(primary_row, "class_suspect_reports_n"):
            raise ValueError(f"{prefix}: class denominator does not match the baseline contrast")
        close(
            as_float(row, "all_same_class_agent_report_share_proxy"),
            agent_all / class_n,
            f"{prefix} all-report share",
        )


def validate_country(
    rows: list[dict[str, str]], context: list[dict[str, str]]
) -> list[dict[str, object]]:
    context_index = {
        (row["drug_class"], row["agent"]): row
        for row in context
        if row["database"] == "FAERS"
    }
    grouped: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in rows:
        grouped.setdefault((row["drug_class"], row["agent"]), []).append(row)
    locked: list[dict[str, object]] = []
    expected_strata = {"GB", "Non-GB observed", "Missing OCCR_COUNTRY"}
    for key, group in grouped.items():
        if {row["country_stratum"] for row in group} != expected_strata:
            raise ValueError(f"country strata are incomplete for {' '.join(key)}")
        context_row = context_index[key]
        sums = {
            cell: sum(as_int(row, cell) for row in group)
            for cell in ("a", "b", "c", "d")
        }
        expected = {
            "a": as_int(context_row, "a_agent_core"),
            "b": as_int(context_row, "b_other_agent_status_core"),
            "c": as_int(context_row, "c_agent_noncore"),
            "d": as_int(context_row, "d_other_agent_status_noncore"),
        }
        if sums != expected:
            raise ValueError(f"country-to-context 2 x 2 sum-back failed for {' '.join(key)}: {sums} != {expected}")
        for row in group:
            a, b, c, d = (as_int(row, cell) for cell in ("a", "b", "c", "d"))
            locked.append({
                "drug_class": row["drug_class"],
                "agent": row["agent"],
                "country_stratum": row["country_stratum"],
                "a": a,
                "b": b,
                "c": c,
                "d": d,
                "difference_pp": 100 * (a / (a + b) - c / (c + d)),
                "ror": ror(a, b, c, d),
                "sparse_display": int(min(a, c) < 5),
            })
    return locked


def validate_burden(rows: list[dict[str, str]], class_context: list[dict[str, str]]) -> None:
    faers_context = {
        row["drug_class"]: row
        for row in class_context
        if row["database"] == "FAERS"
    }
    by_class: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        by_class.setdefault(row["drug_class"], []).append(row)
    for drug_class, group in by_class.items():
        sums = {
            "a_target_suspect_core": sum(as_int(row, "a_target_suspect_core") for row in group),
            "b_target_suspect_noncore": sum(as_int(row, "b_target_suspect_noncore") for row in group),
            "c_non_target_core": sum(as_int(row, "c_non_target_core") for row in group),
            "d_non_target_noncore": sum(as_int(row, "d_non_target_noncore") for row in group),
        }
        base = faers_context[drug_class]
        expected = {
            "a_target_suspect_core": as_int(base, "a_target_suspect_core"),
            "b_target_suspect_noncore": as_int(base, "b_target_suspect_noncore"),
            "c_non_target_core": as_int(base, "c_non_target_core"),
            "d_non_target_noncore": as_int(base, "d_non_target_noncore"),
        }
        if sums != expected:
            raise ValueError(f"drug-burden sum-back failed for {drug_class}: {sums} != {expected}")


def validate_role_profiles(rows: list[dict[str, str]], context: list[dict[str, str]]) -> None:
    expected: dict[tuple[str, str], int] = {}
    for row in context:
        if row["database"] != "FAERS":
            continue
        expected[(row["drug_class"], "core")] = as_int(
            row, "context_core_class_suspect_reports_n"
        )
        expected[(row["drug_class"], "non-core")] = as_int(
            row, "context_noncore_class_suspect_reports_n"
        )
    observed: dict[tuple[str, str], int] = {}
    for row in rows:
        key = (row["drug_class"], row["outcome_group"])
        observed[key] = observed.get(key, 0) + as_int(row, "report_n")
    for key, total in expected.items():
        if observed.get(key) != total:
            raise ValueError(f"role-profile sum-back failed for {key}: {observed.get(key)} != {total}")


def validate_jader_structure(rows: list[dict[str, str]], class_context: list[dict[str, str]]) -> None:
    row = next(
        item
        for item in rows
        if item["scenario"] == "baseline" and item["scope"] == "nonculprit_medication_marker_context"
    )
    base = next(
        item
        for item in class_context
        if item["database"] == "JADER" and item["drug_class"] == "chlorhexidine"
    )
    checks = (
        ("a_chlorhexidine_suspect_core", "a_target_suspect_core"),
        ("b_chlorhexidine_suspect_noncore", "b_target_suspect_noncore"),
        ("c_nonchlorhexidine_core", "c_non_target_core"),
        ("d_nonchlorhexidine_noncore", "d_non_target_noncore"),
    )
    for observed_key, expected_key in checks:
        if as_int(row, observed_key) != as_int(base, expected_key):
            raise ValueError(f"JADER structure audit does not reconcile {observed_key}")


def validate_smq(rows: list[dict[str, str]]) -> None:
    for row in rows:
        a, b, c, d = (as_int(row, cell) for cell in ("a", "b", "c", "d"))
        prefix = f"official SMQ: {row['database']} {row['smq_scope']} {row['agent']}"
        if a + b != as_int(row, "class_core_n") or c + d != as_int(row, "class_noncore_n"):
            raise ValueError(f"{prefix}: 2 x 2 counts do not sum to reported class totals")
        sparse = as_int(row, "sparse_display")
        expected = ror(a, b, c, d)
        if sparse:
            if as_float(row, "ror") is not None:
                raise ValueError(f"{prefix}: sparse/zero cell has a displayed ROR")
        elif expected is not None:
            close(as_float(row, "ror"), expected, f"{prefix} ROR")


def validate_breslow_day(rows: list[dict[str, str]]) -> None:
    expected_agents = {"rocuronium", "succinylcholine", "cefazolin", "teicoplanin"}
    observed_agents = {row["agent"] for row in rows}
    if observed_agents != expected_agents:
        raise ValueError(
            "Independent Breslow-Day validation does not cover the locked country contrasts"
        )
    required = (
        "conditional_mle_bd_p",
        "desctools_mh_bd_p",
        "desctools_tarone_p",
    )
    for row in rows:
        for field in required:
            value = as_float(row, field)
            if value is None or not 0 <= value <= 1:
                raise ValueError(f"Invalid independent Breslow-Day p value for {row['agent']}: {field}")


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def database_coverage() -> dict[str, dict[str, int]]:
    faers = sqlite3.connect(f"file:{FAERS_FLAGS.as_posix()}?mode=ro", uri=True)
    latest = sqlite3.connect(f"file:{FAERS_LATEST.as_posix()}?mode=ro", uri=True)
    canada_jader = sqlite3.connect(f"file:{CJ_FLAGS.as_posix()}?mode=ro", uri=True)
    try:
        faers_context_sql = (
            "any_sedative_hypnotic=1 or any_opioid=1 or any_volatile_anesthetic=1 "
            "or any_local_anesthetic=1 or any_vasoactive_rescue=1"
        )
        coverage = {
            "FAERS": {
                "retained_reports": int(latest.execute("select count(*) from latest_demo").fetchone()[0]),
                "core_anaphylaxis_reports": int(faers.execute("select sum(event_anaphylaxis_core) from event_agg").fetchone()[0]),
                "medication_context_reports": int(faers.execute(f"select count(*) from context_flags where {faers_context_sql}").fetchone()[0]),
            }
        }
        for database in ("Canada", "JADER"):
            total, core = canada_jader.execute(
                "select count(*), sum(event_anaphylaxis_core) from report_level_flags where database=?",
                (database,),
            ).fetchone()
            context = canada_jader.execute(
                "select count(*) from report_level_flags where database=? and "
                "(any_sedative_hypnotic=1 or any_opioid=1 or any_volatile_anesthetic=1 "
                "or any_local_anesthetic=1 or any_vasoactive_rescue=1)",
                (database,),
            ).fetchone()[0]
            coverage[database] = {
                "retained_reports": int(total),
                "core_anaphylaxis_reports": int(core),
                "medication_context_reports": int(context),
            }
        return coverage
    finally:
        faers.close()
        latest.close()
        canada_jader.close()


def main() -> None:
    primary = read_csv("agent_primary")
    context = read_csv("agent_context")
    volume = read_csv("agent_volume_proxy")
    class_context = read_csv("class_context")
    country = read_csv("country")
    burden = read_csv("burden")
    role_profile = read_csv("role_profile")
    jader_structure = read_csv("jader_structure")
    smq = read_csv("smq")
    breslow_day = read_csv("breslow_day")

    locked_primary = validate_agent_rows(primary, "role-symmetric baseline", context=False)
    locked_context = validate_agent_rows(context, "medication context", context=True)
    validate_volume_proxy(volume, primary)
    locked_country = validate_country(country, context)
    validate_burden(burden, class_context)
    validate_role_profiles(role_profile, context)
    validate_jader_structure(jader_structure, class_context)
    validate_smq(smq)
    validate_breslow_day(breslow_day)

    write_csv(OUT / "canonical_agent_primary_v0_26.csv", locked_primary)
    write_csv(OUT / "canonical_agent_context_v0_26.csv", locked_context)
    write_csv(OUT / "canonical_agent_country_v0_26.csv", locked_country)
    coverage = database_coverage()
    input_paths = {name: OUT / filename for name, filename in FILES.items()}
    input_paths.update({
        "faers_flags_db": FAERS_FLAGS,
        "canada_jader_flags_db": CJ_FLAGS,
        "faers_latest_cache": FAERS_LATEST,
    })
    manifest = {
        "build": "v0.26 schema-fixed canonical artifact lock",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "analysis_unit": "database-specific retained report unit",
        "coverage": coverage,
        "validation_rules": {
            "agent_tables": "2 x 2 cells sum to printed class denominators; displayed shares, percentage-point contrasts, and nonsparse RORs recompute from cells",
            "country": "FAERS GB/non-GB/missing-country 2 x 2 cells sum to the corresponding medication-context agent contrast",
            "burden": "FAERS medication-burden strata sum to the corresponding class-level medication-context contrast",
            "role_profiles": "FAERS PS/SS role profiles sum to the medication-context class denominators",
            "jader_structure": "baseline chlorhexidine medication-context counts reconcile between the primary context output and the structural audit",
            "official_smq": "SMQ 2 x 2 cells and nonsparse RORs recompute from the licensed-term aggregate output",
            "breslow_day": "R/DescTools independently recomputes the Breslow-Day and Tarone variants from the locked country-strata 2 x 2 tables",
        },
        "inputs": {
            name: {
                "file": path.name,
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for name, path in input_paths.items()
        },
        "outputs": {
            name: {
                "file": name,
                "bytes": (OUT / name).stat().st_size,
                "sha256": sha256(OUT / name),
            }
            for name in (
                "canonical_agent_primary_v0_26.csv",
                "canonical_agent_context_v0_26.csv",
                "canonical_agent_country_v0_26.csv",
            )
        },
        "result": "passed",
    }
    (OUT / "canonical_artifact_lock_v0_26.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({"result": "passed", "coverage": coverage}, indent=2))


if __name__ == "__main__":
    main()
