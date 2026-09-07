from __future__ import annotations

"""Lock v0.20 manuscript tables to one report-level analysis branch.

This builder deliberately does not recalculate aliases from raw text.  Its
purpose is narrower and auditable: consolidate the already rerun report-level
analyses that share the same FAERS flags database, verify every sum-back, and
write the only CSVs permitted to feed v0.20 tables, figures, and manuscript
numbers.  A failure stops the build rather than silently mixing derivations.
"""

import csv
import hashlib
import json
import math
import os
from datetime import datetime, timezone
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
WORK = REPO_ROOT / "work"
OUT = REPO_ROOT / "aggregate_outputs"
V14 = WORK / "outputs"
V16 = WORK / "outputs"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def public_path(path: Path) -> str:
    """Expose repository-relative paths without leaking local source locations."""
    try:
        relative = path.resolve().relative_to(REPO_ROOT.resolve())
        if relative.parts and relative.parts[0] == "work":
            return f"intermediate_not_distributed/{path.name}"
        return relative.as_posix()
    except ValueError:
        return f"external_input/{path.name}"


def as_int(row: dict[str, str], name: str) -> int:
    return int(row[name])


def ror_ci(a: int, b: int, c: int, d: int) -> tuple[float | None, float | None, float | None]:
    if min(a, b, c, d) == 0:
        return None, None, None
    estimate = a * d / (b * c)
    se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return estimate, math.exp(math.log(estimate) - 1.96 * se), math.exp(math.log(estimate) + 1.96 * se)


def wilson(x: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    p = x / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt((p * (1 - p) + z * z / (4 * n)) / n) / denom
    return centre - half, centre + half


def mover_pp_ci(a: int, n1: int, c: int, n0: int) -> tuple[float, float]:
    p1, p0 = a / n1, c / n0
    l1, u1 = wilson(a, n1)
    l0, u0 = wilson(c, n0)
    diff = p1 - p0
    lower = diff - math.sqrt((p1 - l1) ** 2 + (u0 - p0) ** 2)
    upper = diff + math.sqrt((u1 - p1) ** 2 + (p0 - l0) ** 2)
    return lower * 100, upper * 100


def gamma_q_shape_one(x: float) -> float:
    """Chi-square survival function for df=2, the country-table design."""
    return math.exp(-x / 2)


def bd_moments(a: int, b: int, c: int, d: int, theta: float) -> tuple[float, float]:
    """Exact conditional mean/variance for one 2x2 stratum at common OR theta."""
    row1, row2, col1 = a + b, c + d, a + c
    low, high = max(0, col1 - row2), min(row1, col1)
    values = []
    for x in range(low, high + 1):
        logw = math.lgamma(row1 + 1) - math.lgamma(x + 1) - math.lgamma(row1 - x + 1)
        logw += math.lgamma(row2 + 1) - math.lgamma(col1 - x + 1) - math.lgamma(row2 - col1 + x + 1)
        values.append((x, logw + x * math.log(theta)))
    maximum = max(weight for _, weight in values)
    scaled = [(x, math.exp(weight - maximum)) for x, weight in values]
    total = sum(weight for _, weight in scaled)
    mean = sum(x * weight for x, weight in scaled) / total
    variance = sum((x - mean) ** 2 * weight for x, weight in scaled) / total
    return mean, variance


def conditional_mle_common_or(tables: list[tuple[int, int, int, int]]) -> float:
    """Conditional MLE used by the standard Breslow-Day statistic."""
    observed = sum(a for a, _, _, _ in tables)

    def score(log_theta: float) -> float:
        theta = math.exp(log_theta)
        return sum(bd_moments(a, b, c, d, theta)[0] for a, b, c, d in tables) - observed

    low, high = -30.0, 30.0
    if score(low) > 0 or score(high) < 0:
        raise RuntimeError("Unable to bracket conditional common odds ratio.")
    for _ in range(100):
        mid = (low + high) / 2
        if score(mid) < 0:
            low = mid
        else:
            high = mid
    return math.exp((low + high) / 2)


def mh_bd(tables: list[tuple[int, int, int, int]]) -> tuple[float, float, float, int, float]:
    mh_num = sum(a * d / (a + b + c + d) for a, b, c, d in tables)
    mh_den = sum(b * c / (a + b + c + d) for a, b, c, d in tables)
    mh = mh_num / mh_den
    conditional_mle = conditional_mle_common_or(tables)
    statistic = 0.0
    for a, b, c, d in tables:
        expected, variance = bd_moments(a, b, c, d, conditional_mle)
        statistic += (a - expected) ** 2 / variance
    df = len(tables) - 1
    if df != 2:
        raise ValueError("This implementation is intentionally limited to the three country strata.")
    return mh, conditional_mle, statistic, df, gamma_q_shape_one(statistic)


def main() -> None:
    context_path = V14 / "agent_within_nonculprit_context_v0_14.csv"
    burden_path = V16 / "faers_context_class_ror_by_drug_burden_v0_16.csv"
    country_path = V16 / "faers_agent_occr_country_stratified_v0_16.csv"
    role_path = V16 / "faers_agent_ps_only_role_contrast_v0_16.csv"
    role_profile_path = V16 / "faers_context_suspect_role_profile_v0_16.csv"
    inputs = (context_path, burden_path, country_path, role_path, role_profile_path)
    missing = [str(path) for path in inputs if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing canonical input(s): " + "; ".join(missing))

    context = read_csv(context_path)
    if len(context) != 12:
        raise ValueError(f"Expected 12 agent-context rows; found {len(context)}")
    context_out: list[dict[str, object]] = []
    for row in context:
        a, b, c, d = (as_int(row, key) for key in ("a_agent_core", "b_other_agent_status_core", "c_agent_noncore", "d_other_agent_status_noncore"))
        core_n, noncore_n = a + b, c + d
        if core_n != as_int(row, "context_core_class_suspect_reports_n") or noncore_n != as_int(row, "context_noncore_class_suspect_reports_n"):
            raise ValueError(f"Context sum-back failed for {row['database']} {row['agent']}")
        estimate, low, high = ror_ci(a, b, c, d)
        context_out.append({
            "database": row["database"], "drug_class": row["drug_class"], "agent": row["agent"],
            "core_class_suspect_n": core_n, "noncore_class_suspect_n": noncore_n,
            "a_agent_core": a, "b_other_agent_status_core": b,
            "c_agent_noncore": c, "d_other_agent_status_noncore": d,
            "difference_pp": (a / core_n - c / noncore_n) * 100,
            "difference_ci_low_pp": mover_pp_ci(a, core_n, c, noncore_n)[0],
            "difference_ci_high_pp": mover_pp_ci(a, core_n, c, noncore_n)[1],
            "ror": estimate, "ror_ci_low": low, "ror_ci_high": high,
            "sparse_agent_contrast": row["sparse_agent_contrast"],
            "zero_cell_handling": row["zero_cell_correction"],
        })
    write_csv(OUT / "canonical_agent_context_v0_20.csv", context_out)

    country = read_csv(country_path)
    country_out: list[dict[str, object]] = []
    grouped: dict[tuple[str, str], list[dict[str, object]]] = {}
    for row in country:
        a, b, c, d = (as_int(row, key) for key in ("a", "b", "c", "d"))
        estimate, low, high = ror_ci(a, b, c, d)
        output = {
            "database": row["database"], "drug_class": row["drug_class"], "agent": row["agent"],
            "country_field": "OCCR_COUNTRY", "country_stratum": row["country_stratum"],
            "class_core_n": a + b, "class_noncore_n": c + d,
            "a": a, "b": b, "c": c, "d": d,
            "difference_pp": (a / (a + b) - c / (c + d)) * 100,
            "difference_ci_low_pp": mover_pp_ci(a, a + b, c, c + d)[0],
            "difference_ci_high_pp": mover_pp_ci(a, a + b, c, c + d)[1],
            "ror": estimate, "ror_ci_low": low, "ror_ci_high": high,
            "sparse_interpretation": "No meaningful contrast support" if min(a, c) < 3 else "Interpret with country-composition context",
        }
        country_out.append(output)
        grouped.setdefault((row["drug_class"], row["agent"]), []).append(output)
    write_csv(OUT / "canonical_agent_country_strata_v0_20.csv", country_out)

    summary_out: list[dict[str, object]] = []
    context_index = {(row["database"], row["drug_class"], row["agent"]): row for row in context_out}
    for (drug_class, agent), rows in grouped.items():
        if len(rows) != 3:
            raise ValueError(f"Expected three country strata for {agent}")
        a, b, c, d = (sum(int(row[key]) for row in rows) for key in ("a", "b", "c", "d"))
        faers = context_index[("FAERS", drug_class, agent)]
        if (a, b, c, d) != (int(faers["a_agent_core"]), int(faers["b_other_agent_status_core"]), int(faers["c_agent_noncore"]), int(faers["d_other_agent_status_noncore"])):
            raise ValueError(f"Country-to-context sum-back failed for {agent}")
        mh, conditional_mle, bd, df, p = mh_bd([(int(row["a"]), int(row["b"]), int(row["c"]), int(row["d"])) for row in rows])
        summary_out.append({
            "database": "FAERS", "drug_class": drug_class, "agent": agent,
            "crude_ror": faers["ror"], "mh_ror": mh, "conditional_mle_common_or": conditional_mle, "breslow_day_statistic": bd,
            "breslow_day_df": df, "breslow_day_p_value": p,
            "heterogeneity_interpretation": "Rejected at alpha=0.05" if p < 0.05 else "Not rejected at alpha=0.05",
        })
    write_csv(OUT / "canonical_country_mh_breslow_day_v0_20.csv", summary_out)

    composition_rows: list[dict[str, object]] = []
    for outcome, n_key, agent_key in (
        ("Core antibiotic reports", "class_core_n", "a"),
        ("Non-core antibiotic reports", "class_noncore_n", "c"),
    ):
        for stratum in ("GB", "Non-GB observed"):
            subset = [row for row in country_out if row["drug_class"] == "antibiotic_anchor" and row["country_stratum"] == stratum]
            for row in subset:
                count, total = int(row[agent_key]), int(row[n_key])
                composition_rows.append({
                    "outcome_group": outcome, "country_stratum": stratum,
                    "agent": row["agent"], "agent_reports_n": count,
                    "class_reports_n": total, "agent_share_percent": count / total * 100,
                })
    write_csv(OUT / "canonical_antibiotic_country_composition_v0_20.csv", composition_rows)

    burden = read_csv(burden_path)
    chlorhexidine = [row for row in burden if row["drug_class"] == "chlorhexidine"]
    totals = {key: sum(as_int(row, key) for row in chlorhexidine) for key in ("a_target_suspect_core", "b_target_suspect_noncore", "c_non_target_core", "d_non_target_noncore")}
    if tuple(totals.values()) != (248, 1359, 7547, 357731):
        raise ValueError(f"Chlorhexidine burden sum-back failed: {totals}")
    write_csv(OUT / "canonical_chlorhexidine_burden_v0_20.csv", chlorhexidine)

    role = read_csv(role_path)
    role_out: list[dict[str, object]] = []
    for row in role:
        a, b, c, d = (as_int(row, key) for key in ("a_agent_ps_core", "b_other_ps_agent_status_core", "c_agent_ps_noncore", "d_other_ps_agent_status_noncore"))
        estimate, low, high = ror_ci(a, b, c, d)
        base = context_index[("FAERS", row["drug_class"], row["agent"])]
        role_out.append({
            "database": "FAERS", "drug_class": row["drug_class"], "agent": row["agent"],
            "ps_agent_core": a, "ps_class_core_n": a + b,
            "ps_agent_noncore": c, "ps_class_noncore_n": c + d,
            "ps_share_among_suspect_agent_core": a / int(base["a_agent_core"]),
            "ps_only_difference_pp": (a / (a + b) - c / (c + d)) * 100,
            "ps_only_difference_ci_low_pp": mover_pp_ci(a, a + b, c, c + d)[0],
            "ps_only_difference_ci_high_pp": mover_pp_ci(a, a + b, c, c + d)[1],
            "ps_only_ror": estimate, "ps_only_ror_ci_low": low, "ps_only_ror_ci_high": high,
            "role_definition": row["role_definition"],
        })
    write_csv(OUT / "canonical_ps_only_agent_contrasts_v0_20.csv", role_out)
    write_csv(OUT / "canonical_role_profiles_v0_20.csv", read_csv(role_profile_path))

    tracked = (
        Path(os.environ.get("PV_FAERS_FLAGS_DB", WORK / "faers_periop_hypersensitivity_prelim.sqlite")),
        Path(os.environ.get("PV_CANADA_JADER_FLAGS_DB", WORK / "report_level_flags_canada_jader.sqlite")),
        REPO_ROOT / "config" / "drug_class_aliases.csv",
        REPO_ROOT / "config" / "event_pt_groups.csv",
        REPO_ROOT / "code" / "run_estimand_context_sensitivities.py",
        REPO_ROOT / "code" / "run_reference_setting_audits.py",
        REPO_ROOT / "code" / "run_country_role_probes.py",
        REPO_ROOT / "code" / "build_canonical_artifacts.py",
        REPO_ROOT / "code" / "validate_breslow_day.R",
        REPO_ROOT / "code" / "run_broad_pt_sensitivity.py",
        REPO_ROOT / "code" / "run_teicoplanin_time_audit.py",
        *inputs,
        OUT / "faers_context_broad_pt_sensitivity_v0_20.csv",
        OUT / "faers_teicoplanin_context_core_year_v0_20.csv",
        OUT / "breslow_day_variant_validation.csv",
    )
    manifest = {
        "build": "v0.20 canonical manuscript artifact lock",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "analysis_unit": "report-level flags; FAERS context analysis uses the case-level branch stored in faers_periop_hypersensitivity_prelim.sqlite",
        "inputs": [{"path": public_path(path), "bytes": path.stat().st_size, "sha256": sha256(path)} for path in tracked],
        "invariants": {
            "context_rows": 12,
            "faers_country_sumback": "country a/b/c/d totals equal the corresponding FAERS agent-context 2x2",
            "chlorhexidine_burden_sumback": totals,
            "country_strata": ["GB", "Non-GB observed", "Missing OCCR_COUNTRY"],
            "antibiotic_composition_scope": "GB and observed non-GB strata with the same antibiotic class denominators used in Table 4",
        },
        "scope_boundary": "This lock does not reconstruct historical raw ZIP provenance, DELETE application, or a licensed SMQ. Those items must not be represented as completed analyses.",
    }
    (REPO_ROOT / "manifests" / "canonical_build_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    (REPO_ROOT / "docs" / "README_CANONICAL_BUILD.md").write_text(
        "# v0.20 Canonical Build\n\n"
        "These files are the sole numeric source for the v0.20 manuscript. Every agent table derives from the same report-level flags branch. "
        "The build halts if context, country, or burden totals do not reconcile. Historical source provenance and licensed terminology products remain separate unresolved reproducibility items and are not inferred by this artifact lock.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
