"""Build the non-licensed, aggregate-only v0.26 supplementary package."""

from __future__ import annotations

import csv
import hashlib
import json
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "supplementary" / "v0_26"
AGG = ROOT / "aggregate_outputs" / "schema_fixed"
CONFIG = ROOT / "config"
MANIFESTS = ROOT / "manifests"
FAERS_BUILD = ROOT / "work" / "faers_raw_2004q1_2025q4_schema_fixed" / "faers_source_build_manifest_v0_26.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def copy(source: Path, target_name: str) -> None:
    target = OUT / target_name
    shutil.copy2(source, target)


def canonical_source_name(name: str) -> str:
    for suffix in (".baiduyun.p.downloading", ".downloading"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def write_public_source_manifest(source: Path, target_name: str) -> None:
    with source.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])
    for row in rows:
        if "local_name" in row:
            row["local_name"] = canonical_source_name(row["local_name"])
    write_csv(target_name, rows, fieldnames)


def write_public_smq_manifest(source: Path, target_name: str) -> None:
    payload = read_json(source)
    for entry in payload.get("licensed_source_files", {}).values():
        if "file_name" in entry:
            entry["file_name"] = canonical_source_name(entry["file_name"])
    (OUT / target_name).write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def write_csv(name: str, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with (OUT / name).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    hash_manifest = OUT / "Supplementary_File_S4_manifest_v0_26.csv"
    if hash_manifest.exists():
        hash_manifest.unlink()
    faers = read_json(FAERS_BUILD)
    source_audit = read_json(MANIFESTS / "canada_jader_source_audit_v0_26.json")
    canada_official = read_json(MANIFESTS / "canada_official_archive_audit_v0_26.json")

    write_csv(
        "Supplementary_Table_S1_data_sources_and_retention_v0_26.csv",
        [
            {
                "database": "FAERS",
                "analysis_source": "Official FDA quarterly ASCII archives",
                "source_window": "2004 Q1-2025 Q4",
                "retained_unit": "One retained report per case family after DELETE and latest-version procedure",
                "retained_report_n": faers["cache_audit"]["retained_latest_demo_rows"],
                "deduplication_audit": faers["deduplication"],
                "provenance_note": "All 88 archives were re-downloaded, ZIP/CRC/hash audited, and parsed from raw DEMO/DRUG/REAC files.",
            },
            {
                "database": "Canada Vigilance",
                "analysis_source": "Frozen local extract labeled 31 Jul 2025",
                "source_window": "Local extract date 31 Jul 2025",
                "retained_unit": "Unique supplied report identifier",
                "retained_report_n": 1204167,
                "deduplication_audit": source_audit["source_audits"][0]["deduplication_decision"],
                "provenance_note": "The current official Open Canada archive was separately hash audited but ended 30 Nov 2024 and was not substituted into the prespecified analysis.",
            },
            {
                "database": "JADER",
                "analysis_source": "Frozen local extract labeled Nov 2025",
                "source_window": "Local extract date Nov 2025",
                "retained_unit": "Unique supplied report identifier; latest-round rule implemented",
                "retained_report_n": 998397,
                "deduplication_audit": source_audit["source_audits"][1]["deduplication_decision"],
                "provenance_note": "Local source hashes and linkages were audited. PMDA download terms require user acceptance; no automated replacement download was performed.",
            },
        ],
        [
            "database",
            "analysis_source",
            "source_window",
            "retained_unit",
            "retained_report_n",
            "deduplication_audit",
            "provenance_note",
        ],
    )

    copies = {
        CONFIG / "drug_class_aliases.csv": "Supplementary_Table_S2_exposure_aliases_v0_26.csv",
        CONFIG / "event_pt_groups.csv": "Supplementary_Table_S3_event_PT_mapping_v0_26.csv",
        AGG / "agent_primary_suspect_suspect_v0_26.csv": "Supplementary_Table_S4_role_symmetric_agent_contrasts_v0_26.csv",
        AGG / "agent_within_nonculprit_context_v0_26.csv": "Supplementary_Table_S5_medication_context_agent_contrasts_v0_26.csv",
        AGG / "agent_same_class_report_volume_proxy_v0_26.csv": "Supplementary_Table_S6_same_class_report_volume_proxy_v0_26.csv",
        AGG / "faers_coexposure_counts.csv": "Supplementary_Table_S7a_FAERS_co_suspected_classes_v0_26.csv",
        AGG / "coexposure_counts_canada_jader.csv": "Supplementary_Table_S7b_Canada_JADER_co_suspected_classes_v0_26.csv",
        AGG / "faers_context_medication_burden_by_outcome_v0_16.csv": "Supplementary_Table_S8a_FAERS_medication_burden_audit_v0_26.csv",
        AGG / "faers_context_class_ror_by_drug_burden_v0_16.csv": "Supplementary_Table_S8b_FAERS_burden_stratified_contrasts_v0_26.csv",
        AGG / "faers_context_suspect_role_profile_v0_16.csv": "Supplementary_Table_S8c_FAERS_role_profile_audit_v0_26.csv",
        AGG / "faers_agent_ps_only_role_contrast_v0_16.csv": "Supplementary_Table_S8d_FAERS_PS_only_contrasts_v0_26.csv",
        AGG / "canonical_agent_country_v0_26.csv": "Supplementary_Table_S9a_FAERS_country_strata_v0_26.csv",
        AGG / "faers_agent_occr_country_mh_standardization_v0_16.csv": "Supplementary_Table_S9b_FAERS_country_standardization_v0_26.csv",
        AGG / "breslow_day_variant_validation_v0_26.csv": "Supplementary_Table_S9c_Breslow_Day_validation_v0_26.csv",
        AGG / "faers_context_broad_pt_sensitivity_v0_20.csv": "Supplementary_Table_S10a_broad_PT_and_bronchospasm_sensitivity_v0_26.csv",
        AGG / "official_anaphylactic_reaction_smq_v28_pilot_v0_26.csv": "Supplementary_Table_S10b_official_SMQ_sensitivity_v0_26.csv",
        AGG / "official_anaphylactic_reaction_smq_v28_country_strata_v0_26.csv": "Supplementary_Table_S10c_official_SMQ_country_sensitivity_v0_26.csv",
        AGG / "official_anaphylactic_reaction_smq_v28_country_summary_v0_26.csv": "Supplementary_Table_S10d_official_SMQ_country_summary_v0_26.csv",
        AGG / "jader_structure_and_vaccine_audit_v0_26.csv": "Supplementary_Table_S11_JADER_structural_audit_v0_26.csv",
        MANIFESTS / "canada_official_archive_audit_v0_26.json": "Supplementary_File_S1_Canada_official_archive_audit_v0_26.json",
        AGG / "canonical_artifact_lock_v0_26.json": "Supplementary_File_S3_canonical_artifact_lock_v0_26.json",
    }
    for source, target in copies.items():
        copy(source, target)
    write_public_source_manifest(
        MANIFESTS / "canada_jader_source_file_manifest_v0_26.csv",
        "Supplementary_Table_S12_Canada_JADER_source_file_manifest_v0_26.csv",
    )
    write_public_smq_manifest(
        AGG / "official_smq_manifest_v0_26.json",
        "Supplementary_File_S2_official_SMQ_manifest_v0_26.json",
    )

    lines = [
        "# Supplementary material: v0.26 raw rebuild",
        "",
        "This package contains aggregate derived outputs, bounded project mappings, and source-provenance documentation. It does not contain raw regulatory reports or licensed MedDRA/SMQ membership content.",
        "",
        "## Contents",
        "",
        "- S1: data sources, retained reporting units, and database-specific deduplication rules.",
        "- S2: complete bounded exposure-alias mapping.",
        "- S3: outcome preferred-term mapping and group assignment.",
        "- S4-S6: role-symmetric contrasts, medication-context contrasts, and same-class report-volume proxies.",
        "- S7: co-suspected culprit-class counts and pairs.",
        "- S8: FAERS medication-burden and PS/SS role diagnostics.",
        "- S9: FAERS country-stratified, standardized, and independently validated heterogeneity results.",
        "- S10: broad PT, bronchospasm, and official SMQ terminology sensitivities.",
        "- S11: JADER vaccine-like and 2021 structural sensitivity audit.",
        "- S12: Canada and JADER frozen-source file hashes and linkage audit.",
        "- Files S1-S4: Canada official archive audit, official SMQ metadata, canonical artifact lock, and supplementary-file hashes.",
        "",
        "The official Open Canada archive audit documents a 30 November 2024 package that was not substituted for the prespecified 31 July 2025 Canada extract. JADER source extraction remains subject to PMDA access terms. These boundaries are stated in the main manuscript.",
    ]
    (OUT / "SUPPLEMENTARY_MANIFEST_v0_26.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    manifest_rows = []
    for path in sorted(OUT.iterdir()):
        if path.is_file() and path != hash_manifest:
            manifest_rows.append(
                {
                    "file": path.name,
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                    "content_boundary": "Aggregate or bounded mapping only; no raw reports or licensed MedDRA/SMQ membership terms.",
                }
            )
    write_csv(
        hash_manifest.name,
        manifest_rows,
        ["file", "bytes", "sha256", "content_boundary"],
    )


if __name__ == "__main__":
    main()
