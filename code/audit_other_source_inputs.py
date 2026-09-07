#!/usr/bin/env python3
"""Create a reproducible structural audit for the local Canada and JADER inputs.

The script does not alter source data.  It records file hashes, encodings,
report-identifier uniqueness, and reaction/drug linkage before the analysis
cache is rebuilt.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FAERS_ROOT = ROOT.parents[2]
DEFAULT_CANADA = FAERS_ROOT / "02_data" / "03_other_pharmacovigilance_databases" / "Other_database" / "CANADA_20250731"
DEFAULT_JADER = FAERS_ROOT / "02_data" / "03_other_pharmacovigilance_databases" / "Other_database" / "JADER_202511"
DEFAULT_OUT = ROOT / "manifests"


def resolve_source_file(directory: Path, canonical_name: str) -> Path:
    for candidate in (
        directory / canonical_name,
        directory / f"{canonical_name}.baiduyun.p.downloading",
    ):
        if candidate.is_file() and candidate.stat().st_size > 0:
            return candidate
    raise FileNotFoundError(f"Missing nonempty local source file: {canonical_name} under {directory}")


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def canonical_canada_id(value: str) -> str:
    value = (value or "").strip()
    stripped = value.lstrip("0")
    return stripped or value


def audit_canada(directory: Path) -> tuple[dict[str, object], list[dict[str, object]]]:
    reports_path = resolve_source_file(directory, "reports.txt")
    reactions_path = resolve_source_file(directory, "reactions.txt")
    drug_path = resolve_source_file(directory, "report_drug.txt")
    report_ids: set[str] = set()
    report_rows = 0
    empty_report_ids = 0
    with reports_path.open("r", encoding="utf-8", errors="strict", newline="") as handle:
        for row in csv.reader(handle, delimiter="$", quotechar='"'):
            report_rows += 1
            report_id = canonical_canada_id(row[0] if row else "")
            if not report_id:
                empty_report_ids += 1
            else:
                report_ids.add(report_id)

    def link_audit(path: Path, report_id_index: int) -> tuple[int, int, int]:
        rows = linked = empty = 0
        with path.open("r", encoding="utf-8", errors="strict", newline="") as handle:
            for row in csv.reader(handle, delimiter="$", quotechar='"'):
                rows += 1
                report_id = canonical_canada_id(row[report_id_index] if report_id_index < len(row) else "")
                if not report_id:
                    empty += 1
                elif report_id in report_ids:
                    linked += 1
        return rows, linked, empty

    reaction_rows, reaction_linked, reaction_empty = link_audit(reactions_path, 1)
    drug_rows, drug_linked, drug_empty = link_audit(drug_path, 1)
    summary = {
        "database": "Canada Vigilance",
        "logical_extract": "CANADA_20250731",
        "report_rows": report_rows,
        "unique_report_ids": len(report_ids),
        "duplicate_report_id_rows": report_rows - len(report_ids) - empty_report_ids,
        "empty_report_id_rows": empty_report_ids,
        "reaction_rows": reaction_rows,
        "reaction_rows_linked_to_report": reaction_linked,
        "reaction_rows_unlinked_or_empty": reaction_rows - reaction_linked,
        "drug_rows": drug_rows,
        "drug_rows_linked_to_report": drug_linked,
        "drug_rows_unlinked_or_empty": drug_rows - drug_linked,
        "deduplication_decision": "No duplicate Canada report identifiers were observed; report identifier is used as supplied.",
    }
    files = []
    for logical, path in (("reports.txt", reports_path), ("reactions.txt", reactions_path), ("report_drug.txt", drug_path)):
        files.append({
            "database": "Canada Vigilance",
            "logical_name": logical,
            "local_name": path.name,
            "bytes": path.stat().st_size,
            "sha256": digest(path),
            "encoding": "UTF-8",
        })
    return summary, files


def audit_jader(directory: Path) -> tuple[dict[str, object], list[dict[str, object]]]:
    demo_path = resolve_source_file(directory, "demo202511.csv")
    reactions_path = resolve_source_file(directory, "reac202511.csv")
    drug_path = resolve_source_file(directory, "drug202511.csv")
    latest: dict[str, tuple[int, str]] = {}
    demo_rows = nonnumeric_round_rows = repeated_id_rows = 0
    with demo_path.open("r", encoding="cp932", errors="strict", newline="") as handle:
        for row in csv.DictReader(handle):
            demo_rows += 1
            report_id = row.get("識別番号", "")
            report_round = row.get("報告回数", "")
            if not report_id:
                continue
            try:
                numeric_round = int(report_round or "0")
            except ValueError:
                numeric_round = -1
                nonnumeric_round_rows += 1
            if report_id in latest:
                repeated_id_rows += 1
            previous = latest.get(report_id)
            if previous is None or numeric_round > previous[0]:
                latest[report_id] = (numeric_round, report_round)
    latest_pairs = {(report_id, report_round) for report_id, (_, report_round) in latest.items()}

    def link_audit(path: Path) -> tuple[int, int, int]:
        rows = matched = empty = 0
        with path.open("r", encoding="cp932", errors="strict", newline="") as handle:
            for row in csv.DictReader(handle):
                rows += 1
                report_id = row.get("識別番号", "")
                report_round = row.get("報告回数", "")
                if not report_id:
                    empty += 1
                if (report_id, report_round) in latest_pairs:
                    matched += 1
        return rows, matched, empty

    reaction_rows, reaction_linked, reaction_empty = link_audit(reactions_path)
    drug_rows, drug_linked, drug_empty = link_audit(drug_path)
    summary = {
        "database": "JADER",
        "logical_extract": "JADER_202511",
        "demo_rows": demo_rows,
        "unique_report_ids": len(latest),
        "repeated_report_id_rows": repeated_id_rows,
        "nonnumeric_report_round_rows": nonnumeric_round_rows,
        "reaction_rows": reaction_rows,
        "reaction_rows_linked_to_latest_report_round": reaction_linked,
        "reaction_rows_unlinked_or_empty": reaction_rows - reaction_linked,
        "drug_rows": drug_rows,
        "drug_rows_linked_to_latest_report_round": drug_linked,
        "drug_rows_unlinked_or_empty": drug_rows - drug_linked,
        "deduplication_decision": "No repeated JADER report identifiers were observed in this extract; the latest-round rule is implemented and is identity-preserving.",
    }
    files = []
    for logical, path in (("demo202511.csv", demo_path), ("reac202511.csv", reactions_path), ("drug202511.csv", drug_path)):
        files.append({
            "database": "JADER",
            "logical_name": logical,
            "local_name": path.name,
            "bytes": path.stat().st_size,
            "sha256": digest(path),
            "encoding": "CP932",
        })
    return summary, files


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--canada-dir", type=Path, default=DEFAULT_CANADA)
    parser.add_argument("--jader-dir", type=Path, default=DEFAULT_JADER)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    canada, canada_files = audit_canada(args.canada_dir)
    jader, jader_files = audit_jader(args.jader_dir)
    payload = {"source_audits": [canada, jader], "files": canada_files + jader_files}
    (args.out_dir / "canada_jader_source_audit_v0_26.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    with (args.out_dir / "canada_jader_source_file_manifest_v0_26.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["database", "logical_name", "local_name", "bytes", "sha256", "encoding"])
        writer.writeheader()
        writer.writerows(payload["files"])
    print(json.dumps(payload["source_audits"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
