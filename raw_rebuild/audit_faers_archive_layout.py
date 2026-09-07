#!/usr/bin/env python3
"""Audit table-member layouts across all downloaded FAERS ZIP archives."""

from __future__ import annotations

import argparse
import csv
import json
import zipfile
from pathlib import Path

from build_faers_raw_branch import MEMBER_RE, member_table


ROOT = Path(__file__).resolve().parent
DEFAULT_MANIFEST = ROOT / "faers_ascii_manifest_2004q1_2025q4.csv"
DEFAULT_ZIPS = ROOT.parent / "raw_inputs" / "faers_ascii_2004q1_2025q4"
DEFAULT_OUT = ROOT.parent / "run_logs" / "faers_archive_layout_audit_v0_26.csv"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--zip-dir", type=Path, default=DEFAULT_ZIPS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    with args.manifest.open(encoding="utf-8-sig", newline="") as handle:
        sources = list(csv.DictReader(handle))
    rows: list[dict[str, str]] = []
    for source in sources:
        path = args.zip_dir / source["filename"]
        with zipfile.ZipFile(path) as archive:
            members = [item.filename for item in archive.infolist() if not item.is_dir()]
        found = {table: [] for table in ("demo", "drug", "reac")}
        for member in members:
            table = member_table(member)
            if table:
                found[table].append(member)
        missing = [table for table, names in found.items() if not names]
        rows.append({
            "year": source["year"], "quarter": source["quarter"], "filename": source["filename"],
            "demo_members": " | ".join(found["demo"]), "drug_members": " | ".join(found["drug"]),
            "reac_members": " | ".join(found["reac"]), "status": "ok" if not missing else "missing:" + ",".join(missing),
        })
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    failures = [row for row in rows if row["status"] != "ok"]
    print(json.dumps({"archives": len(rows), "valid_layouts": len(rows) - len(failures), "failures": failures}, indent=2))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
