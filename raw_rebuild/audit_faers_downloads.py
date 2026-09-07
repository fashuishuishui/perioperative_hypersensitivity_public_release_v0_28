#!/usr/bin/env python3
"""Independently audit the completed FAERS source ZIP collection."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DEFAULT_MANIFEST = ROOT / "faers_ascii_manifest_2004q1_2025q4.csv"
DEFAULT_ZIPS = ROOT.parent / "raw_inputs" / "faers_ascii_2004q1_2025q4"
DEFAULT_OUT = ROOT.parent / "run_logs" / "faers_independent_integrity_audit_v0_26.csv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate(path: Path) -> tuple[bool, str]:
    if not path.exists():
        return False, "missing"
    try:
        with zipfile.ZipFile(path) as archive:
            invalid_member = archive.testzip()
            if invalid_member:
                return False, f"crc_error:{invalid_member}"
            names = [member.filename.upper() for member in archive.infolist() if not member.is_dir()]
            missing = [name for name in ("DEMO", "DRUG", "REAC") if not any(name in item and item.endswith(".TXT") for item in names)]
            if missing:
                return False, "missing_required:" + ",".join(missing)
    except (OSError, zipfile.BadZipFile) as exc:
        return False, f"invalid_zip:{exc.__class__.__name__}"
    return True, ""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--zip-dir", type=Path, default=DEFAULT_ZIPS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    with args.manifest.open(encoding="utf-8-sig", newline="") as handle:
        manifest = list(csv.DictReader(handle))
    rows: list[dict[str, str]] = []
    for source in manifest:
        path = args.zip_dir / source["filename"]
        valid, issue = validate(path)
        rows.append({
            "year": source["year"], "quarter": source["quarter"], "filename": source["filename"],
            "url": source["url"], "bytes": str(path.stat().st_size) if path.exists() else "0",
            "sha256": sha256(path) if valid else "", "status": "valid" if valid else "invalid", "issue": issue,
        })
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    invalid = [row for row in rows if row["status"] != "valid"]
    summary = {
        "expected_quarters": len(rows),
        "valid_quarters": len(rows) - len(invalid),
        "invalid_quarters": invalid,
        "valid_source_gb": round(sum(int(row["bytes"]) for row in rows if row["status"] == "valid") / 1024**3, 2),
        "audit_file": str(args.output),
        "audit_sha256": sha256(args.output),
    }
    print(json.dumps(summary, indent=2))
    return 0 if not invalid else 1


if __name__ == "__main__":
    raise SystemExit(main())
