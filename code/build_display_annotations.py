"""Build display-only sparse annotations from the locked v0.20 country table."""

from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SOURCE = REPO_ROOT / "aggregate_outputs" / "canonical_agent_country_strata_v0_20.csv"
OUTPUT = REPO_ROOT / "aggregate_outputs" / "agent_country_strata_display_v0_23.csv"
MANIFEST = REPO_ROOT / "manifests" / "display_manifest.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    with SOURCE.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))

    sparse_keys = set()
    for row in rows:
        a = int(row["a"])
        c = int(row["c"])
        sparse = a < 5 or c < 5
        row["sparse_display"] = "1" if sparse else "0"
        row["display_note"] = (
            "Sparse/descriptive: agent-present core or non-core cell <5"
            if sparse
            else "Country-composition diagnostic"
        )
        if sparse:
            sparse_keys.add((row["agent"], row["country_stratum"], a, c))

    expected = {
        ("cefazolin", "GB", 2, 2),
        ("teicoplanin", "Non-GB observed", 3, 108),
    }
    if sparse_keys != expected:
        raise AssertionError(f"Unexpected sparse display rows: {sparse_keys}")

    fieldnames = list(rows[0])
    with OUTPUT.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    manifest = {
        "build": "v0.23 display-only unit and sparse hardening",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "numeric_source": SOURCE.name,
        "numeric_source_sha256": sha256(SOURCE),
        "output": OUTPUT.name,
        "output_sha256": sha256(OUTPUT),
        "rule": "sparse_display=1 when agent-present core a<5 or agent-present non-core c<5",
        "sparse_rows": [list(item) for item in sorted(sparse_keys)],
        "numeric_change": False,
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"PASS: wrote {OUTPUT.name}; sparse rows={sorted(sparse_keys)}")


if __name__ == "__main__":
    main()
