#!/usr/bin/env python3
"""Download, validate, and fingerprint official FAERS quarterly ASCII ZIPs.

The downloader is intentionally sequential and resumable.  It never treats a
partial file as an input: a ZIP must pass its CRC check and contain the three
tables used by the v0.26 raw build before it is logged as valid.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import sys
import time
import urllib.error
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent
DEFAULT_MANIFEST = ROOT / "faers_ascii_manifest_2004q1_2025q4.csv"
DEFAULT_OUT = ROOT.parent / "raw_inputs" / "faers_ascii_2004q1_2025q4"
DEFAULT_LOG = ROOT.parent / "run_logs" / "faers_download_manifest_v0_26.csv"
REQUIRED_MARKERS = ("DEMO", "DRUG", "REAC")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_zip(path: Path) -> tuple[bool, str]:
    if not path.exists() or path.stat().st_size == 0:
        return False, "missing_or_empty"
    try:
        with zipfile.ZipFile(path) as archive:
            invalid_member = archive.testzip()
            if invalid_member:
                return False, f"crc_error:{invalid_member}"
            names = [member.filename.upper() for member in archive.infolist() if not member.is_dir()]
            missing = [marker for marker in REQUIRED_MARKERS if not any(marker in name and name.endswith(".TXT") for name in names)]
            if missing:
                return False, "missing_required:" + ",".join(missing)
    except (OSError, zipfile.BadZipFile) as exc:
        return False, f"invalid_zip:{exc.__class__.__name__}"
    return True, "valid"


def fetch(url: str, destination: Path, retries: int, timeout: int) -> None:
    partial = destination.with_suffix(destination.suffix + ".part")
    for attempt in range(1, retries + 1):
        try:
            start = partial.stat().st_size if partial.exists() else 0
            headers = {"User-Agent": "periop-anaphylaxis-reanalysis/0.26"}
            if start:
                headers["Range"] = f"bytes={start}-"
            request = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(request, timeout=timeout) as response:
                status = getattr(response, "status", response.getcode())
                append = start > 0 and status == 206
                if start and not append:
                    partial.unlink(missing_ok=True)
                    start = 0
                mode = "ab" if append else "wb"
                expected = response.headers.get("Content-Length")
                expected_n = int(expected) if expected and expected.isdigit() else None
                written = 0
                with partial.open(mode) as output:
                    while chunk := response.read(1024 * 1024):
                        output.write(chunk)
                        written += len(chunk)
                if expected_n is not None and written != expected_n:
                    raise OSError(f"incomplete_http_body: expected={expected_n} received={written}")
            partial.replace(destination)
            valid, issue = validate_zip(destination)
            if valid:
                return
            destination.unlink(missing_ok=True)
            raise OSError(f"post_download_validation_failed:{issue}")
        except (OSError, urllib.error.URLError, urllib.error.HTTPError, TimeoutError):
            if attempt == retries:
                raise
            time.sleep(min(120, 10 * attempt))


def append_log(path: Path, row: dict[str, str]) -> None:
    fields = [
        "timestamp_utc", "year", "quarter", "filename", "url", "status",
        "bytes", "sha256", "validation", "message",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    exists = path.exists()
    with path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        if not exists:
            writer.writeheader()
        writer.writerow(row)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--log", type=Path, default=DEFAULT_LOG)
    parser.add_argument("--quarter", action="append", help="Limit to a quarter such as 2025q4; repeatable.")
    parser.add_argument("--retries", type=int, default=5)
    parser.add_argument("--timeout", type=int, default=180)
    args = parser.parse_args()

    with args.manifest.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    desired = {item.casefold() for item in args.quarter or []}
    if desired:
        rows = [row for row in rows if f"{row['year']}{row['quarter']}".casefold() in desired]
    if not rows:
        raise SystemExit("No source quarters selected.")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    failed = 0
    for row in rows:
        target = args.out_dir / row["filename"]
        valid, validation = validate_zip(target)
        status = "existing_valid" if valid else "downloaded_valid"
        message = ""
        if not valid:
            # An incomplete ZIP from a prior interrupted run cannot be trusted
            # as a resume base.  A retained .part file, in contrast, is an
            # explicitly interrupted HTTP body and is eligible for range resume.
            target.unlink(missing_ok=True)
            try:
                print(f"Retrieving {row['year']}{row['quarter']} ({row['filename']})", flush=True)
                fetch(row["url"], target, args.retries, args.timeout)
                valid, validation = validate_zip(target)
                if not valid:
                    status = "downloaded_invalid"
            except Exception as exc:  # Manifest records the exception, then continues to the next quarter.
                valid = False
                validation = "download_failed"
                status = "download_failed"
                message = f"{exc.__class__.__name__}: {exc}"
        if not valid:
            failed += 1
            print(f"FAILED {row['year']}{row['quarter']}: {validation} {message}", file=sys.stderr, flush=True)
        size = str(target.stat().st_size) if target.exists() else "0"
        digest = sha256(target) if valid else ""
        append_log(args.log, {
            "timestamp_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
            "year": row["year"], "quarter": row["quarter"], "filename": row["filename"],
            "url": row["url"], "status": status, "bytes": size, "sha256": digest,
            "validation": validation, "message": message,
        })

    summary = {
        "quarters_selected": len(rows), "quarters_failed": failed,
        "output_dir": str(args.out_dir), "log": str(args.log),
        "completed_at_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
    }
    print(json.dumps(summary, indent=2))
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
