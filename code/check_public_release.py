from __future__ import annotations

"""Check the public candidate before a GitHub or Zenodo upload."""

import ast
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FORBIDDEN_SUFFIXES = {".sqlite", ".db", ".zip", ".7z", ".rar", ".asc"}
FORBIDDEN_TEXT = ("C:\\", "C:/", "Users/", "OneDrive")
FORBIDDEN_DIRS = {
    ".r_libs",
    "__pycache__",
    ".pytest_cache",
    ".mypy_cache",
    "external_licensed",
    "raw_inputs",
    "work",
}
SECRET_PATTERN = re.compile(
    r"(?i)(api[_-]?key|access[_-]?token|password|client[_-]?secret)\s*[:=]\s*['\"][^'\"]+"
)
TOKEN_PATTERN = re.compile(r"(?:gh[pousr]_[A-Za-z0-9]{20,}|sk-[A-Za-z0-9_-]{20,}|BEGIN [A-Z ]*PRIVATE KEY)")
REQUIRED_DIRS = ("code", "config", "aggregate_outputs", "figures", "manifests")
TEXT_SUFFIXES = {".py", ".md", ".json", ".cff", ".r", ".csv", ".tsv", ".ps1", ".txt", ".yml", ".yaml", ".toml", ".svg"}


def main() -> None:
    errors: list[str] = []
    for name in REQUIRED_DIRS:
        if not (ROOT / name).is_dir():
            errors.append(f"missing directory: {name}")

    for path in ROOT.rglob("*"):
        if any(part in FORBIDDEN_DIRS for part in path.parts):
            errors.append(f"forbidden cache/data directory: {path.relative_to(ROOT)}")
            continue
        if not path.is_file() or ".git" in path.parts:
            continue
        if path.resolve() == Path(__file__).resolve():
            continue
        if path.suffix.lower() in FORBIDDEN_SUFFIXES or path.name.endswith(".downloading"):
            errors.append(f"forbidden data/archive file: {path.relative_to(ROOT)}")
        if path.stat().st_size > 20 * 1024 * 1024:
            errors.append(f"unexpectedly large file: {path.relative_to(ROOT)}")
        if path.suffix.lower() in TEXT_SUFFIXES:
            text = path.read_text(encoding="utf-8", errors="replace")
            if any(token in text for token in FORBIDDEN_TEXT):
                errors.append(f"local-path or download marker in: {path.relative_to(ROOT)}")
            if SECRET_PATTERN.search(text) or TOKEN_PATTERN.search(text):
                errors.append(f"credential-like assignment in: {path.relative_to(ROOT)}")

    for path in ROOT.rglob("*.py"):
        if any(part in FORBIDDEN_DIRS for part in path.parts):
            continue
        try:
            ast.parse(path.read_text(encoding="utf-8"))
        except SyntaxError as exc:
            errors.append(f"Python syntax error in {path.relative_to(ROOT)}: {exc}")

    for path in ROOT.rglob("*.json"):
        if any(part in FORBIDDEN_DIRS for part in path.parts):
            continue
        try:
            json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            errors.append(f"JSON error in {path.relative_to(ROOT)}: {exc}")

    if errors:
        raise SystemExit("PUBLIC RELEASE CHECK FAILED\n" + "\n".join(f"- {e}" for e in errors))
    print("PASS public release boundary, Python syntax, and manifest checks")


if __name__ == "__main__":
    main()
