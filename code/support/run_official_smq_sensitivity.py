"""Existing MedDRA dollar-file and SMQ parsers. Licensed membership inputs are external."""
from pathlib import Path
from collections import defaultdict
import csv

def read_dollar(path: Path):
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        yield from csv.reader(handle, delimiter="$", quotechar='"')

def load_smq():
    smq_metadata = None
    for row in read_dollar(SMQ_LIST):
        if row and row[0] == SMQ_CODE:
            smq_metadata = {
                "code": row[0],
                "name": row[1],
                "level": row[2],
                "version": row[6],
                "status": row[7],
                "algorithm": row[8],
            }
            break
    if smq_metadata is None:
        raise RuntimeError("Anaphylactic reaction SMQ not found")
    if smq_metadata["name"] != SMQ_NAME:
        raise AssertionError(smq_metadata)
    if smq_metadata["version"] != "28.0" or smq_metadata["status"] != "A":
        raise AssertionError(smq_metadata)
    if smq_metadata["algorithm"] != EXPECTED_ALGORITHM:
        raise AssertionError(smq_metadata)

    pt_names = {}
    for row in read_dollar(PT_ASC):
        if len(row) >= 2:
            pt_names[row[0]] = row[1]
    llt_to_pt = {}
    for row in read_dollar(LLT_ASC):
        if len(row) >= 3:
            llt_to_pt[row[0]] = row[2]

    pt_membership: dict[str, dict[str, set[str]]] = defaultdict(
        lambda: {"categories": set(), "scopes": set()}
    )
    missing_codes = set()
    content_rows = 0
    for row in read_dollar(SMQ_CONTENT):
        if not row or row[0] != SMQ_CODE or len(row) < 7 or row[6] != "A":
            continue
        content_rows += 1
        term_code, scope, category = row[1], row[3], row[4]
        if term_code in pt_names:
            pt_code = term_code
        elif term_code in llt_to_pt:
            pt_code = llt_to_pt[term_code]
        else:
            missing_codes.add(term_code)
            continue
        pt_name = pt_names.get(pt_code)
        if pt_name is None:
            missing_codes.add(term_code)
            continue
        pt_membership[pt_name]["categories"].add(category)
        pt_membership[pt_name]["scopes"].add(scope)

    if missing_codes:
        raise RuntimeError(f"Unmapped licensed SMQ term codes: {len(missing_codes)}")
    categories = sorted(
        {category for item in pt_membership.values() for category in item["categories"]}
    )
    scopes = sorted({scope for item in pt_membership.values() for scope in item["scopes"]})
    if categories != ["A", "B", "C", "D"] or scopes != ["1", "2"]:
        raise AssertionError({"categories": categories, "scopes": scopes})
    return smq_metadata, pt_membership, content_rows
