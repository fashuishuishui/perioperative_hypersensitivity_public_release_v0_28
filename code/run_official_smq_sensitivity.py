"""Pilot the licensed MedDRA v28.0 Anaphylactic reaction SMQ.

The script keeps the licensed term membership local. Outputs contain only
aggregate membership counts, hashes, event counts, and derived 2x2 estimates.
FAERS and Canada can be evaluated with English PT strings. JADER cannot be
called an official SMQ analysis without a licensed MedDRA/J PT-code crosswalk.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import sqlite3
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
WORK = REPO_ROOT / "work"
OUT = Path(os.environ.get("PV_ANALYSIS_OUTPUT_DIR", REPO_ROOT / "aggregate_outputs"))
OUT.mkdir(parents=True, exist_ok=True)


def resolve_meddra_ascii(configured: Path) -> Path:
    """Accept either the MedAscii directory or its distribution parent."""
    for candidate in (configured, configured / "MedAscii"):
        if all((candidate / name).is_file() for name in ("smq_list.asc", "smq_content.asc", "pt.asc")):
            return candidate
    raise FileNotFoundError(
        "Could not locate a MedDRA MedAscii directory containing smq_list.asc, "
        "smq_content.asc, and pt.asc under "
        f"{configured}"
    )


def resolve_licensed_file(directory: Path, canonical_name: str) -> Path:
    """Resolve a local licensed source without renaming or copying it.

    Some local archives retain a downloader suffix.  The suffix is recorded in
    the aggregate manifest so provenance remains explicit; licensed content is
    never redistributed.
    """
    for candidate in (
        directory / canonical_name,
        directory / f"{canonical_name}.baiduyun.p.downloading",
    ):
        if candidate.is_file() and candidate.stat().st_size > 0:
            return candidate
    raise FileNotFoundError(
        f"Missing nonempty licensed source {canonical_name} under {directory}"
    )


MEDDRA_CONFIGURED = Path(
    os.environ.get("PV_MEDDRA_ASCII", REPO_ROOT / "external_licensed" / "meddra_v28_0")
)
MEDDRA_ASCII = resolve_meddra_ascii(MEDDRA_CONFIGURED)
MEDDRA_REFERENCE_ASCII = resolve_meddra_ascii(
    Path(os.environ.get("PV_MEDDRA_REFERENCE_ASCII", MEDDRA_ASCII))
)
SMQ_LIST = resolve_licensed_file(MEDDRA_ASCII, "smq_list.asc")
SMQ_CONTENT = resolve_licensed_file(MEDDRA_ASCII, "smq_content.asc")
PT_ASC = resolve_licensed_file(MEDDRA_ASCII, "pt.asc")
LLT_ASC = resolve_licensed_file(MEDDRA_REFERENCE_ASCII, "llt.asc")

RAW_FAERS = Path(os.environ.get("PV_FAERS_RAW_DB", REPO_ROOT / "external_data" / "faers.sqlite"))
LATEST_CACHE = Path(os.environ.get("PV_FAERS_CACHE_DB", REPO_ROOT / "external_data" / "faers_cache.sqlite"))
FAERS_FLAGS = Path(os.environ.get("PV_FAERS_FLAGS_DB", WORK / "faers_periop_hypersensitivity_prelim.sqlite"))
CJ_FLAGS = Path(os.environ.get("PV_CANADA_JADER_FLAGS_DB", WORK / "report_level_flags_canada_jader.sqlite"))
CANADA_REACTIONS = Path(os.environ.get("PV_CANADA_REACTIONS", REPO_ROOT / "external_data" / "canada_20250731" / "reactions.txt"))

SMQ_CODE = "20000021"
SMQ_NAME = "Anaphylactic reaction (SMQ)"
EXPECTED_ALGORITHM = "A or (B and C) or (D and (B or C))"
CONTEXT_CLASSES = [
    "sedative_hypnotic",
    "opioid",
    "volatile_anesthetic",
    "local_anesthetic",
    "vasoactive_rescue",
]
AGENTS = [
    ("NMBA", "rocuronium"),
    ("NMBA", "succinylcholine"),
    ("antibiotic_anchor", "cefazolin"),
    ("antibiotic_anchor", "teicoplanin"),
]
COUNTRY_STRATA = ("GB", "Non-GB observed", "Missing OCCR_COUNTRY")
Z = 1.96


def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_dollar(path: Path):
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        yield from csv.reader(handle, delimiter="$", quotechar='"')


def canonical_canada_report_id(value: str) -> str:
    value = (value or "").strip()
    stripped = value.lstrip("0")
    return stripped or value


def wilson(x: int, n: int) -> tuple[float, float]:
    p = x / n
    z2 = 1.959963984540054**2
    denominator = 1 + z2 / n
    center = (p + z2 / (2 * n)) / denominator
    half = 1.959963984540054 * math.sqrt(
        p * (1 - p) / n + z2 / (4 * n * n)
    ) / denominator
    return center - half, center + half


def newcombe(a: int, n1: int, c: int, n0: int) -> tuple[float, float]:
    p1, p0 = a / n1, c / n0
    l1, u1 = wilson(a, n1)
    l0, u0 = wilson(c, n0)
    diff = p1 - p0
    return (
        diff - math.sqrt((p1 - l1) ** 2 + (u0 - p0) ** 2),
        diff + math.sqrt((u1 - p1) ** 2 + (p0 - l0) ** 2),
    )


def ror_ci(a: int, b: int, c: int, d: int):
    if min(a, b, c, d) == 0:
        return None, None, None
    estimate = a * d / (b * c)
    se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return (
        estimate,
        math.exp(math.log(estimate) - Z * se),
        math.exp(math.log(estimate) + Z * se),
    )


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


def create_term_table(connection: sqlite3.Connection, membership) -> None:
    connection.executescript(
        "drop table if exists temp.smq_pt_category;"
        "create temp table smq_pt_category("
        "pt_name text not null, category text not null, scope text not null,"
        "primary key(pt_name, category, scope));"
    )
    rows = []
    for pt_name, item in membership.items():
        for category in item["categories"]:
            for scope in item["scopes"]:
                rows.append((pt_name, category, scope))
    connection.executemany("insert into smq_pt_category values (?,?,?)", rows)


def create_faers_smq_events(connection: sqlite3.Connection) -> None:
    log("FAERS: attaching read-only raw and latest-version databases")
    connection.execute("attach database ? as raw", (str(RAW_FAERS),))
    connection.execute("attach database ? as latest", (str(LATEST_CACHE),))
    plan = connection.execute(
        "explain query plan select r.primaryid from raw.reac r "
        "join smq_pt_category t on t.pt_name=r.pt_name "
        "join latest.latest_demo l on l.primaryid=r.primaryid and l.source_year=r.source_year and l.source_quarter=r.source_quarter limit 1"
    ).fetchall()
    log("FAERS query plan: " + " | ".join(str(row) for row in plan))
    log("FAERS: constructing report-category hits from licensed SMQ membership")
    connection.executescript(
        "drop table if exists temp.smq_category_hits;"
        "create temp table smq_category_hits as "
        "select r.primaryid, "
        "max(case when t.scope='2' then 1 else 0 end) as narrow_hit, "
        "max(case when t.category='A' then 1 else 0 end) as has_a, "
        "max(case when t.category='B' then 1 else 0 end) as has_b, "
        "max(case when t.category='C' then 1 else 0 end) as has_c, "
        "max(case when t.category='D' then 1 else 0 end) as has_d "
        "from raw.reac r "
        "join smq_pt_category t on t.pt_name=r.pt_name "
        "join latest.latest_demo l on l.primaryid=r.primaryid and l.source_year=r.source_year and l.source_quarter=r.source_quarter "
        "group by r.primaryid;"
        "create unique index temp.idx_smq_category_hits_pid on smq_category_hits(primaryid);"
        "drop table if exists temp.smq_events;"
        "create temp table smq_events as "
        "select primaryid, narrow_hit, "
        "case when has_a=1 or (has_b=1 and has_c=1) "
        "or (has_d=1 and (has_b=1 or has_c=1)) then 1 else 0 end as broad_hit "
        "from smq_category_hits;"
        "create unique index temp.idx_smq_events_pid on smq_events(primaryid);"
    )
    log("FAERS: SMQ event table complete")


def create_context(connection: sqlite3.Connection, database: str | None = None) -> None:
    predicate = " OR ".join(f"any_{name}=1" for name in CONTEXT_CLASSES)
    connection.execute("drop table if exists temp.context_units")
    if database is None:
        connection.execute(
            f"create temp table context_units as select primaryid as unit_id "
            f"from context_flags where {predicate}"
        )
    else:
        connection.execute(
            f"create temp table context_units as select report_id as unit_id "
            f"from report_level_flags where database=? and ({predicate})",
            (database,),
        )
    connection.execute(
        "create unique index temp.idx_context_units_id on context_units(unit_id)"
    )


def estimate_rows(
    connection: sqlite3.Connection,
    database: str,
    id_column: str,
    database_where: str,
    database_parameters: tuple,
) -> list[dict[str, object]]:
    rows = []
    for scope, outcome_column in (("narrow", "narrow_hit"), ("broad_algorithm", "broad_hit")):
        for drug_class, agent in AGENTS:
            class_total = connection.execute(
                f"select count(distinct h.{id_column}) from class_hits h "
                f"join context_units x on x.unit_id=h.{id_column} "
                f"where {database_where} h.drug_class=? and h.suspect_role=1",
                database_parameters + (drug_class,),
            ).fetchone()[0]
            class_core = connection.execute(
                f"select count(distinct h.{id_column}) from class_hits h "
                f"join context_units x on x.unit_id=h.{id_column} "
                f"join smq_events e on e.{id_column}=h.{id_column} "
                f"where {database_where} h.drug_class=? and h.suspect_role=1 "
                f"and e.{outcome_column}=1",
                database_parameters + (drug_class,),
            ).fetchone()[0]
            agent_total = connection.execute(
                f"select count(distinct h.{id_column}) from class_hits h "
                f"join context_units x on x.unit_id=h.{id_column} "
                f"where {database_where} h.drug_class=? and h.generic_or_group=? "
                f"and h.suspect_role=1",
                database_parameters + (drug_class, agent),
            ).fetchone()[0]
            agent_core = connection.execute(
                f"select count(distinct h.{id_column}) from class_hits h "
                f"join context_units x on x.unit_id=h.{id_column} "
                f"join smq_events e on e.{id_column}=h.{id_column} "
                f"where {database_where} h.drug_class=? and h.generic_or_group=? "
                f"and h.suspect_role=1 and e.{outcome_column}=1",
                database_parameters + (drug_class, agent),
            ).fetchone()[0]
            class_noncore = class_total - class_core
            agent_noncore = agent_total - agent_core
            a, b = agent_core, class_core - agent_core
            c, d = agent_noncore, class_noncore - agent_noncore
            p1, p0 = a / (a + b), c / (c + d)
            low, high = newcombe(a, a + b, c, c + d)
            estimate, ror_low, ror_high = ror_ci(a, b, c, d)
            rows.append(
                {
                    "database": database,
                    "smq_scope": scope,
                    "drug_class": drug_class,
                    "agent": agent,
                    "class_core_n": a + b,
                    "class_noncore_n": c + d,
                    "a": a,
                    "b": b,
                    "c": c,
                    "d": d,
                    "core_share": p1,
                    "noncore_share": p0,
                    "difference_pp": 100 * (p1 - p0),
                    "difference_ci_low_pp": 100 * low,
                    "difference_ci_high_pp": 100 * high,
                    "ror": estimate,
                    "ror_ci_low": ror_low,
                    "ror_ci_high": ror_high,
                    "sparse_display": int(a < 5 or c < 5),
                }
            )
    return rows


def create_faers_country_context(connection: sqlite3.Connection) -> None:
    connection.executescript(
        "drop table if exists temp.faers_context_country;"
        "create temp table faers_context_country as "
        "select x.unit_id as primaryid, "
        "case when trim(coalesce(d.occr_country,''))='GB' then 'GB' "
        "when trim(coalesce(d.occr_country,''))='' then 'Missing OCCR_COUNTRY' "
        "else 'Non-GB observed' end as country_stratum "
        "from context_units x join latest.latest_demo d on d.primaryid=x.unit_id;"
        "create unique index temp.idx_faers_context_country_pid "
        "on faers_context_country(primaryid);"
    )


def estimate_faers_country_rows(
    connection: sqlite3.Connection,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    rows: list[dict[str, object]] = []
    for scope, outcome_column in (("narrow", "narrow_hit"), ("broad_algorithm", "broad_hit")):
        for drug_class, agent in AGENTS:
            for country in COUNTRY_STRATA:
                class_total = connection.execute(
                    "select count(distinct h.primaryid) from class_hits h "
                    "join faers_context_country x on x.primaryid=h.primaryid "
                    "where h.drug_class=? and h.suspect_role=1 and x.country_stratum=?",
                    (drug_class, country),
                ).fetchone()[0]
                class_core = connection.execute(
                    f"select count(distinct h.primaryid) from class_hits h "
                    f"join faers_context_country x on x.primaryid=h.primaryid "
                    f"join smq_events e on e.primaryid=h.primaryid "
                    f"where h.drug_class=? and h.suspect_role=1 "
                    f"and x.country_stratum=? and e.{outcome_column}=1",
                    (drug_class, country),
                ).fetchone()[0]
                agent_total = connection.execute(
                    "select count(distinct h.primaryid) from class_hits h "
                    "join faers_context_country x on x.primaryid=h.primaryid "
                    "where h.drug_class=? and h.generic_or_group=? "
                    "and h.suspect_role=1 and x.country_stratum=?",
                    (drug_class, agent, country),
                ).fetchone()[0]
                agent_core = connection.execute(
                    f"select count(distinct h.primaryid) from class_hits h "
                    f"join faers_context_country x on x.primaryid=h.primaryid "
                    f"join smq_events e on e.primaryid=h.primaryid "
                    f"where h.drug_class=? and h.generic_or_group=? "
                    f"and h.suspect_role=1 and x.country_stratum=? "
                    f"and e.{outcome_column}=1",
                    (drug_class, agent, country),
                ).fetchone()[0]
                class_noncore = class_total - class_core
                agent_noncore = agent_total - agent_core
                a, b = agent_core, class_core - agent_core
                c, d = agent_noncore, class_noncore - agent_noncore
                p1, p0 = a / (a + b), c / (c + d)
                low, high = newcombe(a, a + b, c, c + d)
                estimate, ror_low, ror_high = ror_ci(a, b, c, d)
                rows.append(
                    {
                        "database": "FAERS",
                        "smq_scope": scope,
                        "drug_class": drug_class,
                        "agent": agent,
                        "country_field": "OCCR_COUNTRY",
                        "country_stratum": country,
                        "class_core_n": a + b,
                        "class_noncore_n": c + d,
                        "a": a,
                        "b": b,
                        "c": c,
                        "d": d,
                        "core_share": p1,
                        "noncore_share": p0,
                        "difference_pp": 100 * (p1 - p0),
                        "difference_ci_low_pp": 100 * low,
                        "difference_ci_high_pp": 100 * high,
                        "ror": estimate,
                        "ror_ci_low": ror_low,
                        "ror_ci_high": ror_high,
                        "sparse_display": int(a < 5 or c < 5),
                    }
                )

    summaries: list[dict[str, object]] = []
    for scope, _ in (("narrow", "narrow_hit"), ("broad_algorithm", "broad_hit")):
        for drug_class, agent in AGENTS:
            agent_rows = [
                row
                for row in rows
                if row["smq_scope"] == scope and row["agent"] == agent
            ]
            core_total = sum(int(row["class_core_n"]) for row in agent_rows)
            noncore_total = sum(int(row["class_noncore_n"]) for row in agent_rows)
            core_agent_total = sum(int(row["a"]) for row in agent_rows)
            noncore_agent_total = sum(int(row["c"]) for row in agent_rows)
            mh_num = sum(
                int(row["a"]) * int(row["d"])
                / (int(row["a"]) + int(row["b"]) + int(row["c"]) + int(row["d"]))
                for row in agent_rows
            )
            mh_den = sum(
                int(row["b"]) * int(row["c"])
                / (int(row["a"]) + int(row["b"]) + int(row["c"]) + int(row["d"]))
                for row in agent_rows
            )
            core_to_noncore = 100 * sum(
                int(row["class_noncore_n"]) / noncore_total
                * (float(row["core_share"]) - float(row["noncore_share"]))
                for row in agent_rows
            )
            noncore_to_core = 100 * sum(
                int(row["class_core_n"]) / core_total
                * (float(row["core_share"]) - float(row["noncore_share"]))
                for row in agent_rows
            )
            by_country = {str(row["country_stratum"]): row for row in agent_rows}
            gb_difference = float(by_country["GB"]["difference_pp"])
            non_gb_difference = float(by_country["Non-GB observed"]["difference_pp"])
            summaries.append(
                {
                    "database": "FAERS",
                    "smq_scope": scope,
                    "drug_class": drug_class,
                    "agent": agent,
                    "crude_difference_pp": 100
                    * (core_agent_total / core_total - noncore_agent_total / noncore_total),
                    "mh_ror": mh_num / mh_den if mh_den else None,
                    "core_standardized_to_noncore_country_mix_pp": core_to_noncore,
                    "noncore_standardized_to_core_country_mix_pp": noncore_to_core,
                    "gb_difference_pp": gb_difference,
                    "non_gb_observed_difference_pp": non_gb_difference,
                    "gb_vs_non_gb_sign_reversal": int(gb_difference * non_gb_difference < 0),
                }
            )
    return rows, summaries


def canada_event_units(membership):
    hits: dict[str, dict[str, object]] = defaultdict(
        lambda: {"narrow": False, "categories": set()}
    )
    for row in read_dollar(CANADA_REACTIONS):
        if len(row) <= 5:
            continue
        report_id = canonical_canada_report_id(row[1])
        item = membership.get(row[5])
        if not report_id or item is None:
            continue
        if "2" in item["scopes"]:
            hits[report_id]["narrow"] = True
        hits[report_id]["categories"].update(item["categories"])
    events = []
    for report_id, item in hits.items():
        categories = item["categories"]
        broad = (
            "A" in categories
            or ("B" in categories and "C" in categories)
            or ("D" in categories and ("B" in categories or "C" in categories))
        )
        events.append((report_id, int(item["narrow"]), int(broad)))
    return events


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    log("Loading licensed MedDRA v28.0 SMQ membership")
    metadata, membership, content_rows = load_smq()
    category_pt_counts = {
        category: sum(category in item["categories"] for item in membership.values())
        for category in ("A", "B", "C", "D")
    }
    # MedDRA distribution format v28.0: scope 1=broad, scope 2=narrow.
    narrow_pt_count = sum("2" in item["scopes"] for item in membership.values())

    faers = sqlite3.connect(f"file:{FAERS_FLAGS.as_posix()}?mode=ro", uri=True)
    create_term_table(faers, membership)
    create_faers_smq_events(faers)
    create_context(faers)
    faers_narrow_units = faers.execute(
        "select count(*) from smq_events where narrow_hit=1"
    ).fetchone()[0]
    faers_broad_units = faers.execute(
        "select count(*) from smq_events where broad_hit=1"
    ).fetchone()[0]
    result_rows = estimate_rows(faers, "FAERS", "primaryid", "1=1 and", ())
    create_faers_country_context(faers)
    country_rows, country_summaries = estimate_faers_country_rows(faers)
    faers.close()

    log("Canada: streaming reaction rows and applying the official SMQ algorithm")
    canada_events = canada_event_units(membership)
    canada = sqlite3.connect(f"file:{CJ_FLAGS.as_posix()}?mode=ro", uri=True)
    canada.execute(
        "create temp table smq_events(report_id text primary key, narrow_hit integer, broad_hit integer)"
    )
    canada.executemany("insert into smq_events values (?,?,?)", canada_events)
    create_context(canada, "Canada")
    result_rows.extend(
        estimate_rows(canada, "Canada", "report_id", "h.database=? and", ("Canada",))
    )
    canada.close()

    output_csv = OUT / "official_anaphylactic_reaction_smq_v28_pilot_v0_26.csv"
    write_csv(output_csv, result_rows)
    country_csv = OUT / "official_anaphylactic_reaction_smq_v28_country_strata_v0_26.csv"
    country_summary_csv = OUT / "official_anaphylactic_reaction_smq_v28_country_summary_v0_26.csv"
    write_csv(country_csv, country_rows)
    write_csv(country_summary_csv, country_summaries)
    manifest = {
        "build": "v0.26 licensed official Anaphylactic reaction SMQ sensitivity",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "smq_metadata": metadata,
        "licensed_membership_not_redistributed": True,
        "smq_content_active_rows": content_rows,
        "unique_pt_members": len(membership),
        "narrow_unique_pt_members": narrow_pt_count,
        "category_unique_pt_counts": category_pt_counts,
        "licensed_source_files": {
            "smq_list": {"file_name": SMQ_LIST.name, "sha256": sha256(SMQ_LIST)},
            "smq_content": {"file_name": SMQ_CONTENT.name, "sha256": sha256(SMQ_CONTENT)},
            "pt": {"file_name": PT_ASC.name, "sha256": sha256(PT_ASC)},
            "llt": {"file_name": LLT_ASC.name, "sha256": sha256(LLT_ASC)},
        },
        "faers_retained_event_units": {
            "narrow": faers_narrow_units,
            "broad_algorithm": faers_broad_units,
        },
        "canada_event_units": {
            "narrow": sum(row[1] for row in canada_events),
            "broad_algorithm": sum(row[2] for row in canada_events),
        },
        "jader_status": (
            "not run as official SMQ: the local JADER extract stores Japanese PT strings "
            "and no licensed MedDRA/J v28.1 PT-code crosswalk was found"
        ),
        "output_csv": output_csv.name,
        "output_sha256": sha256(output_csv),
        "country_strata_output_csv": country_csv.name,
        "country_strata_output_sha256": sha256(country_csv),
        "country_summary_output_csv": country_summary_csv.name,
        "country_summary_output_sha256": sha256(country_summary_csv),
    }
    manifest_path = OUT / "official_smq_manifest_v0_26.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    log(
        f"PASS: FAERS narrow={faers_narrow_units:,}, broad={faers_broad_units:,}; "
        f"Canada narrow={manifest['canada_event_units']['narrow']:,}, "
        f"broad={manifest['canada_event_units']['broad_algorithm']:,}"
    )
    teicoplanin_summary = {
        row["smq_scope"]: row
        for row in country_summaries
        if row["agent"] == "teicoplanin"
    }
    log(
        "FAERS teicoplanin GB/non-GB differences: "
        f"narrow={teicoplanin_summary['narrow']['gb_difference_pp']:.2f}/"
        f"{teicoplanin_summary['narrow']['non_gb_observed_difference_pp']:.2f} pp; "
        f"broad={teicoplanin_summary['broad_algorithm']['gb_difference_pp']:.2f}/"
        f"{teicoplanin_summary['broad_algorithm']['non_gb_observed_difference_pp']:.2f} pp"
    )


if __name__ == "__main__":
    main()
