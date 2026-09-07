from __future__ import annotations

"""FAERS outcome-concept sensitivity using the existing broad PT group.

This is intentionally a construct-expansion sensitivity, not an SMQ and not a
replacement for the narrow three-PT anaphylaxis outcome.
"""

import csv
import math
import os
import sqlite3
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
DB = Path(os.environ.get("PV_FAERS_FLAGS_DB", REPO_ROOT / "work" / "faers_periop_hypersensitivity_prelim.sqlite"))
RAW_DB = Path(os.environ.get("PV_FAERS_RAW_DB", REPO_ROOT / "external_data" / "faers.sqlite"))
LATEST_CACHE = Path(os.environ.get("PV_FAERS_CACHE_DB", REPO_ROOT / "external_data" / "faers_cache.sqlite"))
OUT = Path(os.environ.get("PV_ANALYSIS_OUTPUT_DIR", REPO_ROOT / "aggregate_outputs"))
PAIRS = (("NMBA", "rocuronium"), ("NMBA", "succinylcholine"), ("antibiotic_anchor", "cefazolin"), ("antibiotic_anchor", "teicoplanin"))
EVENTS = (
    ("anaphylaxis_core", "event_anaphylaxis_core", "Narrow three-PT anaphylaxis outcome"),
    ("periop_hypersensitivity_broad", "event_periop_hypersensitivity_broad", "Predefined broad PT expansion; construct sensitivity only"),
    ("bronchospasm_only", None, "Bronchospasm PT only; phenotype sensitivity for the succinylcholine context"),
)


def ror_ci(a: int, b: int, c: int, d: int) -> tuple[float | None, float | None, float | None]:
    if min(a, b, c, d) == 0:
        return None, None, None
    value = a * d / (b * c)
    se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return value, math.exp(math.log(value) - 1.96 * se), math.exp(math.log(value) + 1.96 * se)


def main() -> None:
    con = sqlite3.connect(DB)
    con.execute("attach database ? as raw", (str(RAW_DB),))
    con.execute("attach database ? as latest", (str(LATEST_CACHE),))
    predicate = " OR ".join(f"any_{name}=1" for name in (
        "sedative_hypnotic", "opioid", "volatile_anesthetic", "local_anesthetic", "vasoactive_rescue"
    ))
    con.execute("drop table if exists temp.context")
    con.execute(f"create temp table context as select primaryid from context_flags where {predicate}")
    con.execute("create unique index temp.idx_context_primaryid on context(primaryid)")
    con.execute("drop table if exists temp.bronchospasm_events")
    con.execute(
        """create temp table bronchospasm_events as
            select distinct r.primaryid
            from raw.reac r indexed by idx_reac_pt_name
            join latest.latest_demo l on l.primaryid=r.primaryid
                                   and l.source_year=r.source_year
                                   and l.source_quarter=r.source_quarter
            where r.pt_name='Bronchospasm'"""
    )
    con.execute("create unique index temp.idx_bronchospasm_events_primaryid on bronchospasm_events(primaryid)")
    rows: list[dict[str, object]] = []
    for event_group, event_column, note in EVENTS:
        for drug_class, agent in PAIRS:
            class_total = con.execute("select count(distinct h.primaryid) from class_hits h join context x on x.primaryid=h.primaryid where h.drug_class=? and h.suspect_role=1", (drug_class,)).fetchone()[0]
            if event_column is None:
                event_join = "join bronchospasm_events e on e.primaryid=h.primaryid"
                event_predicate = "1=1"
            else:
                event_join = "join event_agg e on e.primaryid=h.primaryid"
                event_predicate = f"e.{event_column}=1"
            class_event = con.execute(f"select count(distinct h.primaryid) from class_hits h join context x on x.primaryid=h.primaryid {event_join} where h.drug_class=? and h.suspect_role=1 and {event_predicate}", (drug_class,)).fetchone()[0]
            agent_total = con.execute("select count(distinct h.primaryid) from class_hits h join context x on x.primaryid=h.primaryid where h.drug_class=? and h.generic_or_group=? and h.suspect_role=1", (drug_class, agent)).fetchone()[0]
            agent_event = con.execute(f"select count(distinct h.primaryid) from class_hits h join context x on x.primaryid=h.primaryid {event_join} where h.drug_class=? and h.generic_or_group=? and h.suspect_role=1 and {event_predicate}", (drug_class, agent)).fetchone()[0]
            a, b, c, d = agent_event, class_event - agent_event, agent_total - agent_event, (class_total - class_event) - (agent_total - agent_event)
            sparse = min(a, c) < 5
            estimate, low, high = ror_ci(a, b, c, d) if not sparse else (None, None, None)
            rows.append({
                "database": "FAERS", "event_group": event_group, "event_note": note,
                "drug_class": drug_class, "agent": agent,
                "context_class_suspect_n": class_total, "event_class_suspect_n": class_event,
                "a_agent_event": a, "b_other_agent_status_event": b,
                "c_agent_nonevent": c, "d_other_agent_status_nonevent": d,
                "difference_pp": (a / class_event - c / (class_total - class_event)) * 100,
                "ror": estimate, "ror_ci_low": low, "ror_ci_high": high,
                "sparse": int(sparse),
                "sparse_display_rule": "ROR withheld when the agent-event or agent-nonevent cell is <5",
            })
    con.close()
    OUT.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    with (OUT / "faers_context_broad_pt_sensitivity_v0_20.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)


if __name__ == "__main__":
    main()
