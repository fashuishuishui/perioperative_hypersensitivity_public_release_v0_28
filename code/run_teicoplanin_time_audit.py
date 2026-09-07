from __future__ import annotations

"""Report, rather than presume, temporal concentration of the teicoplanin example."""

import csv
import os
import sqlite3
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
FLAG_DB = Path(os.environ.get("PV_FAERS_FLAGS_DB", REPO_ROOT / "work" / "faers_periop_hypersensitivity_prelim.sqlite"))
RAW_DB = Path(os.environ.get("PV_FAERS_RAW_DB", REPO_ROOT / "external_data" / "faers.sqlite"))
LATEST_CACHE = Path(os.environ.get("PV_FAERS_CACHE_DB", REPO_ROOT / "external_data" / "faers_cache.sqlite"))
OUT = Path(os.environ.get("PV_ANALYSIS_OUTPUT_DIR", REPO_ROOT / "aggregate_outputs"))


def main() -> None:
    con = sqlite3.connect(FLAG_DB)
    con.execute("attach database ? as raw", (str(RAW_DB),))
    con.execute("attach database ? as latest", (str(LATEST_CACHE),))
    predicate = " OR ".join(f"c.any_{name}=1" for name in (
        "sedative_hypnotic", "opioid", "volatile_anesthetic", "local_anesthetic", "vasoactive_rescue"
    ))
    rows = con.execute(
        f"""select d.source_year, count(distinct h.primaryid) as report_n
        from class_hits h
        join event_agg e on e.primaryid=h.primaryid
        join context_flags c on c.primaryid=h.primaryid
        join latest.latest_demo d on d.primaryid=h.primaryid
        where h.drug_class='antibiotic_anchor' and h.generic_or_group='teicoplanin'
          and h.suspect_role=1 and e.event_anaphylaxis_core=1 and ({predicate})
        group by d.source_year order by cast(d.source_year as integer)"""
    ).fetchall()
    con.close()
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "faers_teicoplanin_context_core_year_v0_20.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("year", "context_core_teicoplanin_reports_n", "scope"))
        writer.writeheader()
        writer.writerows({"year": year, "context_core_teicoplanin_reports_n": count, "scope": "FAERS target-independent medication context; core anaphylaxis; suspect role"} for year, count in rows)


if __name__ == "__main__":
    main()
