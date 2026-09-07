from __future__ import annotations

"""Country-stratified and PS/SS role probes for the v0.16 audit fork."""

import csv
import math
import os
import re
import sqlite3
import unicodedata
from collections import defaultdict
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
WORK = REPO_ROOT / "work"
FLAG_DB = Path(os.environ.get("PV_FAERS_FLAGS_DB", WORK / "faers_periop_hypersensitivity_prelim.sqlite"))
RAW_DB = Path(os.environ.get("PV_FAERS_RAW_DB", REPO_ROOT / "external_data" / "faers.sqlite"))
LATEST_CACHE = Path(os.environ.get("PV_FAERS_CACHE_DB", REPO_ROOT / "external_data" / "faers_cache.sqlite"))
CONFIG = REPO_ROOT / "config" / "drug_class_aliases.csv"
OUT = Path(os.environ.get("PV_ANALYSIS_OUTPUT_DIR", REPO_ROOT / "aggregate_outputs"))
PAIRS = (("NMBA", "rocuronium"), ("NMBA", "succinylcholine"), ("antibiotic_anchor", "cefazolin"), ("antibiotic_anchor", "teicoplanin"))
CLASSES = ("NMBA", "antibiotic_anchor", "chlorhexidine", "dye_anchor")


def write_csv(name: str, rows: list[dict[str, object]]) -> None:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / name).open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)


def ror(a: int, b: int, c: int, d: int) -> float | None:
    if min(a, b, c, d) == 0:
        return None
    return a * d / (b * c)


def gammaincc(a: float, x: float) -> float:
    if x <= 0:
        return 1.0
    if x < a + 1:
        term = total = 1.0 / a
        n = 1
        while abs(term) > abs(total) * 1e-14:
            term *= x / (a + n); total += term; n += 1
        return 1.0 - total * math.exp(-x + a * math.log(x) - math.lgamma(a))
    b, c, d, h = x + 1 - a, 1e300, 1 / (x + 1 - a), 1 / (x + 1 - a)
    for i in range(1, 10000):
        an = -i * (i - a); b += 2; d = an * d + b
        if abs(d) < 1e-300: d = 1e-300
        c = b + an / c
        if abs(c) < 1e-300: c = 1e-300
        d = 1 / d; delta = d * c; h *= delta
        if abs(delta - 1) < 1e-14: break
    return math.exp(-x + a * math.log(x) - math.lgamma(a)) * h


def bd_moments(a: int, b: int, c: int, d: int, theta: float) -> tuple[float, float]:
    r1, r0, c1 = a + b, c + d, a + c
    low, high = max(0, c1 - r0), min(r1, c1)
    values = []
    for x in range(low, high + 1):
        logw = math.lgamma(r1 + 1) - math.lgamma(x + 1) - math.lgamma(r1 - x + 1)
        logw += math.lgamma(r0 + 1) - math.lgamma(c1 - x + 1) - math.lgamma(r0 - c1 + x + 1) + x * math.log(theta)
        values.append((x, logw))
    maximum = max(weight for _, weight in values)
    scaled = [(x, math.exp(weight - maximum)) for x, weight in values]
    total = sum(weight for _, weight in scaled)
    mean = sum(x * weight for x, weight in scaled) / total
    return mean, sum((x - mean) ** 2 * weight for x, weight in scaled) / total


def norm(value: str | None) -> str:
    return unicodedata.normalize("NFKC", value or "").casefold().strip()


def matches(value: str | None, alias: str) -> bool:
    return re.search(rf"(?<![a-z0-9]){re.escape(norm(alias))}(?![a-z0-9])", norm(value)) is not None


def create_context(con: sqlite3.Connection) -> None:
    predicate = " OR ".join(f"c.any_{name}=1" for name in ("sedative_hypnotic", "opioid", "volatile_anesthetic", "local_anesthetic", "vasoactive_rescue"))
    con.execute("drop table if exists temp.ctx")
    con.execute(f"""create temp table ctx as
        select c.primaryid, coalesce(e.event_anaphylaxis_core,0) as core,
               case when trim(coalesce(d.occr_country,''))='GB' then 'GB'
                    when trim(coalesce(d.occr_country,''))='' then 'Missing OCCR_COUNTRY'
                    else 'Non-GB observed' end as country_stratum
        from context_flags c left join event_agg e on e.primaryid=c.primaryid
        join latest.latest_demo d indexed by idx_latest_demo_primaryid on d.primaryid=c.primaryid
        where {predicate}""")
    con.execute("create unique index temp.idx_ctx_primaryid on ctx(primaryid)")


def country_agent_audit(con: sqlite3.Connection) -> None:
    rows: list[dict[str, object]] = []
    for drug_class, agent in PAIRS:
        for country in ("GB", "Non-GB observed", "Missing OCCR_COUNTRY"):
            class_core, class_noncore, agent_core, agent_noncore = [int(value) for value in con.execute(
                """select
                    count(distinct case when x.core=1 then h.primaryid end),
                    count(distinct case when x.core=0 then h.primaryid end),
                    count(distinct case when x.core=1 and h.generic_or_group=? then h.primaryid end),
                    count(distinct case when x.core=0 and h.generic_or_group=? then h.primaryid end)
                from class_hits h join ctx x on x.primaryid=h.primaryid
                where h.drug_class=? and h.suspect_role=1 and x.country_stratum=?""",
                (agent, agent, drug_class, country),
            ).fetchone()]
            a, b, c, d = agent_core, class_core - agent_core, agent_noncore, class_noncore - agent_noncore
            rows.append({
                "database": "FAERS", "country_field": "OCCR_COUNTRY", "country_stratum": country,
                "drug_class": drug_class, "agent": agent,
                "class_core_n": class_core, "class_noncore_n": class_noncore,
                "agent_core_n": agent_core, "agent_noncore_n": agent_noncore,
                "core_agent_share": agent_core / class_core if class_core else None,
                "noncore_agent_share": agent_noncore / class_noncore if class_noncore else None,
                "percentage_point_difference": (agent_core / class_core - agent_noncore / class_noncore) * 100 if class_core and class_noncore else None,
                "a": a, "b": b, "c": c, "d": d, "ror_no_correction": ror(a, c, b, d),
                "note": "Country-stratified descriptive contrast; country is not an exposure denominator.",
            })
    write_csv("faers_agent_occr_country_stratified_v0_16.csv", rows)
    summaries: list[dict[str, object]] = []
    for drug_class, agent in PAIRS:
        agent_rows = [row for row in rows if row["agent"] == agent]
        tables = [(int(row["a"]), int(row["b"]), int(row["c"]), int(row["d"])) for row in agent_rows]
        mh_num = sum(a * d / (a + b + c + d) for a, b, c, d in tables)
        mh_den = sum(b * c / (a + b + c + d) for a, b, c, d in tables)
        mh = mh_num / mh_den if mh_den else None
        bd = 0.0
        if mh:
            for a, b, c, d in tables:
                expected, variance = bd_moments(a, b, c, d, mh)
                if variance > 0: bd += (a - expected) ** 2 / variance
        core_total, noncore_total = sum(int(row["class_core_n"]) for row in agent_rows), sum(int(row["class_noncore_n"]) for row in agent_rows)
        core_to_noncore = sum((int(row["class_noncore_n"]) / noncore_total) * (float(row["core_agent_share"]) - float(row["noncore_agent_share"])) for row in agent_rows) * 100
        noncore_to_core = sum((int(row["class_core_n"]) / core_total) * (float(row["core_agent_share"]) - float(row["noncore_agent_share"])) for row in agent_rows) * 100
        summaries.append({"database": "FAERS", "country_field": "OCCR_COUNTRY", "drug_class": drug_class, "agent": agent, "crude_percentage_point_difference": sum(int(row["agent_core_n"]) for row in agent_rows) / core_total * 100 - sum(int(row["agent_noncore_n"]) for row in agent_rows) / noncore_total * 100, "mh_ror": mh, "breslow_day_statistic": bd, "breslow_day_df": len(tables)-1, "breslow_day_p_value": gammaincc((len(tables)-1)/2, bd/2), "core_standardized_to_noncore_country_mix_pp": core_to_noncore, "noncore_standardized_to_core_country_mix_pp": noncore_to_core, "missing_occr_core_n": next(int(row["class_core_n"]) for row in agent_rows if row["country_stratum"] == "Missing OCCR_COUNTRY"), "missing_occr_noncore_n": next(int(row["class_noncore_n"]) for row in agent_rows if row["country_stratum"] == "Missing OCCR_COUNTRY")})
    write_csv("faers_agent_occr_country_mh_standardization_v0_16.csv", summaries)


def ps_role_audit(con: sqlite3.Connection) -> None:
    aliases: list[dict[str, str]] = []
    with CONFIG.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            if row["database_scope"] in {"FAERS", "FAERS_Canada"} and row["drug_class"] in CLASSES:
                aliases.append(row)
    candidate_ids = [row[0] for row in con.execute("""select distinct h.primaryid from class_hits h
        join ctx x on x.primaryid=h.primaryid where h.suspect_role=1 and h.drug_class in ('NMBA','antibiotic_anchor','chlorhexidine','dye_anchor')""")]
    con.execute("drop table if exists temp.ps_flags")
    con.execute("create temp table ps_flags(primaryid text, drug_class text, generic_or_group text, has_ps integer, has_ss integer, primary key(primaryid,drug_class,generic_or_group))")
    batch = 5000
    flags: dict[tuple[str, str, str], list[int]] = defaultdict(lambda: [0, 0])
    for start in range(0, len(candidate_ids), batch):
        chunk = candidate_ids[start:start + batch]
        placeholders = ",".join("?" for _ in chunk)
        for primaryid, role, drugname, prod_ai in con.execute(
            f"""select d.primaryid, d.role_cod, d.drugname, d.prod_ai
                from raw.drug d indexed by idx_drug_primaryid
                join latest.latest_demo l on l.primaryid=d.primaryid
                                        and l.source_year=d.source_year
                                        and l.source_quarter=d.source_quarter
                where d.primaryid in ({placeholders})""",
            chunk,
        ):
            for cfg in aliases:
                if matches(drugname, cfg["alias"]) or matches(prod_ai, cfg["alias"]):
                    key = (primaryid, cfg["drug_class"], cfg["generic_or_group"])
                    if role == "PS": flags[key][0] = 1
                    if role == "SS": flags[key][1] = 1
    con.executemany("insert into ps_flags values (?,?,?,?,?)", [(primaryid, cls, generic, ps, ss) for (primaryid, cls, generic), (ps, ss) in flags.items()])
    con.execute("drop table if exists temp.ps_class_flags")
    con.execute("create temp table ps_class_flags as select primaryid, drug_class, max(has_ps) as has_ps, max(has_ss) as has_ss from ps_flags group by primaryid, drug_class")
    role_rows = []
    for drug_class in CLASSES:
        for core in (0, 1):
            for profile in ("PS present", "SS present without PS", "suspect class hit without reproduced PS/SS alias row"):
                if profile == "PS present": condition = "p.has_ps=1"
                elif profile == "SS present without PS": condition = "p.has_ss=1 and p.has_ps=0"
                else: condition = "coalesce(p.has_ps,0)=0 and coalesce(p.has_ss,0)=0"
                count = con.execute(f"""select count(distinct h.primaryid) from class_hits h join ctx x on x.primaryid=h.primaryid
                    left join ps_class_flags p on p.primaryid=h.primaryid and p.drug_class=h.drug_class
                    where h.drug_class=? and h.suspect_role=1 and x.core=? and ({condition})""", (drug_class, core)).fetchone()[0]
                role_rows.append({"database": "FAERS", "drug_class": drug_class, "outcome_group": "core" if core else "non-core", "role_profile": profile, "report_n": count})
    write_csv("faers_context_suspect_role_profile_v0_16.csv", role_rows)
    class_rows = []
    for drug_class in CLASSES:
        a, b, c, d = [int(value) for value in con.execute("""select
            count(distinct case when x.core=1 and p.has_ps=1 then x.primaryid end),
            count(distinct case when x.core=0 and p.has_ps=1 then x.primaryid end),
            count(distinct case when x.core=1 and coalesce(p.has_ps,0)=0 then x.primaryid end),
            count(distinct case when x.core=0 and coalesce(p.has_ps,0)=0 then x.primaryid end)
            from ctx x left join ps_class_flags p on p.primaryid=x.primaryid and p.drug_class=?""", (drug_class,)).fetchone()]
        class_rows.append({"database": "FAERS", "drug_class": drug_class, "role_definition": "primary suspect (PS) configured-alias reconstruction", "a": a, "b": b, "c": c, "d": d, "ror_no_correction": ror(a,b,c,d), "note": "Exploratory role-specific reconstruction; external dictionary mapping remains unresolved."})
    write_csv("faers_context_ps_only_class_ror_v0_16.csv", class_rows)
    agent_rows = []
    for drug_class, agent in PAIRS:
        a, b, c, d = [int(value) for value in con.execute("""select
            count(distinct case when x.core=1 and p.has_ps=1 then x.primaryid end),
            count(distinct case when x.core=1 and coalesce(p.has_ps,0)=0 then x.primaryid end),
            count(distinct case when x.core=0 and p.has_ps=1 then x.primaryid end),
            count(distinct case when x.core=0 and coalesce(p.has_ps,0)=0 then x.primaryid end)
            from ctx x
            join ps_class_flags cls on cls.primaryid=x.primaryid and cls.drug_class=? and cls.has_ps=1
            left join ps_flags p on p.primaryid=x.primaryid and p.drug_class=? and p.generic_or_group=?""", (drug_class, drug_class, agent)).fetchone()]
        agent_rows.append({"database": "FAERS", "drug_class": drug_class, "agent": agent, "role_definition": "primary suspect (PS) configured-alias reconstruction within PS-class reports", "a_agent_ps_core": a, "b_other_ps_agent_status_core": b, "c_agent_ps_noncore": c, "d_other_ps_agent_status_noncore": d, "ror_no_correction": ror(a,b,c,d)})
    write_csv("faers_agent_ps_only_role_contrast_v0_16.csv", agent_rows)


def main() -> None:
    con = sqlite3.connect(FLAG_DB)
    try:
        con.execute("attach database ? as raw", (str(RAW_DB),))
        con.execute("attach database ? as latest", (str(LATEST_CACHE),))
        create_context(con)
        country_agent_audit(con)
        ps_role_audit(con)
    finally:
        con.close()


if __name__ == "__main__":
    main()
