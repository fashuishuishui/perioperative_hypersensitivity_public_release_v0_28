"""Burden and country slices of the same locked duplicate-candidate rules."""
import pandas as pd
from settings import ROOT
from run_targeted_diagnostics import metric,PAIRS,write
def main():
 f=pd.read_pickle(ROOT/'work/fingerprint_features_private.pkl');h=pd.read_pickle(ROOT/'work/normalized_role_hits_private.pkl');h=h[h.has_ps.eq(1)|h.has_ss.eq(1)];core=set(f.loc[f.core.eq(1),'primaryid']);ctx=set(f.loc[f.in_context,'primaryid']);out=[];geo=[]
 for cl in ['NMBA','antibiotic_anchor']:
  p=set(h.loc[h.drug_class.eq(cl),'primaryid'])&ctx;g=f[f.primaryid.isin(p)]
  for rule,keys in [('unmerged',None),('strict_demographic_date_PT',['strict_fp']),('strict_plus_role_signature',['strict_fp','role_signature'])]:
   kept=p if keys is None else set(g.loc[~g.complete_fingerprint,'primaryid'])|set(g[g.complete_fingerprint].drop_duplicates(keys,keep='last').primaryid)
   for ccl,ag in PAIRS:
    if ccl!=cl:continue
    t=set(h.loc[h.drug_class.eq(cl)&h.agent.eq(ag),'primaryid'])
    for band in ['1–5','6–10','11–20','21+']:
     e=kept&set(f.loc[f.burden.eq(band),'primaryid']);out.append(dict(agent=ag,fingerprint_rule=rule,burden_stratum=band,**metric(len(e&core&t),len((e&core)-t),len((e-core)&t),len((e-core)-t))))
   for country in ['CA','GB','JP']:
    e=kept&set(f.loc[f.country.eq(country),'primaryid']);geo.append(dict(drug_class=cl,country=country,fingerprint_rule=rule,core_n=len(e&core),noncore_n=len(e-core),retained_n=len(e)))
 write('fingerprint_burden_sensitivity.csv',out);write('fingerprint_country_counts.csv',geo)
 print(pd.DataFrame(out).query("agent=='rocuronium' and burden_stratum=='21+'")[["agent","fingerprint_rule","a","b","c","d","ror"]].to_string(index=False))
if __name__=='__main__':main()
