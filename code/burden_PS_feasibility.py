"""PS-only class-parent feasibility within the same medication-burden strata."""
import pandas as pd
from settings import ROOT
from run_targeted_diagnostics import metric,PAIRS,write
def main():
 f=pd.read_pickle(ROOT/'work/fingerprint_features_private.pkl');h=pd.read_pickle(ROOT/'work/normalized_role_hits_private.pkl');h=h[h.has_ps.eq(1)];core=set(f.loc[f.core.eq(1),'primaryid']);ctx=set(f.loc[f.in_context,'primaryid']);rows=[]
 for cl,ag in PAIRS:
  parent=set(h.loc[h.drug_class.eq(cl),'primaryid'])&ctx;t=set(h.loc[h.drug_class.eq(cl)&h.agent.eq(ag),'primaryid'])
  for band in ['1–5','6–10','11–20','21+']:
   p=parent&set(f.loc[f.burden.eq(band),'primaryid']);rows.append(dict(agent=ag,role_scope='PS_only',burden_stratum=band,**metric(len(p&core&t),len((p&core)-t),len((p-core)&t),len((p-core)-t))))
 write('burden_PS_only_feasibility.csv',rows);print(pd.DataFrame(rows).query("burden_stratum=='21+'")[["agent","a","b","c","d","ror"]].to_string(index=False))
if __name__=='__main__':main()
