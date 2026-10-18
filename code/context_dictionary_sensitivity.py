"""Bounded expanded-context context coverage probe; primary context is unchanged."""
import sqlite3,json
import pandas as pd
from settings import ROOT,RAW
from run_targeted_diagnostics import metric,PAIRS,write
from build_drug_event_flags import prefixes
def main():
 f=pd.read_pickle(ROOT/'work/normalized_report_features_private.pkl');h=pd.read_pickle(ROOT/'work/normalized_role_hits_private.pkl');hs=h[h.has_ps.eq(1)|h.has_ss.eq(1)];ids=set(f.primaryid);core=set(f.loc[f.core.eq(1),'primaryid']);ctx=set(f.loc[f.in_context,'primaryid'])
 c=sqlite3.connect(':memory:',uri=True);c.execute('attach database ? as raw',(RAW.as_uri()+'?mode=ro',));f[['primaryid','source_year','source_quarter']].to_sql('candidate',c,index=False);c.execute('create unique index candidate_id on candidate(primaryid)')
 groups={'two_additional_brands':['DIPRIVAN','ULTIVA'],'expanded_generic_markers':['DIPRIVAN','ULTIVA','ETOMIDATE','THIOPENTAL','MORPHINE','SUFENTANIL','SUFENTA']};rows=[];counts=[]
 found={}
 for alias in groups['expanded_generic_markers']:
  hit=set()
  for field,index in [('drugname','idx_drug_drugname'),('prod_ai','idx_drug_prod_ai')]:
   for lo in prefixes(alias):
    hi=lo[:-1]+chr(ord(lo[-1])+1)
    hit.update(r[0] for r in c.execute(f'select d.primaryid from raw.drug d indexed by {index} join candidate k on k.primaryid=d.primaryid and k.source_year=d.source_year and k.source_quarter=d.source_quarter where d.{field}>=? and d.{field}<? and upper(d.{field}) like ?',(lo,hi,alias+'%')))
  found[alias]=hit;counts.append(dict(alias=alias,candidate_reports=len(hit),additional_to_primary_context=len(hit-ctx)))
 for label,aliases in groups.items():
  e=ctx|set().union(*(found[a] for a in aliases))
  for cl,ag in PAIRS:
   p=set(hs.loc[hs.drug_class.eq(cl),'primaryid'])&e;t=set(hs.loc[hs.drug_class.eq(cl)&hs.agent.eq(ag),'primaryid']);rows.append(dict(database='FAERS',agent=ag,context_rule=label,additional_class_reports=len(p-ctx),**metric(len(p&core&t),len((p&core)-t),len((p-core)&t),len((p-core)-t))))
 write('context_dictionary_coverage_sensitivity.csv',rows);write('context_dictionary_probe_hits.csv',counts);print(pd.DataFrame(rows)[['agent','context_rule','additional_class_reports','ror']].to_string(index=False))
if __name__=='__main__':main()
