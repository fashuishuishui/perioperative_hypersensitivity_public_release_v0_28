"""Supporting report-role, context, geography and burden diagnostics."""
from pathlib import Path
from itertools import combinations
import json,sqlite3
import pandas as pd
from run_targeted_diagnostics import ROOT,RAW,LATEST,MARKERS,PAIRS,ro,metric,write

def main():
 c=ro(ROOT/'work/analysis_flags_private.sqlite')
 hits=pd.read_sql_query('select * from class_hits',c)
 f=pd.read_pickle(ROOT/'work/normalized_report_features_private.pkl')
 demo=pd.read_pickle(ROOT/'work/normalized_demo_fields_private.pkl')
 f=f.merge(demo[['primaryid','occp_cod']],on='primaryid',validate='one_to_one');ids=set(f.primaryid);core=set(r[0] for r in c.execute("select primaryid from event_hits where event_group='anaphylaxis_core'"))
 ctx=set(hits.loc[hits.drug_class.isin(MARKERS),'primaryid']);strict=set(hits.loc[hits.drug_class.isin(MARKERS[:-1]),'primaryid'])
 subsets={'PS_only':hits.has_ps.eq(1),'PS_SS':hits.has_ps.eq(1)|hits.has_ss.eq(1),'PS_SS_C':hits[['has_ps','has_ss','has_c']].max(axis=1).eq(1),'any_role':pd.Series(True,index=hits.index)}
 roles={role:{(cl,ag):set(g.primaryid) for (cl,ag),g in hits[mask].groupby(['drug_class','agent'])} for role,mask in subsets.items()}
 def parent(role,cl):return set().union(*(v for (cc,aa),v in roles[role].items() if cc==cl))
 def row(cl,ag,role='PS_SS',setting='medication_context',eligible=None,**extra):
  p=parent(role,cl)&(ctx if setting=='medication_context' else ids)
  if eligible is not None:p&=eligible
  t=roles[role].get((cl,ag),set());return dict(database='FAERS',drug_class=cl,agent=ag,role_scope=role,setting=setting,**extra,**metric(len(p&core&t),len((p&core)-t),len((p-core)&t),len((p-core)-t)))
 occ=f.occp_cod.fillna('').str.strip().str.upper();f['reporter_group']=occ.map(lambda x:'HCP_MD_PH_OT' if x in ['MD','PH','OT'] else ('HP_separate' if x=='HP' else ('nonHCP_CN_LW' if x in ['CN','LW'] else ('missing' if not x else 'other_code'))))
 reporters={'all_reporters':ids,'HCP_MD_PH_OT':set(f.loc[f.reporter_group.eq('HCP_MD_PH_OT'),'primaryid']),'HCP_with_HP':set(f.loc[f.reporter_group.isin(['HCP_MD_PH_OT','HP_separate']),'primaryid'])}
 for v in sorted(f.reporter_group.unique()):reporters['partition_'+v]=set(f.loc[f.reporter_group.eq(v),'primaryid'])
 rr=[];profiles=[]
 for role in roles:
  for setting in ['baseline','medication_context']:
   for cl,ag in PAIRS:
    p=row(cl,ag,role,setting);parts=[]
    for label,v in reporters.items():
     r=row(cl,ag,role,setting,v,reporter=label);rr.append(r)
     if label.startswith('partition_'):parts.append(r)
    assert all(sum(r[k] for r in parts)==p[k] for k in 'abcd')
 write('Table_S14_exact_roles_and_reporter.csv',rr)
 for setting in ['baseline','medication_context']:
  for cl in ['NMBA','antibiotic_anchor']:
   p=parent('PS_SS',cl)&(ctx if setting=='medication_context' else ids)
   for outcome,eligible in [('core',p&core),('noncore',p-core)]:
    g=f[f.primaryid.isin(eligible)]
    for (code,group),v in g.groupby(['occp_cod','reporter_group'],dropna=False):profiles.append(dict(setting=setting,drug_class=cl,outcome=outcome,occp_cod=code if code else 'MISSING',reporter_group=group,n=len(v),parent_n=len(g)))
 write('Table_S14b_reporter_composition.csv',profiles)
 deletion=[]
 for setting in ['baseline','medication_context']:
  for cl,ag in roles['PS_SS']:
   if cl=='NMBA' and ag!='succinylcholine':deletion.append(row('NMBA','succinylcholine',setting=setting,eligible=ids-roles['PS_SS'][cl,ag],exclusion='exclude_same_class_agent',excluded_agent=ag))
  for cl,ag in PAIRS:
   other=set().union(*(v for (cc,aa),v in roles['PS_SS'].items() if cc in [('antibiotic_anchor' if cl=='NMBA' else 'NMBA'),'chlorhexidine','dye_anchor']))
   deletion.append(row(cl,ag,setting=setting,eligible=ids-other,exclusion='exclude_other_culprit_classes',excluded_agent='other_three_classes'))
 write('Table_S15_comparator_and_cosuspect_exclusions.csv',deletion)
 country=f.occr_country.fillna('').str.strip();f['country_stratum']=country.map(lambda x:'GB' if x=='GB' else ('observed_non_GB' if x else 'missing'))
 countries=[];fixed=[]
 for setting in ['baseline','medication_context']:
  for cl,ag in PAIRS:
   pp=row(cl,ag,setting=setting);parts=[row(cl,ag,setting=setting,eligible=set(f.loc[f.country_stratum.eq(v),'primaryid']),country_stratum=v) for v in ['GB','observed_non_GB','missing']]
   assert all(sum(r[k] for r in parts)==pp[k] for k in 'abcd');N=sum(r['class_core_n']+r['class_noncore_n'] for r in parts);p1=p0=0.;ok=True
   for r in parts:
    w=(r['class_core_n']+r['class_noncore_n'])/N;r['common_weight']=w
    if r['core_share'] is None or r['noncore_share'] is None:ok=False
    else:p1+=w*r['core_share'];p0+=w*r['noncore_share']
   countries+=parts
   fixed.append(dict(agent=ag,setting=setting,crude_ror=pp['ror'],fixed_mix_core_share=p1 if ok else None,fixed_mix_noncore_share=p0 if ok else None,fixed_mix_ror=p1/(1-p1)/(p0/(1-p0)) if ok and 0<p0<1 and 0<p1<1 else None,uncertainty='descriptive; no CI; no causal adjustment'))
 write('Table_S16_country_strata.csv',countries);write('Table_S16b_common_country_mix.csv',fixed)
 write('Table_S17_context_marker_sensitivity.csv',[row(cl,ag,eligible=e,context_sensitivity=label) for cl,ag in PAIRS for label,e in [('primary_five_marker_context',ctx),('four_marker_context_no_rescue_only',strict),('rescue_only_context',ctx-strict)]])
 # Exclusive pairwise NMBA comparison: co-listed A+B reports are counted and excluded.
 pairwise=[];nmbas={ag:v for (cl,ag),v in roles['PS_SS'].items() if cl=='NMBA'}
 for a,b in combinations(sorted(nmbas),2):
  both=nmbas[a]&nmbas[b]&ctx;aa=(nmbas[a]-nmbas[b])&ctx;bb=(nmbas[b]-nmbas[a])&ctx
  pairwise.append(dict(agent=a,comparator=b,excluded_dual_listed_core_n=len(both&core),excluded_dual_listed_noncore_n=len(both-core),**metric(len(aa&core),len(bb&core),len(aa-core),len(bb-core))))
 write('Table_S29_exclusive_pairwise_NMBA.csv',pairwise)
 # Medication burden is reconstructed from raw rows for the same retained source periods.
 bc=sqlite3.connect(':memory:',uri=True);bc.execute('attach database ? as raw',(RAW.as_uri()+'?mode=ro',));f[['primaryid','source_year','source_quarter']].to_sql('candidate',bc,index=False)
 burden=pd.read_sql_query("select k.primaryid,count(distinct d.drug_seq) as drug_sequence_count,count(distinct case when upper(d.role_cod) in ('PS','SS') then d.drug_seq end) as suspect_drug_sequence_count from candidate k cross join raw.drug d indexed by idx_drug_primaryid on d.primaryid=k.primaryid and d.source_year=k.source_year and d.source_quarter=k.source_quarter group by k.primaryid",bc);bc.close();assert set(burden.primaryid)==ids
 f=f.merge(burden,on='primaryid',validate='one_to_one');f['burden_stratum']=pd.cut(f.drug_sequence_count,[0,5,10,20,float('inf')],labels=['1-5','6-10','11-20','21+']).astype(str)
 bsum=[];brows=[]
 for cl in ['NMBA','antibiotic_anchor']:
  for inside in [False,True]:
   p=parent('PS_SS',cl)&(ctx if inside else ids-ctx)
   for outcome,e in [('core',p&core),('noncore',p-core)]:
    g=f[f.primaryid.isin(e)].drug_sequence_count
    bsum.append(dict(drug_class=cl,in_context=inside,outcome=outcome,n=len(g),median=g.median(),q1=g.quantile(.25),q3=g.quantile(.75)))
 for cl,ag in PAIRS:
  parts=[row(cl,ag,eligible=set(f.loc[f.burden_stratum.eq(v),'primaryid']),burden_stratum=v) for v in ['1-5','6-10','11-20','21+']];pp=row(cl,ag);assert all(sum(r[k] for r in parts)==pp[k] for k in 'abcd');brows+=parts
 write('Table_S8a_medication_burden.csv',bsum);write('Table_S8b_burden_stratified_contrasts.csv',brows)
 f.to_pickle(ROOT/'work/normalized_support_features_private.pkl')
 flow=dict(FAERS_retained=19961236,FAERS_core=len(core),FAERS_context=len(ctx),FAERS_context_core=len(ctx&core))
 (ROOT/'metadata/flow_counts.json').write_text(json.dumps(flow,indent=2),encoding='utf-8')
 # Current Canada/JADER source branch is built separately from full raw files.
 faers=pd.read_csv(ROOT/'aggregate_data/FAERS_primary.csv');cj=pd.read_csv(ROOT/'aggregate_data/Canada_JADER_primary.csv');write('Table_S4_S5_canonical_primary.csv',pd.concat([faers,cj],ignore_index=True))
 # Count class co-suspicion using normalized FAERS flags.
 co=[]
 for a,b in combinations(['NMBA','antibiotic_anchor','chlorhexidine','dye_anchor'],2):
  aa=parent('PS_SS',a)&core;bb=parent('PS_SS',b)&core;co.append(dict(database='FAERS',class_A=a,class_B=b,class_A_core_n=len(aa),class_B_core_n=len(bb),both_core_n=len(aa&bb)))
 write('Table_S7_FAERS_cosuspected_classes.csv',co)
 print(json.dumps(flow));print('Supporting diagnostics complete')
if __name__=='__main__':main()
