"""Post hoc duplicate-candidate and burden diagnostics; aggregates only exported."""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from decimal import Decimal,InvalidOperation
from collections import Counter
import json,re,zipfile,sqlite3
import pandas as pd
import numpy as np
from scipy.stats import chi2
from settings import ROOT,RAW,path
from run_targeted_diagnostics import PAIRS,MARKERS,metric,write,ro,log

def demographics(f):
 p=ROOT/'work/fingerprint_demo_private.pkl'
 if p.exists():return pd.read_pickle(p)
 registry=pd.read_csv(ROOT/'config/faers_DEMO_inventory.csv',dtype=str)
 groups={(str(y),str(q)):set(g.primaryid) for (y,q),g in f.groupby(['source_year','source_quarter'])}
 def job(r):
  y,q=r['year'],r['quarter'];ids=groups.get((y,q),set());out=[]
  with zipfile.ZipFile(path('faers_zip_dir')/r['filename']) as z,z.open(r['demo_member']) as src:
   header=src.readline().decode('utf-8-sig',errors='strict').strip().lower().split('$');ix={s:i for i,s in enumerate(header)}
   aliases={'age':['age'],'age_unit':['age_cod'],'sex':['sex','gndr_cod'],'event_dt':['event_dt'],'sender':['mfr_sndr'],'lit_ref':['lit_ref'],'rpt_cod':['rpt_cod','rept_cod']}
   columns={k:next((ix[s] for s in v if s in ix),None) for k,v in aliases.items()}
   for line in src:
    pid=line.split(b'$',1)[0].decode('ascii').strip()
    if pid not in ids:continue
    vals=line.decode('utf-8',errors='replace').rstrip('\r\n').split('$');out.append(dict(primaryid=pid,**{k:vals[i].strip() if i is not None and i<len(vals) else '' for k,i in columns.items()}))
  assert len(out)==len(ids),(y,q,len(out),len(ids));return pd.DataFrame(out)
 log('Recovering complete fingerprint fields from 88 source quarters')
 with ThreadPoolExecutor(max_workers=4) as pool:d=pd.concat(list(pool.map(job,registry.to_dict('records'))),ignore_index=True)
 assert set(d.primaryid)==set(f.primaryid) and not d.primaryid.duplicated().any();d.to_pickle(p);return d

def main():
 f=pd.read_pickle(ROOT/'work/normalized_report_features_private.pkl');rx=pd.read_pickle(ROOT/'work/normalized_reactions_private.pkl');demo=demographics(f);f=f.merge(demo,on='primaryid',validate='one_to_one');f=f.set_index('primaryid',drop=False)
 h=pd.read_pickle(ROOT/'work/normalized_role_hits_private.pkl');hs=h[h.has_ps.eq(1)|h.has_ss.eq(1)]
 c=ro(ROOT/'work/analysis_flags_private.sqlite');allhits=pd.read_sql_query('select * from class_hits',c);c.close()
 ctx=set(allhits.loc[allhits.drug_class.isin(MARKERS),'primaryid']);core=set(f.loc[f.core.eq(1),'primaryid'])
 ptsets=rx.groupby('primaryid').pt_name.agg(lambda x:tuple(sorted(set(s.strip().upper() for s in x if s.strip()))));f['PT_set']=f.primaryid.map(ptsets)
 # Full class/agent role signature, independent of lexical/source row duplication.
 signatures={rid:tuple(sorted((r.drug_class,r.agent,int(r.has_ps),int(r.has_ss)) for r in g.itertuples())) for rid,g in hs.groupby('primaryid')}
 f['role_signature']=f.primaryid.map(signatures)
 def age(s):
  try:
   v=Decimal(s);return str(v.normalize()) if v.is_finite() and v>=0 else ''
  except InvalidOperation:return ''
 f['age_norm']=f.age.map(age);f['age_unit']=f.age_unit.str.upper();f['sex']=f.sex.str.upper();f['country']=f.occr_country.fillna('').str.strip().str.upper()
 def date(s):
  if not re.fullmatch(r'\d{8}',s):return ''
  try:datetime.strptime(s,'%Y%m%d');return s
  except ValueError:return ''
 f['event_full_date']=f.event_dt.map(date)
 f['complete_fingerprint']=f.age_norm.ne('')&f.age_unit.ne('')&f.sex.isin(['M','F'])&f.event_full_date.ne('')&f.country.ne('')&f.PT_set.notna()
 keys=['age_norm','age_unit','sex','event_full_date','country','PT_set'];f['strict_fp']=list(zip(*(f[k] for k in keys)))
 keysb=['age_norm','age_unit','sex','event_dt','country','PT_set'];f['loose_fp']=list(zip(*(f[k] for k in keysb)))
 # Latest source receipt wins; tie-break numerically by PRIMARYID, not agent status.
 f['pid_order']=pd.to_numeric(f.primaryid,errors='raise');f=f.sort_values(['fda_dt','pid_order'],kind='stable')
 diagnostics=[];coverage=[];clusters=[];profiles=[];burdenrows=[];roleprofiles=[]
 for cl in ['NMBA','antibiotic_anchor']:
  classids=set(hs.loc[hs.drug_class.eq(cl),'primaryid']);g=f[f.primaryid.isin(classids)&f.primaryid.isin(ctx)].copy()
  for outcome,iscore in [('core',True),('noncore',False)]:
   arm=g[g.core.eq(int(iscore))];valid=arm[arm.complete_fingerprint]
   coverage.append(dict(drug_class=cl,outcome=outcome,reports=len(arm),complete_fingerprint_n=len(valid),incomplete_n=len(arm)-len(valid),strict_unique_fingerprint_n=valid.strict_fp.nunique(),strict_removed_n=len(valid)-valid.strict_fp.nunique(),loose_unique_pattern_n=arm.loose_fp.nunique(),loose_collision_excess_n=len(arm)-arm.loose_fp.nunique()))
  for fp,grp in g[g.complete_fingerprint].groupby('strict_fp',sort=False):
   if len(grp)<2:continue
   clusters.append(dict(drug_class=cl,outcome='core' if grp.core.iloc[0] else 'noncore',cluster_reports=len(grp),distinct_caseids=grp.caseid.nunique(),distinct_senders=grp.sender[grp.sender.ne('')].nunique(),year=grp.event_full_date.iloc[0][:4],country=grp.country.iloc[0],role_signatures=grp.role_signature.nunique(),PT_count=len(fp[-1])))
  # Main strict sensitivity and more conservative addition of identical class-role signatures.
  for rule,dupcols in [('strict_demographic_date_PT', ['strict_fp']),('strict_plus_role_signature',['strict_fp','role_signature'])]:
   complete=g[g.complete_fingerprint];kept=set(g.loc[~g.complete_fingerprint,'primaryid'])|set(complete.drop_duplicates(dupcols,keep='last').primaryid)
   for ccl,ag in PAIRS:
    if ccl!=cl:continue
    t=set(hs.loc[hs.drug_class.eq(cl)&hs.agent.eq(ag),'primaryid']);p=set(g.primaryid)
    for label,e in [('unmerged',p),(rule,kept)]:
     diagnostics.append(dict(database='FAERS',drug_class=cl,agent=ag,fingerprint_rule=label,removed_reports=len(p)-len(e),**metric(len(e&core&t),len((e&core)-t),len((e-core)&t),len((e-core)-t))))
  # Agent-specific completeness; no public patient-level fingerprints or sender-linked profiles.
  for ccl,ag in PAIRS:
   if ccl!=cl:continue
   t=set(hs.loc[hs.drug_class.eq(cl)&hs.agent.eq(ag),'primaryid']);a=g[g.primaryid.isin(t)&g.core.eq(1)];valid=a[a.complete_fingerprint]
   profiles.append(dict(agent=ag,context_core_n=len(a),complete_fingerprint_n=len(valid),incomplete_n=len(a)-len(valid),strict_unique_n=valid.strict_fp.nunique(),loose_unique_patterns_n=a.loose_fp.nunique(),distinct_senders_n=a.sender[a.sender.ne('')].nunique(),distinct_caseids_n=a.caseid.nunique(),full_event_date_n=a.event_full_date.ne('').sum(),age_sex_country_complete_n=(a.age_norm.ne('')&a.sex.isin(['M','F'])&a.country.ne('')).sum()))
 # Current burden table is reconstructed directly from source DRUG_SEQ to avoid stale aliases.
 bc=sqlite3.connect(':memory:',uri=True);bc.execute('attach database ? as raw',(RAW.as_uri()+'?mode=ro',));f[['primaryid','source_year','source_quarter']].to_sql('candidate',bc,index=False)
 burden=pd.read_sql_query('select k.primaryid,count(distinct d.drug_seq) as drug_count from candidate k cross join raw.drug d indexed by idx_drug_primaryid on d.primaryid=k.primaryid and d.source_year=k.source_year and d.source_quarter=k.source_quarter group by k.primaryid',bc);bc.close()
 f['drug_count']=f.primaryid.map(burden.set_index('primaryid').drug_count);assert f.drug_count.notna().all();f['burden']=pd.cut(f.drug_count,[0,5,10,20,float('inf')],labels=['1–5','6–10','11–20','21+'])
 hetero=[]
 for cl,ag in PAIRS:
  p=set(hs.loc[hs.drug_class.eq(cl),'primaryid'])&ctx;t=set(hs.loc[hs.drug_class.eq(cl)&hs.agent.eq(ag),'primaryid']);total=metric(len(p&core&t),len((p&core)-t),len((p-core)&t),len((p-core)-t));parts=[]
  hit=hs[hs.drug_class.eq(cl)&hs.agent.eq(ag)].set_index('primaryid');ps=set(hit.loc[hit.has_ps.eq(1)].index);ss=set(hit.loc[hit.has_ps.eq(0)&hit.has_ss.eq(1)].index)
  for label in ['1–5','6–10','11–20','21+']:
   e=set(f.loc[f.burden.eq(label),'primaryid'])&p;r=dict(database='FAERS',drug_class=cl,agent=ag,burden_stratum=label,**metric(len(e&core&t),len((e&core)-t),len((e-core)&t),len((e-core)-t)));parts.append(r);burdenrows.append(r)
   for outcome,ee in [('core',e&core&t),('noncore',(e-core)&t)]:roleprofiles.append(dict(agent=ag,burden_stratum=label,outcome=outcome,target_n=len(ee),PS_present_n=len(ee&ps),SS_only_n=len(ee&ss),sender_n=f.loc[f.primaryid.isin(ee),'sender'].nunique()))
  assert all(sum(r[k] for r in parts)==total[k] for k in 'abcd')
  mh_num=sum(r['a']*r['d']/sum(r[k] for k in 'abcd') for r in parts);mh_den=sum(r['b']*r['c']/sum(r[k] for k in 'abcd') for r in parts)
  eligible=[r for r in parts if all(r[k]>0 for k in 'abcd')];weight=[1/sum(1/r[k] for k in 'abcd') for r in eligible];logs=[__import__('math').log(r['a']*r['d']/(r['b']*r['c'])) for r in eligible];mean=sum(w*x for w,x in zip(weight,logs))/sum(weight);Q=sum(w*(x-mean)**2 for w,x in zip(weight,logs));hetero.append(dict(agent=ag,crude_ror=total['ror'],MH_descriptive_ror=mh_num/mh_den,Woolf_Q=Q,df=len(eligible)-1,p=float(chi2.sf(Q,len(eligible)-1)),nonzero_strata_n=len(eligible),caution='exploratory heterogeneity; MH not causal adjustment; sparse strata limit asymptotics'))
 write('fingerprint_coverage.csv',coverage);write('fingerprint_cluster_size_summary.csv',pd.DataFrame(clusters).value_counts(dropna=False).reset_index(name='cluster_n') if clusters else [{'cluster_n':0}]);write('fingerprint_primary_sensitivity.csv',pd.DataFrame(diagnostics).drop_duplicates(['agent','fingerprint_rule']));write('fingerprint_agent_core_profile.csv',profiles);write('burden_current_strata.csv',burdenrows);write('burden_MH_Woolf.csv',hetero);write('burden_role_sender_profile.csv',roleprofiles)
 # Lexical audit: unrecognized co-medications are not an unmapped target-exposure rate.
 auditids=set(hs.primaryid)&ctx;cc=sqlite3.connect(':memory:',uri=True);cc.execute('attach database ? as raw',(RAW.as_uri()+'?mode=ro',));f.loc[f.primaryid.isin(auditids),['primaryid','source_year','source_quarter']].to_sql('candidate',cc,index=False)
 names=pd.read_sql_query("select d.drugname,count(*) as source_row_n,count(distinct d.primaryid) as report_n from candidate k cross join raw.drug d indexed by idx_drug_primaryid on d.primaryid=k.primaryid and d.source_year=k.source_year and d.source_quarter=k.source_quarter where trim(coalesce(d.prod_ai,''))='' group by d.drugname order by report_n desc",cc);cc.close()
 cfg=list(__import__('csv').DictReader((ROOT/'config/drug_class_aliases.csv').open(encoding='utf-8-sig')));aliases=[r['alias'].upper() for r in cfg if r['database_scope'] in ['FAERS','FAERS_Canada','All']];names['dictionary_matched']=names.drugname.fillna('').str.upper().map(lambda s:any(s.startswith(a) for a in aliases));write('blank_ingredient_DRUGNAME_top100.csv',names[~names.dictionary_matched].head(100))
 f.to_pickle(ROOT/'work/fingerprint_features_private.pkl');log('Fingerprint, burden and lexical diagnostics complete');print(pd.DataFrame(profiles).to_string(index=False));print(pd.DataFrame(diagnostics).drop_duplicates(['agent','fingerprint_rule'])[['agent','fingerprint_rule','a','c','ror']].to_string(index=False))
if __name__=='__main__':main()
