"""Source-qualified FAERS reporting comparisons and diagnostics."""
from pathlib import Path
from collections import defaultdict, Counter
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import csv, io, json, math, re, sqlite3, zipfile, platform
import scipy
import pandas as pd
from scipy.stats import norm

from settings import ROOT,RAW,LATEST,LICENSE,REFERENCE,path
PAIRS=[('NMBA','rocuronium'),('NMBA','succinylcholine'),('antibiotic_anchor','cefazolin'),('antibiotic_anchor','teicoplanin')]
MARKERS=['sedative_hypnotic','opioid','volatile_anesthetic','local_anesthetic','vasoactive_rescue']
PRIMARY_COMP={'MALIGNANT HYPERTHERMIA','HYPERKALAEMIA','BLOOD POTASSIUM INCREASED','NEUROMUSCULAR BLOCK PROLONGED','PROLONGED NEUROMUSCULAR BLOCKADE'}
BROAD_COMP=PRIMARY_COMP|{'PROLONGED APNOEA','APNOEA','RESPIRATORY PARALYSIS'}
ADJACENT={'HYPOTENSION','BRONCHOSPASM','URTICARIA'}
ANAPHYLACTOID={'ANAPHYLACTOID REACTION','ANAPHYLACTOID SHOCK'}
CURRENT_SOURCE_BUILD=True
TAG='normalized_'


def log(x):print(datetime.now().strftime('%H:%M:%S'),x,flush=True)
def write(name,rows):
 d=rows if isinstance(rows,pd.DataFrame) else pd.DataFrame(rows)
 d.to_csv(ROOT/'aggregate_data'/name,index=False,encoding='utf-8-sig',float_format='%.12g');return d
def ro(p):
 assert p.exists() and p.stat().st_size>1000,p
 return sqlite3.connect(p.as_uri()+'?mode=ro',uri=True)
def metric(a,b,c,d,z=1.959963984540054):
 a,b,c,d=map(int,(a,b,c,d));assert min(a,b,c,d)>=0
 out=dict(a=a,b=b,c=c,d=d,class_core_n=a+b,class_noncore_n=c+d,sparse=min(a,c)<5 or min(a,b,c,d)==0)
 p=a/(a+b) if a+b else None;q=c/(c+d) if c+d else None
 out.update(core_share=p,noncore_share=q,difference_pp=100*(p-q) if p is not None and q is not None else None,ror=None,ror_ci_low=None,ror_ci_high=None)
 if not out['sparse']:
  r=a*d/(b*c);se=math.sqrt(sum(1/n for n in (a,b,c,d)))
  out.update(ror=r,ror_ci_low=math.exp(math.log(r)-z*se),ror_ci_high=math.exp(math.log(r)+z*se))
 if p is not None and q is not None:
  def wilson(x,n):
   den=1+z*z/n;mid=(x/n+z*z/(2*n))/den;rad=z*math.sqrt((x/n)*(1-x/n)/n+z*z/(4*n*n))/den;return mid-rad,mid+rad
  l1,u1=wilson(a,a+b);l0,u0=wilson(c,c+d)
  out.update(difference_ci_low=100*(p-q-math.sqrt((p-l1)**2+(u0-q)**2)),difference_ci_high=100*(p-q+math.sqrt((u1-p)**2+(q-l0)**2)))
 return out

def source_frames():
 c=ro(ROOT/'work/analysis_flags_private.sqlite')
 allhits=pd.read_sql_query("select * from class_hits where drug_class in ('NMBA','antibiotic_anchor')",c)
 h=allhits[allhits.has_ps.eq(1)|allhits.has_ss.eq(1)][['primaryid','drug_class','agent']].copy()
 l=ro(LATEST);l.execute('create temp table candidate_pid(primaryid text primary key)');l.executemany('insert into candidate_pid values(?)',((x,) for x in allhits.primaryid.unique()))
 f=pd.read_sql_query('select d.* from candidate_pid k cross join latest_demo d indexed by idx_latest_demo_primaryid on d.primaryid=k.primaryid',l)
 core=set(r[0] for r in c.execute("select primaryid from event_hits where event_group='anaphylaxis_core'"));broad=set(r[0] for r in c.execute("select primaryid from event_hits where event_group='periop_hypersensitivity_broad'"))
 f['core']=f.primaryid.isin(core).astype(int);f['broad']=f.primaryid.isin(broad).astype(int)
 for marker in MARKERS:
  hit=set(r[0] for r in c.execute('select distinct primaryid from class_hits where drug_class=?',(marker,)));f['any_'+marker]=f.primaryid.isin(hit).astype(int)
 f['in_context']=f[['any_'+m for m in MARKERS]].max(axis=1).eq(1)
 total=l.execute('select count(*) from latest_demo').fetchone()[0]
 f.to_pickle(ROOT/'work/normalized_report_features_private.pkl');h.to_pickle(ROOT/'work/normalized_agent_membership_private.pkl');allhits.to_pickle(ROOT/'work/normalized_role_hits_private.pkl')
 c.close();l.close();return f,h,core,total

def reactions(f):
 cache=ROOT/('work/'+TAG+'reactions_private.pkl')
 if cache.exists():
  cached=pd.read_pickle(cache)
  if 'source_term' in cached:return cached
 c=sqlite3.connect(':memory:',uri=True);c.execute('attach database ? as raw',(RAW.as_uri()+'?mode=ro',))
 f[['primaryid','source_year','source_quarter']].to_sql('candidate',c,index=False)
 c.execute('create unique index candidate_pid on candidate(primaryid)')
 log('Extracting source-qualified reaction rows for retained same-class reports')
 d=pd.read_sql_query('select distinct k.primaryid,r.pt_name from candidate k cross join raw.reac r indexed by idx_reac_primaryid on r.primaryid=k.primaryid and r.source_year=k.source_year and r.source_quarter=k.source_quarter',c)
 d['source_term']=d.pt_name.fillna('').str.strip();d.pt_name=d.source_term.str.upper();d=d.drop_duplicates();d.to_pickle(cache);c.close();log(f'Reaction rows {len(d):,}');return d

def parent_terms():
 def read(p):
  with p.open(encoding='utf-8-sig') as f:return list(csv.reader(f,delimiter='$'))
 pts={r[0]:r[1].strip().upper() for r in read(LICENSE/'pt.asc') if len(r)>1}
 llt={r[1].strip().upper():pts.get(r[2],'') for r in read(REFERENCE/'llt.asc') if len(r)>2}
 return pts,llt

def sugammadex(f):
 c=ro(ROOT/'work/analysis_flags_private.sqlite');d=pd.read_sql_query("select * from class_hits where drug_class='reversal_probe'",c);c.close();d=d[d.primaryid.isin(f.primaryid)]
 write('Table_S21a_sugammadex_match_audit.csv',[dict(match='case-insensitive SUGAMMADEX/BRIDION prefix in DRUGNAME or PROD_AI',any_role_report_n=len(d),PS_SS_report_n=int((d.has_ps.eq(1)|d.has_ss.eq(1)).sum()))])
 return set(d.primaryid),set(d.loc[d.has_ps.eq(1)|d.has_ss.eq(1),'primaryid'])

def demo_fields(f):
 target=ROOT/('work/'+TAG+'demo_fields_private.pkl')
 if target.exists():return pd.read_pickle(target)
 groups={ (str(y),str(q)):set(g.primaryid.astype(str)) for (y,q),g in f.groupby(['source_year','source_quarter'])}
 registry=pd.read_csv(ROOT/'config/faers_DEMO_inventory.csv',dtype=str)
 def job(r):
  y,q=r['year'],r['quarter'];ids=groups.get((y,q),set())
  p=path('faers_zip_dir')/r['filename'];rows=[]
  with zipfile.ZipFile(p) as z, z.open(r['demo_member']) as bf:
   header=bf.readline().decode('utf-8-sig',errors='replace').strip().lower().split('$');ix={v:i for i,v in enumerate(header)};key='primaryid' if 'primaryid' in ix else 'isr';ki=ix[key]
   fields={k:next((ix[v] for v in aliases if v in ix),None) for k,aliases in {'rpt_cod':['rpt_cod','rept_cod'],'mfr_sndr':['mfr_sndr'],'lit_ref':['lit_ref'],'fda_dt':['fda_dt'],'event_dt':['event_dt'],'occp_cod':['occp_cod']}.items()}
   for raw in bf:
    # ISR/PRIMARYID is the first column in every audited DEMO source.
    assert ki==0
    pid=raw.split(b'$',1)[0].decode('ascii',errors='replace').strip()
    if pid not in ids:continue
    vals=raw.decode('utf-8',errors='replace').rstrip('\r\n').split('$')
    rows.append(dict(primaryid=pid,source_year=y,source_quarter=q,**{k:vals[i].strip() if i is not None and i<len(vals) else '' for k,i in fields.items()}))
  d=pd.DataFrame(rows);assert len(d)==len(ids) and not d.primaryid.duplicated().any(),(y,q,len(d),len(ids))
  return d,dict(year=y,quarter=q,candidate_n=len(ids),matched_n=len(d),**{k+'_field_present':v is not None for k,v in fields.items()})
 log('Recovering report-channel and literature fields from 88 original DEMO members')
 with ThreadPoolExecutor(max_workers=4) as ex:result=list(ex.map(job,registry.to_dict('records')))
 d=pd.concat([a for a,b in result],ignore_index=True);assert set(d.primaryid)==set(f.primaryid)
 d.to_pickle(target);write('Table_S25a_DEMO_field_coverage.csv',[b for a,b in result]);return d

def smq_flags(rx):
 # Reuse the existing, audited membership loader, without invoking its hashing or analysis runner.
 import ast
 # The licensed PT/SMQ files and the LLT file are stored in separate original
 # folders. Load only the existing parsing functions; no copying or relocation.
 src=ROOT/'code/support/run_official_smq_sensitivity.py'
 tree=ast.parse(src.read_text(encoding='utf-8'))
 selected=ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ['read_dollar','load_smq']],type_ignores=[])
 ns=dict(defaultdict=defaultdict,Path=Path,csv=csv,SMQ_CODE='20000021',SMQ_NAME='Anaphylactic reaction (SMQ)',EXPECTED_ALGORITHM='A or (B and C) or (D and (B or C))',SMQ_LIST=LICENSE/'smq_list.asc',SMQ_CONTENT=LICENSE/'smq_content.asc',PT_ASC=LICENSE/'pt.asc',LLT_ASC=REFERENCE/'llt.asc')
 exec(compile(selected,str(src),'exec'),ns)
 meta,membership,n=ns['load_smq']()
 pts={r[0]:r[1].upper() for r in ns['read_dollar'](ns['PT_ASC']) if len(r)>1}
 name_parent={r[1].upper():pts.get(r[2],'') for r in ns['read_dollar'](ns['LLT_ASC']) if len(r)>2}
 mem={k.upper():v for k,v in membership.items()}
 records=[]
 for source in sorted(set(rx.pt_name)):
  exact=mem.get(source);parent=name_parent.get(source,source);mapped=mem.get(parent)
  records.append(dict(source_term=source,parent_term=parent,llt_parent_used=parent!=source,exact_narrow=bool(exact and '2' in exact['scopes']),normalized_narrow=bool(mapped and '2' in mapped['scopes']),exact_categories=''.join(sorted(exact['categories'])) if exact else '',normalized_categories=''.join(sorted(mapped['categories'])) if mapped else ''))
 mapping=pd.DataFrame(records);mapping.to_pickle(ROOT/'work/licensed_term_audit_private.pkl')
 rows=[]
 # Exact PT and source-roll-up scopes.
 for matching in ['exact','normalized']:
  cats={cat:set(mapping.loc[mapping[matching+'_categories'].str.contains(cat),'source_term']) for cat in 'ABCD'}
  sets={cat:set(rx.loc[rx.pt_name.isin(terms),'primaryid']) for cat,terms in cats.items()}
  narrow=set(rx.loc[rx.pt_name.isin(set(mapping.loc[mapping[matching+'_narrow'],'source_term'])),'primaryid'])
  algo=sets['A']|(sets['B']&sets['C'])|(sets['D']&(sets['B']|sets['C']))
  rows.append((matching,narrow,algo))
 # A small source-term audit, not a redistributed SMQ membership list.
 public=[]
 for term in ['ANAPHYLAXIS','ANAPHYLACTIC REACTION','ANAPHYLACTIC SHOCK','ANAPHYLACTOID REACTION','ANAPHYLACTOID SHOCK']:
  v=next((r for r in records if r['source_term']==term),None)
  if v is None:
   parent=name_parent.get(term,term);mapped=mem.get(parent);v=dict(llt_parent_used=parent!=term,exact_narrow=bool(mem.get(term) and '2' in mem[term]['scopes']),normalized_narrow=bool(mapped and '2' in mapped['scopes']))
  public.append(dict(source_term=term,source_report_n=int(rx.loc[rx.pt_name.eq(term),'primaryid'].nunique()),denominator_report_n=int(rx.primaryid.nunique()),denominator='all-role NMBA or antibiotic candidate reports with source reactions',is_current_pt=term in set(pts.values()),parent_mapping_used=v['llt_parent_used'],narrow_exact=v['exact_narrow'],narrow_after_source_rollup=v['normalized_narrow']))
 write('Table_S23a_source_terminology_audit.csv',public)
 return rows

def main():
 f,h,core_all,total=source_frames();ids=set(f.primaryid);core=core_all&ids;ctx=set(f.loc[f.in_context,'primaryid'])
 parents={cl:set(g.primaryid) for cl,g in h.groupby('drug_class')};targets={(cl,ag):set(g.primaryid) for (cl,ag),g in h.groupby(['drug_class','agent'])}
 def row(cl,ag,setting='medication_context',eligible=None,outcome=None,**extra):
  p=parents[cl]&(ctx if setting=='medication_context' else ids)
  if eligible is not None:p&=eligible
  e=core if outcome is None else outcome;t=targets[cl,ag];a=len(p&e&t);b=len((p&e)-t);c=len((p-e)&t);d=len((p-e)-t)
  return dict(database='FAERS',drug_class=cl,agent=ag,setting=setting,**extra,**metric(a,b,c,d))
 write('FAERS_primary.csv',[row(cl,ag,setting) for setting in ['baseline','medication_context'] for cl,ag in PAIRS])
 rx=reactions(f);sugany,sugsusp=sugammadex(f)
 co=[];sens=[]
 for cl,ag in PAIRS:
  for setting in ['baseline','medication_context']:
   p=parents[cl]&(ctx if setting=='medication_context' else ids);t=targets[cl,ag]
   for scope,sg in [('any_role',sugany),('PS_SS',sugsusp)]:
    for arm,pp in [('core',p&core),('noncore',p-core)]:
     co.append(dict(database='FAERS',agent=ag,setting=setting,sugammadex_role=scope,outcome=arm,class_n=len(pp),target_n=len(pp&t),class_sugammadex_n=len(pp&sg),target_sugammadex_n=len(pp&t&sg)))
    sens.append(row(cl,ag,setting,ids-sg,diagnostic='exclude_sugammadex_'+scope))
 write('Table_S21b_sugammadex_coreporting.csv',co);write('Table_S21c_sugammadex_exclusion.csv',sens)
 tops=[]
 for cl,ag in PAIRS:
  for setting in ['baseline','medication_context']:
   p=targets[cl,ag]&(ctx if setting=='medication_context' else ids)-core
   counts=rx[rx.primaryid.isin(p)].groupby('pt_name').primaryid.nunique().sort_values(ascending=False)
   for rank,(pt,n) in enumerate(counts.head(10).items(),1):tops.append(dict(database='FAERS',agent=ag,setting=setting,rank=rank,source_term=pt,report_n=int(n),noncore_target_report_n=len(p),fraction=n/len(p)))
 write('Table_S22a_noncore_top10_PT.csv',tops)
 pts,llt=parent_terms()
 events={}
 matchrows=[]
 for label,terms in [('pharmacologic_competitors',PRIMARY_COMP),('pharmacologic_plus_respiratory',BROAD_COMP)]:
  expanded_terms=terms|{llt.get(t,t) for t in terms}
  source_terms={t for t in set(rx.pt_name) if t in expanded_terms or llt.get(t,t) in expanded_terms}
  events[label]=set(rx.loc[rx.pt_name.isin(source_terms),'primaryid'])
  for t in sorted(source_terms):matchrows.append(dict(family=label,source_term=t,report_n=int(rx.loc[rx.pt_name.eq(t),'primaryid'].nunique())))
 write('Table_S22b_event_exclusions.csv',[row(cl,ag,setting,ids-e,diagnostic='exclude_'+label) for setting in ['baseline','medication_context'] for cl,ag in PAIRS if cl=='NMBA' for label,e in events.items()])
 write('Table_S22c_event_family_matches.csv',matchrows)
 expanded=core|set(rx.loc[rx.pt_name.isin(ANAPHYLACTOID),'primaryid']);outcomes=[('core_plus_anaphylactoid',expanded)]
 for matching,narrow,algo in smq_flags(rx):outcomes.extend([('smq_narrow_'+matching,narrow),('smq_algorithm_'+matching,algo)])
 outr=[]
 for label,e in outcomes:
  for cl,ag in PAIRS:outr.append(row(cl,ag,outcome=e,diagnostic=label))
 write('Table_S23b_outcome_sensitivity.csv',outr)
 countries=f.occr_country.fillna('').str.strip();parts={'JP':set(f.loc[countries.eq('JP'),'primaryid']),'CA':set(f.loc[countries.eq('CA'),'primaryid']),'missing':set(f.loc[countries.eq(''),'primaryid']),'other_known':set(f.loc[~countries.isin(['','JP','CA']),'primaryid'])}
 geo=[]
 for setting in ['baseline','medication_context']:
  for cl,ag in PAIRS:
   rr=[row(cl,ag,setting,v,country_partition=k) for k,v in parts.items()];parent=row(cl,ag,setting)
   assert all(sum(r[k] for r in rr)==parent[k] for k in 'abcd');geo+=rr
   geo.append(row(cl,ag,setting,ids-parts['JP']-parts['CA'],country_partition='exclude_JP_CA'))
 write('Table_S24_geography_overlap_sensitivity.csv',geo)
 demo=demo_fields(f);demo.rpt_cod=demo.rpt_cod.str.upper().str.strip();d=f.merge(demo,on=['primaryid','source_year','source_quarter'],validate='one_to_one',suffixes=('','_source'))
 assert d.fda_dt.fillna('').eq(d.fda_dt_source.fillna('')).all()
 channels={k:set(d.loc[d.rpt_cod.eq(k),'primaryid']) for k in ['EXP','PER','DIR']};channels['other_missing']=ids-set().union(*channels.values())
 ch=[]
 for setting in ['baseline','medication_context']:
  for cl,ag in PAIRS:
   rr=[row(cl,ag,setting,v,report_channel=k) for k,v in channels.items()];parent=row(cl,ag,setting);assert all(sum(r[k] for r in rr)==parent[k] for k in 'abcd');ch+=rr
 write('Table_S25b_report_channel_sensitivity.csv',ch)
 # PS-only channel contrasts must reconstruct PS-specific class parents, not restrict only target drugs.
 rh=pd.read_pickle(ROOT/'work/normalized_role_hits_private.pkl');rh=rh[rh.has_ps.eq(1)]
 psparents={cl:set(g.primaryid)&ids for cl,g in rh.groupby('drug_class')};pstarget={(cl,ag):set(g.primaryid)&ids for (cl,ag),g in rh.groupby(['drug_class','agent'])}
 ps=[]
 for cl,ag in PAIRS:
  for label,eligible in [('all',ids),*channels.items()]:
   p=psparents[cl]&ctx&eligible;t=pstarget.get((cl,ag),set());ps.append(dict(database='FAERS',agent=ag,role_scope='PS_only',report_channel=label,**metric(len(p&core&t),len((p&core)-t),len((p-core)&t),len((p-core)-t))))
 write('Table_S25c_PS_only_report_channel.csv',ps)
 # Teicoplanin report provenance: counts, never IDs or full literature strings.
 teico=targets['antibiotic_anchor','teicoplanin']&ctx&core;td=d[d.primaryid.isin(teico)].copy();td['receipt_year']=td.fda_dt.str[:4];td['literature_reference_present']=td.lit_ref.ne('');td['country']=td.occr_country.fillna('').replace('','MISSING');td.to_pickle(ROOT/'work/teicoplanin_provenance_private.pkl')
 profile=[]
 for field in ['receipt_year','country','rpt_cod','mfr_sndr','literature_reference_present']:
  for val,n in td[field].value_counts(dropna=False).items():profile.append(dict(dimension=field,value=str(val),report_n=int(n),parent_n=len(td)))
 write('Table_S25d_teicoplanin_provenance.csv',profile)
 # Annual numerator AND class denominator to avoid counts-only notoriety claims.
 annual=[]
 for cl,ag in PAIRS:
  for year,gg in d.groupby(d.fda_dt.str[:4]):annual.append(row(cl,ag,eligible=set(gg.primaryid),receipt_year=year))
 write('Table_S25e_annual_context_diagnostics.csv',annual)
 # Conventional whole-database background, same role and endpoint; separate estimand panel.
 whole=[]
 for cl,ag in PAIRS:
  t=targets[cl,ag];a=len(t&core_all);c=len(t-core_all);whole.append(dict(database='FAERS',agent=ag,background='whole_database',**metric(a,len(core_all)-a,c,total-len(core_all)-c)))
 write('Table_S26_whole_database_background.csv',whole)
 bonf=[]
 for cl,ag in PAIRS:
  r=row(cl,ag);bonf.append(dict(agent=ag,family_n=4,interval='Bonferroni 98.75% individual interval',**metric(*(r[k] for k in 'abcd'),z=float(norm.ppf(1-.05/8)))))
 write('Table_S27_primary_multiplicity_sensitivity.csv',bonf)
 manifest={'analysis_lock':'2026-10-18','FAERS_class_report_n':len(f),'FAERS_reaction_rows':len(rx),'software':{'python':platform.python_version(),'pandas':pd.__version__,'scipy':scipy.__version__,'sqlite':sqlite3.sqlite_version},'private_exports':'report identifiers and licensed membership stay in work'}
 (ROOT/'metadata/analysis_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
 log('COMPLETE');print(pd.DataFrame(sens).query("setting=='medication_context'")[['agent','diagnostic','a','c','ror']].to_string(index=False));print(pd.DataFrame(outr)[['agent','diagnostic','a','ror']].to_string(index=False))

if __name__=='__main__':main()
