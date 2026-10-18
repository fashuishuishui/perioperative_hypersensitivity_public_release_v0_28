"""Rebuild current complete source extracts; no use of the old report cache."""
from pathlib import Path
from functools import lru_cache
from collections import defaultdict,Counter
import csv,json,re,unicodedata
import pandas as pd
from settings import ROOT,INPUTS,path
from run_targeted_diagnostics import metric,PAIRS,MARKERS,write,log

def norm(x):return unicodedata.normalize('NFKC',str(x or '')).casefold().strip()
def chunks(p,**kw):return pd.read_csv(p,dtype=str,keep_default_na=False,chunksize=150000,**kw)
def main_db(db):
 target=ROOT/'work'/f'{db}_current_private.pkl'
 if target.exists():return pd.read_pickle(target)
 aliases=list(csv.DictReader((ROOT/'config/drug_class_aliases.csv').open(encoding='utf-8-sig')))
 cfg=[r for r in aliases if r['database_scope'] in [db,'All'] or (db=='Canada' and r['database_scope']=='FAERS_Canada')]
 for al in (['sugammadex','bridion'] if db=='Canada' else ['スガマデクス','ブリディオン']):cfg.append(dict(alias=al,drug_class='reversal_probe',generic_or_group='sugammadex'))
 lookup=defaultdict(set)
 for r in cfg:lookup[norm(r['alias'])].add((r['drug_class'],r['generic_or_group']))
 regex=re.compile(('(?<![a-z0-9])(?:' if db=='Canada' else '(?:')+'|'.join(re.escape(s) for s in sorted(lookup,key=len,reverse=True))+(')(?![a-z0-9])' if db=='Canada' else ')'))
 @lru_cache(maxsize=300000)
 def match(text):return set().union(*(lookup[m.group()] for m in regex.finditer(norm(text))))
 folder=path('canada_raw' if db=='Canada' else 'jader_raw')
 if db=='Canada':
  parts=[]
  for g in chunks(folder/'reports.txt',sep='$',header=None,usecols=[0,2,3,4,10,12,14,26,34]):
   g.columns=['primaryid','round','first_received','last_received','sex','age','age_unit','serious','reporter'];parts.append(g)
  reports=pd.concat(parts,ignore_index=True);seen=len(reports);reports['round_num']=pd.to_numeric(reports['round'],errors='coerce').fillna(0)
  reports=reports.sort_values('round_num',kind='stable').drop_duplicates('primaryid',keep='last');latest=None
 else:
  reports=pd.read_csv(folder/'demo202609.csv',encoding='cp932',dtype=str,keep_default_na=False);seen=len(reports);reports['round_num']=pd.to_numeric(reports['報告回数'],errors='raise');reports=reports.sort_values('round_num',kind='stable').drop_duplicates('識別番号',keep='last').rename(columns={'識別番号':'primaryid','報告回数':'round','報告年度・四半期':'report_quarter','性別':'sex','年齢':'age'})
  latest=dict(zip(reports.primaryid,reports['round']))
 ids=set(reports.primaryid);log(f'{db}: {seen:,} source metadata rows -> {len(ids):,} retained')
 hit={};drugrows=unlinked=nonlatest=0
 if db=='Canada':iterator=chunks(folder/'report_drug.txt',sep='$',header=None,usecols=[1,3,4])
 else:iterator=chunks(folder/'drug202609.csv',encoding='cp932',usecols=['識別番号','報告回数','医薬品（一般名）','医薬品（販売名）','医薬品の関与'])
 for g in iterator:
  if db=='Canada':iterator2=((r[0],None,r[1],r[2]) for r in g.itertuples(index=False,name=None))
  else:iterator2=((r[0],r[1],r[3]+' '+r[4],r[2]) for r in g[['識別番号','報告回数','医薬品の関与','医薬品（一般名）','医薬品（販売名）']].itertuples(index=False,name=None))
  for rid,ver,text,role in iterator2:
   drugrows+=1
   if rid not in ids:unlinked+=1;continue
   if latest is not None and latest[rid]!=ver:nonlatest+=1;continue
   suspect=('suspect' in norm(role)) if db=='Canada' else ('被疑薬' in role)
   for cl,ag in match(text):hit[rid,cl,ag]=bool(hit.get((rid,cl,ag),False) or suspect)
 log(f'{db}: drug matching complete; {len(hit):,} report-agent hits')
 h=pd.DataFrame([(rid,cl,ag,sus) for (rid,cl,ag),sus in hit.items()],columns=['primaryid','drug_class','agent','suspect'])
 candidate=set(h.loc[h.drug_class.isin(['NMBA','antibiotic_anchor']),'primaryid']);ctx=set(h.loc[h.drug_class.isin(MARKERS),'primaryid'])
 eventcfg=list(csv.DictReader((ROOT/'config/event_pt_groups.csv').open(encoding='utf-8-sig')))
 col='faers_canada_pt' if db=='Canada' else 'jader_pt_ja';coreterms={norm(r[col]) for r in eventcfg if r['event_group']=='anaphylaxis_core'}
 core=set();rxrows=[];reactionrows=reactunlinked=reactnonlatest=0;versions=Counter()
 if db=='Canada':iterator=chunks(folder/'reactions.txt',sep='$',header=None,usecols=[1,5,9])
 else:iterator=chunks(folder/'reac202609.csv',encoding='cp932',usecols=['識別番号','報告回数','有害事象','有害事象の発現日'])
 for g in iterator:
  if db=='Canada':iterator2=((r[0],None,r[1],r[2]) for r in g.itertuples(index=False,name=None))
  else:iterator2=((r[0],r[1],r[2],'') for r in g[['識別番号','報告回数','有害事象']].itertuples(index=False,name=None))
  for rid,ver,term,version in iterator2:
   reactionrows+=1
   if rid not in ids:reactunlinked+=1;continue
   if latest is not None and latest[rid]!=ver:reactnonlatest+=1;continue
   versions[version]+=1
   if norm(term) in coreterms:core.add(rid)
   if rid in candidate:rxrows.append((rid,term))
 rx=pd.DataFrame(rxrows,columns=['primaryid','source_term']).drop_duplicates();f=reports[reports.primaryid.isin(candidate)].copy();f['core']=f.primaryid.isin(core);f['in_context']=f.primaryid.isin(ctx)
 meta=dict(database=db,metadata_rows=seen,retained_reports=len(ids),metadata_duplicate_rows=seen-len(ids),core_reports=len(core),context_reports=len(ctx),context_core_reports=len(ctx&core),drug_rows=drugrows,unlinked_drug_rows=unlinked,nonlatest_round_drug_rows=nonlatest,reaction_rows=reactionrows,unlinked_reaction_rows=reactunlinked,nonlatest_round_reaction_rows=reactnonlatest,reaction_version_values=dict(versions),source_label='coverage through 2026-05-31' if db=='Canada' else 'September 2026 release',analysis_lock='2026-10-18')
 if db=='Canada':
  def date(s):
   m=re.fullmatch(r'(\d{2})-([A-Z]{3})-(\d{2}|\d{4})',s)
   if not m:return ''
   d,mo,y=m.groups();y=int(y);y=y+(1900 if y>=65 else 2000) if y<100 else y
   return f'{y:04}-{list(__import__("calendar").month_abbr).index(mo.title()):02}-{d}'
  for field in ['first_received','last_received']:
   values=reports[field].map(date);valid=values[values.ne('')];meta[field+'_range']=[valid.min(),valid.max()];meta[field+'_missing_n']=int(values.eq('').sum())
 else:meta['report_quarter_values']=sorted(reports.report_quarter.unique().tolist())
 obj=dict(f=f,h=h,rx=rx,meta=meta);pd.to_pickle(obj,target);log(f'{db}: source extraction complete');return obj

def run():
 primary=[];sug=[];co=[];tops=[];events=[];flows=[];termaudit=[]
 for db in ['Canada','JADER']:
  o=main_db(db);f,h,rx,meta=[o[k] for k in ['f','h','rx','meta']];ids=set(f.primaryid);core=set(f.loc[f.core,'primaryid']);ctx=set(f.loc[f.in_context,'primaryid']);flows.append(meta)
  hs=h[h.suspect];parents={cl:set(g.primaryid) for cl,g in hs.groupby('drug_class')};targets={(cl,ag):set(g.primaryid) for (cl,ag),g in hs.groupby(['drug_class','agent'])}
  def row(cl,ag,setting='medication_context',eligible=None,outcome=None,**kw):
   p=parents.get(cl,set())&(ctx if setting=='medication_context' else ids)
   if eligible is not None:p&=eligible
   e=core if outcome is None else outcome;t=targets.get((cl,ag),set())
   return dict(database=db,drug_class=cl,agent=ag,setting=setting,**kw,**metric(len(p&e&t),len((p&e)-t),len((p-e)&t),len((p-e)-t)))
  for cl,ag in PAIRS:
   for setting in ['baseline','medication_context']:primary.append(row(cl,ag,setting))
   for label,mask in [('any_role',h.agent.eq('sugammadex')),('suspect',h.agent.eq('sugammadex')&h.suspect)]:
    sg=set(h.loc[mask,'primaryid']);sug.append(row(cl,ag,eligible=ids-sg,diagnostic='exclude_sugammadex_'+label))
    t=targets.get((cl,ag),set())&ctx
    for outcome,p in [('core',t&core),('noncore',t-core)]:co.append(dict(database=db,agent=ag,scope=label,outcome=outcome,target_n=len(p),sugammadex_n=len(p&sg)))
   p=(targets.get((cl,ag),set())&ctx)-core;freq=rx[rx.primaryid.isin(p)].source_term.value_counts()
   for rank,(term,n) in enumerate(freq.head(10).items(),1):tops.append(dict(database=db,agent=ag,rank=rank,source_term=term,report_n=int(n),noncore_target_report_n=len(p)))
  # Actual source terms, kept separate from an official cross-language SMQ.
  pharmac={'malignant hyperthermia','hyperthermia malignant','hyperkalaemia','blood potassium increased','neuromuscular block prolonged','prolonged neuromuscular blockade'} if db=='Canada' else {'悪性高熱','悪性高熱症','高カリウム血症','血中カリウム増加','神経筋ブロック遷延','神経筋ブロックの遷延','筋弛緩作用遷延'}
  matched=rx.source_term.map(norm).isin({norm(t) for t in pharmac});ee=set(rx.loc[matched,'primaryid'])
  for term,g in rx[matched].groupby('source_term'):termaudit.append(dict(database=db,family='pharmacologic_competitors',source_term=term,report_n=g.primaryid.nunique()))
  for cl,ag in PAIRS:
   if cl=='NMBA':events.append(row(cl,ag,eligible=ids-ee,diagnostic='exclude_pharmacologic_competitors'))
 write('Canada_JADER_primary.csv',primary);write('current_Canada_JADER_sugammadex_exclusion.csv',sug);write('current_Canada_JADER_sugammadex_coreporting.csv',co);write('current_Canada_JADER_noncore_top10.csv',tops);write('current_Canada_JADER_competing_event_exclusion.csv',events);write('current_Canada_JADER_event_matches.csv',termaudit)
 (ROOT/'metadata/current_Canada_JADER_source_audit.json').write_text(json.dumps(flows,ensure_ascii=False,indent=2),encoding='utf-8')
 print(pd.DataFrame(primary).query("setting=='medication_context'")[['database','agent','a','c','ror']].to_string(index=False));print(pd.DataFrame(sug).to_string(index=False))
if __name__=='__main__':run()
