"""Build case-insensitive lexical flags from the retained source branch.

New private database only. Input source databases are read-only. Three-character indexed ranges avoid repeated full raw-table scans.
"""
from pathlib import Path
from itertools import product
from collections import defaultdict
from datetime import datetime
import csv,json,sqlite3

from settings import ROOT as R,RAW,LATEST


def prefixes(s):
 head=s[:3];return sorted({''.join(t) for t in product(*[(x.lower(),x.upper()) if x.isalpha() and x.isascii() else (x,) for x in head])})
def log(s):print(datetime.now().strftime('%H:%M:%S'),s,flush=True)
def main():
 c=sqlite3.connect((R/'work/analysis_flags_private.sqlite').as_uri(),uri=True)
 c.execute('pragma journal_mode=wal');c.execute('pragma synchronous=normal')
 c.execute('attach database ? as raw',(RAW.as_uri()+'?mode=ro',));c.execute('attach database ? as latest',(LATEST.as_uri()+'?mode=ro',))
 c.executescript('''create table if not exists class_hits(primaryid text,drug_class text,agent text,has_ps integer,has_ss integer,has_c integer,has_other integer,primary key(primaryid,drug_class,agent));
 create table if not exists event_hits(primaryid text,event_group text,primary key(primaryid,event_group));
 create table if not exists completed(kind text,key text,primary key(kind,key));''')
 rows=list(csv.DictReader((R/'config/drug_class_aliases.csv').open(encoding='utf-8-sig')))
 rows=[r for r in rows if r['database_scope'] in ['FAERS','FAERS_Canada','All']]
 for n,r in enumerate(rows,1):
  alias=r['alias'].strip().upper();key=r['drug_class']+'|'+r['generic_or_group']+'|'+alias
  if c.execute("select 1 from completed where kind='drug' and key=?",(key,)).fetchone():continue
  for field,ix in [('drugname','idx_drug_drugname'),('prod_ai','idx_drug_prod_ai')]:
   for p in prefixes(alias):
    high=p[:-1]+chr(ord(p[-1])+1)
    sql=f'''insert into class_hits select d.primaryid,?,?,max(upper(coalesce(d.role_cod,''))='PS'),max(upper(coalesce(d.role_cod,''))='SS'),max(upper(coalesce(d.role_cod,''))='C'),max(upper(coalesce(d.role_cod,'')) not in ('PS','SS','C'))
      from raw.drug d indexed by {ix} join latest.latest_demo l on l.primaryid=d.primaryid and l.source_year=d.source_year and l.source_quarter=d.source_quarter
      where d.{field}>=? and d.{field}<? and upper(d.{field}) like ? group by d.primaryid
      on conflict(primaryid,drug_class,agent) do update set has_ps=max(has_ps,excluded.has_ps),has_ss=max(has_ss,excluded.has_ss),has_c=max(has_c,excluded.has_c),has_other=max(has_other,excluded.has_other)'''
    c.execute(sql,(r['drug_class'],r['generic_or_group'],p,high,alias+'%'))
  c.execute("insert into completed values('drug',?)",(key,));c.commit();log(f"drug {n}/{len(rows)} {alias}")
 # Sugammadex is a diagnostic co-medication, never an eligibility or culprit class.
 for alias in ['SUGAMMADEX','BRIDION']:
  key='reversal_probe|sugammadex|'+alias
  if c.execute("select 1 from completed where kind='drug' and key=?",(key,)).fetchone():continue
  for field,ix in [('drugname','idx_drug_drugname'),('prod_ai','idx_drug_prod_ai')]:
   for p in prefixes(alias):
    hi=p[:-1]+chr(ord(p[-1])+1)
    c.execute(f'''insert into class_hits select d.primaryid,'reversal_probe','sugammadex',max(upper(coalesce(d.role_cod,''))='PS'),max(upper(coalesce(d.role_cod,''))='SS'),max(upper(coalesce(d.role_cod,''))='C'),max(upper(coalesce(d.role_cod,'')) not in ('PS','SS','C')) from raw.drug d indexed by {ix} join latest.latest_demo l on l.primaryid=d.primaryid and l.source_year=d.source_year and l.source_quarter=d.source_quarter where d.{field}>=? and d.{field}<? and upper(d.{field}) like ? group by d.primaryid on conflict(primaryid,drug_class,agent) do update set has_ps=max(has_ps,excluded.has_ps),has_ss=max(has_ss,excluded.has_ss),has_c=max(has_c,excluded.has_c),has_other=max(has_other,excluded.has_other)''',(p,hi,alias+'%'))
  c.execute("insert into completed values('drug',?)",(key,));c.commit();log(key)
 events=list(csv.DictReader((R/'config/event_pt_groups.csv').open(encoding='utf-8-sig')))
 for r in events:
  term=r['faers_canada_pt'].strip().upper()
  if not term:continue
  key=r['event_group']+'|'+term
  if c.execute("select 1 from completed where kind='event' and key=?",(key,)).fetchone():continue
  for p in prefixes(term):
   hi=p[:-1]+chr(ord(p[-1])+1)
   c.execute('''insert or ignore into event_hits select r.primaryid,? from raw.reac r indexed by idx_reac_pt_name join latest.latest_demo l on l.primaryid=r.primaryid and l.source_year=r.source_year and l.source_quarter=r.source_quarter where r.pt_name>=? and r.pt_name<? and upper(r.pt_name)=?''',(r['event_group'],p,hi,term))
  c.execute("insert into completed values('event',?)",(key,));c.commit();log('event '+term)
 c.execute('create index if not exists idx_hits_class on class_hits(drug_class,agent)');c.commit()
 counts={k:v for k,v in c.execute('select event_group,count(*) from event_hits group by event_group')}
 meta={'source_raw':str(RAW),'retained_latest':str(LATEST),'match':'case-insensitive original prefix aliases; case-insensitive exact original event source terms','event_counts':counts,'class_report_counts':dict(c.execute('select drug_class,count(distinct primaryid) from class_hits group by drug_class')),'source_files_unchanged':True}
 (R/'metadata/analysis_build.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2),encoding='utf-8')
 c.execute('pragma wal_checkpoint(truncate)');c.close();log(json.dumps(meta['event_counts']))
if __name__=='__main__':main()
