"""Recalculate all exported four-cell estimates without production helpers."""
from pathlib import Path
from statistics import NormalDist
import csv,json,math,sys
folder=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parents[1]
checked=0
for p in sorted(folder.glob('Table_*.csv')):
 rows=list(csv.DictReader(p.open(encoding='utf-8-sig')))
 if not rows or not all(k in rows[0] for k in 'abcd'):continue
 for r in rows:
  a,b,c,d=[int(r[k]) for k in 'abcd'];assert min(a,b,c,d)>=0
  sparse=min(a,c)<5 or min(a,b,c,d)==0
  assert r.get('sparse',str(sparse))==str(sparse),(p.name,r)
  z=NormalDist().inv_cdf(1-.05/(2*int(r.get('family_n') or 1)))
  def close(k,x):
   if r.get(k,'')!='':assert math.isclose(float(r[k]),x,rel_tol=2e-10,abs_tol=2e-9),(p.name,k,r[k],x)
  if r.get('ror','')!='':
   assert not sparse;rr=a*d/(b*c);se=math.sqrt(1/a+1/b+1/c+1/d);close('ror',rr);close('ror_ci_low',rr*math.exp(-z*se));close('ror_ci_high',rr*math.exp(z*se))
  elif 'ror' in r:assert sparse
  if a+b and c+d:
   p1=a/(a+b);p0=c/(c+d);delta=p1-p0
   def wilson(x,n):
    mid=(x+z*z/2)/(n+z*z);rad=z*math.sqrt(x*(1-x/n)+z*z/4)/(n+z*z);return mid-rad,mid+rad
   l1,u1=wilson(a,a+b);l0,u0=wilson(c,c+d);close('difference_pp',delta*100);close('difference_ci_low',(delta-math.hypot(p1-l1,u0-p0))*100);close('difference_ci_high',(delta+math.hypot(u1-p1,p0-l0))*100)
  checked+=1
print(json.dumps(dict(four_cell_rows_recomputed=checked,status='passed'),indent=2))
