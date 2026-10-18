"""Exploratory JADER source-drug burden and PT-relaxed FAERS merging."""
import json, math
import pandas as pd
from scipy.stats import chi2
from settings import ROOT, path
from run_targeted_diagnostics import metric, PAIRS, write

LABELS=['1–5','6–10','11–20','21+']

def heterogeneity(parts):
    positive=[r for r in parts if min(r[k] for k in 'abcd')>0]
    weights=[1/sum(1/r[k] for k in 'abcd') for r in positive]
    logs=[math.log(r['a']*r['d']/(r['b']*r['c'])) for r in positive]
    k=len(positive)
    mean=sum(w*x for w,x in zip(weights,logs))/sum(weights) if k else float('nan')
    q=sum(w*(x-mean)**2 for w,x in zip(weights,logs)) if k else float('nan')
    n=lambda r:sum(r[s] for s in 'abcd')
    den=sum(r['b']*r['c']/n(r) for r in parts if n(r))
    return dict(MH_descriptive_ror=sum(r['a']*r['d']/n(r) for r in parts if n(r))/den if den else float('nan'),Woolf_Q=q,df=k-1,p=chi2.sf(q,k-1) if k>1 else float('nan'),nonzero_strata_n=k,caution='exploratory; native source drug-sequence count; sparse support limits asymptotics')

def jader_burden():
    source=pd.read_pickle(ROOT/'work/JADER_current_private.pkl')
    f,h=source['f'],source['h'];latest=dict(zip(f.primaryid,f['round']))
    hs=h[h.suspect];core=set(f.loc[f.core,'primaryid']);ctx=set(f.loc[f.in_context,'primaryid'])
    pairs=[]
    for g in pd.read_csv(path('jader_raw')/'drug202609.csv',encoding='cp932',dtype=str,keep_default_na=False,chunksize=150000,usecols=['識別番号','報告回数','医薬品連番']):
        keep=g['識別番号'].map(latest).eq(g['報告回数']);g=g[keep]
        assert g['医薬品連番'].ne('').all()
        pairs.append(g[['識別番号','医薬品連番']])
    seq=pd.concat(pairs,ignore_index=True).drop_duplicates()
    count=seq.groupby('識別番号')['医薬品連番'].nunique()
    f=f.copy();f['drug_count']=f.primaryid.map(count);assert f.drug_count.notna().all()
    f['burden']=pd.cut(f.drug_count,[0,5,10,20,float('inf')],labels=LABELS)
    canonical=pd.read_csv(ROOT/'aggregate_data/retained_Canada_JADER_primary.csv')
    rows=[];stats=[]
    for cl,ag in PAIRS:
        parent=set(hs.loc[hs.drug_class.eq(cl),'primaryid'])&ctx
        target=set(hs.loc[hs.drug_class.eq(cl)&hs.agent.eq(ag),'primaryid'])
        parts=[]
        for label in LABELS:
            p=parent&set(f.loc[f.burden.eq(label),'primaryid'])
            row=dict(database='JADER',drug_class=cl,agent=ag,burden_stratum=label,**metric(len(p&core&target),len((p&core)-target),len((p-core)&target),len((p-core)-target)))
            parts.append(row);rows.append(row)
        v=canonical[canonical.database.eq('JADER')&canonical.setting.eq('medication_context')&canonical.agent.eq(ag)].iloc[0]
        assert all(sum(r[k] for r in parts)==int(v[k]) for k in 'abcd')
        summary=heterogeneity(parts)
        if bool(v['sparse']):summary['MH_descriptive_ror']=float('nan')
        stats.append(dict(database='JADER',agent=ag,crude_ror=v.ror,**summary))
    write('JADER_burden_strata.csv',rows);write('JADER_burden_MH_Woolf.csv',stats)
    print(pd.DataFrame(rows)[['agent','burden_stratum','a','b','c','d','ror','ror_ci_low','ror_ci_high']].to_string(index=False))

def relaxed_fingerprint():
    f=pd.read_pickle(ROOT/'work/fingerprint_features_private.pkl')
    h=pd.read_pickle(ROOT/'work/normalized_role_hits_private.pkl');hs=h[h.has_ps.eq(1)|h.has_ss.eq(1)]
    f=f.sort_values(['fda_dt','pid_order'],kind='stable')
    keys=['age_norm','age_unit','sex','event_full_date','country','core']
    rows=[];coverage=[]
    for cl in ['NMBA','antibiotic_anchor']:
        parent=set(hs.loc[hs.drug_class.eq(cl),'primaryid'])
        g=f[f.primaryid.isin(parent)&f.in_context].copy()
        complete=g[g.complete_fingerprint]
        kept=set(g.loc[~g.complete_fingerprint,'primaryid'])|set(complete.drop_duplicates(keys,keep='last').primaryid)
        core=set(g.loc[g.core.eq(1),'primaryid'])
        for outcome,flag in [('core',1),('noncore',0)]:
            arm=g[g.core.eq(flag)];valid=arm[arm.complete_fingerprint]
            coverage.append(dict(drug_class=cl,outcome=outcome,reports=len(arm),complete_fingerprint_n=len(valid),incomplete_n=len(arm)-len(valid),retained_n=len(set(arm.primaryid)&kept),removed_n=len(arm)-len(set(arm.primaryid)&kept)))
        for ccl,ag in PAIRS:
            if ccl!=cl:continue
            target=set(hs.loc[hs.drug_class.eq(cl)&hs.agent.eq(ag),'primaryid'])
            rows.append(dict(database='FAERS',drug_class=cl,agent=ag,fingerprint_rule='PT_relaxed_complete_demographics_date_country_outcome',**metric(len(kept&core&target),len((kept&core)-target),len((kept-core)&target),len((kept-core)-target))))
    write('fingerprint_PT_relaxed_sensitivity.csv',rows);write('fingerprint_PT_relaxed_coverage.csv',coverage)
    print(pd.DataFrame(rows)[['agent','a','b','c','d','ror','ror_ci_low','ror_ci_high']].to_string(index=False))

if __name__=='__main__':
    jader_burden();relaxed_fingerprint()
    (ROOT/'metadata/additional_composition_checks.json').write_text(json.dumps(dict(JADER_burden_all_four_agents_sum_to_primary=True,JADER_unit='distinct native source drug sequence at retained report round; not standardized active ingredient count',FAERS_PT_relaxed='complete age/unit, sex, full date, country, same core status; PT identities omitted; missing fields remain singletons',interpretation='scenario sensitivities, not upper/lower bounds on duplicates or patients'),indent=2),encoding='utf-8')
