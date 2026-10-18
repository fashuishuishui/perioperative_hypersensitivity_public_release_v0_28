"""Regenerate the three current manuscript figures from released aggregates."""
from pathlib import Path
import argparse,os
import pandas as pd
R=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(R/'work/matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
AGENTS=['rocuronium','succinylcholine','cefazolin','teicoplanin']
NAMES={a:a.capitalize() for a in AGENTS}
COLORS={'FAERS':'#287CB3','Canada':'#D87C27','JADER':'#398853'}
MAPPING={'Table_S4_S5_canonical_primary.csv': 'Table_S4a_primary_four_cells.csv', 'Table_S7_FAERS_cosuspected_classes.csv': 'Table_S4b_FAERS_cosuspected_classes.csv', 'Table_S8a_medication_burden.csv': 'Table_S5a_burden_summary.csv', 'burden_current_strata.csv': 'Table_S5b_burden_strata.csv', 'burden_MH_Woolf.csv': 'Table_S5c_burden_MH_Woolf.csv', 'burden_role_sender_profile.csv': 'Table_S5d_burden_roles_and_senders.csv', 'Table_S14_exact_roles_and_reporter.csv': 'Table_S6a_roles_and_reporters.csv', 'Table_S14b_reporter_composition.csv': 'Table_S6b_reporter_composition.csv', 'Table_S15_comparator_and_cosuspect_exclusions.csv': 'Table_S7a_comparator_exclusions.csv', 'Table_S29_exclusive_pairwise_NMBA.csv': 'Table_S7b_exclusive_NMBA_pairs.csv', 'Table_S17_context_marker_sensitivity.csv': 'Table_S8a_context_marker_sensitivity.csv', 'context_dictionary_coverage_sensitivity.csv': 'Table_S8b_context_coverage_probe.csv', 'context_dictionary_probe_hits.csv': 'Table_S8c_probe_alias_hits.csv', 'Table_S21b_sugammadex_coreporting.csv': 'Table_S9a_FAERS_sugammadex_coreports.csv', 'Table_S21c_sugammadex_exclusion.csv': 'Table_S9b_FAERS_sugammadex_exclusion.csv', 'current_Canada_JADER_sugammadex_coreporting.csv': 'Table_S9c_Canada_JADER_sugammadex_coreports.csv', 'current_Canada_JADER_sugammadex_exclusion.csv': 'Table_S9d_Canada_JADER_sugammadex_exclusion.csv', 'Table_S22a_noncore_top10_PT.csv': 'Table_S10a_FAERS_noncore_PT.csv', 'current_Canada_JADER_noncore_top10.csv': 'Table_S10b_Canada_JADER_noncore_PT.csv', 'Table_S22b_event_exclusions.csv': 'Table_S11a_FAERS_pharmacologic_exclusions.csv', 'current_Canada_JADER_competing_event_exclusion.csv': 'Table_S11b_Canada_JADER_pharmacologic_exclusions.csv', 'Table_S22c_event_family_matches.csv': 'Table_S11c_FAERS_family_matches.csv', 'current_Canada_JADER_event_matches.csv': 'Table_S11d_Canada_JADER_family_matches.csv', 'Table_S23b_outcome_sensitivity.csv': 'Table_S12_FAERS_outcome_sensitivity.csv', 'Table_S16_country_strata.csv': 'Table_S13a_country_strata.csv', 'Table_S16b_common_country_mix.csv': 'Table_S13b_fixed_country_mix.csv', 'Table_S24_geography_overlap_sensitivity.csv': 'Table_S13c_geographic_overlap.csv', 'Table_S25a_DEMO_field_coverage.csv': 'Table_S14a_quarterly_field_coverage.csv', 'Table_S25b_report_channel_sensitivity.csv': 'Table_S14b_report_channels.csv', 'Table_S25c_PS_only_report_channel.csv': 'Table_S14c_PS_only_channels.csv', 'Table_S25d_teicoplanin_provenance.csv': 'Table_S14d_teicoplanin_provenance.csv', 'Table_S25e_annual_context_diagnostics.csv': 'Table_S14e_annual_context.csv', 'fingerprint_coverage.csv': 'Table_S15a_fingerprint_coverage.csv', 'fingerprint_cluster_size_summary.csv': 'Table_S15b_cluster_summary.csv', 'fingerprint_primary_sensitivity.csv': 'Table_S15c_fingerprint_contrasts.csv', 'fingerprint_agent_core_profile.csv': 'Table_S15d_agent_fingerprint_coverage.csv', 'fingerprint_burden_sensitivity.csv': 'Table_S15e_fingerprint_burden.csv', 'fingerprint_country_counts.csv': 'Table_S15f_fingerprint_country_counts.csv', 'Table_S26_whole_database_background.csv': 'Table_S16_whole_database.csv', 'Table_S27_primary_multiplicity_sensitivity.csv': 'Table_S17_multiplicity.csv', 'blank_ingredient_DRUGNAME_top100.csv': 'Table_S2c_blank_ingredient_name_audit.csv', 'Table_S23a_source_terminology_audit.csv': 'Table_S3b_terminology_audit.csv', 'burden_PS_only_feasibility.csv': 'Table_S5e_PS_only_burden_feasibility.csv', 'JADER_burden_strata.csv': 'Table_S5f_JADER_burden_strata.csv', 'JADER_burden_MH_Woolf.csv': 'Table_S5g_JADER_burden_MH_Woolf.csv', 'fingerprint_PT_relaxed_sensitivity.csv': 'Table_S15g_PT_relaxed_fingerprint_contrasts.csv', 'fingerprint_PT_relaxed_coverage.csv': 'Table_S15h_PT_relaxed_fingerprint_coverage.csv'}

def read(name):return pd.read_csv(R/MAPPING[name])
data=read
def num(n):return f'{int(n):,}'
def export(fig,name):
    for ext in ['png','pdf','svg']:fig.savefig(OUT/f'{name}.{ext}',dpi=300,bbox_inches='tight',facecolor='white')
    plt.close(fig)

def pick(df,**kw):
 for k,v in kw.items():df=df[df[k].eq(v)]
 assert len(df)==1,(kw,len(df));return df.iloc[0]

def plots(p,flow,cj):
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6,'pdf.fonttype':42,'svg.fonttype':'none'})
 fig,ax=plt.subplots(figsize=(7.1,6.7));ax.axis('off');ax.set(xlim=(0,3),ylim=(0,7))
 for x,db in enumerate(['FAERS','Canada','JADER']):
  m=flow if db=='FAERS' else next(v for v in cj if v['database']==db)
  total=m['FAERS_retained'] if db=='FAERS' else m['retained_reports'];core=m['FAERS_core'] if db=='FAERS' else m['core_reports'];ctx=m['FAERS_context'] if db=='FAERS' else m['context_reports'];cx=m['FAERS_context_core'] if db=='FAERS' else m['context_core_reports']
  strings=([db+'\n2004 Q1–2025 Q4 releases','DEMO input rows\n23,992,742','After DELETE / case-version\nretention: 19,963,024','After 1,788 identifier-collision\nexclusions: 19,961,236'] if db=='FAERS' else [db+('\nThrough May 2026' if db=='Canada' else '\nSeptember 2026 release'),'Source metadata rows\n'+num(m['metadata_rows']),'After latest-version retention\n'+num(total),'Retained reports\n'+num(total)+'\nUnlinked drug rows: 0'])
  strings=[s.replace('After latest-version retention','After latest-version\nretention') for s in strings]
  strings+=['Outcome partition\nCore '+num(core)+'\nNon-core '+num(total-core),'Medication-marker context\nTotal '+num(ctx)+'\nCore '+num(cx)+'; non-core '+num(ctx-cx)]
  ys=[6.5,5.45,4.4,3.35,2.2,1.05]
  for y,s in zip(ys,strings):ax.text(x+.5,y,s,ha='center',va='center',fontsize=7.6,bbox=dict(boxstyle='square,pad=.55',fc='white',ec=COLORS[db],lw=.8))
  for a,b in zip(ys,ys[1:]):ax.annotate('',xy=(x+.5,b+.39),xytext=(x+.5,a-.32),arrowprops=dict(arrowstyle='->',lw=.6,color='#777'))
 ax.text(1.5,.1,'Separate reporting systems; identical medication-context eligibility in core and non-core arms.',ha='center',fontsize=7.5);export(fig,'Fig1_workflow')
 fig,axes=plt.subplots(1,2,figsize=(7.1,3.8))
 for j,ax in enumerate(axes):
  ax.set_xscale('log');ax.set_xlim(.06,30);ax.set_ylim(3.6,-.6);ax.set_yticks(range(4),[NAMES[a] for a in AGENTS] if j==0 else ['']*4);ax.axvline(1,color='#888',ls='--',lw=.8);ax.grid(axis='x',color='#eee',lw=.5);ax.set_xlabel('ROR (95% CI)');ax.text(-.16,1.04,chr(97+j),transform=ax.transAxes,fontweight='bold',fontsize=10)
  for i,ag in enumerate(AGENTS):
   specs=[('FAERS','baseline',-.09,'#777','o'),('FAERS','medication_context',.09,COLORS['FAERS'],'s')] if j==0 else [(db,'medication_context',(k-1)*.19,COLORS[db],m) for k,(db,m) in enumerate(zip(COLORS,['o','s','^']))]
   for db,setting,off,col,m in specs:
    v=pick(p,database=db,agent=ag,setting=setting);y=i+off
    if pd.notna(v.ror):
     ax.errorbar(v.ror,y,xerr=[[v.ror-v.ror_ci_low],[v.ror_ci_high-v.ror]],fmt=m,color=col,ms=4,lw=.8)
     if j:ax.annotate(f'n={int(v.a)}',(v.ror_ci_high,y),xytext=(3,0),textcoords='offset points',va='center',fontsize=6.3,color=col)
    else:ax.text(.07,y,f'{db}: n={int(v.a)}; sparse',fontsize=6.3,color=col,va='center')
 for label,col,m in [('Unrestricted same-class','#777','o'),('Medication context',COLORS['FAERS'],'s')]:axes[0].plot([],[],color=col,marker=m,label=label)
 for db,m in zip(COLORS,['o','s','^']):axes[1].plot([],[],color=COLORS[db],marker=m,label=db)
 for ax in axes:ax.legend(frameon=False,fontsize=6.5,loc='lower left',bbox_to_anchor=(0,1.01))
 fig.tight_layout(w_pad=2);export(fig,'Fig2_primary_contrasts')
 burden=data('burden_current_strata.csv');fp=data('fingerprint_primary_sensitivity.csv');role=data('Table_S14_exact_roles_and_reporter.csv');sug=data('Table_S21c_sugammadex_exclusion.csv');event=data('Table_S22b_event_exclusions.csv');out=data('Table_S23b_outcome_sensitivity.csv');geo=data('Table_S24_geography_overlap_sensitivity.csv');ctx=data('Table_S17_context_marker_sensitivity.csv');remove=data('Table_S15_comparator_and_cosuspect_exclusions.csv')
 labels=['Primary context','Drug count 1–5','Drug count 6–10','Drug count 11–20','Drug count 21+','Strict fingerprint merge','Fingerprint + role signature','No sugammadex (any role)','No rescue-only eligibility','No JP/CA reports','HCP-restricted PS+SS','PS only','Core + anaphylactoid','Official narrow SMQ','No rocuronium reports','No pharmacologic competitors']
 fig,axes=plt.subplots(1,4,figsize=(7.2,6.7),sharey=True)
 for j,(ax,ag) in enumerate(zip(axes,AGENTS)):
  rr=[pick(p,database='FAERS',agent=ag,setting='medication_context')]+[pick(burden,agent=ag,burden_stratum=s) for s in ['1–5','6–10','11–20','21+']]+[pick(fp,agent=ag,fingerprint_rule=s) for s in ['strict_demographic_date_PT','strict_plus_role_signature']]+[pick(sug,agent=ag,setting='medication_context',diagnostic='exclude_sugammadex_any_role'),pick(ctx,agent=ag,context_sensitivity='four_marker_context_no_rescue_only'),pick(geo,agent=ag,setting='medication_context',country_partition='exclude_JP_CA'),pick(role,agent=ag,role_scope='PS_SS',reporter='HCP_MD_PH_OT',setting='medication_context'),pick(role,agent=ag,role_scope='PS_only',reporter='all_reporters',setting='medication_context'),pick(out,agent=ag,diagnostic='core_plus_anaphylactoid'),pick(out,agent=ag,diagnostic='smq_narrow_normalized'),pick(remove,agent=ag,setting='medication_context',excluded_agent='rocuronium') if ag=='succinylcholine' else None,pick(event,agent=ag,setting='medication_context',diagnostic='exclude_pharmacologic_competitors') if ag in AGENTS[:2] else None]
  ax.set_xscale('log');ax.set_xlim(.035,40);ax.set_ylim(15.6,-.6);ax.set_yticks(range(16),labels);ax.tick_params(axis='y',labelsize=7);ax.tick_params(axis='x',labelsize=7);ax.axvline(1,color='#888',ls='--',lw=.7);ax.grid(axis='x',color='#eee',lw=.5);ax.set_title(NAMES[ag],fontsize=8);ax.set_xlabel('ROR (95% CI)',fontsize=7)
  for y,v in enumerate(rr):
   if v is None:ax.text(.06,y,'Not applicable',fontsize=5.6,color='#777',va='center')
   elif pd.notna(v.ror):ax.errorbar(v.ror,y,xerr=[[v.ror-v.ror_ci_low],[v.ror_ci_high-v.ror]],fmt='o',color=COLORS['FAERS'] if y==0 else '#555',ms=3,lw=.7)
   else:ax.text(.06,y,f'n={int(v.a)}; sparse',fontsize=5.6,color='#555',va='center')
  for y in [.5,4.5,6.5,11.5,13.5]:ax.axhline(y,color='#ddd',lw=.5)
 fig.subplots_adjust(left=.31,right=.99,wspace=.22,bottom=.09,top=.95);export(fig,'Fig3_diagnostics')

def grouped_figure():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    p=read('Table_S4_S5_canonical_primary.csv');burden=read('burden_current_strata.csv');fp=read('fingerprint_primary_sensitivity.csv');role=read('Table_S14_exact_roles_and_reporter.csv');sug=read('Table_S21c_sugammadex_exclusion.csv');event=read('Table_S22b_event_exclusions.csv');out=read('Table_S23b_outcome_sensitivity.csv');geo=read('Table_S24_geography_overlap_sensitivity.csv');ctx=read('Table_S17_context_marker_sensitivity.csv');remove=read('Table_S15_comparator_and_cosuspect_exclusions.csv')
    labels=['Primary context','Drug count 1–5','Drug count 6–10','Drug count 11–20','Drug count 21+','Strict fingerprint merge','Fingerprint + role signature','No sugammadex (any role)','No rescue-only eligibility','No JP/CA reports','HCP-restricted PS+SS','PS only','Core + anaphylactoid','Official narrow SMQ','No rocuronium reports','No pharmacologic competitors']
    fig,axes=plt.subplots(1,4,figsize=(7.2,6.9),sharey=True)
    for ax,ag in zip(axes,AGENTS):
        rr=[pick(p,database='FAERS',agent=ag,setting='medication_context')]+[pick(burden,agent=ag,burden_stratum=s) for s in ['1–5','6–10','11–20','21+']]+[pick(fp,agent=ag,fingerprint_rule=s) for s in ['strict_demographic_date_PT','strict_plus_role_signature']]+[pick(sug,agent=ag,setting='medication_context',diagnostic='exclude_sugammadex_any_role'),pick(ctx,agent=ag,context_sensitivity='four_marker_context_no_rescue_only'),pick(geo,agent=ag,setting='medication_context',country_partition='exclude_JP_CA'),pick(role,agent=ag,role_scope='PS_SS',reporter='HCP_MD_PH_OT',setting='medication_context'),pick(role,agent=ag,role_scope='PS_only',reporter='all_reporters',setting='medication_context'),pick(out,agent=ag,diagnostic='core_plus_anaphylactoid'),pick(out,agent=ag,diagnostic='smq_narrow_normalized'),pick(remove,agent=ag,setting='medication_context',excluded_agent='rocuronium') if ag=='succinylcholine' else None,pick(event,agent=ag,setting='medication_context',diagnostic='exclude_pharmacologic_competitors') if ag in AGENTS[:2] else None]
        ax.set_xscale('log');ax.set_xlim(.035,40);ax.set_ylim(15.6,-.6);ax.set_yticks(range(16),labels);ax.tick_params(axis='y',labelsize=6.8);ax.tick_params(axis='x',labelsize=7);ax.axvline(1,color='#888',ls='--',lw=.7);ax.grid(axis='x',color='#eee',lw=.5);ax.set_title(NAMES[ag],fontsize=8);ax.set_xlabel('ROR (95% CI)',fontsize=7)
        for y,v in enumerate(rr):
            if v is None:ax.text(.045,y,'Not applicable',fontsize=5.5,color='#777',va='center')
            elif pd.notna(v.ror):ax.errorbar(v.ror,y,xerr=[[v.ror-v.ror_ci_low],[v.ror_ci_high-v.ror]],fmt='o',color=COLORS['FAERS'] if y==0 else '#555',ms=3,lw=.7)
            else:ax.text(.045,y,f'n={int(v.a)}; sparse',fontsize=5.5,color='#555',va='center')
        for y in [.5,4.5,6.5,11.5,13.5]:ax.axhline(y,color='#ddd',lw=.5)
    fig.subplots_adjust(left=.44,right=.99,wspace=.22,bottom=.08,top=.95)
    for label,y in [('Reference',0),('Medication\nburden',2.5),('Duplicate\nmerging',5.5),('Population\n& role',9.5),('Outcome',12.5),('Comparator/\nevents',14.5)]:
        yy=fig.transFigure.inverted().transform(axes[0].transData.transform((1,y)))[1]
        fig.text(.012,yy,label,fontsize=6.8,fontweight='bold',va='center',ha='left',color='#333')
    export(fig,'Fig3_diagnostics')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,default=R/'figures')
    args=parser.parse_args();OUT=args.output_dir;OUT.mkdir(parents=True,exist_ok=True)
    src=pd.read_csv(R/'Table_S1a_sources_and_retention.csv');fa=src[src.database.eq('FAERS')].iloc[0]
    flow={'FAERS_retained':int(fa.retained_reports),'FAERS_core':int(fa.core_reports),'FAERS_context':int(fa.context_reports),'FAERS_context_core':int(fa.context_core_reports)}
    cj=src[src.database.ne('FAERS')].to_dict('records')
    plots(read('Table_S4_S5_canonical_primary.csv'),flow,cj);grouped_figure()
    print('Generated Fig1–Fig3 from released aggregate tables.')
