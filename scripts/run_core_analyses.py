#!/usr/bin/env python3
"""Core statistics from the frozen 7,071-genome source.
Usage: python scripts/run_core_analyses.py --input data/PA_7071.xlsx --output results
For inferential residuals, see paper; individual Pearson residuals are descriptive.
"""
import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, kendalltau

PERIODS = {'2005-2009':(2005,2009),'2010-2014':(2010,2014),
           '2015-2019':(2015,2019),'2020-2024':(2020,2024)}
SEROTYPES = ['O11','O6','O1','O4','O2','O3','O5','O12','O10','O7','O9','O13','O15','O17']
INVALID_SEROTYPES = {'','-','NT','MULTIPLE','NAN','NA','N/A','NONE','<NA>'}


def norm_st(x):
    if pd.isna(x) or str(x).strip().upper() in {'','-','UNKNOWN','NAN','NONE','NA','N/A','<NA>'}:
        return 'Unknown'
    s = str(x).strip()
    if s.upper().startswith('ST'):
        s = s[2:].strip()
    if re.fullmatch(r'\d+(?:\.0+)?',s):
        return f'ST{int(float(s))}'
    return 'ST'+s


def norm_serotype(x):
    if pd.isna(x):return pd.NA
    s = str(x).strip().upper()
    return pd.NA if s in INVALID_SEROTYPES else s


def load(path):
    raw = pd.read_excel(path)
    if 'Serotyope' in raw and 'Serotype' not in raw:
        raw = raw.rename(columns={'Serotyope':'Serotype'})
    for name in ('ST','Serotype','Year','country','Sample'):
        if name not in raw:raise ValueError(f'Missing input column: {name}')
    d = raw[['Sample','ST','Serotype','Year','country']].copy()
    d['ST'] = d['ST'].map(norm_st)
    d['Serotype'] = d['Serotype'].map(norm_serotype).astype('string')
    d['Year'] = pd.to_numeric(d['Year'],errors='coerce')
    d['country'] = d['country'].astype('string').str.strip()
    d['country'] = d['country'].mask(d['country'].str.upper().isin(['','NAN','NA','-','NONE','NOT COLLECTED','RESTRICTED ACCESS']), pd.NA)
    d['Sample'] = d['Sample'].astype(str)
    if d['Sample'].duplicated().any():raise ValueError('Duplicate Sample accession; deduplicate before analysis')
    return d


def contingency(d):
    table = pd.crosstab(d['ST'],d['Serotype'])
    chi2, p, dof, expected = chi2_contingency(table.values, correction=False)
    n=int(table.values.sum())
    v=float(np.sqrt(chi2/(n*(min(table.shape)-1))))
    exp=pd.DataFrame(expected,index=table.index,columns=table.columns)
    residuals=(table-exp)/np.sqrt(exp)
    return table,exp,residuals,{'n':n,'chi2':float(chi2),'dof':int(dof),
         'chi_square_pvalue':float(p),'cramers_v':v,'expected_lt5_fraction':float((expected<5).mean())}


def shift_share(early, late, serotype, w0=None,w1=None):
    """Symmetric exact two-component identity. Missing-period p(s|ST)=0 by convention."""
    early=early.reset_index(drop=True);late=late.reset_index(drop=True)
    w0=np.ones(len(early),dtype=float) if w0 is None else np.asarray(w0,dtype=float)
    w1=np.ones(len(late),dtype=float) if w1 is None else np.asarray(w1,dtype=float)
    if min(w0.sum(),w1.sum())<=0: raise ValueError('Empty/zero-weight period')
    w0=w0/w0.sum();w1=w1/w1.sum()
    d0=early[['ST','Serotype']].copy();d1=late[['ST','Serotype']].copy()
    d0['wt']=w0;d1['wt']=w1
    f0=d0.groupby('ST').wt.sum();f1=d1.groupby('ST').wt.sum()
    g0=d0[d0.Serotype.eq(serotype)].groupby('ST').wt.sum()
    g1=d1[d1.Serotype.eq(serotype)].groupby('ST').wt.sum()
    union=sorted(set(f0.index)|set(f1.index))
    details=[]
    for st in union:
        a=float(f0.get(st,0));b=float(f1.get(st,0))
        pa=float(g0.get(st,0))/a if a else 0
        pb=float(g1.get(st,0))/b if b else 0
        C=(b-a)*(pa+pb)/2;W=(pb-pa)*(a+b)/2
        details.append((st,C,W,a,b,pa,pb))
    C=sum(v[1] for v in details);W=sum(v[2] for v in details)
    obs=float(d1.loc[d1.Serotype.eq(serotype),'wt'].sum()-d0.loc[d0.Serotype.eq(serotype),'wt'].sum())
    if not np.isclose(C+W,obs,atol=1e-11):
        raise AssertionError(f'Shift-share reconstruction failed for {serotype}: {C+W} vs {obs}')
    result={'Serotype':serotype,'Observed_change_pp':obs*100,
            'ST_frequency_component_pp':C*100,'Within_ST_component_pp':W*100,
            'ST_component_percent_of_change':(C/obs*100 if abs(obs)>1e-14 else np.nan),
            'reconstruction_error_pp':(C+W-obs)*100}
    detail=pd.DataFrame(details,columns=['ST','C','W','f_early','f_late','p_early','p_late'])
    detail['C_pp']=detail.pop('C')*100
    detail['W_pp']=detail.pop('W')*100
    return result, detail


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--input',default='data/PA_7071.xlsx')
    ap.add_argument('--output',default='results')
    args=ap.parse_args()
    out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    d=load(args.input)
    d.to_csv(out/'cleaned_metadata.csv',index=False)
    typed=d[d.Serotype.notna()].copy()
    known=typed[typed.ST.ne('Unknown')].copy()
    checks={'source_genomes':len(d),'resolved_serotype_all_years':len(typed),'known_ST_and_serotype':len(known),
            'year_present_all_years':int(d.Year.notna().sum()),'country_present_all_years':int(d.country.notna().sum())}
    a,exp,res,stats=contingency(known)
    _,_,_,stats_unknown=contingency(typed)
    checks.update({'known_ST':stats,'unknown_ST_included':stats_unknown})
    a.to_csv(out/'ST_serotype_observed_knownST.csv')
    exp.to_csv(out/'ST_serotype_expected_knownST.csv')
    res.to_csv(out/'ST_serotype_Pearson_residual_knownST.csv')
    purity=a.div(a.sum(axis=1),axis=0)
    purity.to_csv(out/'ST_serotype_fidelity_proportions.csv')
    modern=typed[typed.Year.between(2005,2024)].copy()
    all_modern=d[d.Year.between(2005,2024)].copy()
    per=[];ss=[]
    for p,(start,end) in PERIODS.items():
        sub=modern[modern.Year.between(start,end)]
        checks['typed_'+p]=len(sub)
        for sero,n in sub.Serotype.value_counts().items():
            per.append({'period':p,'Serotype':sero,'n':int(n),'total_n':len(sub),'prevalence_pct':n/len(sub)*100})
    period_sero=pd.DataFrame(per)
    period_sero.to_csv(out/'serotype_period_prevalence.csv',index=False)
    early=modern[modern.Year.between(2005,2009)]
    late=modern[modern.Year.between(2020,2024)]
    for sero in sorted(modern.Serotype.dropna().unique()):
        r,detail=shift_share(early,late,sero)
        ss.append(r)
        if sero=='O11': detail.to_csv(out/'O11_shift_share_ST_contributions.csv',index=False)
    ss=pd.DataFrame(ss).sort_values('Serotype')
    ss.to_csv(out/'shift_share_all_serotypes.csv',index=False)
    annual=[]
    for denom,dataset in [('typed_only',modern),('all_year_assigned',all_modern)]:
        for year in range(2005,2025):
            x=dataset[dataset.Year.eq(year)]
            n=len(x);o=int(x.Serotype.eq('O11').sum())
            if n: annual.append({'denominator':denom,'Year':year,'n':n,'O11_n':o,'O11_percent':o/n*100})
    annual=pd.DataFrame(annual)
    annual.to_csv(out/'annual_O11_two_denominators.csv',index=False)
    for denom in annual.denominator.unique():
        z=annual[annual.denominator.eq(denom)]
        kt=kendalltau(z.Year,z.O11_percent)
        checks['kendall_'+denom]={'tau':float(kt.statistic),'p':float(kt.pvalue),'n_years':len(z)}
        z=z[z.n>=50]
        kt=kendalltau(z.Year,z.O11_percent)
        checks['kendall_'+denom+'_n_ge_50']={'tau':float(kt.statistic),'p':float(kt.pvalue),'n_years':len(z)}
    st_rows=[]
    for st in ['ST463','ST773','ST235','ST111','ST253','ST395']:
        all_st=all_modern[all_modern.ST.eq(st)]
        a0=all_modern[all_modern.Year.between(2005,2009)]
        a1=all_modern[all_modern.Year.between(2020,2024)]
        c0=int(a0.ST.eq(st).sum());c1=int(a1.ST.eq(st).sum())
        p0=c0/len(a0)*100;p1=c1/len(a1)*100
        st_rows.append({'ST':st,'n_2005_2024':len(all_st),'early_n':c0,'early_total_n':len(a0),'early_percent':p0,
                        'late_n':c1,'late_total_n':len(a1),'late_percent':p1,
                        'fold_change':p1/p0 if p0 else np.nan,
                        'first_year_2005_2024':int(all_st.Year.min()) if len(all_st) else None})
    pd.DataFrame(st_rows).to_csv(out/'selected_ST_early_late_changes.csv',index=False)
    checks['ST_early_raw_n']=len(all_modern[all_modern.Year.between(2005,2009)])
    checks['ST_late_raw_n']=len(all_modern[all_modern.Year.between(2020,2024)])
    checks['shift_O11']=ss.loc[ss.Serotype.eq('O11')].iloc[0].to_dict()
    checks['shift_O6']=ss.loc[ss.Serotype.eq('O6')].iloc[0].to_dict()
    # Invariant checks are appropriate for this specific 7,071-genome source collection.
    assert checks['source_genomes']==7071
    assert checks['known_ST_and_serotype']==4475
    assert [checks['typed_'+p] for p in PERIODS]==[361,955,2157,2196]
    assert abs(checks['known_ST']['cramers_v']-0.8873)<0.0001
    assert abs(checks['shift_O11']['ST_component_percent_of_change']-53.6956)<0.005
    with open(out/'audit_summary.json','w',encoding='utf-8') as f:json.dump(checks,f,indent=2,ensure_ascii=False)
    print('PASS: verified major calculations; written to',out)
    print(json.dumps(checks,indent=2,ensure_ascii=False))

if __name__=='__main__': main()
