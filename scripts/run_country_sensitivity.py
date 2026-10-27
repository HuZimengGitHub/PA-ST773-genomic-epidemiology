#!/usr/bin/env python3
"""Country-stratified descriptive sensitivity analysis (NO BioProject metadata needed).
Recomputes O11/O6 shift-share ST-frequency and within-ST components under sampling choices.
"""
import argparse, sys
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parent))
from run_core_analyses import load,shift_share


def summarize(label,base,weighting='isolate',min_n=1):
    early=base[base.Year.between(2005,2009)].reset_index(drop=True)
    late=base[base.Year.between(2020,2024)].reset_index(drop=True)
    if not len(early) or not len(late):return []
    if weighting=='equal_country':
        weights=[]
        for sub in (early,late):
            size=sub.groupby('country').size()
            eligible=size[size>=min_n].index
            w=np.array([1/size.at[c] if c in eligible else 0.0 for c in sub.country],dtype=float)
            if w.sum()==0:return []
            weights.append(w)
        w0,w1=weights
        c0=int((w0>0).sum());c1=int((w1>0).sum())
        n0=int((w0>0).sum());n1=int((w1>0).sum())
    else:
        w0=w1=None
        n0=len(early);n1=len(late)
    output=[]
    for sero in ('O11','O6'):
        result,_=shift_share(early,late,sero,w0,w1)
        result.update({'analysis':label,'weighting':weighting,'min_n_country_period':min_n,'early_n':n0,'late_n':n1,
                       'early_countries':early.country.nunique(),'late_countries':late.country.nunique()})
        output.append(result)
    return output


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--input',default='data/PA_7071.xlsx')
    ap.add_argument('--output',default='results')
    args=ap.parse_args()
    out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    d=load(args.input)
    d=d[d.Serotype.notna() & d.Year.between(2005,2024)].copy()
    known=d[d.country.notna()].copy()
    rows=[]
    rows+=summarize('all genomes (primary)',d)
    rows+=summarize('country metadata available, isolate pooled',known)
    for threshold in [1,5,10,20]:
        rows+=summarize('equal country weighting',known,'equal_country',threshold)
        early_ct=known[known.Year.between(2005,2009)].groupby('country').size()
        late_ct=known[known.Year.between(2020,2024)].groupby('country').size()
        stable=set(early_ct[early_ct>=threshold].index)&set(late_ct[late_ct>=threshold].index)
        if stable:
            sub=known[known.country.isin(stable)]
            rows+=summarize(f'stable countries (both periods n≥{threshold})',sub)
            rows+=summarize(f'stable countries equal weighted (both periods n≥{threshold})',sub,'equal_country',threshold)
    top_countries=known.country.value_counts().head(10).index
    for country in top_countries:
        rows+=summarize('leave country out: '+str(country),known[known.country.ne(country)])
    result=pd.DataFrame(rows)
    result.to_csv(out/'country_shiftshare_sensitivity.csv',index=False)
    # Report country-period sample sizes to disclose which countries can contribute to comparisons.
    n=known.assign(period=np.select([known.Year.between(2005,2009),known.Year.between(2020,2024)],
                                     ['early','late'],default='middle'))
    n.groupby(['country','period']).size().unstack(fill_value=0).to_csv(out/'country_period_sample_sizes.csv')
    print('Saved',out,'with',len(result),'O11/O6 scenario rows')
    print(result.loc[result.Serotype.eq('O11'),['analysis','weighting','min_n_country_period','early_n','late_n',
               'Observed_change_pp','ST_frequency_component_pp','ST_component_percent_of_change']].round(3).to_string(index=False))

if __name__=='__main__':main()
