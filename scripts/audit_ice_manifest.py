#!/usr/bin/env python3
"""Audit 17 explicitly labeled ICE tips; this does NOT infer a missing 12th ICE-A."""
import csv
import json
from pathlib import Path
from collections import Counter
from Bio import Phylo
ROOT=Path(__file__).resolve().parents[1]

def main():
    tree=Phylo.read(ROOT/'data/ST773/core.aln.treefile','newick')
    tips=list(tree.get_terminals()); names=[t.name for t in tips]
    with (ROOT/'data/ST773/ICE_assignments_from_legacy_plot.csv').open(newline='') as f:
        labels=list(csv.DictReader(f))
    count=Counter(x['Variant'] for x in labels)
    absent=[x['Assembly'] for x in labels if not any(x['Assembly'] in str(t) for t in names)]
    summary={'total_tree_tips':len(tips),'ICE_A_in_legacy_plot':count['ICE-A'],
             'ICE_B_in_legacy_plot':count['ICE-B'], 'ICE_A_expected_in_manuscript':12,
             'missing_labeled_tips':absent,
             'unresolved':'No sequence-level ICE classification supplied for the manuscript 12th ICE-A.'}
    dest=ROOT/'results/ice_manifest_audit.json';dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
