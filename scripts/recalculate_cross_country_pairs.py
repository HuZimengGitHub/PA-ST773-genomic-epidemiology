#!/usr/bin/env python3
"""Recalculate auditable cross-country ST773 near-neighbor pairs from frozen inputs.

Independent GCA accessions only; excludes 'Reference' and 'snippy_reference.fna'.
Unknown / restricted locations are excluded from pair denominators. Each pair
is unordered and counted exactly once. A distance strictly LESS than the
threshold qualifies (default: <5 SNP). No geographic information is inferred
from phylogeny or SNP distances. Use frozen version for paper QA; do not silently
substitute this for the existing Figure 4 illustration.
"""
import argparse
import csv
import json
import re
from pathlib import Path
from collections import Counter

ROOT = Path(__file__).resolve().parents[1]
ACC = 'Assembly Accession'
LOC = 'Assembly BioSample Geographic location'
INVALID = {'', 'unknown', 'missing', 'not available', 'not provided', 'not collected',
           'not applicable', 'restricted access', 'unreported', 'n/a', 'na', 'none', '-'}
ALIASES = {'czechia': 'Czech Republic', 'czech republic': 'Czech Republic',
           'korea': 'South Korea', 'republic of korea': 'South Korea',
           'south korea': 'South Korea', 'uk': 'United Kingdom',
           'united kingdom': 'United Kingdom', 'united states': 'USA',
           'united states of america': 'USA', 'usa': 'USA'}


def country(geo):
    value = str(geo or '').strip()
    if value.lower() in INVALID or value.lower().startswith('missing:'):
        return None
    value = value.split(':', 1)[0].strip()
    if value.lower() in INVALID:
        return None
    return ALIASES.get(value.lower(), value)


def read_tsv(path):
    with open(path, encoding='utf-8-sig', newline='') as fp:
        rd = csv.DictReader(fp, delimiter='\t')
        return list(rd), rd.fieldnames


def write_tsv(path, rows, columns):
    with open(path, 'w', encoding='utf-8', newline='') as fp:
        w = csv.DictWriter(fp, delimiter='\t', fieldnames=columns)
        w.writeheader()
        w.writerows(rows)


def accession(text):
    match = re.search(r'GCA_\d+\.\d+', str(text))
    return match.group(0) if match else None


def recalculate(metadata, matrix, resfinder, output, threshold=5):
    output.mkdir(parents=True, exist_ok=True)
    meta, headers = read_tsv(metadata)
    if ACC not in headers or LOC not in headers:
        raise ValueError('Expected assembly accession and geographic location columns')
    ids = [r[ACC] for r in meta]
    if len(ids) != 626 or len(set(ids)) != 626:
        raise ValueError(f'Expected 626 unique assembly IDs; got {len(ids)} rows / {len(set(ids))} unique')
    geocoded = {r[ACC]: country(r[LOC]) for r in meta}

    with open(matrix, newline='', encoding='utf-8') as fp:
        rd = csv.reader(fp, delimiter='\t')
        head = next(rd)
        colids = [accession(x) for x in head[1:]]
        mat = []; rowids = []
        for row in rd:
            rowids.append(accession(row[0])); mat.append(row[1:])
    if len(mat) != 628 or len(colids) != 628 or any(len(r) != 628 for r in mat):
        raise ValueError('Unexpected SNP matrix dimensions (expected 628 x 628)')
    if colids != rowids:
        raise ValueError('Matrix row and column orders differ')
    actual = [x for x in rowids if x]
    if len(actual) != 626 or len(set(actual)) != 626 or set(actual) != set(ids):
        raise ValueError('SNP matrix GCA IDs do not match the 626 metadata accessions')

    genes, fields = read_tsv(resfinder)
    if '#FILE' not in fields:
        raise ValueError('ResFinder input missing #FILE')
    ndm_cols = [x for x in fields if 'NDM' in x.upper()]
    if not ndm_cols:
        raise ValueError('ResFinder input lacks NDM gene columns')
    ndm = {}
    for row in genes:
        a = accession(row['#FILE'])
        if a:
            if a in ndm:
                raise ValueError(f'Duplicate ResFinder accession {a}')
            ndm[a] = any(str(row[x] or '').strip().lower() not in {'', '.', '-', 'nan', 'none', '0'}
                         for x in ndm_cols)
    if set(ndm) != set(ids):
        raise ValueError('ResFinder and metadata GCA accession sets differ')

    indices = [i for i, name in enumerate(rowids) if name and geocoded[name]]
    eligible = 0; close = []; cnt=Counter()
    for pos, i in enumerate(indices):
        ai = rowids[i]; a_country = geocoded[ai]
        for j in indices[pos+1:]:
            aj = rowids[j]; b_country = geocoded[aj]
            if a_country == b_country: continue
            eligible += 1
            distance = int(mat[i][j])
            if distance < threshold:
                first, second = sorted([(ai, a_country), (aj, b_country)])
                both = int(ndm[ai] and ndm[aj])
                close.append({'GCA_1':first[0],'Country_1':first[1],
                              'GCA_2':second[0],'Country_2':second[1],
                              'SNP_distance':distance,'Both_NDM_positive':both})
                cnt['both_ndm'] += both
    close.sort(key=lambda x:(x['GCA_1'],x['GCA_2']))
    write_tsv(output/'cross_country_pairs_lt5.tsv', close,
              ['GCA_1','Country_1','GCA_2','Country_2','SNP_distance','Both_NDM_positive'])
    def show_path(p):
        try: return str(Path(p).resolve().relative_to(ROOT))
        except ValueError: return str(p)
    sums = {'assemblies':len(ids), 'matrix_tip_count':len(mat),
            'geocoded_assemblies':len(indices), 'unknown_country':len(ids)-len(indices),
            'cross_country_pairs_tested':eligible,'snp_threshold_strictly_less_than':threshold,
            'cross_country_pairs_below_threshold':len(close),
            'both_blaNDM_positive':cnt['both_ndm'],
            'at_least_one_blaNDM_negative':len(close)-cnt['both_ndm'],
            'reference_tips_excluded':2,
            'source_metadata':show_path(metadata),'source_snp':show_path(matrix),'source_resfinder':show_path(resfinder)}
    (output/'ST773_pair_audit.json').write_text(json.dumps(sums,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps(sums,indent=2,ensure_ascii=False))
    return sums


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--metadata',type=Path,default=ROOT/'data/ST773/metadata_626.tsv')
    p.add_argument('--snp',type=Path,default=ROOT/'data/ST773/snp_distance_matrix.tsv')
    p.add_argument('--resfinder',type=Path,default=ROOT/'data/ST773/resfinder_summary.tsv')
    p.add_argument('--output',type=Path,default=ROOT/'results')
    p.add_argument('--threshold',type=int,default=5,help='Use strict distance < threshold; default 5')
    a=p.parse_args()
    if a.threshold <= 0: p.error('threshold must be positive')
    recalculate(a.metadata,a.snp,a.resfinder,a.output,a.threshold)

if __name__=='__main__':main()
