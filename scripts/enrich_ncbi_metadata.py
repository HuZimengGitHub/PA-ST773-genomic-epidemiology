#!/usr/bin/env python3
"""Enrich an auditable 626-assembly ST773 table from NCBI BioSample, then recalculate Figure 4.

No extra Python packages required. Never assigns country from SNP similarity.
Usage: python scripts/enrich_ncbi_metadata.py [--online] [--truststore]
Default: offline; --online sends actual NCBI BioSample requests.
"""
import argparse
import csv
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
START = BASE / 'data/ST773/metadata_626.tsv'
SNP = BASE / 'data/ST773/snp_distance_matrix.tsv'
RES = BASE / 'data/ST773/resfinder_summary.tsv'
OUTPUT = BASE / 'results/NCBI_enriched_metadata.tsv'
REPORT = BASE / 'results/NCBI_metadata_enrichment_report.tsv'
PAIRS = BASE / 'results/NCBI_recalculated_close_pairs.tsv'
SUMMARY = BASE / 'results/NCBI_recalculated_summary.txt'
ACC = 'Assembly Accession'
BS = 'Assembly BioSample Accession'
DATE = 'Assembly BioSample Collection date'
LOC = 'Assembly BioSample Geographic location'


def read_tsv(path):
    with open(path, newline='', encoding='utf-8-sig') as f:
        rd = csv.DictReader(f, delimiter='\t')
        return list(rd), rd.fieldnames


def write_tsv(path, rows, fields):
    with open(path, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def country(geo):
    x = (geo or '').strip()
    if not x or x.lower() in {
        'unknown', 'missing', 'not available', 'not provided', 'not collected',
        'not applicable', 'restricted access', 'unreported', 'n/a', 'na', 'none', '-',
    } or x.lower().startswith('missing:'):
        return None
    x = x.split(':', 1)[0].strip()
    return {'Czechia': 'Czech Republic', 'Korea': 'South Korea',
            'Republic of Korea': 'South Korea', 'UK': 'United Kingdom',
            'United States': 'USA', 'United States of America': 'USA'}.get(x, x)


def need_lookup(row):
    # BioSamples with explicit 'restricted access' / 'not collected' are respected.
    raw = (row[LOC] or '').strip().lower()
    return raw in {'', 'unknown', 'not provided', 'not available', 'missing', 'n/a', 'na'}


def get_xml(endpoint, args, tries=2):
    params = dict(args)
    params['tool'] = 'ST773MetadataVerification'
    email = os.environ.get('NCBI_EMAIL', '').strip()
    if email:
        params['email'] = email
    url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/{endpoint}?{urllib.parse.urlencode(params)}'
    for trial in range(tries):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'ST773MetadataVerification/1.0'})
            with urllib.request.urlopen(req, timeout=35) as res:
                return ET.fromstring(res.read())
        except Exception as e:
            if trial + 1 == tries:
                raise
            time.sleep(3.0)
    return None


def assembly_biosample(accession):
    # Resolve assembly accession to its *linked* BioSample. No approximate matches.
    query = f'{accession}[Assembly Accession]'
    root = get_xml('esearch.fcgi', {'db': 'assembly', 'term': query, 'retmax': '10'})
    uids = [e.text for e in root.findall('.//IdList/Id') if e.text]
    if not uids:
        return None
    summary = get_xml('esummary.fcgi', {'db': 'assembly', 'id': ','.join(uids), 'retmode': 'xml'})
    for item in summary.findall('.//DocumentSummary'):
        found = (item.findtext('AssemblyAccession') or '').strip()
        if found == accession:
            bs = (item.findtext('BioSampleAccn') or '').strip()
            return bs if bs.startswith('SAM') else None
    return None


def fetch_biosample(biosample):
    root = get_xml('efetch.fcgi', {'db': 'biosample', 'id': biosample, 'retmode': 'xml'})
    for sample in root.findall('.//BioSample'):
        if sample.get('accession', '') != biosample:
            continue
        fields = {}
        for a in sample.findall('.//Attribute'):
            key = (a.get('harmonized_name') or a.get('attribute_name') or '').lower().strip()
            val = (a.text or '').strip()
            if key and val:
                fields[key] = val
        loc = fields.get('geo_loc_name') or fields.get('geographic location') or fields.get('geographic_location') or fields.get('country') or ''
        date = fields.get('collection_date') or fields.get('collection date') or ''
        return loc, date
    return '', ''


def recalculate(rows):
    by_id = {x[ACC]: x for x in rows}
    with open(SNP, newline='', encoding='utf-8') as f:
        rd = csv.reader(f, delimiter='\t')
        head = next(rd)
        genomes = []
        matrix = []
        for line in rd:
            name = line[0]
            m = re.search(r'GCA_\d+\.\d+', name)
            genomes.append(m.group(0) if m else None)
            matrix.append(line[1:])
    if len(genomes) != 628 or len({g for g in genomes if g}) != 626:
        raise ValueError('Unexpected SNP matrix dimensions or assembly IDs')
    if set(x for x in genomes if x) != set(by_id):
        raise ValueError('Metadata and SNP matrix differ in assembly IDs')
    geocoded = {a: country(r[LOC]) for a, r in by_id.items()}
    ndm = {}
    res, _ = read_tsv(RES)
    for r in res:
        m = re.search(r'GCA_\d+\.\d+', r['#FILE'])
        if not m:
            continue
        ndm[m.group(0)] = any((val or '').strip() not in ('', '.', '-')
                                 for key, val in r.items() if 'NDM' in key.upper())
    if set(ndm) != set(by_id):
        raise ValueError('ResFinder and metadata accession IDs do not match')
    valid_indices = [i for i, a in enumerate(genomes) if a and geocoded[a]]
    eligible = 0
    close = []
    for a_idx, i in enumerate(valid_indices):
        a = genomes[i]
        for j in valid_indices[a_idx+1:]:
            b = genomes[j]
            if geocoded[a] == geocoded[b]:
                continue
            eligible += 1
            d = int(matrix[i][j])
            if d < 5:
                close.append({'GCA_1': a, 'Country_1': geocoded[a], 'GCA_2': b,
                              'Country_2': geocoded[b], 'SNP_distance': d,
                              'Both_NDM_positive': int(ndm[a] and ndm[b])})
    cols = ['GCA_1', 'Country_1', 'GCA_2', 'Country_2', 'SNP_distance', 'Both_NDM_positive']
    write_tsv(PAIRS, close, cols)
    both = sum(x['Both_NDM_positive'] for x in close)
    report = (f'Independent assemblies: {len(by_id)}\n'
              f'Country-resolved assemblies: {len(valid_indices)}\n'
              f'Cross-country pairs eligible for SNP comparison: {eligible}\n'
              f'Cross-country pairs with SNP < 5: {len(close)}\n'
              f'Both blaNDM-positive: {both}\n'
              f'At least one blaNDM-negative: {len(close)-both}\n'
              f'Unresolved geographic origin: {len(by_id)-len(valid_indices)}\n'
              f'Reference and snippy_reference.fna excluded from genome pair counts.\n')
    SUMMARY.write_text(report, encoding='utf-8')
    print('\n' + report)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--online', action='store_true', help='Opt into NCBI queries (offline by default)')
    parser.add_argument('--truststore', action='store_true', help='Use OS certificate store (requires truststore package)')
    args = parser.parse_args()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    if args.truststore:
        import truststore
        truststore.inject_into_ssl()
    rows, headers = read_tsv(START)
    if len(rows) != 626 or len({r[ACC] for r in rows}) != 626:
        raise ValueError('Input is not 626 unique assemblies')
    targets = [r for r in rows if need_lookup(r)]
    print(f'626 assembly records; {len(targets)} blanks/unknown to re-check via BioSample')
    changes = []
    for n, row in enumerate(targets, 1):
        acc = row[ACC]
        if not args.online:
            changes.append({'Assembly': acc, 'BioSample': row.get(BS, ''), 'New_location': '', 'New_date': '', 'Status': 'offline'})
            continue
        try:
            bs = (row.get(BS) or '').strip()
            if not bs:
                bs = assembly_biosample(acc) or ''
                time.sleep(0.45)
            if bs and not row.get(BS):
                row[BS] = bs
            loc, date = fetch_biosample(bs) if bs else ('', '')
            if loc and country(loc):
                row[LOC] = loc
                row['Country_source'] = 'NCBI BioSample ' + bs
            if date and (not row[DATE].strip()):
                row[DATE] = date
                row['Date_source'] = 'NCBI BioSample ' + bs
            status = 'country filled' if loc and country(loc) else ('country still unknown' if bs else 'no linked BioSample')
            changes.append({'Assembly': acc, 'BioSample': bs, 'New_location': loc,
                            'New_date': date, 'Status': status})
            print(f'[{n}/{len(targets)}] {acc}: {bs} | {loc or "no country"}')
            time.sleep(0.45)
        except Exception as e:
            print(f'[{n}/{len(targets)}] {acc}: request failed ({e})', file=sys.stderr)
            changes.append({'Assembly': acc, 'BioSample': row.get(BS, ''), 'New_location': '',
                            'New_date': '', 'Status': 'error: ' + str(e)[:150]})
    write_tsv(OUTPUT, rows, headers)
    write_tsv(REPORT, changes, ['Assembly', 'BioSample', 'New_location', 'New_date', 'Status'])
    recalculate(rows)
    print('Output:', OUTPUT, REPORT, PAIRS, SUMMARY, sep='\n  ')


if __name__ == '__main__':
    main()
