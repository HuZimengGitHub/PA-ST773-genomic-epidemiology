#!/usr/bin/env python3
"""
Figure 4B - ST773 cross-country genomic links stratified by blaNDM status
Improvements: European inset, node color legend, cleaner labels
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from collections import defaultdict
import re
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

import cartopy.crs as ccrs
import cartopy.feature as cfeature

# ============ PATHS ============
ROOT = Path(__file__).resolve().parents[1]
SNP_MATRIX = str(ROOT / "data/ST773/snp_distance_matrix.tsv")
METADATA = str(ROOT / "data/ST773/metadata_626.tsv")
RESFINDER = str(ROOT / "data/ST773/resfinder_summary.tsv")
OUTDIR = str(ROOT / "results")
Path(OUTDIR).mkdir(parents=True, exist_ok=True)
SNP_CUTOFF = 5

# ============ COORDS ============
COORDS = {
    "China": (104.0, 35.0), "India": (78.0, 22.0), "USA": (-98.0, 38.0),
    "Nigeria": (8.0, 10.0), "Ghana": (-1.0, 8.0), "Cote d'Ivoire": (-5.5, 7.5),
    "South Korea": (127.8, 36.0), "Ukraine": (31.0, 49.0),
    "Netherlands": (5.3, 52.1), "Spain": (-3.7, 40.4), "Germany": (10.4, 51.2),
    "Czech Republic": (15.5, 49.8), "Switzerland": (8.2, 46.8),
    "United Kingdom": (-1.2, 52.5), "Norway": (8.5, 61.0),
    "Fiji": (178.0, -18.0), "South Africa": (25.0, -30.0),
    "Japan": (138.3, 36.2), "France": (2.2, 46.6), "Italy": (12.6, 42.5),
    "Brazil": (-51.9, -14.2), "Australia": (134.0, -25.0),
    "Canada": (-106.0, 56.0), "Morocco": (-7.1, 31.8),
    "Turkey": (35.2, 39.9), "Egypt": (30.8, 26.8),
    "Saudi Arabia": (45.1, 23.9), "Singapore": (103.8, 1.4),
    "Bangladesh": (90.4, 23.7),
}

MAIN_LABEL_OFFSETS = {
    "USA": (-15, -18), "Canada": (-10, 10), "Nigeria": (10, -12),
    "Ghana": (-30, -12), "Cote d'Ivoire": (-45, -15),
    "China": (10, -12), "India": (10, -12), "South Korea": (10, 8),
    "Japan": (10, 5), "Fiji": (-10, -15), "South Africa": (10, -10),
    "Brazil": (-10, -12), "Australia": (10, -10),
    "Morocco": (-10, -15), "Egypt": (10, -10),
    "Saudi Arabia": (10, -10), "Bangladesh": (10, -10),
    "Singapore": (10, -10), "Turkey": (12, 8),
    "Italy": (10, 8),
}

EUROPE_COUNTRIES = {"Netherlands", "Germany", "Spain", "Ukraine", "Czech Republic",
                     "Switzerland", "Norway", "France", "United Kingdom", "Italy"}

INSET_LABEL_OFFSETS = {
    "Netherlands": (12, 8), "Germany": (12, -5), "Spain": (-15, -15),
    "Ukraine": (12, -5), "Czech Republic": (12, -12), "Switzerland": (-15, -12),
    "Norway": (12, 8), "France": (-15, 8), "United Kingdom": (-40, 8),
    "Italy": (12, -8),
}

def extract_gca(name):
    m = re.search(r'(GCA_\d+\.\d+)', str(name))
    return m.group(1) if m else None

def parse_country(loc_string):
    """Normalize country values without assuming unknown countries from SNP links."""
    if pd.isna(loc_string):
        return 'Unknown'
    value = str(loc_string).split(':', 1)[0].strip()
    if value.lower() in {'', 'nan', 'unknown', 'not collected', 'not applicable',
                         'restricted access', 'missing', 'n/a', 'none', '-'}:
        return 'Unknown'
    aliases = {'republic of korea': 'South Korea', 'korea': 'South Korea',
               'south korea': 'South Korea', 'czechia': 'Czech Republic',
               'uk': 'United Kingdom', 'united states': 'USA',
               'united states of america': 'USA'}
    return aliases.get(value.lower(), value)


# ============ DATA PROCESSING ============
print("Reading data...")
meta = pd.read_csv(METADATA, sep='\t')
meta.columns = [c.strip() for c in meta.columns]
acc_col = [c for c in meta.columns if 'accession' in c.lower()][0]
loc_col = [c for c in meta.columns if 'location' in c.lower() or 'geographic' in c.lower()][0]

gca_to_country = {}
for _, row in meta.iterrows():
    gca = extract_gca(row[acc_col])
    if gca:
        gca_to_country[gca] = parse_country(row[loc_col])

resf = pd.read_csv(RESFINDER, sep='\t')
ndm_cols = [c for c in resf.columns if 'NDM' in c.upper()]
if not ndm_cols:
    raise ValueError('No blaNDM columns found in ResFinder table')
gca_to_ndm = {}
for _, row in resf.iterrows():
    gca = extract_gca(str(row.iloc[0]))
    if gca:
        gca_to_ndm[gca] = any(str(row[col]).strip().lower() not in ['.', '', 'nan', '0', 'none', '-']
                              for col in ndm_cols)

snp = pd.read_csv(SNP_MATRIX, sep='\t', index_col=0)
samples = list(snp.columns)
snp_to_gca = {}
for s in samples:
    gca = extract_gca(s)
    if gca:
        snp_to_gca[s] = gca
    # Ignore 'Reference' and 'snippy_reference.fna' (non-independent records)

ndm_plus_pairs = defaultdict(lambda: {"count": 0, "min_snp": 999999})
ndm_minus_pairs = defaultdict(lambda: {"count": 0, "min_snp": 999999})
skip_labels = {"Unknown"}
n = len(samples)
for i in range(n):
    for j in range(i+1, n):
        try:
            dist = int(snp.iloc[i, j])
        except:
            continue
        if dist >= SNP_CUTOFF:
            continue
        gca1, gca2 = snp_to_gca.get(samples[i]), snp_to_gca.get(samples[j])
        if not gca1 or not gca2 or gca1 == gca2:
            continue
        c1, c2 = gca_to_country.get(gca1, "Unknown"), gca_to_country.get(gca2, "Unknown")
        if c1 in skip_labels or c2 in skip_labels or c1 == c2:
            continue
        pair_key = tuple(sorted([c1, c2]))
        ndm1, ndm2 = gca_to_ndm.get(gca1, False), gca_to_ndm.get(gca2, False)
        target = ndm_plus_pairs if (ndm1 and ndm2) else ndm_minus_pairs
        target[pair_key]["count"] += 1
        target[pair_key]["min_snp"] = min(target[pair_key]["min_snp"], dist)

country_counts = defaultdict(int)
country_ndm_counts = defaultdict(int)
country_total = defaultdict(int)
for s in samples:
    gca = snp_to_gca.get(s)
    if gca:
        c = gca_to_country.get(gca, "Unknown")
        if c not in skip_labels:
            country_counts[c] += 1
            country_total[c] += 1
            if gca_to_ndm.get(gca, False):
                country_ndm_counts[c] += 1

total_np = sum(v["count"] for v in ndm_plus_pairs.values())
total_nm = sum(v["count"] for v in ndm_minus_pairs.values())
print(f"  NDM+ pairs: {total_np}, NDM- pairs: {total_nm}, fraction: {total_np/max(total_np+total_nm, 1)*100:.1f}%")

# ============ DRAW FUNCTION ============
def draw_map(ax, transform, countries_to_draw=None, show_labels=None,
             label_offsets=None, fontsize=6, node_scale=1.0):
    if countries_to_draw is None:
        countries_to_draw = set(country_counts.keys())
    if show_labels is None:
        show_labels = countries_to_draw
    if label_offsets is None:
        label_offsets = {}

    for pair, info in ndm_minus_pairs.items():
        c1, c2 = pair
        if c1 in countries_to_draw and c2 in countries_to_draw:
            if c1 in COORDS and c2 in COORDS:
                x1, y1 = COORDS[c1]
                x2, y2 = COORDS[c2]
                lw = min(0.5 + info["count"] * 0.15, 2.5) * node_scale
                ax.plot([x1, x2], [y1, y2], transform=ccrs.Geodetic(),
                        color='#4393C3', alpha=0.5, linewidth=lw, zorder=2)

    for pair, info in ndm_plus_pairs.items():
        c1, c2 = pair
        if c1 in countries_to_draw and c2 in countries_to_draw:
            if c1 in COORDS and c2 in COORDS:
                x1, y1 = COORDS[c1]
                x2, y2 = COORDS[c2]
                lw = min(0.7 + info["count"] * 0.2, 4.0) * node_scale
                ax.plot([x1, x2], [y1, y2], transform=ccrs.Geodetic(),
                        color='#D6604D', alpha=0.6, linewidth=lw, zorder=3)

    for country in countries_to_draw:
        if country not in COORDS or country not in country_counts:
            continue
        x, y = COORDS[country]
        count = country_counts[country]
        size = min(30 + count * 3.5, 220) * (node_scale ** 2)

        ndm_frac = country_ndm_counts.get(country, 0) / max(country_total.get(country, 1), 1)
        if ndm_frac > 0.7:
            color = '#D6604D'
        elif ndm_frac > 0.3:
            color = '#FC8D59'
        else:
            color = '#4393C3'

        ax.scatter(x, y, s=size, c=color, edgecolors='black',
                   linewidths=0.6, alpha=0.85, zorder=5, transform=transform)

        if country in show_labels:
            offset = label_offsets.get(country, (10, 8))
            ax.annotate(
                f"{country}\n(n={count})",
                xy=(x, y),
                xycoords=transform._as_mpl_transform(ax),
                textcoords='offset points',
                xytext=offset,
                fontsize=fontsize,
                color='#222222',
                fontweight='bold',
                zorder=6,
            )

# ============ PLOT ============
print("Plotting...")

fig = plt.figure(figsize=(16, 10))

# Main map
ax_main = fig.add_axes([0.02, 0.05, 0.96, 0.88], projection=ccrs.Robinson())
ax_main.set_global()
ax_main.add_feature(cfeature.LAND, facecolor='#f5f5f5', edgecolor='none')
ax_main.add_feature(cfeature.OCEAN, facecolor='white')
ax_main.add_feature(cfeature.BORDERS, linewidth=0.2, edgecolor='#dddddd')
ax_main.add_feature(cfeature.COASTLINE, linewidth=0.3, edgecolor='#aaaaaa')
transform = ccrs.PlateCarree()

all_countries = set(country_counts.keys())
non_europe_labels = all_countries - EUROPE_COUNTRIES
draw_map(ax_main, transform,
         countries_to_draw=all_countries,
         show_labels=non_europe_labels,
         label_offsets=MAIN_LABEL_OFFSETS,
         fontsize=5.5, node_scale=1.0)

ax_main.annotate("see inset →", xy=(15, 47),
                 xycoords=transform._as_mpl_transform(ax_main),
                 fontsize=7, color='#666666', fontstyle='italic', zorder=10)

# European inset
ax_inset = fig.add_axes([0.58, 0.52, 0.40, 0.44],
                         projection=ccrs.LambertConformal(central_longitude=12, central_latitude=50))
ax_inset.set_extent([-12, 42, 36, 66], crs=transform)
ax_inset.add_feature(cfeature.LAND, facecolor='#f5f5f5', edgecolor='none')
ax_inset.add_feature(cfeature.OCEAN, facecolor='white')
ax_inset.add_feature(cfeature.BORDERS, linewidth=0.3, edgecolor='#cccccc')
ax_inset.add_feature(cfeature.COASTLINE, linewidth=0.4, edgecolor='#999999')

draw_map(ax_inset, transform,
         countries_to_draw=EUROPE_COUNTRIES,
         show_labels=EUROPE_COUNTRIES,
         label_offsets=INSET_LABEL_OFFSETS,
         fontsize=7, node_scale=1.3)

for spine in ax_inset.spines.values():
    spine.set_edgecolor('#888888')
    spine.set_linewidth(1.2)

ax_inset.set_title("European genomic relatedness", fontsize=9,
                    fontweight='bold', pad=6, color='#333333')

# Legend
legend_elements = [
    Line2D([0], [0], color='#D6604D', linewidth=2.5, alpha=0.7,
           label=f'Both bla$_{{NDM}}$+ ({total_np} pairs, {total_np/max(total_np+total_nm, 1)*100:.0f}%)'),
    Line2D([0], [0], color='#4393C3', linewidth=1.5, alpha=0.5,
           label=f'≥1 bla$_{{NDM}}$− ({total_nm} pairs, {total_nm/max(total_np+total_nm, 1)*100:.0f}%)'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#D6604D',
           markersize=8, label='Node: >70% NDM+'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#FC8D59',
           markersize=8, label='Node: 30–70% NDM+'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#4393C3',
           markersize=8, label='Node: <30% NDM+'),
]

ax_main.legend(handles=legend_elements, loc='lower left', fontsize=8,
               frameon=True, fancybox=False, edgecolor='#cccccc',
               ncol=1, handlelength=2)

ax_main.set_title(
    f"ST773 cross-country genomic links (<{SNP_CUTOFF} SNPs) stratified by bla$_{{NDM}}$ status",
    fontsize=13, fontweight='bold', pad=12)

plt.savefig(f"{OUTDIR}/figure4B_final.pdf", dpi=300, bbox_inches='tight')
plt.savefig(f"{OUTDIR}/figure4B_final.png", dpi=300, bbox_inches='tight')
print(f"Saved to {OUTDIR}/figure4B_final.pdf and .png")
