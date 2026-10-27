#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 3C: ST773 phylogeny with resistome architecture
========================================================
Integrated figure showing:
  - Maximum-likelihood phylogeny of 628 ST773 genomes (left panel)
  - Annotation tracks (right of tree):
      Track 1: blaNDM carriage (red = NDM+, blue = NDM-)
      Track 2: ICE6660 variant (orange = A, purple = B, yellow = NDM+ unresolved)
      Track 3: Country (top 15 colored)
      Track 4: Year (gradient)
      Track 5: Integron cassette (sul1 + qnrVC1 + aadA11 co-occurrence)
      Track 6+: Key resistance genes heatmap

Key findings this figure supports:
  - ICE variant annotations from the existing figure source
  - blaNDM+ subpopulation concentrated in specific clades
  - Resistome shows clone-level coherence

Input files (update paths below if needed):
  - TREE_FILE:       ML tree (.treefile from IQ-TREE)
  - RESFINDER_FILE:  ABRicate resfinder summary TSV
  - METADATA_FILE:   合并结果.xlsx or ST773_merged_626_filled.tsv
  - ICE alignments already processed; groups hardcoded below

Output:
  - figure3C_tree.pdf
  - figure3C_tree.png

Usage:
  cd <repository>
  python scripts/plot_figure3C.py

Author: Grace Garcia / Claude (Anthropic)
"""

import os
from pathlib import Path
import sys
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from Bio import Phylo

# ============================================================================
# PATHS — update if your files are elsewhere
# ============================================================================
ROOT = Path(__file__).resolve().parents[1]
TREE_FILE = str(ROOT / "data/ST773/core.aln.treefile")
RESFINDER_FILE = str(ROOT / "data/ST773/resfinder_summary.tsv")
METADATA_FILE = str(ROOT / "data/ST773/metadata_626.tsv")
OUTPUT_PREFIX = str(ROOT / "results/figure3C_tree")

# ============================================================================
# ICE VARIANT GROUPS (11 A and 6 B IDs in existing plot annotation)
# ============================================================================
ICE_A_GCA = {
    'GCA_009664165.1',  # reference CP041945.1
    'GCA_018448985.1',  # CP075176.1
    'GCA_033221155.1',  # CP137493.1
    'GCA_033225445.1',  # CP137504.1
    'GCA_036418295.1',  # CP143908.1
    'GCA_041931775.1',  # CP142443.1
    'GCA_041931785.1',  # CP142446.1
    'GCA_041931795.1',  # CP142445.1
    'GCA_041933115.1',  # CP142449.1
    'GCA_050240275.1',  # CP163153.1
    'GCA_050908855.1',  # CP194329.1
}

ICE_B_GCA = {
    'GCA_013255565.1',  # CP053917.1
    'GCA_036872625.1',  # CP145538.1
    'GCA_037223945.1',  # CP147725.1
    'GCA_041931805.1',  # CP142447.1
    'GCA_041931815.1',  # CP142448.1
    'GCA_050239745.1',  # CP166022.1
}

# ============================================================================
# KEY RESISTANCE GENES TO VISUALIZE (order matters — left to right in heatmap)
# ============================================================================
KEY_GENES = [
    'blaNDM', 'blaVIM-2', 'qnrVC1', 'aadA11', 'sul1',
    'rmtB', 'tet(G)', 'CmlA9',
]

# ============================================================================
# HELPERS
# ============================================================================
def normalize_gca(name):
    """Extract GCA_XXXXXXXXX.X from any string containing it."""
    m = re.search(r'GCA_\d+\.\d+', str(name))
    return m.group(0) if m else str(name).strip()


def gca_short(name):
    """GCA_XXXXXXXXX.X → GCA_XXXXXXXXX (strip version)."""
    name = normalize_gca(name)
    return '.'.join(name.split('.')[:1]) if '.' in name else name


# ============================================================================
# 1. LOAD TREE
# ============================================================================
print("=" * 70)
print("FIGURE 3C: ST773 phylogeny with resistome architecture")
print("=" * 70)

if not os.path.exists(TREE_FILE):
    print(f"ERROR: Tree not found: {TREE_FILE}")
    sys.exit(1)

print(f"\nLoading tree: {TREE_FILE}")
tree = Phylo.read(TREE_FILE, "newick")
tips = [t for t in tree.get_terminals()]
tip_names = [t.name for t in tips]
n_tips = len(tips)
print(f"  {n_tips} tips")

# Ladderize for cleaner display
tree.ladderize(reverse=False)
tips = [t for t in tree.get_terminals()]
tip_names = [t.name for t in tips]

# Build GCA lookup for each tip
tip_to_gca = {name: normalize_gca(name) for name in tip_names}

# ============================================================================
# 2. LOAD RESISTOME DATA
# ============================================================================
print(f"\nLoading resistome: {RESFINDER_FILE}")
if not os.path.exists(RESFINDER_FILE):
    print(f"ERROR: Resfinder summary not found: {RESFINDER_FILE}")
    sys.exit(1)

res = pd.read_csv(RESFINDER_FILE, sep='\t')
print(f"  {res.shape[0]} rows × {res.shape[1]} cols")

# First column is filename; rest are genes
res_genes = list(res.columns[1:])
# Check how ABRicate encodes presence — often a "." for absent, "ID:COV" for present
# We'll define "present" as any non-"." value

# Build a per-genome dict: {gca -> {gene -> 0/1}}
filename_col = res.columns[0]
res_by_gca = {}
for _, row in res.iterrows():
    fn = str(row[filename_col])
    gca = normalize_gca(fn)
    gene_hits = {}
    for g in res_genes:
        val = str(row[g]).strip()
        gene_hits[g] = 0 if val in ('.', '', 'nan', 'NaN') else 1
    res_by_gca[gca] = gene_hits

# Helper: check if a gene is present, accepting partial gene-name matches
def gene_present(gca, gene_query):
    """Check if any gene matching gene_query is present in this genome."""
    hits = res_by_gca.get(gca, {})
    query_lower = gene_query.lower()
    for g, present in hits.items():
        if present and query_lower in g.lower():
            return 1
    return 0

# blaNDM status
def is_ndm_positive(gca):
    return gene_present(gca, 'blaNDM') or gene_present(gca, 'NDM')

# Integron cassette (sul1 + qnrVC1 + aadA11)
def has_cassette(gca):
    return all([
        gene_present(gca, 'sul1'),
        gene_present(gca, 'qnrVC1'),
        gene_present(gca, 'aadA11'),
    ])

# ============================================================================
# 3. LOAD METADATA (country, year)
# ============================================================================
print(f"\nLoading metadata: {METADATA_FILE}")
if not os.path.exists(METADATA_FILE):
    print(f"ERROR: Metadata not found: {METADATA_FILE}")
    sys.exit(1)

if METADATA_FILE.endswith('.xlsx'):
    meta_df = pd.read_excel(METADATA_FILE)
else:
    meta_df = pd.read_csv(METADATA_FILE, sep='\t')

# Find relevant columns (handle different naming)
gca_col = None
for c in meta_df.columns:
    if 'genome' in c.lower() or 'assembly' in c.lower() or 'name' in c.lower():
        gca_col = c
        break
if gca_col is None:
    gca_col = meta_df.columns[0]

country_col = None
for c in meta_df.columns:
    if 'country' in c.lower() or 'geo' in c.lower():
        country_col = c
        break

year_col = None
for c in meta_df.columns:
    if c.lower() == 'year' or 'date' in c.lower():
        year_col = c
        break

print(f"  Columns: GCA={gca_col!r}, Country={country_col!r}, Year={year_col!r}")

# Build lookup: GCA_short -> (country, year)
meta_lookup = {}
for _, row in meta_df.iterrows():
    gca = normalize_gca(row[gca_col])
    country = str(row[country_col]).strip() if country_col and pd.notna(row[country_col]) else 'Unknown'
    # Clean country (take before colon if present)
    if ':' in country:
        country = country.split(':')[0].strip()
    if country.lower() in ('nan', '', 'missing', 'not collected', 'restricted access', 'unknown', 'not applicable', 'n/a'):
        country = 'Unknown'
    
    year = None
    if year_col and pd.notna(row[year_col]):
        try:
            match_year = re.search(r"(?:19|20)\d{2}", str(row[year_col]))
            year = int(match_year.group(0)) if match_year else None
        except (ValueError, TypeError):
            year = None
    
    meta_lookup[gca] = {'country': country, 'year': year}

# ============================================================================
# 4. BUILD PER-TIP ANNOTATION TABLE
# ============================================================================
print("\nBuilding per-tip annotation table...")

ndm_arr = []
ice_arr = []
country_arr = []
year_arr = []
cassette_arr = []
resistome_mat = np.zeros((n_tips, len(KEY_GENES)), dtype=int)

for i, tip_name in enumerate(tip_names):
    gca = tip_to_gca[tip_name]
    gca_s = gca_short(gca)
    
    ndm = is_ndm_positive(gca)
    ndm_arr.append(ndm)
    
    # ICE classification: A, B, NDM+ unresolved, or NDM-
    if gca in ICE_A_GCA:
        ice = 'ICE-A'
    elif gca in ICE_B_GCA:
        ice = 'ICE-B'
    elif ndm:
        ice = 'NDM+ unresolved'
    else:
        ice = 'NDM-'
    ice_arr.append(ice)
    
    cassette_arr.append(has_cassette(gca))
    
    meta = meta_lookup.get(gca, {'country': 'Unknown', 'year': None})
    country_arr.append(meta['country'])
    year_arr.append(meta['year'])
    
    for j, gene in enumerate(KEY_GENES):
        resistome_mat[i, j] = gene_present(gca, gene)

print(f"  NDM+: {sum(ndm_arr)}/{n_tips}")
print(f"  ICE-A: {ice_arr.count('ICE-A')}, ICE-B: {ice_arr.count('ICE-B')}, "
      f"NDM+ unresolved: {ice_arr.count('NDM+ unresolved')}, NDM-: {ice_arr.count('NDM-')}")
print(f"  Integron cassette (sul1+qnrVC1+aadA11): {sum(cassette_arr)}/{n_tips}")

# Country top 15
from collections import Counter
country_counter = Counter(c for c in country_arr if c != 'Unknown')
top_countries = [c for c, _ in country_counter.most_common(15)]
print(f"  Top 15 countries: {top_countries}")

# Country color map (using tab20 palette)
country_color_map = {}
for i, c in enumerate(top_countries):
    country_color_map[c] = plt.cm.tab20(i / 15.0)
country_color_map['Unknown'] = '#EEEEEE'
country_color_map['Other'] = '#BBBBBB'

# Year range
years_valid = [y for y in year_arr if y is not None]
year_min = min(years_valid) if years_valid else 2005
year_max = max(years_valid) if years_valid else 2024
year_cmap = plt.cm.viridis
year_norm = plt.Normalize(vmin=year_min, vmax=year_max)

# ============================================================================
# 5. BUILD TREE COORDINATES
# ============================================================================
print("\nComputing tree layout...")

tip_y = {t.name: i for i, t in enumerate(tips)}
max_x = max(tree.distance(tree.root, t) for t in tips)

def get_y(clade):
    if clade.is_terminal():
        return tip_y[clade.name]
    return np.mean([get_y(c) for c in clade.clades])

def get_x(clade):
    return tree.distance(tree.root, clade)

# ============================================================================
# 6. PLOT
# ============================================================================
print("\nBuilding figure...")

fig = plt.figure(figsize=(14, max(12, n_tips * 0.025)))

# Layout: [tree] [NDM] [ICE] [country] [year] [cassette] [resistome heatmap]
# Widths proportional
widths = {
    'tree': 0.40,
    'ndm': 0.025,
    'ice': 0.025,
    'country': 0.025,
    'year': 0.025,
    'cassette': 0.025,
    'resistome': 0.035 * len(KEY_GENES),
    'gap': 0.004,
}

left = 0.06
bottom = 0.14
height = 0.80

pos = left
ax_tree = fig.add_axes([pos, bottom, widths['tree'], height])
pos += widths['tree'] + widths['gap']
ax_ndm = fig.add_axes([pos, bottom, widths['ndm'], height])
pos += widths['ndm'] + widths['gap']
ax_ice = fig.add_axes([pos, bottom, widths['ice'], height])
pos += widths['ice'] + widths['gap']
ax_country = fig.add_axes([pos, bottom, widths['country'], height])
pos += widths['country'] + widths['gap']
ax_year = fig.add_axes([pos, bottom, widths['year'], height])
pos += widths['year'] + widths['gap']
ax_cassette = fig.add_axes([pos, bottom, widths['cassette'], height])
pos += widths['cassette'] + widths['gap']
ax_resistome = fig.add_axes([pos, bottom, widths['resistome'], height])

# --- TREE ---
print("  Drawing tree...")
# Draw branches manually for speed
lines = []
for clade in tree.find_clades():
    if clade is tree.root:
        continue
    parent = None
    # Find parent
    for c in tree.find_clades():
        if clade in c.clades:
            parent = c
            break
    if parent is None:
        continue
    x_parent = get_x(parent)
    x_child = get_x(clade)
    y_child = get_y(clade)
    # Horizontal
    lines.append([(x_parent, y_child), (x_child, y_child)])

# Vertical connectors
for clade in tree.find_clades():
    if clade.is_terminal():
        continue
    if not clade.clades:
        continue
    x = get_x(clade)
    ys = [get_y(c) for c in clade.clades]
    lines.append([(x, min(ys)), (x, max(ys))])

from matplotlib.collections import LineCollection
lc = LineCollection(lines, colors='black', linewidths=0.3)
ax_tree.add_collection(lc)

ax_tree.set_xlim(-max_x * 0.02, max_x * 1.05)
ax_tree.set_ylim(-1, n_tips)
ax_tree.invert_yaxis()  # first tip at top
ax_tree.set_xlabel('Substitutions per site', fontsize=10)
ax_tree.set_yticks([])
ax_tree.spines['top'].set_visible(False)
ax_tree.spines['right'].set_visible(False)
ax_tree.spines['left'].set_visible(False)
ax_tree.set_title(f'ST773 maximum-likelihood phylogeny (n={n_tips})',
                  fontsize=11, loc='left', pad=8)

# --- TRACK 1: blaNDM ---
def draw_track(ax, values, color_map, title, cmap_for_none=None):
    """Draw a single-column track."""
    for i, v in enumerate(values):
        c = color_map.get(v, '#FFFFFF')
        ax.barh(i, 1, height=1, color=c, edgecolor='none')
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.5, n_tips - 0.5)
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(title, fontsize=8, rotation=90, x=0.5, y=1.005, ha='center', va='bottom')

ndm_colors = {1: '#C0392B', 0: '#D5E8F2'}
draw_track(ax_ndm, ndm_arr, ndm_colors, 'blaNDM')

# --- TRACK 2: ICE variant ---
ice_colors = {
    'ICE-A': '#E67E22',
    'ICE-B': '#8E44AD',
    'NDM+ unresolved': '#F7DC6F',
    'NDM-': '#ECF0F1',
}
draw_track(ax_ice, ice_arr, ice_colors, 'ICE variant')

# --- TRACK 3: Country ---
country_colors_per_tip = []
for c in country_arr:
    if c in country_color_map:
        country_colors_per_tip.append(country_color_map[c])
    elif c == 'Unknown':
        country_colors_per_tip.append(country_color_map['Unknown'])
    else:
        country_colors_per_tip.append(country_color_map['Other'])
for i, col in enumerate(country_colors_per_tip):
    ax_country.barh(i, 1, height=1, color=col, edgecolor='none')
ax_country.set_xlim(0, 1); ax_country.set_ylim(-0.5, n_tips - 0.5)
ax_country.invert_yaxis()
ax_country.set_xticks([]); ax_country.set_yticks([])
for s in ax_country.spines.values():
    s.set_visible(False)
ax_country.set_title('Country', fontsize=8, rotation=90, x=0.5, y=1.005, ha='center', va='bottom')

# --- TRACK 4: Year ---
for i, y in enumerate(year_arr):
    if y is None:
        col = '#EEEEEE'
    else:
        col = year_cmap(year_norm(y))
    ax_year.barh(i, 1, height=1, color=col, edgecolor='none')
ax_year.set_xlim(0, 1); ax_year.set_ylim(-0.5, n_tips - 0.5)
ax_year.invert_yaxis()
ax_year.set_xticks([]); ax_year.set_yticks([])
for s in ax_year.spines.values():
    s.set_visible(False)
ax_year.set_title('Year', fontsize=8, rotation=90, x=0.5, y=1.005, ha='center', va='bottom')

# --- TRACK 5: Cassette ---
cassette_colors = {1: '#2C3E50', 0: '#FFFFFF'}
draw_track(ax_cassette, cassette_arr, cassette_colors, 'MDR cassette')

# --- TRACK 6+: Resistome heatmap ---
gene_cmap = ListedColormap(['#FFFFFF', '#1F4E79'])
ax_resistome.imshow(resistome_mat, aspect='auto', cmap=gene_cmap,
                    interpolation='nearest')
ax_resistome.set_xticks(range(len(KEY_GENES)))
ax_resistome.set_xticklabels(KEY_GENES, rotation=90, fontsize=8)
ax_resistome.set_yticks([])
for s in ax_resistome.spines.values():
    s.set_visible(False)
ax_resistome.set_title('Resistance genes', fontsize=8, pad=8)

# ============================================================================
# 7. LEGEND PANEL (below the plot)
# ============================================================================
# Add legends as text boxes at the bottom
legend_y = 0.02
legend_items = []

# blaNDM legend
legend_items.append(mpatches.Patch(color='#C0392B', label='blaNDM+'))
legend_items.append(mpatches.Patch(color='#D5E8F2', label='blaNDM−'))

# ICE legend
legend_items.append(mpatches.Patch(color='#E67E22',
                                    label=f'ICE-A (n={ice_arr.count("ICE-A")})'))
legend_items.append(mpatches.Patch(color='#8E44AD',
                                    label=f'ICE-B (n={ice_arr.count("ICE-B")})'))
legend_items.append(mpatches.Patch(color='#F7DC6F',
                                    label=f'NDM+ unresolved (n={ice_arr.count("NDM+ unresolved")})'))
legend_items.append(mpatches.Patch(color='#ECF0F1',
                                    label=f'NDM− (n={ice_arr.count("NDM-")})'))

fig.legend(handles=legend_items, loc='lower left', ncol=6,
           bbox_to_anchor=(0.06, 0.005), fontsize=7, frameon=False)

# Country legend (far right, vertically centered — NOT overlapping heatmap)
country_legend_handles = [
    mpatches.Patch(color=country_color_map[c], label=f'{c} ({country_counter[c]})')
    for c in top_countries
]
fig.legend(handles=country_legend_handles, loc='lower right',
           bbox_to_anchor=(0.98, 0.005), fontsize=6, frameon=False,
           title='Top 15 countries', title_fontsize=7, ncol=5)

# Year colorbar (middle bottom, between the two legends)
cbar_ax = fig.add_axes([0.45, 0.075, 0.10, 0.008])
sm = plt.cm.ScalarMappable(cmap=year_cmap, norm=year_norm)
sm.set_array([])
cb = fig.colorbar(sm, cax=cbar_ax, orientation='horizontal')
cb.set_label('Year', fontsize=7)
cb.ax.tick_params(labelsize=6)

# ============================================================================
# 8. SAVE
# ============================================================================
Path(OUTPUT_PREFIX).parent.mkdir(parents=True, exist_ok=True)
plt.savefig(f'{OUTPUT_PREFIX}.pdf', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(f'{OUTPUT_PREFIX}.png', dpi=200, bbox_inches='tight', facecolor='white')

print(f"\n[Saved] {OUTPUT_PREFIX}.pdf")
print(f"[Saved] {OUTPUT_PREFIX}.png")

# ============================================================================
# 9. SUMMARY PRINTED
# ============================================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"Total tips: {n_tips}")
print(f"NDM+: {sum(ndm_arr)} ({100*sum(ndm_arr)/n_tips:.1f}%)")
print(f"MDR cassette (sul1+qnrVC1+aadA11): {sum(cassette_arr)} ({100*sum(cassette_arr)/n_tips:.1f}%)")
print(f"ICE-A characterized: {ice_arr.count('ICE-A')}")
print(f"ICE-B characterized: {ice_arr.count('ICE-B')}")
print(f"NDM+ unresolved (likely ICE but assembly fragmented): {ice_arr.count('NDM+ unresolved')}")
print(f"\nTop 5 countries:")
for c, n in country_counter.most_common(5):
    print(f"  {c}: {n}")
print(f"\nYear range: {year_min}-{year_max}")
print("\nDONE")