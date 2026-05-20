#!/usr/bin/env python3
"""
Publication-quality figures for Gene Target Identification Analysis
==================================================================

Generates 10 figures illustrating gene-level engineering targets for
lysine biofortification in maize endosperm, informed by o2 mutant
metabolic reprogramming.

All figures exclude non-GPR (transport/model-defined) reactions.

Author: Computational analysis pipeline
Date: March 2026
"""

import os
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch, Patch
from matplotlib.ticker import MaxNLocator
import matplotlib.patches as mpatches
from collections import defaultdict

# ===========================================================================
# Configuration
# ===========================================================================
BASE_DIR = '/mnt/nrdstor/ssbio/nampimah/ENDOSPERM/2026/Endosperm_Model_MARY'
DATA_DIR = os.path.join(BASE_DIR, 'gene_target_analysis')
ANALYSIS_DIR = os.path.join(BASE_DIR, 'analysis_results')
FIG_DIR = os.path.join(DATA_DIR, 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

# ---------- Color palette (colorblind-friendly) ----------
WT_COLOR   = '#2166ac'
O2_COLOR   = '#b2182b'
KO_COLOR   = '#d73027'
OE_COLOR   = '#4575b4'
MIM_COLOR  = '#74add1'
COMBO_COLOR = '#fdae61'
VERIFIED_COLOR   = '#1a9850'
UNVERIFIED_COLOR = '#bdbdbd'
HIGHLIGHT_COLOR  = '#fc8d59'

# ---------- Publication rcParams ----------
plt.rcParams.update({
    'font.size':        11,
    'font.family':      'sans-serif',
    'font.sans-serif':  ['DejaVu Sans', 'Arial', 'Helvetica'],
    'axes.titlesize':   13,
    'axes.titleweight': 'bold',
    'axes.labelsize':   11,
    'axes.labelweight': 'regular',
    'axes.linewidth':   0.8,
    'xtick.labelsize':  9,
    'ytick.labelsize':  9,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'legend.fontsize':  9,
    'legend.framealpha': 0.9,
    'figure.dpi':       150,
    'savefig.dpi':      300,
    'savefig.bbox':     'tight',
    'savefig.pad_inches': 0.1,
    'axes.spines.top':   False,
    'axes.spines.right': False,
    'lines.linewidth':  2,
    'lines.markersize': 7,
    'patch.linewidth':  0.6,
    'grid.linewidth':   0.4,
    'grid.alpha':       0.3,
})

DAPS = [6, 8, 10, 12, 15, 18, 22, 30]

# ===========================================================================
# GPR filter definitions
# ===========================================================================
NO_GPR_REACTIONS = {
    'Exe2[K]', 'Trans2[K]', 'cpTransport_C00007[K]',
    'MR00558[K,p]', 'R01308[K,p]', 'cpTransport_C00712[K]',
    'cpTransport_C00712[K,p]',
}
_NO_GPR_BASES = {'Exe2', 'Trans2', 'cpTransport_C00007', 'cpTransport_C00712',
                 'MR00558', 'R01308'}


def has_gpr(rxn_name):
    """Return True if the reaction has a GPR association."""
    name = rxn_name.strip()
    if name in NO_GPR_REACTIONS:
        return False
    base = name.split('[')[0]
    return base not in _NO_GPR_BASES


def combo_has_gpr(combo_str):
    """Return True if ALL reactions in a combination string have GPR."""
    for base in _NO_GPR_BASES:
        if base in combo_str:
            return False
    return True


# ===========================================================================
# Gene-ID lookup for key targets (from endosperm_GSM GPR mapping)
# ===========================================================================
GENE_LOOKUP = {
    'R00216[K,p]': ('Pyruvate DH',           ['GRMZM2G085019', '+4']),
    'R00282[K,c]': ('Acetyl-CoA carboxylase', ['GRMZM2G000236']),
    'R00883[K,c]': ('Phosphofructokinase',    ['GRMZM2G119300', 'GRMZM2G149265']),
    'R00512[K,c]': ('Pyruvate kinase',        ['GRMZM2G067453', '+3']),
    'R01015[K,c]': ('GAPDH',                  ['GRMZM2G002807', 'GRMZM2G017257']),
    'R00667[K,m]': ('TCA cycle',              ['GRMZM2G080828', 'GRMZM2G119583']),
    'R01280[K,p]': ('Chorismate synthase',    ['GRMZM2G002614', '+12']),
    'R00342[K,p]': ('Transketolase',          ['GRMZM2G035767', '+10']),
    'MR00422[K,p]':('Lipid metab.',           ['GRMZM2G002656', '+14']),
    'MR00732[K,m]':('Biotin metab.',          ['GRMZM2G377341', 'GRMZM5G858094']),
    'R02292[K,p]': ('DHDPS',                  ['GRMZM2G027835']),
    'R00451[K,p]': ('DAP decarboxylase',      ['GRMZM2G020446']),
}

COMPARTMENT_LABELS = {
    'K,p': 'plastid', 'K,c': 'cytosol', 'K,m': 'mito.',
    'K,x': 'perox.', 'K': '',
}


def format_rxn(rxn, with_compartment=True):
    """Pretty-print a reaction name: R00216[K,p] -> R00216 (plastid)."""
    name = rxn.strip()
    for code, label in COMPARTMENT_LABELS.items():
        suffix = '[' + code + ']'
        if name.endswith(suffix):
            base = name[:-len(suffix)]
            if with_compartment and label:
                return base + ' (' + label + ')'
            return base
    return name


def format_rxn_gene(rxn):
    """Reaction name with enzyme annotation on second line."""
    base = format_rxn(rxn)
    key = rxn.strip()
    if key in GENE_LOOKUP:
        return base + '\n' + GENE_LOOKUP[key][0]
    return base


# ===========================================================================
# Robust CSV loader (handles commas inside brackets like [K,p])
# ===========================================================================
def load_csv(filepath):
    """Load standard CSV with Python csv module (handles quoting)."""
    rows = []
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    return rows


def load_csv_robust(filepath):
    """Load CSV where reaction names contain unquoted commas (e.g. R00216[K,p]).

    Strategy: after splitting on comma, rejoin any field whose opening bracket
    '[' is not closed by ']' with the next field(s).
    """
    with open(filepath, 'r') as f:
        header_line = f.readline().strip()
        headers = header_line.split(',')
        n_cols = len(headers)

        rows = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(',')
            # Rejoin fields split inside brackets
            rejoined = []
            i = 0
            while i < len(parts):
                field = parts[i]
                while field.count('[') > field.count(']') and i + 1 < len(parts):
                    i += 1
                    field = field + ',' + parts[i]
                rejoined.append(field)
                i += 1
            if len(rejoined) == n_cols:
                row = dict(zip(headers, rejoined))
                rows.append(row)
    return rows


# ===========================================================================
# Helper: save figure in both PNG and PDF
# ===========================================================================
def savefig(fig, name):
    fig.savefig(os.path.join(FIG_DIR, name + '.png'))
    fig.savefig(os.path.join(FIG_DIR, name + '.pdf'))
    plt.close(fig)


def screen_color(st):
    """Return color for a Screen_Types string."""
    if 'KO' in st and ('OE' in st or 'MIM' in st):
        return COMBO_COLOR
    if 'MIM' in st and 'OE' in st:
        return COMBO_COLOR
    if 'KO' in st:
        return KO_COLOR
    if 'OE' in st:
        return OE_COLOR
    if 'MIM' in st:
        return MIM_COLOR
    return '#999999'


# ===========================================================================
# Load all datasets
# ===========================================================================
print("Loading data...")

baseline  = load_csv(os.path.join(DATA_DIR, '01_baseline_fba.csv'))

ranking_raw = load_csv_robust(os.path.join(DATA_DIR, '07_master_target_ranking.csv'))
ranking = [r for r in ranking_raw if has_gpr(r['Reaction'])]

overexp   = load_csv_robust(os.path.join(DATA_DIR, '03_overexpression_screen.csv'))
knockout  = load_csv_robust(os.path.join(DATA_DIR, '04_knockout_screen.csv'))
mimicry   = load_csv_robust(os.path.join(DATA_DIR, '05_o2_bound_mimicry.csv'))
combos    = load_csv(os.path.join(DATA_DIR, '06_combination_targets.csv'))
verification = load_csv_robust(os.path.join(DATA_DIR, '08_target_verification.csv'))
lys_pathway  = load_csv_robust(os.path.join(ANALYSIS_DIR, '02_lysine_pathway_fluxes.csv'))
biomass_comp = load_csv(os.path.join(ANALYSIS_DIR, '01_biomass_comparison.csv'))
pathway_summary = load_csv(os.path.join(ANALYSIS_DIR, '07_pathway_flux_summary.csv'))
diff_flux = load_csv_robust(os.path.join(ANALYSIS_DIR, '06_differential_flux_all.csv'))

print("Data loaded successfully.\n")
# Quick sanity check
print("  GPR-filtered ranking: %d targets" % len(ranking))
if ranking:
    print("  Top target: %s  (Score=%.2f, FC=%.1f)" % (
        ranking[0]['Reaction'],
        float(ranking[0]['Composite_Score']),
        float(ranking[0]['Max_Fold_Change'])))


# ===================================================================
# FIGURE 1: O2 vs WT Lysine Biosynthesis Context
# ===================================================================
print("\nGenerating Figure 1: O2 vs WT Lysine Context...")

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# --- Panel A: R00451 flux across DAPs ---
ax = axes[0, 0]
wt_r00451, o2_r00451 = [], []
for row in lys_pathway:
    if row.get('Reaction', '').strip() == 'R00451[K,p]':
        for dap in DAPS:
            wt_r00451.append(float(row.get('WT_DAP%d' % dap, 0)))
            o2_r00451.append(float(row.get('O2_DAP%d' % dap, 0)))
        break
if not wt_r00451:
    for row in lys_pathway:
        if 'R00451' in row.get('Reaction', ''):
            for dap in DAPS:
                wt_r00451.append(float(row.get('WT_DAP%d' % dap, 0)))
                o2_r00451.append(float(row.get('O2_DAP%d' % dap, 0)))
            break
if not wt_r00451:
    print("  WARNING: R00451 not found, using hardcoded fallback")
    wt_r00451 = [3e-6, 2e-6, 3e-6, 3e-6, 7e-6, 5e-6, 2e-6, 0]
    o2_r00451 = [4e-6, 4e-6, 3e-6, 3e-6, 5.8e-3, 7e-6, 1.25, 0]

ax.plot(DAPS, wt_r00451, 'o-', color=WT_COLOR, label='Wild Type', zorder=3)
ax.plot(DAPS, o2_r00451, 's-', color=O2_COLOR, label='o2 Mutant', zorder=3)
ax.set_xlabel('Days After Pollination (DAP)')
ax.set_ylabel('R00451 Flux (DAP Decarboxylase)')
ax.set_title('A. Lysine Biosynthesis Flux (R00451)')
ax.legend(frameon=True, fancybox=True)
ax.set_xticks(DAPS)
ax.grid(True, axis='both')
if len(o2_r00451) > 6 and o2_r00451[6] > 0 and wt_r00451[6] > 0:
    fc_22 = o2_r00451[6] / wt_r00451[6]
    ax.annotate('{:,.0f}x at DAP 22'.format(fc_22),
                xy=(22, o2_r00451[6]),
                xytext=(24, o2_r00451[6] * 0.8),
                fontsize=9, color=O2_COLOR, fontweight='bold',
                arrowprops=dict(arrowstyle='->', color=O2_COLOR, lw=1.5))

# --- Panel B: Fold-change bar chart ---
ax = axes[0, 1]
fc_vals = []
for i, dap in enumerate(DAPS):
    if wt_r00451[i] > 1e-15:
        fc_vals.append(o2_r00451[i] / wt_r00451[i])
    else:
        fc_vals.append(1.0)
log_fc = [np.log10(max(f, 0.01)) for f in fc_vals]
colors_fc = [O2_COLOR if f > 1 else WT_COLOR for f in fc_vals]
ax.bar(range(len(DAPS)), log_fc, color=colors_fc, alpha=0.85,
       edgecolor='black', linewidth=0.5)
ax.axhline(y=0, color='black', ls='--', lw=0.8, alpha=0.5)
ax.set_xticks(range(len(DAPS)))
ax.set_xticklabels(['DAP %d' % d for d in DAPS], rotation=45, ha='right')
ax.set_ylabel('log10(Fold Change, o2/WT)')
ax.set_title('B. Lysine Flux Fold Change')
ax.grid(True, axis='y')
if len(log_fc) > 6:
    ax.annotate('FC = {:,.0f}'.format(fc_vals[6]),
                xy=(6, log_fc[6]), xytext=(4.5, log_fc[6] - 0.5),
                fontsize=9, fontweight='bold', color=O2_COLOR,
                arrowprops=dict(arrowstyle='->', color=O2_COLOR, lw=1.5))

# --- Panel C: Biomass comparison ---
ax = axes[1, 0]
wt_bm = [float(r['WT_Biomass']) for r in biomass_comp]
o2_bm = [float(r['O2_Biomass']) for r in biomass_comp]
ax.plot(DAPS, wt_bm, 'o-', color=WT_COLOR, label='Wild Type')
ax.plot(DAPS, o2_bm, 's-', color=O2_COLOR, label='o2 Mutant')
ax.set_xlabel('Days After Pollination (DAP)')
ax.set_ylabel('Seed Biomass Flux')
ax.set_title('C. Seed Biomass Across Development')
ax.legend(frameon=True, fancybox=True)
ax.set_xticks(DAPS)
ax.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))
ax.grid(True)

# --- Panel D: Pathway flux at DAP 22 ---
ax = axes[1, 1]
pathways_of_interest = ['Lysine_Biosynthesis', 'Aspartate_Pathway', 'TCA_Cycle', 'Glycolysis']
pathway_labels = ['Lysine\nBiosynthesis', 'Aspartate\nPathway', 'TCA\nCycle', 'Glycolysis']
wt_pf, o2_pf = [], []
for pw in pathways_of_interest:
    for row in pathway_summary:
        if row['Pathway'] == pw and int(row['DAP']) == 22:
            wt_pf.append(float(row['WT_TotalAbsFlux']))
            o2_pf.append(float(row['O2_TotalAbsFlux']))
            break

x = np.arange(len(pathway_labels))
w = 0.35
ax.bar(x - w / 2, wt_pf, w, label='Wild Type', color=WT_COLOR, alpha=0.85,
       edgecolor='black', linewidth=0.5)
ax.bar(x + w / 2, o2_pf, w, label='o2 Mutant', color=O2_COLOR, alpha=0.85,
       edgecolor='black', linewidth=0.5)
ax.set_xticks(x)
ax.set_xticklabels(pathway_labels, fontsize=9)
ax.set_ylabel('Total Absolute Flux')
ax.set_title('D. Pathway Activity at DAP 22')
ax.legend(frameon=True, fancybox=True)
ax.set_yscale('log')
ax.grid(True, axis='y')

fig.tight_layout(h_pad=3, w_pad=3)
savefig(fig, 'Fig1_O2_vs_WT_Lysine_Context')
print("  Saved Fig1_O2_vs_WT_Lysine_Context")


# ===================================================================
# FIGURE 2: Master Target Ranking - Top 25 GPR targets
# ===================================================================
print("Generating Figure 2: Master Target Ranking...")

top_n = min(25, len(ranking))

rxns         = [r['Reaction'] for r in ranking[:top_n]]
scores       = [float(r['Composite_Score']) for r in ranking[:top_n]]
fc_vals_t    = [float(r['Max_Fold_Change']) for r in ranking[:top_n]]
screen_types = [r['Screen_Types'] for r in ranking[:top_n]]
bar_colors   = [screen_color(s) for s in screen_types]

# Pretty labels: reaction (compartment) - enzyme name
display_labels = []
for rxn in rxns:
    key = rxn.strip()
    base = format_rxn(rxn, with_compartment=True)
    if key in GENE_LOOKUP:
        enzyme = GENE_LOOKUP[key][0]
        display_labels.append('%s  [%s]' % (base, enzyme))
    else:
        display_labels.append(base)

fig = plt.figure(figsize=(18, 10))
gs = gridspec.GridSpec(1, 2, width_ratios=[1.1, 1], wspace=0.45)

y_pos = np.arange(top_n)

# --- Panel A: Composite Score ---
ax = fig.add_subplot(gs[0])
ax.barh(y_pos, scores, color=bar_colors, alpha=0.85,
        edgecolor='black', linewidth=0.5)
ax.set_yticks(y_pos)
ax.set_yticklabels(display_labels, fontsize=8.5)
ax.invert_yaxis()
ax.set_xlabel('Composite Score (weighted)')
ax.set_title('A. Top %d Gene-Associated Targets' % top_n, pad=10)
ax.grid(True, axis='x')

# Value annotations
for i in range(top_n):
    ax.text(scores[i] + 0.02, i, '%.2f' % scores[i],
            va='center', fontsize=7.5, color='#333333')

legend_elements = [
    Patch(facecolor=KO_COLOR, edgecolor='black', label='Knockout'),
    Patch(facecolor=OE_COLOR, edgecolor='black', label='Overexpression'),
    Patch(facecolor=MIM_COLOR, edgecolor='black', label='o2 Bound Mimicry'),
    Patch(facecolor=COMBO_COLOR, edgecolor='black', label='Multi-screen hit'),
]
ax.legend(handles=legend_elements, loc='lower right', fontsize=8, frameon=True)

# --- Panel B: Fold Change (log10) ---
ax2 = fig.add_subplot(gs[1])
log_fc_t = [np.log10(max(f, 1.0)) for f in fc_vals_t]
ax2.barh(y_pos, log_fc_t, color=bar_colors, alpha=0.85,
         edgecolor='black', linewidth=0.5)
ax2.set_yticks(y_pos)
ax2.set_yticklabels(display_labels, fontsize=8.5)
ax2.invert_yaxis()
ax2.set_xlabel('log10(Max Fold Change in Lysine Flux)')
ax2.set_title('B. Lysine Flux Enhancement', pad=10)
ax2.grid(True, axis='x')

# FC value annotations
for i in range(top_n):
    fc = fc_vals_t[i]
    if fc > 2:
        ax2.text(log_fc_t[i] + 0.05, i, '{:,.0f}x'.format(fc),
                 va='center', fontsize=8, fontweight='bold', color=KO_COLOR)
    elif fc > 1.05:
        ax2.text(log_fc_t[i] + 0.005, i, '%.2fx' % fc,
                 va='center', fontsize=7, color='#555555')

fig.suptitle('Master Target Ranking  -  GPR-Associated Reactions Only',
             fontsize=14, fontweight='bold', y=0.98)
fig.subplots_adjust(left=0.18, right=0.98, top=0.92, bottom=0.06, wspace=0.45)
savefig(fig, 'Fig2_Master_Target_Ranking')
print("  Saved Fig2_Master_Target_Ranking")


# ===================================================================
# FIGURE 3: Screening Results - KO, OE, Mimicry
# ===================================================================
print("Generating Figure 3: Screening Results Overview...")

fig, axes = plt.subplots(1, 3, figsize=(18, 7))

# --- Panel A: Knockout targets ---
ax = axes[0]
ko_data = [(r['Reaction'], float(r['Fold_Change']), float(r['Biomass_Pct_of_WT']))
           for r in knockout if float(r['Fold_Change']) > 1.0 and has_gpr(r['Reaction'])]
ko_data.sort(key=lambda x: x[1], reverse=True)
ko_top = ko_data[:15]

if ko_top:
    ko_labels = [format_rxn(k[0]) for k in ko_top]
    ko_fc = [k[1] for k in ko_top]
    ko_bm = [k[2] for k in ko_top]
    y = np.arange(len(ko_top))
    ax.barh(y, [np.log10(max(f, 1.001)) for f in ko_fc],
            color=KO_COLOR, alpha=0.85, edgecolor='black', linewidth=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(ko_labels, fontsize=8.5)
    ax.invert_yaxis()
    for i, (rxn, fc, bm) in enumerate(ko_top[:5]):
        ax.text(np.log10(max(fc, 1.001)) + 0.05, i,
                '{:,.0f}x (BM:{:.0f}%)'.format(fc, bm),
                va='center', fontsize=7.5, fontweight='bold')
ax.set_xlabel('log10(Fold Change)')
ax.set_title('A. Knockout Targets')
ax.grid(True, axis='x')

# --- Panel B: Overexpression targets ---
ax = axes[1]
oe_data = [(r['Reaction'], float(r['Fold_Change']), float(r['New_Biomass']))
           for r in overexp if float(r['Fold_Change']) > 1.0 and has_gpr(r['Reaction'])]
oe_data.sort(key=lambda x: x[1], reverse=True)
oe_top = oe_data[:15]

if oe_top:
    oe_labels = [format_rxn(k[0]) for k in oe_top]
    oe_fc = [k[1] for k in oe_top]
    y = np.arange(len(oe_top))
    ax.barh(y, [np.log10(max(f, 1.001)) for f in oe_fc],
            color=OE_COLOR, alpha=0.85, edgecolor='black', linewidth=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(oe_labels, fontsize=8.5)
    ax.invert_yaxis()
ax.set_xlabel('log10(Fold Change)')
ax.set_title('B. Overexpression Targets (2x Vmax)')
ax.grid(True, axis='x')

# --- Panel C: o2-bound mimicry targets ---
ax = axes[2]
mim_data = [(r['Reaction'], float(r['Fold_Change']), float(r['Vmax_Change_Pct']))
            for r in mimicry if float(r['Fold_Change']) > 1.0 and has_gpr(r['Reaction'])]
mim_data.sort(key=lambda x: x[1], reverse=True)
mim_top = mim_data[:10] if mim_data else []

if mim_top:
    mim_labels = [format_rxn(k[0]) for k in mim_top]
    mim_fc = [k[1] for k in mim_top]
    mim_vpct = [k[2] for k in mim_top]
    y = np.arange(len(mim_top))
    ax.barh(y, [np.log10(max(f, 1.001)) for f in mim_fc],
            color=MIM_COLOR, alpha=0.85, edgecolor='black', linewidth=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(mim_labels, fontsize=8.5)
    ax.invert_yaxis()
    for i, (rxn, fc, vpct) in enumerate(mim_top):
        arrow = '+' if vpct > 0 else '-'
        ax.text(np.log10(max(fc, 1.001)) + 0.002, i,
                'Vmax %s%.1f%%' % (arrow, abs(vpct)),
                va='center', fontsize=7.5)
ax.set_xlabel('log10(Fold Change)')
ax.set_title('C. o2-Bound Mimicry Targets')
ax.grid(True, axis='x')

fig.tight_layout(w_pad=3)
savefig(fig, 'Fig3_Screening_Results')
print("  Saved Fig3_Screening_Results")


# ===================================================================
# FIGURE 4: Combination Targets - Synergy Analysis
# ===================================================================
print("Generating Figure 4: Combination Targets...")

fig, axes = plt.subplots(1, 2, figsize=(16, 8))

# --- Panel A: Top combination strategies ---
ax = axes[0]
combo_data = [(r['Combination'], float(r['Fold_Change']), float(r['Biomass_Pct']))
              for r in combos if float(r['Biomass_Pct']) > 0
              and combo_has_gpr(r['Combination'])]
combo_data.sort(key=lambda x: x[1], reverse=True)
combo_top = combo_data[:15]

combo_names = []
for c in combo_top:
    name = c[0]
    name = name.replace('[K,p]', '(p)').replace('[K,c]', '(c)')
    name = name.replace('[K,m]', '(m)').replace('[K]', '')
    if len(name) > 42:
        name = name[:39] + '...'
    combo_names.append(name)

combo_fc = [c[1] for c in combo_top]
combo_bm = [c[2] for c in combo_top]
y = np.arange(len(combo_top))
bar_colors_c = [VERIFIED_COLOR if bm > 50 else KO_COLOR for bm in combo_bm]

ax.barh(y, [np.log10(max(f, 1.001)) for f in combo_fc],
        color=bar_colors_c, alpha=0.85, edgecolor='black', linewidth=0.5)
ax.set_yticks(y)
ax.set_yticklabels(combo_names, fontsize=7.5)
ax.invert_yaxis()
ax.set_xlabel('log10(Fold Change in Lysine Flux)')
ax.set_title('A. Top Combination Strategies')
ax.grid(True, axis='x')
ax.legend(handles=[
    Patch(facecolor=VERIFIED_COLOR, edgecolor='black', label='Viable (>50% biomass)'),
    Patch(facecolor=KO_COLOR, edgecolor='black', label='Non-viable'),
], loc='lower right', fontsize=8)

# --- Panel B: Scatter - Biomass vs Fold Change ---
ax = axes[1]
all_fc = [float(r['Fold_Change']) for r in combos
          if float(r['Biomass_Pct']) > 0 and combo_has_gpr(r['Combination'])]
all_bm = [float(r['Biomass_Pct']) for r in combos
          if float(r['Biomass_Pct']) > 0 and combo_has_gpr(r['Combination'])]

viable    = [(fc, bm) for fc, bm in zip(all_fc, all_bm) if bm > 50]
nonviable = [(fc, bm) for fc, bm in zip(all_fc, all_bm) if bm <= 50]

if viable:
    ax.scatter([np.log10(max(fc, 1.001)) for fc, bm in viable],
               [bm for _, bm in viable],
               c=VERIFIED_COLOR, s=55, alpha=0.7, edgecolors='black',
               linewidth=0.4, label='Viable', zorder=3)
if nonviable:
    ax.scatter([np.log10(max(fc, 1.001)) for fc, bm in nonviable],
               [bm for _, bm in nonviable],
               c=KO_COLOR, s=55, alpha=0.7, edgecolors='black',
               linewidth=0.4, label='Non-viable', zorder=3)

ax.axhline(y=50, color='grey', ls='--', lw=0.8, alpha=0.6)
ax.text(0.1, 52, 'Viability threshold (50%)', fontsize=8, color='grey')
ax.set_xlabel('log10(Lysine Flux Fold Change)')
ax.set_ylabel('Biomass Retained (%)')
ax.set_title('B. Lysine Enhancement vs Biomass Viability')
ax.legend(frameon=True, fancybox=True)
ax.grid(True)

fig.tight_layout(w_pad=3)
savefig(fig, 'Fig4_Combination_Targets')
print("  Saved Fig4_Combination_Targets")


# ===================================================================
# FIGURE 5: Verification Summary
# ===================================================================
print("Generating Figure 5: Verification Summary...")

fig, axes = plt.subplots(1, 2, figsize=(14, 7))

verified, unverified = [], []
for r in verification:
    rxn = r['Reaction']
    if not has_gpr(rxn):
        continue
    fc = float(r['Fold_Change'])
    is_v = r['Verified'] == 'True'
    bm = float(r['Biomass_Retained_Pct'])
    (verified if is_v else unverified).append((rxn, fc, bm))

# --- Panel A: Verified targets ---
ax = axes[0]
if verified:
    v_labels = [format_rxn(v[0]) for v in verified]
    v_fc = [v[1] for v in verified]
    y = np.arange(len(verified))
    ax.barh(y, v_fc, color=VERIFIED_COLOR, alpha=0.85,
            edgecolor='black', linewidth=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(v_labels, fontsize=9)
    ax.invert_yaxis()
    for i, (rxn, fc, bm) in enumerate(verified):
        ax.text(fc + 0.001, i, 'BM: %.0f%%' % bm, va='center', fontsize=8)
ax.set_xlabel('Fold Change in Lysine Flux (R00451)')
ax.set_title('A. Verified Gene-Associated Targets', color=VERIFIED_COLOR)
ax.grid(True, axis='x')

# --- Panel B: Pie chart ---
ax = axes[1]
n_v, n_u = len(verified), len(unverified)
total = n_v + n_u
if total > 0:
    wedges, texts, autotexts = ax.pie(
        [n_v, n_u], explode=(0.04, 0),
        labels=['Verified\n(%d/%d)' % (n_v, total),
                'Unconfirmed\n(%d/%d)' % (n_u, total)],
        colors=[VERIFIED_COLOR, UNVERIFIED_COLOR],
        autopct='%1.1f%%', shadow=False, startangle=90,
        textprops={'fontsize': 11})
    for at in autotexts:
        at.set_fontweight('bold')
ax.set_title('B. Target Verification Rate')

fig.tight_layout(w_pad=3)
savefig(fig, 'Fig5_Verification_Summary')
print("  Saved Fig5_Verification_Summary")


# ===================================================================
# FIGURE 6: DAP 22 Metabolic Reprogramming
# ===================================================================
print("Generating Figure 6: DAP 22 Metabolic Reprogramming...")

dap22_diffs = []
for row in diff_flux:
    rxn = row['Reaction']
    wt_22 = float(row.get('WT_DAP22', 0))
    o2_22 = float(row.get('O2_DAP22', 0))
    diff_22 = float(row.get('Diff_DAP22', 0))
    if abs(diff_22) > 0.01:
        dap22_diffs.append((rxn, wt_22, o2_22, diff_22))
dap22_diffs.sort(key=lambda x: abs(x[3]), reverse=True)

fig, axes = plt.subplots(1, 2, figsize=(16, 9))

# --- Panel A: Upregulated ---
ax = axes[0]
upregulated = [(r, w, o, d) for r, w, o, d in dap22_diffs if d > 0.01][:20]
if upregulated:
    up_labels = [format_rxn(r[0])[:32] for r in upregulated]
    up_diffs = [r[3] for r in upregulated]
    y = np.arange(len(upregulated))
    ax.barh(y, up_diffs, color=O2_COLOR, alpha=0.85,
            edgecolor='black', linewidth=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(up_labels, fontsize=8)
    ax.invert_yaxis()
ax.set_xlabel('Flux Difference (o2 - WT)')
ax.set_title('A. Upregulated in o2 at DAP 22', color=O2_COLOR)
ax.grid(True, axis='x')

# --- Panel B: Downregulated ---
ax = axes[1]
downregulated = [(r, w, o, d) for r, w, o, d in dap22_diffs if d < -0.01]
downregulated.sort(key=lambda x: x[3])
downregulated = downregulated[:20]
if downregulated:
    down_labels = [format_rxn(r[0])[:32] for r in downregulated]
    down_diffs = [abs(r[3]) for r in downregulated]
    y = np.arange(len(downregulated))
    ax.barh(y, down_diffs, color=WT_COLOR, alpha=0.85,
            edgecolor='black', linewidth=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(down_labels, fontsize=8)
    ax.invert_yaxis()
ax.set_xlabel('|Flux Difference| (WT - o2)')
ax.set_title('B. Downregulated in o2 at DAP 22', color=WT_COLOR)
ax.grid(True, axis='x')

fig.tight_layout(w_pad=3)
savefig(fig, 'Fig6_DAP22_Metabolic_Reprogramming')
print("  Saved Fig6_DAP22_Metabolic_Reprogramming")


# ===================================================================
# FIGURE 7: Lysine Pathway Heatmap
# ===================================================================
print("Generating Figure 7: Lysine Pathway Heatmap...")

fig, axes = plt.subplots(1, 2, figsize=(16, 8))

lys_rxns_ordered = [
    'R00480[K,p]', 'R00355[K,p]', 'R02291[K,p]', 'R02292[K,p]',
    'R04198[K,p]', 'R02735[K,p]', 'R00451[K,p]',
    'R00715[K,c]', 'cpTransport_C00047[K]', 'ExB_C00047[K]',
    'Exchange_C00047[K,L]', 'Seed_Biomass[K]',
]
lys_labels = [
    'Aspartate kinase\n(R00480)', 'Asp aminotransferase\n(R00355)',
    'Asp-semialdehyde DH\n(R02291)', 'DHDPS\n(R02292)',
    'DHDPR\n(R04198)', 'meso-DAP production\n(R02735)',
    'DAP decarboxylase\n(R00451)',
    'LKR\n(R00715)', 'Lys transport\n(cpTransport)', 'Lys export\n(ExB)',
    'Lys exchange', 'Seed Biomass',
]

wt_matrix = np.zeros((len(lys_rxns_ordered), len(DAPS)))
o2_matrix = np.zeros((len(lys_rxns_ordered), len(DAPS)))

for row in lys_pathway:
    rxn = row.get('Reaction', '').strip()
    if rxn in lys_rxns_ordered:
        idx = lys_rxns_ordered.index(rxn)
        for j, dap in enumerate(DAPS):
            wt_matrix[idx, j] = float(row.get('WT_DAP%d' % dap, 0))
            o2_matrix[idx, j] = float(row.get('O2_DAP%d' % dap, 0))

fc_matrix = np.zeros_like(wt_matrix)
for i in range(fc_matrix.shape[0]):
    for j in range(fc_matrix.shape[1]):
        wv = abs(wt_matrix[i, j])
        ov = abs(o2_matrix[i, j])
        if wv > 1e-15:
            fc_matrix[i, j] = np.log10(max(ov / wv, 0.001))
        elif ov > 1e-15:
            fc_matrix[i, j] = 6
        else:
            fc_matrix[i, j] = 0

# --- Panel A ---
ax = axes[0]
im1 = ax.imshow(np.abs(wt_matrix), aspect='auto', cmap='Blues',
                interpolation='nearest')
ax.set_xticks(range(len(DAPS)))
ax.set_xticklabels(['DAP %d' % d for d in DAPS], fontsize=9,
                    rotation=45, ha='right')
ax.set_yticks(range(len(lys_labels)))
ax.set_yticklabels(lys_labels, fontsize=8)
ax.set_title('A. WT Lysine Pathway Flux (|flux|)')
plt.colorbar(im1, ax=ax, shrink=0.8, label='|Flux|')

# --- Panel B ---
ax = axes[1]
im2 = ax.imshow(fc_matrix, aspect='auto', cmap='RdBu_r',
                interpolation='nearest', vmin=-3, vmax=6)
ax.set_xticks(range(len(DAPS)))
ax.set_xticklabels(['DAP %d' % d for d in DAPS], fontsize=9,
                    rotation=45, ha='right')
ax.set_yticks(range(len(lys_labels)))
ax.set_yticklabels(lys_labels, fontsize=8)
ax.set_title('B. Fold Change (log10, o2/WT)')
plt.colorbar(im2, ax=ax, shrink=0.8, label='log10(o2/WT)')

for i in range(fc_matrix.shape[0]):
    for j in range(fc_matrix.shape[1]):
        val = fc_matrix[i, j]
        if abs(val) > 2:
            color = 'white' if abs(val) > 3 else 'black'
            ax.text(j, i, '%.1f' % val, ha='center', va='center',
                    fontsize=6, color=color, fontweight='bold')

fig.tight_layout(w_pad=3)
savefig(fig, 'Fig7_Lysine_Pathway_Heatmap')
print("  Saved Fig7_Lysine_Pathway_Heatmap")


# ===================================================================
# FIGURE 8: Gene-Level Engineering Strategy Diagram
# ===================================================================
print("Generating Figure 8: Gene-Level Engineering Strategy Diagram...")

fig, ax = plt.subplots(figsize=(16, 12))
ax.set_xlim(0, 16)
ax.set_ylim(0, 12)
ax.axis('off')
ax.set_title('Gene-Level Engineering Strategy for Lysine Biofortification\n'
             'in Wild-Type Maize (GPR-Associated Targets Only)',
             fontsize=14, fontweight='bold', pad=20)

# --- Tier 1 ---
tier1 = FancyBboxPatch((0.5, 8.5), 15, 3, boxstyle="round,pad=0.3",
                        facecolor='#fee0d2', edgecolor=KO_COLOR, linewidth=2)
ax.add_patch(tier1)
ax.text(8, 11.1, 'TIER 1: HIGH-IMPACT KNOCKOUT TARGETS', fontsize=12,
        fontweight='bold', ha='center', color=KO_COLOR)
ax.text(8, 10.55,
        '(Block fatty acid biosynthesis entry  -  100% biomass retained)',
        fontsize=9, ha='center', color='#666666')

tier1_items = [
    ('R00216[K,p]', 'Pyruvate Dehydrogenase (plastid)', '16,234x',
     'GRMZM2G085019 +4 isozymes', '5 genes  -  CRISPR multiplex KO'),
    ('R00282[K,c]', 'Acetyl-CoA Carboxylase (cytosol)', '28x',
     'GRMZM2G000236', '1 gene  -  single CRISPR KO'),
]
for i, (rxn, annot, fc, genes, note) in enumerate(tier1_items):
    yp = 10.0 - i * 0.7
    ax.text(1.0, yp, '*', fontsize=16, color=KO_COLOR, va='center',
            fontweight='bold')
    ax.text(1.5, yp, rxn, fontsize=10, fontweight='bold', va='center')
    ax.text(4.2, yp, annot, fontsize=9, va='center', color='#444444')
    ax.text(8.8, yp, 'FC: ' + fc, fontsize=10, fontweight='bold',
            va='center', color=KO_COLOR)
    ax.text(10.5, yp, genes, fontsize=8, va='center', color='#1a5276',
            style='italic')
    ax.text(1.5, yp - 0.3, '-> ' + note, fontsize=7.5, va='center',
            color='#666666')

# --- Tier 2 ---
tier2 = FancyBboxPatch((0.5, 5.0), 15, 3, boxstyle="round,pad=0.3",
                        facecolor='#deebf7', edgecolor=OE_COLOR, linewidth=2)
ax.add_patch(tier2)
ax.text(8, 7.6, 'TIER 2: OVEREXPRESSION & MIMICRY TARGETS', fontsize=12,
        fontweight='bold', ha='center', color=OE_COLOR)
ax.text(8, 7.1, '(~10% lysine increase  -  enhance carbon supply & capture)',
        fontsize=9, ha='center', color='#666666')

tier2_items = [
    ('R01280[K,p]', 'Chorismate Synthase', 'OE', '1.10x',
     'GRMZM2G002614 (+12 isozymes)'),
    ('R00883[K,c]', 'Phosphofructokinase', 'OE', '1.10x',
     'GRMZM2G119300 or GRMZM2G149265'),
    ('R00512[K,c]', 'Pyruvate Kinase', 'OE', '1.10x',
     'GRMZM2G067453 (+3 isozymes)'),
    ('R00667[K,m]', 'TCA Cycle (mito.)', 'OE', '1.10x',
     'GRMZM2G080828 or GRMZM2G119583'),
    ('R00342[K,p]', 'Transketolase (plastid)', 'MIM', '1.10x',
     'GRMZM2G035767 (+10 isozymes)'),
    ('MR00422[K,p]', 'Lipid Metab. (plastid)', 'MIM', '1.10x',
     'GRMZM2G002656 (+14 isozymes)'),
]
for i, (rxn, annot, strategy, fc, genes) in enumerate(tier2_items):
    yp = 6.6 - i * 0.35
    clr = OE_COLOR if strategy == 'OE' else MIM_COLOR
    ax.text(1.0, yp, '*', fontsize=14, color=clr, va='center',
            fontweight='bold')
    ax.text(1.5, yp, rxn, fontsize=9, fontweight='bold', va='center')
    ax.text(4.0, yp, annot, fontsize=8, va='center', color='#444444')
    ax.text(7.5, yp, strategy, fontsize=8, fontweight='bold',
            va='center', color=clr)
    ax.text(8.3, yp, fc, fontsize=9, fontweight='bold', va='center',
            color=clr)
    ax.text(9.5, yp, genes, fontsize=7.5, va='center', color='#1a5276',
            style='italic')

# --- Tier 3 ---
tier3 = FancyBboxPatch((0.5, 2.0), 15, 2.5, boxstyle="round,pad=0.3",
                        facecolor='#fff7bc', edgecolor='#d95f0e', linewidth=2)
ax.add_patch(tier3)
ax.text(8, 4.1, 'TIER 3: SYNERGISTIC GENE COMBINATIONS', fontsize=12,
        fontweight='bold', ha='center', color='#d95f0e')
ax.text(8, 3.7,
        '(Maximum lysine gain with 100% biomass viability  -  all targets have known genes)',
        fontsize=9, ha='center', color='#666666')

combos_display = [
    ('KO: R00216 + KO: R00282', 'Pyr DH + ACC block', '16,234x',
     'GRMZM2G085019x5 + GRMZM2G000236'),
    ('MIM: MR00422 + KO: R00216', 'Lipid down + Pyr DH block', '16,234x',
     'GRMZM2G002656 + GRMZM2G085019x5'),
    ('OE: R01280 + KO: R00216', 'Chorismate up + Pyr DH block', '16,234x',
     'GRMZM2G002614 + GRMZM2G085019x5'),
    ('OE: R00883 + KO: R00282', 'PFK up + ACC block', '28x',
     'GRMZM2G119300 + GRMZM2G000236'),
]
for i, (combo, annot, fc, genes) in enumerate(combos_display):
    yp = 3.2 - i * 0.32
    ax.text(1.0, yp, '>', fontsize=14, color='#d95f0e', va='center',
            fontweight='bold')
    ax.text(1.5, yp, combo, fontsize=9, fontweight='bold', va='center')
    ax.text(6.5, yp, annot, fontsize=8, va='center', color='#444444')
    ax.text(10.0, yp, 'FC: ' + fc, fontsize=9, fontweight='bold',
            va='center', color='#d95f0e')
    ax.text(12.0, yp, genes, fontsize=7, va='center', color='#1a5276',
            style='italic')

# --- Reference box ---
ref = FancyBboxPatch((0.5, 0.2), 15, 1.5, boxstyle="round,pad=0.3",
                      facecolor='#f0f0f0', edgecolor='#333333', linewidth=2)
ax.add_patch(ref)
ax.text(8, 1.4, 'REFERENCE: Lysine Biosynthetic Pathway Genes',
        fontsize=11, fontweight='bold', ha='center', color='#333333')
ax.text(1.5, 0.9,
        'DHDPS (R02292[K,p]): GRMZM2G027835  --  Committed step of lysine biosynthesis (overexpress)',
        fontsize=8.5, va='center', color='#1a5276')
ax.text(1.5, 0.5,
        'DAP Decarboxylase (R00451[K,p]): GRMZM2G020446  --  Final step producing free L-lysine',
        fontsize=8.5, va='center', color='#1a5276')

savefig(fig, 'Fig8_Engineering_Strategy_Diagram')
print("  Saved Fig8_Engineering_Strategy_Diagram")


# ===================================================================
# FIGURE 9: Temporal Dynamics
# ===================================================================
print("Generating Figure 9: Temporal Pathway Dynamics...")

fig, axes = plt.subplots(2, 3, figsize=(18, 10))

key_reactions = {
    'R00216[K,p]': 'Pyruvate DH (plastid)\nGRMZM2G085019 +4',
    'R00282[K,c]': 'Acetyl-CoA Carboxylase\nGRMZM2G000236',
    'R00342[K,p]': 'Transketolase (plastid)\nGRMZM2G035767 +10',
    'R00883[K,c]': 'Phosphofructokinase\nGRMZM2G119300',
    'R00512[K,c]': 'Pyruvate Kinase\nGRMZM2G067453 +3',
    'R00667[K,m]': 'TCA Cycle (mito.)\nGRMZM2G080828',
}

for idx, (rxn, label) in enumerate(key_reactions.items()):
    ax = axes[idx // 3, idx % 3]
    wt_vals, o2_vals = [], []
    for row in diff_flux:
        if row['Reaction'] == rxn:
            for dap in DAPS:
                wt_vals.append(float(row.get('WT_DAP%d' % dap, 0)))
                o2_vals.append(float(row.get('O2_DAP%d' % dap, 0)))
            break

    if wt_vals and o2_vals:
        ax.plot(DAPS, wt_vals, 'o-', color=WT_COLOR, label='WT')
        ax.plot(DAPS, o2_vals, 's-', color=O2_COLOR, label='o2')
        ax.set_xlabel('DAP')
        ax.set_ylabel('Flux')
        ax.set_title(label, fontsize=10)
        ax.legend(fontsize=8, frameon=True)
        ax.set_xticks(DAPS)
        ax.grid(True)
    else:
        ax.text(0.5, 0.5, '%s\n(not found)' % rxn,
                ha='center', va='center', transform=ax.transAxes,
                fontsize=10, color='#999999')
        ax.set_title(label, fontsize=10)

fig.suptitle('Temporal Dynamics of Key Target Reactions (WT vs o2)',
             fontsize=14, fontweight='bold', y=1.02)
fig.tight_layout()
savefig(fig, 'Fig9_Temporal_Dynamics')
print("  Saved Fig9_Temporal_Dynamics")


# ===================================================================
# FIGURE 10: Screen Type Distribution
# ===================================================================
print("Generating Figure 10: Screen Type Distribution...")

fig, axes = plt.subplots(1, 3, figsize=(16, 5.5))

ko_only  = sum(1 for r in ranking
               if 'KO' in r['Screen_Types']
               and 'OE' not in r['Screen_Types']
               and 'MIM' not in r['Screen_Types'])
oe_only  = sum(1 for r in ranking
               if 'OE' in r['Screen_Types']
               and 'KO' not in r['Screen_Types']
               and 'MIM' not in r['Screen_Types'])
mim_only = sum(1 for r in ranking
               if 'MIM' in r['Screen_Types']
               and 'KO' not in r['Screen_Types']
               and 'OE' not in r['Screen_Types'])
multi    = sum(1 for r in ranking if '+' in r['Screen_Types'])

# --- Panel A: Pie ---
ax = axes[0]
sizes  = [ko_only, oe_only, mim_only, multi]
labels_p = ['KO only\n(%d)' % ko_only, 'OE only\n(%d)' % oe_only,
            'Mimicry only\n(%d)' % mim_only, 'Multi-screen\n(%d)' % multi]
colors_p = [KO_COLOR, OE_COLOR, MIM_COLOR, COMBO_COLOR]
nz = [(s, l, c) for s, l, c in zip(sizes, labels_p, colors_p) if s > 0]
if nz:
    s_, l_, c_ = zip(*nz)
    ax.pie(s_, labels=l_, colors=c_, autopct='%1.0f%%', startangle=90,
           textprops={'fontsize': 9})
ax.set_title('A. Target Distribution by Screen')

# --- Panel B: Boxplot ---
ax = axes[1]
ko_fcs  = [np.log10(max(float(r['Max_Fold_Change']), 1.001))
           for r in ranking if 'KO' in r['Screen_Types']]
oe_fcs  = [np.log10(max(float(r['Max_Fold_Change']), 1.001))
           for r in ranking if 'OE' in r['Screen_Types']]
mim_fcs = [np.log10(max(float(r['Max_Fold_Change']), 1.001))
           for r in ranking if 'MIM' in r['Screen_Types']]

bp = ax.boxplot([ko_fcs, oe_fcs, mim_fcs],
                tick_labels=['Knockout', 'Overexpr.', 'Mimicry'],
                patch_artist=True, widths=0.5)
for patch, color in zip(bp['boxes'], [KO_COLOR, OE_COLOR, MIM_COLOR]):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
ax.set_ylabel('log10(Max Fold Change)')
ax.set_title('B. Fold Change by Screen Type')
ax.grid(True, axis='y')

# --- Panel C: Scatter ---
ax = axes[2]
for r in ranking:
    score = float(r['Composite_Score'])
    fc = float(r['Max_Fold_Change'])
    c = screen_color(r['Screen_Types'])
    ax.scatter(np.log10(max(fc, 1.001)), score, c=c, s=40, alpha=0.7,
               edgecolors='black', linewidth=0.3)

ax.set_xlabel('log10(Max Fold Change)')
ax.set_ylabel('Composite Score')
ax.set_title('C. Composite Score vs Fold Change')
ax.grid(True)
ax.legend(handles=[
    plt.Line2D([0], [0], marker='o', color='w',
               markerfacecolor=KO_COLOR, markersize=8, label='KO'),
    plt.Line2D([0], [0], marker='o', color='w',
               markerfacecolor=OE_COLOR, markersize=8, label='OE'),
    plt.Line2D([0], [0], marker='o', color='w',
               markerfacecolor=MIM_COLOR, markersize=8, label='Mimicry'),
    plt.Line2D([0], [0], marker='o', color='w',
               markerfacecolor=COMBO_COLOR, markersize=8, label='Multi'),
], fontsize=8, frameon=True)

fig.tight_layout(w_pad=3)
savefig(fig, 'Fig10_Screen_Type_Distribution')
print("  Saved Fig10_Screen_Type_Distribution")


# ===================================================================
# Summary
# ===================================================================
print("\n" + "=" * 70)
print("All figures saved to: %s" % FIG_DIR)
print("=" * 70)
print("\nFigure inventory:")
for f in sorted(os.listdir(FIG_DIR)):
    if f.endswith('.png'):
        print("  %s" % f)
