#!/usr/bin/env python3
"""
==============================================================================
Publication-Quality Figures: Flux Sum Analysis
O2 Mutant vs Wild Type Maize Endosperm - Metabolic Reprogramming
==============================================================================

Generates multi-panel figures suitable for journal submission:
  Figure 1: Overview panel (biomass + pathway flux sums + lysine temporal)
  Figure 2: Heatmap of differential flux sums across development
  Figure 3: Metabolite target ranking (lollipop + radar)
  Figure 4: Lysine node flux sum decomposition (reaction contributions)
  Figure 5: TCA/Aspartate pathway rewiring (focused temporal comparison)
  Figure 6: Supplementary — all pathway temporal profiles
==============================================================================
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D


# ============================================================================
# CONFIGURATION
# ============================================================================
BASE_DIR = "/mnt/nrdstor/ssbio/nampimah/ENDOSPERM/2026/Endosperm_Model_MARY"
FSA_DIR = os.path.join(BASE_DIR, "flux_sum_analysis")
ANA_DIR = os.path.join(BASE_DIR, "analysis_results")
FIG_DIR = os.path.join(FSA_DIR, "publication_figures")
os.makedirs(FIG_DIR, exist_ok=True)

DAP_VALUES = [6, 8, 10, 12, 15, 18, 22, 30]

# ---- Publication RC params ----
plt.rcParams.update({
    'font.size': 10,
    'font.family': 'sans-serif',
    'font.sans-serif': ['DejaVu Sans', 'Arial', 'Helvetica'],
    'axes.titlesize': 11,
    'axes.titleweight': 'bold',
    'axes.labelsize': 10,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 8.5,
    'legend.frameon': True,
    'legend.framealpha': 0.9,
    'legend.edgecolor': '0.8',
    'figure.dpi': 150,
    'savefig.dpi': 600,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'lines.linewidth': 1.8,
    'lines.markersize': 5,
    'pdf.fonttype': 42,   # TrueType for editable text in PDF
    'ps.fonttype': 42,
})

# ---- Colorblind-friendly palette ----
WT_COLOR = '#2166ac'
O2_COLOR = '#b2182b'
WT_LIGHT = '#92c5de'
O2_LIGHT = '#f4a582'
NEUTRAL   = '#636363'
INCREASE_COLOR = '#d73027'
DECREASE_COLOR = '#4575b4'
HIGHLIGHT = '#fc8d59'
GREEN = '#1a9850'

# Pathway colors
PATHWAY_COLORS = {
    'Lysine_Biosynthesis': '#d73027',
    'Lysine_Degradation':  '#f46d43',
    'Aspartate_Family':    '#fdae61',
    'TCA_Cycle':           '#66c2a5',
    'Glycolysis':          '#3288bd',
    'Energy_Cofactors':    '#5e4fa2',
    'Amino_Acids':         '#8c510a',
    'Polyamines':          '#c7eae5',
}

PATHWAY_LABELS = {
    'Lysine_Biosynthesis': 'Lysine Biosynthesis',
    'Lysine_Degradation':  'Lysine Degradation',
    'Aspartate_Family':    'Aspartate Family',
    'TCA_Cycle':           'TCA Cycle',
    'Glycolysis':          'Glycolysis',
    'Energy_Cofactors':    'Energy Cofactors',
    'Amino_Acids':         'Amino Acids',
    'Polyamines':          'Polyamines',
}


# ============================================================================
# LOAD DATA
# ============================================================================
print("Loading data ...")

df_biomass = pd.read_csv(os.path.join(ANA_DIR, "01_biomass_comparison.csv"))
df_pathways = pd.read_csv(os.path.join(FSA_DIR, "05_pathway_flux_sum_comparison.csv"))
df_targets = pd.read_csv(os.path.join(FSA_DIR, "06_metabolite_targets_ranked.csv"))
df_diff = pd.read_csv(os.path.join(FSA_DIR, "02_differential_flux_sums_ranked.csv"))
df_lysine_pw = pd.read_csv(os.path.join(FSA_DIR, "03_lysine_pathway_flux_sums.csv"))
df_contrib = pd.read_csv(os.path.join(FSA_DIR, "04_lysine_flux_sum_contributions.csv"))
df_complete = pd.read_csv(os.path.join(FSA_DIR, "01_flux_sums_complete.csv"))

print(f"  Biomass: {len(df_biomass)} rows")
print(f"  Pathways: {len(df_pathways)} rows")
print(f"  Targets: {len(df_targets)} rows")
print(f"  Differential: {len(df_diff)} rows")
print(f"  Contributions: {len(df_contrib)} rows")


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================
def add_panel_label(ax, label, x=-0.12, y=1.08, fontsize=14):
    """Add panel label (A, B, C, ...) to axis."""
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=fontsize, fontweight='bold', va='top', ha='left')


def clean_met_name(name):
    """Shorten metabolite name for figure labels."""
    replacements = {
        'D-Fructose-1,6-bisP': 'Fru-1,6-bisP',
        'D-Glyceraldehyde-3P': 'G3P',
        'alpha-D-Glucose-6P': 'Glc-6P',
        'D-Fructose-6P': 'Fru-6P',
        'Phosphoenolpyruvate': 'PEP',
        '3-Phospho-D-glyceroyl-phosphate': '1,3-BPG',
        '2-Phospho-D-glycerate': '2PG',
        '3-Phospho-D-glycerate': '3PG',
        'Glycerone-phosphate': 'DHAP',
        'L-Aspartate-4-semialdehyde': 'Asp-semialdehyde',
        'meso-Diaminopimelate': 'meso-DAP',
        'Dihydrodipicolinate': 'DHDP',
        'D-Glucose': 'Glucose',
        'D-Fructose': 'Fructose',
        'L-Aspartate': 'Aspartate',
        'L-Glutamate': 'Glutamate',
        'L-Lysine': 'Lysine',
        'L-Leucine': 'Leucine',
        'L-Methionine': 'Methionine',
        'L-Threonine': 'Threonine',
        'L-Asparagine': 'Asparagine',
        'L-Homoserine': 'Homoserine',
        'L-Saccharopine': 'Saccharopine',
        'L-2-Aminoadipate-6-semialdehyde': 'AAS',
        'L-2-Aminoadipate': 'AAA',
        'Oxaloacetate': 'OAA',
        '2-Oxoglutarate': 'α-KG',
        'Orthophosphate': 'Pi',
        'Diphosphate': 'PPi',
        'Acetyl-CoA': 'Acetyl-CoA',
    }
    for old, new in replacements.items():
        name = name.replace(old, new)
    return name


# ============================================================================
# FIGURE 1: Multi-Panel Overview
# 4 panels: (A) Biomass, (B) Lysine flux sum temporal,
#            (C) Key pathway flux sum bar comparison, (D) Aspartate pathway
# ============================================================================
print("\nGenerating Figure 1: Multi-panel overview ...")

fig = plt.figure(figsize=(7.5, 7.5))
gs = gridspec.GridSpec(2, 2, hspace=0.38, wspace=0.35)

# ---- Panel A: Biomass ----
ax_a = fig.add_subplot(gs[0, 0])
x = np.arange(len(DAP_VALUES))
w = 0.35
ax_a.bar(x - w/2, df_biomass['WT_Biomass']*1000, w,
         color=WT_COLOR, edgecolor='white', linewidth=0.5, label='Wild Type', zorder=3)
ax_a.bar(x + w/2, df_biomass['O2_Biomass']*1000, w,
         color=O2_COLOR, edgecolor='white', linewidth=0.5, label='o2 Mutant', zorder=3)
ax_a.set_xticks(x)
ax_a.set_xticklabels([f'{d}' for d in DAP_VALUES])
ax_a.set_xlabel('DAP')
ax_a.set_ylabel('Biomass Flux (×10⁻³)')
ax_a.set_title('Seed Biomass Production')
ax_a.legend(loc='upper right')
ax_a.grid(axis='y', alpha=0.25, zorder=0)
add_panel_label(ax_a, 'A')

# ---- Panel B: Lysine (plastid) Flux Sum Temporal ----
ax_b = fig.add_subplot(gs[0, 1])
# Get lysine plastid data
lys_p = df_complete[df_complete['Metabolite_ID'] == 'C00047[K,p]']
if len(lys_p) > 0:
    wt_lys = [lys_p.iloc[0][f'WT_DAP{d}'] for d in DAP_VALUES]
    o2_lys = [lys_p.iloc[0][f'O2_DAP{d}'] for d in DAP_VALUES]
else:
    wt_lys = [0]*8
    o2_lys = [0]*8

ax_b.plot(DAP_VALUES, wt_lys, 'o-', color=WT_COLOR, label='Wild Type', zorder=5)
ax_b.plot(DAP_VALUES, o2_lys, 's-', color=O2_COLOR, label='o2 Mutant', zorder=5)
ax_b.set_xlabel('DAP')
ax_b.set_ylabel('Flux Sum (Φ)')
ax_b.set_title('Lysine Turnover (Plastid)')
ax_b.set_xticks(DAP_VALUES)
ax_b.legend(loc='upper left')
ax_b.grid(alpha=0.25)
# Annotate the spike at DAP22
o2_max_idx = np.argmax(o2_lys)
if o2_lys[o2_max_idx] > 0.01:
    ax_b.annotate(f'{o2_lys[o2_max_idx]:.2f}',
                  xy=(DAP_VALUES[o2_max_idx], o2_lys[o2_max_idx]),
                  xytext=(DAP_VALUES[o2_max_idx]-3, o2_lys[o2_max_idx]*0.85),
                  fontsize=8, color=O2_COLOR, fontweight='bold',
                  arrowprops=dict(arrowstyle='->', color=O2_COLOR, lw=0.8))
add_panel_label(ax_b, 'B')

# ---- Panel C: Pathway Flux Sum Fold Changes (grouped bar) ----
ax_c = fig.add_subplot(gs[1, 0])
# Pick key pathways (exclude Amino_Acids which is huge, and Polyamines which is tiny)
key_pw = ['Lysine_Biosynthesis', 'Aspartate_Family', 'TCA_Cycle', 'Glycolysis', 'Energy_Cofactors']
pw_fold_changes = {}
for pw in key_pw:
    sub = df_pathways[df_pathways['Pathway'] == pw]
    fc_vals = []
    for d in DAP_VALUES:
        row = sub[sub['DAP'] == d]
        if len(row) > 0:
            fc = row.iloc[0]['Fold_Change']
            if pd.notna(fc) and np.isfinite(fc):
                fc_vals.append(fc)
            else:
                fc_vals.append(1.0)
        else:
            fc_vals.append(1.0)
    pw_fold_changes[pw] = fc_vals

n_pw = len(key_pw)
x_pw = np.arange(len(DAP_VALUES))
total_w = 0.75
bar_w = total_w / n_pw

for i, pw in enumerate(key_pw):
    offset = (i - n_pw/2 + 0.5) * bar_w
    vals = np.array(pw_fold_changes[pw])
    # Cap extreme values for display
    vals_display = np.clip(vals, 0, 20)
    ax_c.bar(x_pw + offset, vals_display, bar_w,
             color=PATHWAY_COLORS[pw], edgecolor='white', linewidth=0.3,
             label=PATHWAY_LABELS[pw], zorder=3)
    # Add text for clipped bars
    for j, (v, vd) in enumerate(zip(vals, vals_display)):
        if v > 20:
            ax_c.text(x_pw[j] + offset, vd + 0.3, f'{v:.0f}×',
                      ha='center', va='bottom', fontsize=5.5, fontweight='bold',
                      color=PATHWAY_COLORS[pw], rotation=90)

ax_c.axhline(y=1, color='black', linestyle='--', linewidth=0.7, alpha=0.5, zorder=2)
ax_c.set_xticks(x_pw)
ax_c.set_xticklabels([f'{d}' for d in DAP_VALUES])
ax_c.set_xlabel('DAP')
ax_c.set_ylabel('Fold Change (o2 / WT)')
ax_c.set_title('Pathway Flux Sum Fold Changes')
ax_c.legend(loc='upper left', fontsize=6.5, ncol=1)
ax_c.set_ylim(0, 22)
ax_c.grid(axis='y', alpha=0.2, zorder=0)
add_panel_label(ax_c, 'C')

# ---- Panel D: Aspartate Family Flux Sum Temporal ----
ax_d = fig.add_subplot(gs[1, 1])
# Aspartate in cytosol
asp_c = df_complete[df_complete['Metabolite_ID'] == 'C00049[K,c]']
oaa_c = df_complete[df_complete['Metabolite_ID'] == 'C00036[K,c]']
akg_c = df_complete[df_complete['Metabolite_ID'] == 'C00026[K,c]']
glu_c = df_complete[df_complete['Metabolite_ID'] == 'C00025[K,c]']

mets_to_plot = [
    ('C00049[K,c]', 'Aspartate [cyt]', '#d73027', 'o'),
    ('C00036[K,c]', 'OAA [cyt]',       '#fc8d59', 's'),
    ('C00026[K,c]', 'α-KG [cyt]',      '#91bfdb', '^'),
    ('C00025[K,c]', 'Glutamate [cyt]',  '#4575b4', 'D'),
]

for met_id, label, color, marker in mets_to_plot:
    row = df_complete[df_complete['Metabolite_ID'] == met_id]
    if len(row) > 0:
        diffs = [row.iloc[0][f'Diff_DAP{d}'] for d in DAP_VALUES]
        ax_d.plot(DAP_VALUES, diffs, f'{marker}-', color=color, label=label,
                  markersize=5, zorder=5)

ax_d.axhline(y=0, color='black', linestyle='-', linewidth=0.5, alpha=0.4)
ax_d.set_xlabel('DAP')
ax_d.set_ylabel('ΔΦ (o2 − WT)')
ax_d.set_title('Key Precursor Flux Sum Changes')
ax_d.set_xticks(DAP_VALUES)
ax_d.legend(loc='upper left', fontsize=7)
ax_d.grid(alpha=0.25)
add_panel_label(ax_d, 'D')

for fmt in ['png', 'pdf']:
    fig.savefig(os.path.join(FIG_DIR, f'Figure1_overview.{fmt}'), dpi=600)
plt.close(fig)
print("  Saved Figure 1.")


# ============================================================================
# FIGURE 2: Heatmap of Differential Flux Sums
# ============================================================================
print("Generating Figure 2: Heatmap ...")

# Select top 40 non-currency metabolites
# Exclude C00001, C00080, C00007, C00011
currency = {'C00001', 'C00080', 'C00007', 'C00011'}
df_filt = df_diff[~df_diff['Base_KEGG_ID'].isin(currency)].head(40).copy()
diff_cols = [f'Diff_DAP{d}' for d in DAP_VALUES]
hm_data = df_filt[diff_cols].values
met_names = [clean_met_name(n) for n in df_filt['Metabolite_Name'].values]

fig, ax = plt.subplots(figsize=(5.5, 9))
vmax = np.percentile(np.abs(hm_data), 95)
norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
im = ax.imshow(hm_data, aspect='auto', cmap='RdBu_r', norm=norm, interpolation='nearest')

ax.set_xticks(range(len(DAP_VALUES)))
ax.set_xticklabels([f'DAP {d}' for d in DAP_VALUES], fontsize=8, rotation=45, ha='right')
ax.set_yticks(range(len(met_names)))
ax.set_yticklabels(met_names, fontsize=7)
ax.set_title('Differential Flux Sum (o2 − WT)\nAcross Endosperm Development', pad=12)

# Add colorbar
cbar = plt.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
cbar.set_label('ΔΦ (o2 − WT)', fontsize=9)
cbar.ax.tick_params(labelsize=8)

# Add thin grid lines
for i in range(len(met_names)):
    ax.axhline(y=i-0.5, color='white', linewidth=0.3)
for j in range(len(DAP_VALUES)):
    ax.axvline(x=j-0.5, color='white', linewidth=0.3)

plt.tight_layout()
for fmt in ['png', 'pdf']:
    fig.savefig(os.path.join(FIG_DIR, f'Figure2_heatmap.{fmt}'), dpi=600)
plt.close(fig)
print("  Saved Figure 2.")


# ============================================================================
# FIGURE 3: Metabolite Target Ranking
# Two panels: (A) Lollipop chart — composite score, (B) Stacked contributions
# ============================================================================
print("Generating Figure 3: Target ranking ...")

fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(10, 6.5),
                                  gridspec_kw={'width_ratios': [1.2, 1], 'wspace': 0.45})

# ---- Panel A: Lollipop chart of top 25 targets ----
top25 = df_targets.head(25).copy()
top25 = top25.iloc[::-1].reset_index(drop=True)
y_pos = np.arange(len(top25))
labels = [clean_met_name(n) for n in top25['Metabolite_Name'].values]
scores = top25['Composite_Score'].values
colors = [INCREASE_COLOR if d == 'INCREASED' else DECREASE_COLOR for d in top25['Direction_in_O2'].values]

ax_a.hlines(y_pos, 0, scores, colors=colors, linewidth=1.5, zorder=3)
ax_a.scatter(scores, y_pos, c=colors, s=45, edgecolors='white', linewidths=0.5, zorder=4)
ax_a.set_yticks(y_pos)
ax_a.set_yticklabels(labels, fontsize=7.5)
ax_a.set_xlabel('Composite Target Score')
ax_a.set_title('Metabolite Targets for\nLysine Enhancement')
ax_a.grid(axis='x', alpha=0.2)

# Legend
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor=INCREASE_COLOR,
           markersize=8, label='↑ Turnover in o2'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor=DECREASE_COLOR,
           markersize=8, label='↓ Turnover in o2'),
]
ax_a.legend(handles=legend_elements, loc='lower right', fontsize=8)
add_panel_label(ax_a, 'A', x=-0.18)

# ---- Panel B: Temporal profile of top biologically relevant targets ----
bio_targets = [
    ('C00049[K,c]', 'Aspartate [cyt]',      '#d73027'),
    ('C00036[K,c]', 'OAA [cyt]',             '#fc8d59'),
    ('C00149[K,c]', 'Malate [cyt]',          '#91bfdb'),
    ('C00026[K,p]', 'α-KG [pla]',            '#4575b4'),
    ('C00025[K,p]', 'Glutamate [pla]',       '#1a9850'),
    ('C00047[K,p]', 'Lysine [pla]',          '#762a83'),
    ('C00022[K,c]', 'Pyruvate [cyt]',        '#b35806'),
]

for met_id, label, color in bio_targets:
    row = df_complete[df_complete['Metabolite_ID'] == met_id]
    if len(row) > 0:
        diffs = [row.iloc[0][f'Diff_DAP{d}'] for d in DAP_VALUES]
        ax_b.plot(DAP_VALUES, diffs, 'o-', color=color, label=label, markersize=4, zorder=5)

ax_b.axhline(y=0, color='black', linestyle='-', linewidth=0.5, alpha=0.4)
ax_b.set_xlabel('DAP')
ax_b.set_ylabel('ΔΦ (o2 − WT)')
ax_b.set_title('Temporal Flux Sum Changes\nof Key Precursors')
ax_b.set_xticks(DAP_VALUES)
ax_b.legend(loc='upper left', fontsize=7, ncol=1)
ax_b.grid(alpha=0.25)
add_panel_label(ax_b, 'B', x=-0.14)

for fmt in ['png', 'pdf']:
    fig.savefig(os.path.join(FIG_DIR, f'Figure3_targets.{fmt}'), dpi=600)
plt.close(fig)
print("  Saved Figure 3.")


# ============================================================================
# FIGURE 4: Lysine Flux Sum Decomposition (Stacked Bar)
# Shows reaction contributions to lysine flux sum at key DAPs
# ============================================================================
print("Generating Figure 4: Lysine flux sum decomposition ...")

# Focus on cytosolic lysine C00047[K,c]
df_lys_c = df_contrib[df_contrib['Metabolite'] == 'C00047[K,c]'].copy()

# Key DAPs
focus_daps = [6, 10, 15, 18, 22, 30]

fig, axes = plt.subplots(2, 3, figsize=(10, 6), sharey=False)
axes = axes.flatten()

for idx, dap in enumerate(focus_daps):
    ax = axes[idx]
    sub = df_lys_c[df_lys_c['DAP'] == dap].copy()

    # Get active reactions (nonzero in either WT or O2)
    sub = sub[(sub['WT_abs_Sv'].abs() > 1e-12) | (sub['O2_abs_Sv'].abs() > 1e-12)].copy()
    sub = sub.sort_values('O2_abs_Sv', ascending=False)

    if len(sub) == 0:
        ax.set_title(f'DAP {dap}', fontsize=9)
        ax.text(0.5, 0.5, 'No active\nreactions', transform=ax.transAxes, ha='center', va='center',
                fontsize=8, color='gray')
        continue

    # Show top 6 reactions
    sub = sub.head(6)
    rxn_labels = []
    for _, r in sub.iterrows():
        rxn = r['Reaction']
        # Shorten reaction names
        short = rxn.replace('[K,', '[').replace(']', ']')
        if len(short) > 22:
            short = short[:20] + '..'
        rxn_labels.append(short)

    y_r = np.arange(len(rxn_labels))
    bw = 0.35
    ax.barh(y_r - bw/2, sub['WT_abs_Sv'].values, bw, color=WT_COLOR, label='WT', edgecolor='white', linewidth=0.3)
    ax.barh(y_r + bw/2, sub['O2_abs_Sv'].values, bw, color=O2_COLOR, label='o2', edgecolor='white', linewidth=0.3)
    ax.set_yticks(y_r)
    ax.set_yticklabels(rxn_labels, fontsize=6.5)
    ax.set_title(f'DAP {dap}', fontsize=9, fontweight='bold')
    ax.set_xlabel('|S·v|', fontsize=8)
    ax.grid(axis='x', alpha=0.2)
    if idx == 0:
        ax.legend(fontsize=7, loc='lower right')

fig.suptitle('Reaction Contributions to Lysine (Cytosol) Flux Sum',
             fontsize=12, fontweight='bold', y=1.01)
plt.tight_layout()
for fmt in ['png', 'pdf']:
    fig.savefig(os.path.join(FIG_DIR, f'Figure4_lysine_decomposition.{fmt}'), dpi=600)
plt.close(fig)
print("  Saved Figure 4.")


# ============================================================================
# FIGURE 5: TCA / Aspartate Pathway Rewiring
# Paired panels showing WT vs O2 for key TCA + aspartate metabolites
# ============================================================================
print("Generating Figure 5: TCA/Aspartate rewiring ...")

fig, axes = plt.subplots(2, 3, figsize=(10, 6.5), sharex=True)

panel_mets = [
    ('C00036[K,c]', 'Oxaloacetate [cyt]'),
    ('C00149[K,c]', 'Malate [cyt]'),
    ('C00026[K,c]', '2-Oxoglutarate [cyt]'),
    ('C00049[K,c]', 'Aspartate [cyt]'),
    ('C00025[K,c]', 'Glutamate [cyt]'),
    ('C00047[K,p]', 'Lysine [pla]'),
]

panel_labels = ['A', 'B', 'C', 'D', 'E', 'F']

for idx, (met_id, title) in enumerate(panel_mets):
    ax = axes.flatten()[idx]
    row = df_complete[df_complete['Metabolite_ID'] == met_id]
    if len(row) > 0:
        wt_vals = [row.iloc[0][f'WT_DAP{d}'] for d in DAP_VALUES]
        o2_vals = [row.iloc[0][f'O2_DAP{d}'] for d in DAP_VALUES]

        ax.fill_between(DAP_VALUES, wt_vals, alpha=0.15, color=WT_COLOR)
        ax.fill_between(DAP_VALUES, o2_vals, alpha=0.15, color=O2_COLOR)
        ax.plot(DAP_VALUES, wt_vals, 'o-', color=WT_COLOR, label='WT', markersize=4, zorder=5)
        ax.plot(DAP_VALUES, o2_vals, 's-', color=O2_COLOR, label='o2', markersize=4, zorder=5)

    ax.set_title(title, fontsize=9.5)
    ax.set_xticks(DAP_VALUES)
    ax.set_ylabel('Φ', fontsize=9)
    ax.grid(alpha=0.2)
    if idx == 0:
        ax.legend(fontsize=7.5, loc='best')
    add_panel_label(ax, panel_labels[idx], x=-0.15, y=1.1, fontsize=12)

for ax in axes[1, :]:
    ax.set_xlabel('DAP')

fig.suptitle('TCA Cycle / Aspartate Pathway Metabolite Turnover',
             fontsize=12, fontweight='bold', y=1.01)
plt.tight_layout()
for fmt in ['png', 'pdf']:
    fig.savefig(os.path.join(FIG_DIR, f'Figure5_TCA_aspartate_rewiring.{fmt}'), dpi=600)
plt.close(fig)
print("  Saved Figure 5.")


# ============================================================================
# FIGURE 6: Pathway Flux Sum Comparison (Grouped Bar Per Pathway)
# Supplementary — full pathway temporal comparison
# ============================================================================
print("Generating Figure 6: Pathway comparisons (supplementary) ...")

all_pathways = ['Lysine_Biosynthesis', 'Lysine_Degradation', 'Aspartate_Family',
                'TCA_Cycle', 'Glycolysis', 'Energy_Cofactors']

fig, axes = plt.subplots(3, 2, figsize=(9, 10), sharex=True)
axes = axes.flatten()

for idx, pw in enumerate(all_pathways):
    ax = axes[idx]
    sub = df_pathways[df_pathways['Pathway'] == pw]

    wt_vals = []
    o2_vals = []
    for d in DAP_VALUES:
        row = sub[sub['DAP'] == d]
        if len(row) > 0:
            wt_vals.append(row.iloc[0]['WT_FluxSum_Total'])
            o2_vals.append(row.iloc[0]['O2_FluxSum_Total'])
        else:
            wt_vals.append(0)
            o2_vals.append(0)

    x_pos = np.arange(len(DAP_VALUES))
    w = 0.35
    ax.bar(x_pos - w/2, wt_vals, w, color=WT_COLOR, edgecolor='white', linewidth=0.3,
           label='Wild Type', zorder=3)
    ax.bar(x_pos + w/2, o2_vals, w, color=O2_COLOR, edgecolor='white', linewidth=0.3,
           label='o2 Mutant', zorder=3)

    ax.set_title(PATHWAY_LABELS[pw], fontsize=10, fontweight='bold',
                 color=PATHWAY_COLORS[pw])
    ax.set_ylabel('Σ Flux Sum')
    ax.set_xticks(x_pos)
    ax.set_xticklabels([f'{d}' for d in DAP_VALUES])
    ax.grid(axis='y', alpha=0.2, zorder=0)

    if idx == 0:
        ax.legend(fontsize=7.5)
    add_panel_label(ax, chr(65+idx), x=-0.12, y=1.08, fontsize=12)

for ax in axes[4:]:
    ax.set_xlabel('DAP')

fig.suptitle('Pathway-Level Flux Sum Comparison Across Development',
             fontsize=13, fontweight='bold', y=1.0)
plt.tight_layout()
for fmt in ['png', 'pdf']:
    fig.savefig(os.path.join(FIG_DIR, f'Figure6_pathway_comparison_supp.{fmt}'), dpi=600)
plt.close(fig)
print("  Saved Figure 6.")


# ============================================================================
# FIGURE 7: Metabolic Strategy Summary (Dot Plot)
# x = Mean |ΔΦ|, y = Consistency, size = Relevance, color = direction
# ============================================================================
print("Generating Figure 7: Strategy dot plot ...")

# Take top 50 targets excluding pure currency
top50 = df_targets[~df_targets['Base_KEGG_ID'].isin(currency)].head(50).copy()

fig, ax = plt.subplots(figsize=(7, 5.5))

sizes = top50['Relevance_Weight'].values * 25
colors = [INCREASE_COLOR if d == 'INCREASED' else DECREASE_COLOR for d in top50['Direction_in_O2'].values]
alphas = 0.5 + 0.5 * (top50['Consistency'].values)

for i in range(len(top50)):
    ax.scatter(top50.iloc[i]['Mean_Abs_Diff'],
               top50.iloc[i]['Consistency'],
               s=sizes[i], c=colors[i], alpha=alphas[i],
               edgecolors='black', linewidths=0.3, zorder=4)

# Label top 15
for i, (_, row) in enumerate(top50.head(15).iterrows()):
    name = clean_met_name(row['Metabolite_Name'])
    if len(name) > 20:
        name = name[:18] + '..'
    ax.annotate(name,
                xy=(row['Mean_Abs_Diff'], row['Consistency']),
                xytext=(5, 4), textcoords='offset points',
                fontsize=6, color='#333333',
                arrowprops=dict(arrowstyle='-', color='#999999', lw=0.3))

ax.set_xlabel('Mean |ΔΦ| (O2 − WT)', fontsize=10)
ax.set_ylabel('Direction Consistency', fontsize=10)
ax.set_title('Metabolite Target Landscape', fontsize=11, fontweight='bold')
ax.grid(alpha=0.2)

# Legend for direction and size
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor=INCREASE_COLOR,
           markersize=10, label='↑ in o2 (enhance in WT)'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor=DECREASE_COLOR,
           markersize=10, label='↓ in o2 (reduce in WT)'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='gray',
           markersize=6, label='Low relevance'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='gray',
           markersize=12, label='High relevance'),
]
ax.legend(handles=legend_elements, loc='upper right', fontsize=7.5, framealpha=0.9)

plt.tight_layout()
for fmt in ['png', 'pdf']:
    fig.savefig(os.path.join(FIG_DIR, f'Figure7_strategy_dotplot.{fmt}'), dpi=600)
plt.close(fig)
print("  Saved Figure 7.")


# ============================================================================
# FIGURE 8: Focused Lysine Biosynthesis Pathway Bar Chart
# ============================================================================
print("Generating Figure 8: Lysine biosynthesis pathway intermediates ...")

lys_bio_mets = [
    ('C00049[K,p]', 'Asp\n[pla]'),
    ('C00441[K,p]', 'Asp-SA\n[pla]'),
    ('C03340[K,p]', 'DHDP\n[pla]'),
    ('C00680[K,p]', 'meso-DAP\n[pla]'),
    ('C00047[K,p]', 'Lys\n[pla]'),
    ('C00047[K,c]', 'Lys\n[cyt]'),
]

fig, axes = plt.subplots(1, 2, figsize=(9, 4.5), gridspec_kw={'width_ratios': [1.5, 1], 'wspace': 0.35})

# Panel A: DAP 22 comparison (where the biggest difference is)
ax = axes[0]
met_labels = []
wt_vals_22 = []
o2_vals_22 = []
for met_id, label in lys_bio_mets:
    row = df_complete[df_complete['Metabolite_ID'] == met_id]
    if len(row) > 0:
        wt_vals_22.append(row.iloc[0]['WT_DAP22'])
        o2_vals_22.append(row.iloc[0]['O2_DAP22'])
    else:
        wt_vals_22.append(0)
        o2_vals_22.append(0)
    met_labels.append(label)

x_pos = np.arange(len(met_labels))
w = 0.35
ax.bar(x_pos - w/2, wt_vals_22, w, color=WT_COLOR, edgecolor='white', linewidth=0.5,
       label='Wild Type', zorder=3)
ax.bar(x_pos + w/2, o2_vals_22, w, color=O2_COLOR, edgecolor='white', linewidth=0.5,
       label='o2 Mutant', zorder=3)
ax.set_xticks(x_pos)
ax.set_xticklabels(met_labels, fontsize=8)
ax.set_ylabel('Flux Sum (Φ)')
ax.set_title('Lysine Biosynthesis Intermediates at DAP 22')
ax.legend(fontsize=8)
ax.grid(axis='y', alpha=0.2, zorder=0)

# Add fold change annotations
for i, (wt, o2) in enumerate(zip(wt_vals_22, o2_vals_22)):
    if wt > 1e-10:
        fc = o2 / wt
        ax.text(i + w/2, o2 + 0.02, f'{fc:.0f}×', ha='center', va='bottom',
                fontsize=7, fontweight='bold', color=O2_COLOR)
add_panel_label(ax, 'A', x=-0.1)

# Panel B: Lysine plastid across all DAPs
ax = axes[1]
row = df_complete[df_complete['Metabolite_ID'] == 'C00047[K,p]']
if len(row) > 0:
    wt_vals = [row.iloc[0][f'WT_DAP{d}'] for d in DAP_VALUES]
    o2_vals = [row.iloc[0][f'O2_DAP{d}'] for d in DAP_VALUES]
    ax.semilogy(DAP_VALUES, [max(v, 1e-8) for v in wt_vals], 'o-', color=WT_COLOR,
                label='Wild Type', markersize=5, zorder=5)
    ax.semilogy(DAP_VALUES, [max(v, 1e-8) for v in o2_vals], 's-', color=O2_COLOR,
                label='o2 Mutant', markersize=5, zorder=5)
ax.set_xlabel('DAP')
ax.set_ylabel('Flux Sum (Φ, log scale)')
ax.set_title('Lysine [pla] Turnover')
ax.set_xticks(DAP_VALUES)
ax.legend(fontsize=8)
ax.grid(alpha=0.25)
add_panel_label(ax, 'B', x=-0.14)

plt.tight_layout()
for fmt in ['png', 'pdf']:
    fig.savefig(os.path.join(FIG_DIR, f'Figure8_lysine_biosynthesis.{fmt}'), dpi=600)
plt.close(fig)
print("  Saved Figure 8.")


# ============================================================================
# FIGURE 9: Glycolysis / Carbon Source Rewiring
# ============================================================================
print("Generating Figure 9: Carbon flux rewiring ...")

carbon_mets = [
    ('C00031[K,c]', 'Glucose [cyt]',   '#e41a1c'),
    ('C00095[K,c]', 'Fructose [cyt]',  '#ff7f00'),
    ('C00668[K,c]', 'Glc-6P [cyt]',    '#984ea3'),
    ('C00354[K,c]', 'Fru-1,6-bisP [cyt]', '#377eb8'),
    ('C00118[K,c]', 'G3P [cyt]',       '#4daf4a'),
    ('C00074[K,c]', 'PEP [cyt]',       '#a65628'),
    ('C00022[K,c]', 'Pyruvate [cyt]',  '#f781bf'),
]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5),
                                gridspec_kw={'width_ratios': [1, 1], 'wspace': 0.3})

# Panel A: WT glycolytic flux sums across DAPs
for met_id, label, color in carbon_mets:
    row = df_complete[df_complete['Metabolite_ID'] == met_id]
    if len(row) > 0:
        wt_vals = [row.iloc[0][f'WT_DAP{d}'] for d in DAP_VALUES]
        ax1.plot(DAP_VALUES, wt_vals, 'o-', color=color, label=label, markersize=4, zorder=5)
ax1.set_xlabel('DAP')
ax1.set_ylabel('Flux Sum (Φ)')
ax1.set_title('Wild Type — Carbon Metabolite Turnover')
ax1.set_xticks(DAP_VALUES)
ax1.legend(fontsize=6.5, ncol=1, loc='upper left')
ax1.grid(alpha=0.25)
add_panel_label(ax1, 'A', x=-0.1)

# Panel B: Difference
for met_id, label, color in carbon_mets:
    row = df_complete[df_complete['Metabolite_ID'] == met_id]
    if len(row) > 0:
        diffs = [row.iloc[0][f'Diff_DAP{d}'] for d in DAP_VALUES]
        ax2.plot(DAP_VALUES, diffs, 'o-', color=color, label=label, markersize=4, zorder=5)
ax2.axhline(y=0, color='black', linestyle='-', linewidth=0.5, alpha=0.4)
ax2.set_xlabel('DAP')
ax2.set_ylabel('ΔΦ (o2 − WT)')
ax2.set_title('Carbon Metabolite Flux Sum Difference')
ax2.set_xticks(DAP_VALUES)
ax2.legend(fontsize=6.5, ncol=1, loc='upper left')
ax2.grid(alpha=0.25)
add_panel_label(ax2, 'B', x=-0.1)

plt.tight_layout()
for fmt in ['png', 'pdf']:
    fig.savefig(os.path.join(FIG_DIR, f'Figure9_carbon_rewiring.{fmt}'), dpi=600)
plt.close(fig)
print("  Saved Figure 9.")


# ============================================================================
# DONE
# ============================================================================
print(f"\n{'='*70}")
print(f"ALL PUBLICATION FIGURES SAVED TO: {FIG_DIR}")
print(f"{'='*70}")

# List all generated files
for f in sorted(os.listdir(FIG_DIR)):
    fpath = os.path.join(FIG_DIR, f)
    size_kb = os.path.getsize(fpath) / 1024
    print(f"  {f:<50s} {size_kb:7.1f} KB")
