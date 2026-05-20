#!/usr/bin/env python3
"""
Visualization script for O2 mutant vs Wild Type metabolic reprogramming analysis.
Generates publication-quality figures.
"""

import os
import csv
import sys

# Check matplotlib availability
try:
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend for HPC
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    from matplotlib.ticker import ScalarFormatter
    HAS_MPL = True
except ImportError:
    HAS_MPL = False
    print("WARNING: matplotlib not available. Install with: pip install matplotlib")
    print("Generating text-based plots instead.")

BASE_DIR = "/mnt/nrdstor/ssbio/nampimah/ENDOSPERM/2026/Endosperm_Model_MARY"
OUTPUT_DIR = os.path.join(BASE_DIR, "analysis_results")
FIG_DIR = os.path.join(OUTPUT_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

DAP_VALUES = [6, 8, 10, 12, 15, 18, 22, 30]

# ============================================================================
# Load data from CSV files
# ============================================================================
def load_csv(filename):
    filepath = os.path.join(OUTPUT_DIR, filename)
    rows = []
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    return rows

def load_dat(filename):
    filepath = os.path.join(OUTPUT_DIR, filename)
    data = []
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            parts = line.strip().split('\t')
            data.append([float(x) for x in parts])
    return data

if not HAS_MPL:
    print("Cannot generate figures without matplotlib. Exiting.")
    sys.exit(1)

# Color scheme
WT_COLOR = '#2166ac'    # Blue
O2_COLOR = '#b2182b'    # Red
DIFF_COLOR = '#4daf4a'  # Green

plt.rcParams.update({
    'font.size': 11,
    'axes.titlesize': 13,
    'axes.labelsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

# ============================================================================
# FIGURE 1: Biomass Comparison
# ============================================================================
print("Generating Figure 1: Biomass Comparison...")
biomass_data = load_dat("plot_biomass.dat")
daps = [row[0] for row in biomass_data]
wt_bm = [row[1] for row in biomass_data]
o2_bm = [row[2] for row in biomass_data]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

ax1.plot(daps, wt_bm, 'o-', color=WT_COLOR, linewidth=2, markersize=8, label='Wild Type (B73)')
ax1.plot(daps, o2_bm, 's-', color=O2_COLOR, linewidth=2, markersize=8, label='O2 Mutant')
ax1.set_xlabel('Days After Pollination (DAP)')
ax1.set_ylabel('Biomass Flux (Seed_Biomass)')
ax1.set_title('A. Seed Biomass Across Development')
ax1.legend()
ax1.set_xticks(daps)
ax1.ticklabel_format(style='scientific', axis='y', scilimits=(0,0))
ax1.grid(True, alpha=0.3)

# Fold change
fc = [o/w if w > 0 else 0 for w, o in zip(wt_bm, o2_bm)]
colors = [O2_COLOR if f > 1 else WT_COLOR for f in fc]
ax2.bar(range(len(daps)), fc, color=colors, alpha=0.7, edgecolor='black', linewidth=0.5)
ax2.axhline(y=1, color='black', linestyle='--', linewidth=1)
ax2.set_xticks(range(len(daps)))
ax2.set_xticklabels([f'DAP {d}' for d in daps], rotation=45)
ax2.set_ylabel('Fold Change (O2/WT)')
ax2.set_title('B. Biomass Fold Change')
ax2.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "Fig1_Biomass_Comparison.png"))
plt.savefig(os.path.join(FIG_DIR, "Fig1_Biomass_Comparison.pdf"))
plt.close()
print("  Saved Fig1_Biomass_Comparison.png/pdf")

# ============================================================================
# FIGURE 2: Lysine Biosynthesis Flux (R00451 - DAP decarboxylase)
# ============================================================================
print("Generating Figure 2: Lysine Biosynthesis...")
bio_data = load_dat("plot_lysine_biosynthesis.dat")
wt_bio = [row[1] for row in bio_data]
o2_bio = [row[2] for row in bio_data]

fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# Full scale
ax = axes[0]
ax.plot(daps, wt_bio, 'o-', color=WT_COLOR, linewidth=2, markersize=8, label='WT')
ax.plot(daps, o2_bio, 's-', color=O2_COLOR, linewidth=2, markersize=8, label='O2')
ax.set_xlabel('DAP')
ax.set_ylabel('Flux (R00451 - DAP Decarboxylase)')
ax.set_title('A. Lysine Biosynthesis (Full Scale)')
ax.legend()
ax.set_xticks(daps)
ax.grid(True, alpha=0.3)

# Zoomed (early DAPs only, where both are comparable)
ax = axes[1]
early_idx = [i for i, d in enumerate(daps) if d <= 12]
ax.plot([daps[i] for i in early_idx], [wt_bio[i] for i in early_idx], 'o-', 
        color=WT_COLOR, linewidth=2, markersize=8, label='WT')
ax.plot([daps[i] for i in early_idx], [o2_bio[i] for i in early_idx], 's-', 
        color=O2_COLOR, linewidth=2, markersize=8, label='O2')
ax.set_xlabel('DAP')
ax.set_ylabel('Flux')
ax.set_title('B. Early Development (DAP 6-12)')
ax.legend()
ax.set_xticks([daps[i] for i in early_idx])
ax.ticklabel_format(style='scientific', axis='y', scilimits=(0,0))
ax.grid(True, alpha=0.3)

# Log scale fold change
ax = axes[2]
import math
fc_bio = []
for w, o in zip(wt_bio, o2_bio):
    if w > 1e-15:
        fc_bio.append(o / w)
    else:
        fc_bio.append(1)
colors = [O2_COLOR if f > 1 else WT_COLOR for f in fc_bio]
ax.bar(range(len(daps)), [math.log10(max(f, 0.001)) for f in fc_bio], 
       color=colors, alpha=0.7, edgecolor='black', linewidth=0.5)
ax.axhline(y=0, color='black', linestyle='--', linewidth=1)
ax.set_xticks(range(len(daps)))
ax.set_xticklabels([f'DAP {d}' for d in daps], rotation=45)
ax.set_ylabel('log10(Fold Change O2/WT)')
ax.set_title('C. Biosynthesis Fold Change (log10)')
ax.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "Fig2_Lysine_Biosynthesis.png"))
plt.savefig(os.path.join(FIG_DIR, "Fig2_Lysine_Biosynthesis.pdf"))
plt.close()
print("  Saved Fig2_Lysine_Biosynthesis.png/pdf")

# ============================================================================
# FIGURE 3: DHDPS (pathway entry) vs DAP decarboxylase (pathway exit)
# ============================================================================
print("Generating Figure 3: Pathway Entry vs Exit...")
dhdps_data = load_dat("plot_DHDPS.dat")
wt_dhdps = [row[1] for row in dhdps_data]
o2_dhdps = [row[2] for row in dhdps_data]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

ax1.plot(daps, wt_dhdps, 'o-', color=WT_COLOR, linewidth=2, markersize=8, label='WT DHDPS')
ax1.plot(daps, o2_dhdps, 's-', color=O2_COLOR, linewidth=2, markersize=8, label='O2 DHDPS')
ax1.set_xlabel('DAP')
ax1.set_ylabel('Flux')
ax1.set_title('A. DHDPS (R02292) - Entry to Lysine Branch')
ax1.legend()
ax1.set_xticks(daps)
ax1.grid(True, alpha=0.3)

# Compare entry and exit (should be 1:1 in steady state)
ax2.scatter(wt_dhdps, wt_bio, color=WT_COLOR, s=80, marker='o', label='WT', zorder=3)
ax2.scatter(o2_dhdps, o2_bio, color=O2_COLOR, s=80, marker='s', label='O2', zorder=3)
max_val = max(max(wt_dhdps + o2_dhdps), max(wt_bio + o2_bio))
ax2.plot([0, max_val], [0, max_val], 'k--', alpha=0.5, label='1:1 line')
ax2.set_xlabel('DHDPS Flux (Pathway Entry)')
ax2.set_ylabel('DAP Decarboxylase Flux (Pathway Exit)')
ax2.set_title('B. Pathway Entry vs Exit (Flux Coupling)')
ax2.legend()
ax2.grid(True, alpha=0.3)
# Add DAP labels
for i, d in enumerate(daps):
    ax2.annotate(f'DAP{d}', (wt_dhdps[i], wt_bio[i]), textcoords="offset points",
                xytext=(5,5), fontsize=7, color=WT_COLOR)
    ax2.annotate(f'DAP{d}', (o2_dhdps[i], o2_bio[i]), textcoords="offset points",
                xytext=(5,5), fontsize=7, color=O2_COLOR)

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "Fig3_Pathway_Entry_Exit.png"))
plt.savefig(os.path.join(FIG_DIR, "Fig3_Pathway_Entry_Exit.pdf"))
plt.close()
print("  Saved Fig3_Pathway_Entry_Exit.png/pdf")

# ============================================================================
# FIGURE 4: Lysine Degradation (LKR)
# ============================================================================
print("Generating Figure 4: Lysine Degradation...")
deg_data = load_dat("plot_lysine_degradation.dat")
wt_lkr_c = [row[1] for row in deg_data]
o2_lkr_c = [row[2] for row in deg_data]
wt_lkr_p = [row[3] for row in deg_data]
o2_lkr_p = [row[4] for row in deg_data]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

x = range(len(daps))
width = 0.35
ax1.bar([i - width/2 for i in x], wt_lkr_c, width, color=WT_COLOR, alpha=0.7, label='WT', edgecolor='black', linewidth=0.5)
ax1.bar([i + width/2 for i in x], o2_lkr_c, width, color=O2_COLOR, alpha=0.7, label='O2', edgecolor='black', linewidth=0.5)
ax1.set_xticks(x)
ax1.set_xticklabels([f'DAP {d}' for d in daps], rotation=45)
ax1.set_ylabel('Flux')
ax1.set_title('A. LKR Cytosolic (R00715[K,c]) - Lys Degradation')
ax1.legend()
ax1.ticklabel_format(style='scientific', axis='y', scilimits=(0,0))
ax1.grid(True, alpha=0.3, axis='y')

# Lysine exchange
exch_data = load_dat("plot_lysine_exchange.dat")
wt_exch = [row[1] for row in exch_data]
o2_exch = [row[2] for row in exch_data]

ax2.plot(daps, wt_exch, 'o-', color=WT_COLOR, linewidth=2, markersize=8, label='WT Exchange')
ax2.plot(daps, o2_exch, 's-', color=O2_COLOR, linewidth=2, markersize=8, label='O2 Exchange')
ax2.set_xlabel('DAP')
ax2.set_ylabel('Exchange Flux')
ax2.set_title('B. Lysine Exchange Reaction')
ax2.legend()
ax2.set_xticks(daps)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "Fig4_Lysine_Degradation_Exchange.png"))
plt.savefig(os.path.join(FIG_DIR, "Fig4_Lysine_Degradation_Exchange.pdf"))
plt.close()
print("  Saved Fig4_Lysine_Degradation_Exchange.png/pdf")

# ============================================================================
# FIGURE 5: Pathway-Level Heatmap
# ============================================================================
print("Generating Figure 5: Pathway Heatmap...")
pathway_data = load_csv("07_pathway_flux_summary.csv")

pathways = []
seen = set()
for row in pathway_data:
    p = row['Pathway']
    if p not in seen:
        pathways.append(p)
        seen.add(p)

# Build percent change matrix
pct_matrix = []
for pathway in pathways:
    row_vals = []
    for dap in DAP_VALUES:
        matching = [r for r in pathway_data if r['Pathway'] == pathway and int(r['DAP']) == dap]
        if matching:
            pct = float(matching[0]['Percent_Change'])
            row_vals.append(pct)
        else:
            row_vals.append(0)
    pct_matrix.append(row_vals)

fig, ax = plt.subplots(figsize=(12, 6))

# Clip extreme values for visualization
import numpy as np
pct_arr = np.array(pct_matrix)
vmax = 200
pct_clipped = np.clip(pct_arr, -vmax, vmax)

im = ax.imshow(pct_clipped, cmap='RdBu_r', aspect='auto', vmin=-vmax, vmax=vmax)
ax.set_xticks(range(len(DAP_VALUES)))
ax.set_xticklabels([f'DAP {d}' for d in DAP_VALUES])
ax.set_yticks(range(len(pathways)))
ax.set_yticklabels([p.replace('_', ' ') for p in pathways])
ax.set_title('Pathway Flux Change: O2 Mutant vs Wild Type (% Change, clipped ±200%)')
ax.set_xlabel('Days After Pollination')

# Add text annotations
for i in range(len(pathways)):
    for j in range(len(DAP_VALUES)):
        val = pct_arr[i, j]
        if abs(val) > 1000:
            text = f'{val/1000:.0f}K'
        else:
            text = f'{val:.0f}'
        color = 'white' if abs(pct_clipped[i, j]) > 100 else 'black'
        ax.text(j, i, text, ha='center', va='center', fontsize=8, color=color)

cbar = plt.colorbar(im, ax=ax, label='% Change (O2 - WT)')
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "Fig5_Pathway_Heatmap.png"))
plt.savefig(os.path.join(FIG_DIR, "Fig5_Pathway_Heatmap.pdf"))
plt.close()
print("  Saved Fig5_Pathway_Heatmap.png/pdf")

# ============================================================================
# FIGURE 6: Lysine Flux Balance (Production vs Consumption)
# ============================================================================
print("Generating Figure 6: Lysine Flux Balance...")
balance = load_csv("03_lysine_flux_balance.csv")

wt_prod = [float(r['WT_Lys_Production']) for r in balance]
o2_prod = [float(r['O2_Lys_Production']) for r in balance]
wt_cons = [abs(float(r['WT_Lys_Consumption'])) for r in balance]
o2_cons = [abs(float(r['O2_Lys_Consumption'])) for r in balance]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

x = np.arange(len(DAP_VALUES))
width = 0.2

ax1.bar(x - 1.5*width, wt_prod, width, color=WT_COLOR, alpha=0.7, label='WT Production', edgecolor='black', linewidth=0.5)
ax1.bar(x - 0.5*width, o2_prod, width, color=O2_COLOR, alpha=0.7, label='O2 Production', edgecolor='black', linewidth=0.5)
ax1.bar(x + 0.5*width, wt_cons, width, color=WT_COLOR, alpha=0.3, hatch='//', label='WT Consumption', edgecolor=WT_COLOR, linewidth=0.5)
ax1.bar(x + 1.5*width, o2_cons, width, color=O2_COLOR, alpha=0.3, hatch='//', label='O2 Consumption', edgecolor=O2_COLOR, linewidth=0.5)
ax1.set_xticks(x)
ax1.set_xticklabels([f'DAP {d}' for d in DAP_VALUES], rotation=45)
ax1.set_ylabel('Absolute Flux')
ax1.set_title('A. Lysine Production & Consumption')
ax1.legend(fontsize=8)
ax1.grid(True, alpha=0.3, axis='y')

# Ratio of production O2/WT 
ratios = [o/w if w > 1e-15 else 1 for w, o in zip(wt_prod, o2_prod)]
colors = [O2_COLOR if r > 1 else WT_COLOR for r in ratios]
ax2.bar(x, [math.log10(max(r, 0.01)) for r in ratios], color=colors, alpha=0.7, 
        edgecolor='black', linewidth=0.5)
ax2.axhline(y=0, color='black', linestyle='--', linewidth=1)
ax2.set_xticks(x)
ax2.set_xticklabels([f'DAP {d}' for d in DAP_VALUES], rotation=45)
ax2.set_ylabel('log10(O2/WT Production Ratio)')
ax2.set_title('B. Lysine Turnover Ratio (log10)')
ax2.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "Fig6_Lysine_Flux_Balance.png"))
plt.savefig(os.path.join(FIG_DIR, "Fig6_Lysine_Flux_Balance.pdf"))
plt.close()
print("  Saved Fig6_Lysine_Flux_Balance.png/pdf")

# ============================================================================
# FIGURE 7: Constraint (Bound) Differences for Key Reactions
# ============================================================================
print("Generating Figure 7: Bound Differences...")
bounds_data = load_csv("08_bound_differences.csv")

# Focus on R00716 (LKR alternative) - shows clear constraint difference
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

key_bound_rxns = [
    ("R00716[K,c]", "LKR Alt (Lys Degradation)"),
    ("R00451[K,p]", "DAP Decarboxylase (Lys Biosynthesis)"),
    ("R02292[K,p]", "DHDPS (Pathway Entry)"),
    ("R03102[K,c]", "Aminoadipate-semialdehyde DH (Degradation)"),
]

for idx, (rxn, title) in enumerate(key_bound_rxns):
    ax = axes[idx // 2][idx % 2]
    rxn_bounds = [r for r in bounds_data if r['Reaction'] == rxn]
    if not rxn_bounds:
        ax.set_title(f'{title}\n(No data)')
        continue
    
    d = [int(r['DAP']) for r in rxn_bounds]
    wt_hi = [float(r['WT_vmax']) for r in rxn_bounds]
    o2_hi = [float(r['O2_vmax']) for r in rxn_bounds]
    
    ax.plot(d, wt_hi, 'o-', color=WT_COLOR, linewidth=2, markersize=7, label='WT v_max')
    ax.plot(d, o2_hi, 's-', color=O2_COLOR, linewidth=2, markersize=7, label='O2 v_max')
    ax.fill_between(d, wt_hi, o2_hi, alpha=0.1, color='gray')
    ax.set_xlabel('DAP')
    ax.set_ylabel('Bound Value')
    ax.set_title(f'{rxn}\n{title}')
    ax.legend(fontsize=9)
    ax.set_xticks(d)
    ax.grid(True, alpha=0.3)

plt.suptitle('Flux Bound (Constraint) Differences: WT vs O2 Mutant', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "Fig7_Bound_Differences.png"))
plt.savefig(os.path.join(FIG_DIR, "Fig7_Bound_Differences.pdf"))
plt.close()
print("  Saved Fig7_Bound_Differences.png/pdf")

# ============================================================================
# FIGURE 8: Comprehensive Summary Panel
# ============================================================================
print("Generating Figure 8: Summary Panel...")
fig = plt.figure(figsize=(18, 14))
gs = gridspec.GridSpec(3, 3, hspace=0.4, wspace=0.35)

# Panel A: Biomass
ax = fig.add_subplot(gs[0, 0])
ax.plot(daps, wt_bm, 'o-', color=WT_COLOR, linewidth=2, markersize=6, label='WT')
ax.plot(daps, o2_bm, 's-', color=O2_COLOR, linewidth=2, markersize=6, label='O2')
ax.set_xlabel('DAP')
ax.set_ylabel('Biomass Flux')
ax.set_title('A. Seed Biomass')
ax.legend(fontsize=8)
ax.set_xticks(daps)
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax.grid(True, alpha=0.3)

# Panel B: Lysine biosynthesis
ax = fig.add_subplot(gs[0, 1])
ax.plot(daps, wt_bio, 'o-', color=WT_COLOR, linewidth=2, markersize=6, label='WT')
ax.plot(daps, o2_bio, 's-', color=O2_COLOR, linewidth=2, markersize=6, label='O2')
ax.set_xlabel('DAP')
ax.set_ylabel('R00451 Flux')
ax.set_title('B. Lysine Biosynthesis (DAP Decarboxylase)')
ax.legend(fontsize=8)
ax.set_xticks(daps)
ax.grid(True, alpha=0.3)

# Panel C: LKR degradation
ax = fig.add_subplot(gs[0, 2])
x = np.arange(len(daps))
width = 0.35
ax.bar(x - width/2, wt_lkr_c, width, color=WT_COLOR, alpha=0.7, label='WT')
ax.bar(x + width/2, o2_lkr_c, width, color=O2_COLOR, alpha=0.7, label='O2')
ax.set_xticks(x)
ax.set_xticklabels([str(d) for d in daps])
ax.set_xlabel('DAP')
ax.set_ylabel('R00715 Flux')
ax.set_title('C. LKR (Lysine Degradation)')
ax.legend(fontsize=8)
ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax.grid(True, alpha=0.3, axis='y')

# Panel D: Lysine flux balance
ax = fig.add_subplot(gs[1, 0])
ax.plot(daps, wt_prod, 'o-', color=WT_COLOR, linewidth=2, markersize=6, label='WT')
ax.plot(daps, o2_prod, 's-', color=O2_COLOR, linewidth=2, markersize=6, label='O2')
ax.set_xlabel('DAP')
ax.set_ylabel('Lysine Total Production')
ax.set_title('D. Lysine Production Flux')
ax.legend(fontsize=8)
ax.set_xticks(daps)
ax.grid(True, alpha=0.3)

# Panel E: Exchange reaction
ax = fig.add_subplot(gs[1, 1])
ax.plot(daps, wt_exch, 'o-', color=WT_COLOR, linewidth=2, markersize=6, label='WT')
ax.plot(daps, o2_exch, 's-', color=O2_COLOR, linewidth=2, markersize=6, label='O2')
ax.set_xlabel('DAP')
ax.set_ylabel('Exchange Flux')
ax.set_title('E. Lysine Exchange')
ax.legend(fontsize=8)
ax.set_xticks(daps)
ax.grid(True, alpha=0.3)

# Panel F: DHDPS
ax = fig.add_subplot(gs[1, 2])
ax.plot(daps, wt_dhdps, 'o-', color=WT_COLOR, linewidth=2, markersize=6, label='WT')
ax.plot(daps, o2_dhdps, 's-', color=O2_COLOR, linewidth=2, markersize=6, label='O2')
ax.set_xlabel('DAP')
ax.set_ylabel('R02292 Flux')
ax.set_title('F. DHDPS (Pathway Entry)')
ax.legend(fontsize=8)
ax.set_xticks(daps)
ax.grid(True, alpha=0.3)

# Panel G: Pathway heatmap (big panel)
ax = fig.add_subplot(gs[2, :])
pct_clipped = np.clip(pct_arr, -200, 200)
im = ax.imshow(pct_clipped, cmap='RdBu_r', aspect='auto', vmin=-200, vmax=200)
ax.set_xticks(range(len(DAP_VALUES)))
ax.set_xticklabels([f'DAP {d}' for d in DAP_VALUES])
ax.set_yticks(range(len(pathways)))
ax.set_yticklabels([p.replace('_', ' ') for p in pathways], fontsize=9)
ax.set_title('G. Pathway-Level Flux Change (% O2 vs WT, clipped ±200%)')
for i in range(len(pathways)):
    for j in range(len(DAP_VALUES)):
        val = pct_arr[i, j]
        text = f'{val/1000:.0f}K' if abs(val) > 1000 else f'{val:.0f}'
        color = 'white' if abs(pct_clipped[i, j]) > 100 else 'black'
        ax.text(j, i, text, ha='center', va='center', fontsize=7, color=color)
plt.colorbar(im, ax=ax, label='% Change', shrink=0.8)

plt.suptitle('Metabolic Reprogramming in O2 Mutant Maize Endosperm:\nLysine Pathway Analysis Across Development',
             fontsize=15, fontweight='bold', y=1.01)
plt.savefig(os.path.join(FIG_DIR, "Fig8_Summary_Panel.png"))
plt.savefig(os.path.join(FIG_DIR, "Fig8_Summary_Panel.pdf"))
plt.close()
print("  Saved Fig8_Summary_Panel.png/pdf")

# ============================================================================
# FIGURE 9: Uniquely Active Reactions Count
# ============================================================================
print("Generating Figure 9: Uniquely Active Reactions...")
unique_data = load_csv("10_uniquely_active_reactions.csv")

wt_unique_count = {d: 0 for d in DAP_VALUES}
o2_unique_count = {d: 0 for d in DAP_VALUES}
for row in unique_data:
    d = int(row['DAP'])
    if row['Active_In'] == 'WT_only':
        wt_unique_count[d] += 1
    else:
        o2_unique_count[d] += 1

fig, ax = plt.subplots(figsize=(10, 5))
x = np.arange(len(DAP_VALUES))
width = 0.35
ax.bar(x - width/2, [wt_unique_count[d] for d in DAP_VALUES], width, 
       color=WT_COLOR, alpha=0.7, label='Unique to WT', edgecolor='black', linewidth=0.5)
ax.bar(x + width/2, [o2_unique_count[d] for d in DAP_VALUES], width, 
       color=O2_COLOR, alpha=0.7, label='Unique to O2', edgecolor='black', linewidth=0.5)
ax.set_xticks(x)
ax.set_xticklabels([f'DAP {d}' for d in DAP_VALUES])
ax.set_ylabel('Number of Uniquely Active Reactions')
ax.set_title('Reactions Active in Only One Genotype')
ax.legend()
ax.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "Fig9_Uniquely_Active_Reactions.png"))
plt.savefig(os.path.join(FIG_DIR, "Fig9_Uniquely_Active_Reactions.pdf"))
plt.close()
print("  Saved Fig9_Uniquely_Active_Reactions.png/pdf")

print("\n" + "=" * 60)
print("ALL FIGURES GENERATED!")
print(f"Location: {FIG_DIR}")
print("=" * 60)
for f in sorted(os.listdir(FIG_DIR)):
    print(f"  {f}")
