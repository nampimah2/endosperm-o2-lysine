#!/usr/bin/env python3
"""
Flux Variability Analysis of WT and o2 mutant maize endosperm.
Reads WILDTYPE_FVA.xlsx and O2_MUTANT_FVA.xlsx (one sheet per DAP).
Outputs:
  fva_results/01_lysine_pathway_fva.csv       -- key pathway reaction FVA ranges
  fva_results/02_flexibility_comparison.csv   -- flux range ratio O2/WT per reaction
  fva_results/03_essential_reactions.csv      -- reactions with zero range (fixed)
  fva_results/04_fva_summary.csv              -- per-DAP summary statistics
  fva_results/05_active_reactions.csv         -- active reaction counts per DAP
  fva_results/figures/                        -- publication figures
"""
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import warnings
warnings.filterwarnings('ignore')

BASE = "/mnt/nrdstor/ssbio/nampimah/ENDOSPERM/2026/Endosperm_Model_MARY"
OUT  = os.path.join(BASE, "fva_results")
FIG  = os.path.join(OUT, "figures")
os.makedirs(FIG, exist_ok=True)

DAP_MAP = {
    'B6': 6, 'B8': 8, 'B10': 10, 'B12': 12,
    'B15': 15, 'B18': 18, 'B22': 22, 'B30': 30
}
DAP_LIST = [6, 8, 10, 12, 15, 18, 22, 30]

# Key reactions for focused analysis
KEY_REACTIONS = {
    'R00451[K,p]':           'DAP decarboxylase (Lys biosyn, final)',
    'R00480[K,p]':           'Aspartate kinase (DAP pathway entry)',
    'R00715[K,c]':           'LKR cytosolic (Lys catabolism)',
    'R00715[K,p]':           'LKR plastidic (Lys catabolism)',
    'R00716[K,c]':           'Alternative LKR (Lys catabolism)',
    'R03102[K,c]':           'SDH (saccharopine dehydrogenase)',
    'cmTransport_C00047[K]': 'Lysine transport (cm)',
    'Sink1[K]':              'Lysine sink/export',
    'Lysine_PoolRxn[K]':     'Lysine pool reaction',
}

# ── Load all sheets ───────────────────────────────────────────────────────────
print("Loading FVA data...")
wt_xl = pd.ExcelFile(os.path.join(BASE, 'WILDTYPE_FVA.xlsx'))
o2_xl = pd.ExcelFile(os.path.join(BASE, 'O2_MUTANT_FVA.xlsx'))

wt_data = {}  # dap -> DataFrame
o2_data = {}

for wt_sh, dap in DAP_MAP.items():
    o2_sh = wt_sh.replace('B', 'O')
    dfw = wt_xl.parse(wt_sh)
    dfw.columns = [c.lower() if c == 'pathway' else c for c in dfw.columns]
    dfw = dfw.rename(columns={'pathway': 'Pathway'})
    dfw['flux_range'] = dfw['max_rate'] - dfw['min_rate']
    wt_data[dap] = dfw.set_index('reaction')

    dfo = o2_xl.parse(o2_sh)
    dfo['flux_range'] = dfo['max_rate'] - dfo['min_rate']
    o2_data[dap] = dfo.set_index('reaction')

# ── 1. Lysine pathway FVA table ───────────────────────────────────────────────
print("Computing lysine pathway FVA...")
rows = []
for rxn, desc in KEY_REACTIONS.items():
    for dap in DAP_LIST:
        wt_row = wt_data[dap].loc[rxn] if rxn in wt_data[dap].index else None
        o2_row = o2_data[dap].loc[rxn] if rxn in o2_data[dap].index else None
        row = {
            'Reaction': rxn,
            'Description': desc,
            'DAP': dap,
            'WT_min': wt_row['min_rate'] if wt_row is not None else np.nan,
            'WT_max': wt_row['max_rate'] if wt_row is not None else np.nan,
            'WT_range': wt_row['flux_range'] if wt_row is not None else np.nan,
            'O2_min': o2_row['min_rate'] if o2_row is not None else np.nan,
            'O2_max': o2_row['max_rate'] if o2_row is not None else np.nan,
            'O2_range': o2_row['flux_range'] if o2_row is not None else np.nan,
        }
        if row['WT_range'] is not None and row['WT_range'] > 1e-10:
            row['range_ratio'] = row['O2_range'] / row['WT_range']
        else:
            row['range_ratio'] = np.nan
        rows.append(row)

df_lys = pd.DataFrame(rows)
df_lys.to_csv(os.path.join(OUT, '01_lysine_pathway_fva.csv'), index=False)
print("  -> Saved 01_lysine_pathway_fva.csv")

# ── 2. Global flux range ratio O2/WT across all reactions ────────────────────
print("Computing global flexibility comparison...")
flex_rows = []
for dap in DAP_LIST:
    wt = wt_data[dap][['min_rate', 'max_rate', 'flux_range']].copy()
    o2 = o2_data[dap][['min_rate', 'max_rate', 'flux_range']].copy()
    merged = wt.join(o2, lsuffix='_wt', rsuffix='_o2')
    merged['DAP'] = dap
    # Ratio: < 1 means more constrained in O2, > 1 means more flexible in O2
    merged['range_ratio'] = np.where(
        merged['flux_range_wt'] > 1e-10,
        merged['flux_range_o2'] / merged['flux_range_wt'],
        np.nan
    )
    flex_rows.append(merged.reset_index())

df_flex = pd.concat(flex_rows, ignore_index=True)
df_flex.to_csv(os.path.join(OUT, '02_flexibility_comparison.csv'), index=False)
print("  -> Saved 02_flexibility_comparison.csv")

# ── 3. Essential reactions (zero range = fixed flux) ─────────────────────────
print("Identifying essential reactions...")
ess_rows = []
for dap in DAP_LIST:
    wt = wt_data[dap]
    o2 = o2_data[dap]
    # Reactions where range drops to zero in O2 but is nonzero in WT (newly constrained)
    wt_flex   = wt[wt['flux_range'] > 1e-10].index
    o2_fixed  = o2[o2['flux_range'] <= 1e-10].index
    newly_fixed = set(wt_flex) & set(o2_fixed)
    # Reactions where range opens in O2 vs zero in WT (newly flexible)
    wt_fixed  = wt[wt['flux_range'] <= 1e-10].index
    o2_flex   = o2[o2['flux_range'] > 1e-10].index
    newly_flex = set(wt_fixed) & set(o2_flex)
    # Reactions fixed in both
    both_fixed = set(wt[wt['flux_range'] <= 1e-10].index) & set(o2[o2['flux_range'] <= 1e-10].index)
    ess_rows.append({
        'DAP': dap,
        'WT_fixed': len(wt[wt['flux_range'] <= 1e-10]),
        'O2_fixed': len(o2[o2['flux_range'] <= 1e-10]),
        'Both_fixed': len(both_fixed),
        'Newly_fixed_in_O2': len(newly_fixed),
        'Newly_flexible_in_O2': len(newly_flex),
    })

df_ess = pd.DataFrame(ess_rows)
df_ess.to_csv(os.path.join(OUT, '03_essential_reactions.csv'), index=False)
print("  -> Saved 03_essential_reactions.csv")

# ── 4. Per-DAP FVA summary statistics ────────────────────────────────────────
print("Computing FVA summary statistics...")
sum_rows = []
for dap in DAP_LIST:
    wt = wt_data[dap]
    o2 = o2_data[dap]
    sum_rows.append({
        'DAP': dap,
        'WT_total_reactions': len(wt),
        'WT_active_max': (wt['max_rate'] > 1e-6).sum(),
        'WT_active_min': (wt['min_rate'] > 1e-6).sum(),
        'WT_mean_range': wt['flux_range'].mean(),
        'WT_sum_range': wt['flux_range'].sum(),
        'O2_total_reactions': len(o2),
        'O2_active_max': (o2['max_rate'] > 1e-6).sum(),
        'O2_active_min': (o2['min_rate'] > 1e-6).sum(),
        'O2_mean_range': o2['flux_range'].mean(),
        'O2_sum_range': o2['flux_range'].sum(),
        # Lysine-specific
        'WT_R00451_max': wt.loc['R00451[K,p]', 'max_rate'] if 'R00451[K,p]' in wt.index else np.nan,
        'O2_R00451_max': o2.loc['R00451[K,p]', 'max_rate'] if 'R00451[K,p]' in o2.index else np.nan,
        'WT_cmLys_max': wt.loc['cmTransport_C00047[K]', 'max_rate'] if 'cmTransport_C00047[K]' in wt.index else np.nan,
        'O2_cmLys_max': o2.loc['cmTransport_C00047[K]', 'max_rate'] if 'cmTransport_C00047[K]' in o2.index else np.nan,
    })

df_sum = pd.DataFrame(sum_rows)
df_sum.to_csv(os.path.join(OUT, '04_fva_summary.csv'), index=False)
print("  -> Saved 04_fva_summary.csv")
print()
print("=== FVA SUMMARY ===")
print(df_sum[['DAP','WT_active_max','O2_active_max','WT_R00451_max','O2_R00451_max','WT_cmLys_max','O2_cmLys_max']].to_string(index=False))

# ── 5. Active reaction counts ─────────────────────────────────────────────────
df_active = df_sum[['DAP','WT_active_max','O2_active_max']].copy()
df_active.to_csv(os.path.join(OUT, '05_active_reactions.csv'), index=False)

# ── FIGURES ──────────────────────────────────────────────────────────────────
print("\nGenerating FVA figures...")

COLORS = {'WT': '#2166AC', 'O2': '#D6604D'}
DAP_TICKS = DAP_LIST
x = np.arange(len(DAP_LIST))
W = 0.38  # bar width

# ── Fig A: Max achievable lysine flux (R00451) ────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

ax = axes[0]
wt_max = df_sum['WT_R00451_max'].values
o2_max = df_sum['O2_R00451_max'].values
ax.bar(x - W/2, wt_max, W, label='Wild Type', color=COLORS['WT'], alpha=0.85, edgecolor='white')
ax.bar(x + W/2, o2_max, W, label='o2 mutant',  color=COLORS['O2'], alpha=0.85, edgecolor='white')
ax.set_xticks(x); ax.set_xticklabels(['DAP %d' % d for d in DAP_LIST], rotation=45, ha='right')
ax.set_xlabel('Developmental Stage'); ax.set_ylabel('Max Achievable Flux\n(mmol g DW⁻¹ h⁻¹)')
ax.set_title('Maximum Achievable Lysine\nBiosynthesis Flux (R00451)', fontsize=11, fontweight='bold')
ax.legend(); ax.grid(axis='y', alpha=0.3)
for i, (w, o) in enumerate(zip(wt_max, o2_max)):
    if o > w * 1.5 and w < 0.1:
        ax.annotate('%.1fx' % (o / w if w > 1e-8 else 0),
                    xy=(i + W/2, o), xytext=(0, 4), textcoords='offset points',
                    ha='center', fontsize=8, color=COLORS['O2'])

ax = axes[1]
wt_cm = df_sum['WT_cmLys_max'].values
o2_cm = df_sum['O2_cmLys_max'].values
ax.bar(x - W/2, wt_cm, W, label='Wild Type', color=COLORS['WT'], alpha=0.85, edgecolor='white')
ax.bar(x + W/2, o2_cm, W, label='o2 mutant',  color=COLORS['O2'], alpha=0.85, edgecolor='white')
ax.set_xticks(x); ax.set_xticklabels(['DAP %d' % d for d in DAP_LIST], rotation=45, ha='right')
ax.set_xlabel('Developmental Stage'); ax.set_ylabel('Max Lysine Transport Flux\n(mmol g DW⁻¹ h⁻¹)')
ax.set_title('Maximum Achievable Lysine\nTransport (cmTransport_C00047)', fontsize=11, fontweight='bold')
ax.legend(); ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(FIG, 'FigFVA1_max_lysine_flux.png'), dpi=200, bbox_inches='tight')
plt.close()
print("  -> FigFVA1_max_lysine_flux.png")

# ── Fig B: Active reaction counts ────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 5))
ax.bar(x - W/2, df_sum['WT_active_max'].values, W, label='Wild Type', color=COLORS['WT'], alpha=0.85, edgecolor='white')
ax.bar(x + W/2, df_sum['O2_active_max'].values, W, label='o2 mutant',  color=COLORS['O2'], alpha=0.85, edgecolor='white')
ax.set_xticks(x); ax.set_xticklabels(['DAP %d' % d for d in DAP_LIST], rotation=45, ha='right')
ax.set_xlabel('Developmental Stage'); ax.set_ylabel('Number of Reactions\nwith max_rate > 10⁻⁶')
ax.set_title('Metabolic Flexibility: Active Reactions per Genotype', fontsize=11, fontweight='bold')
ax.legend(); ax.grid(axis='y', alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG, 'FigFVA2_active_reactions.png'), dpi=200, bbox_inches='tight')
plt.close()
print("  -> FigFVA2_active_reactions.png")

# ── Fig C: Newly fixed vs newly flexible reactions in O2 ─────────────────────
fig, ax = plt.subplots(figsize=(8, 5))
newly_fixed   = df_ess['Newly_fixed_in_O2'].values
newly_flexible = df_ess['Newly_flexible_in_O2'].values
ax.bar(x - W/2, newly_fixed,    W, label='Newly constrained in o2',  color='#C94B4B', alpha=0.85, edgecolor='white')
ax.bar(x + W/2, newly_flexible, W, label='Newly flexible in o2',     color='#4B8EC9', alpha=0.85, edgecolor='white')
ax.set_xticks(x); ax.set_xticklabels(['DAP %d' % d for d in DAP_LIST], rotation=45, ha='right')
ax.set_xlabel('Developmental Stage'); ax.set_ylabel('Number of Reactions')
ax.set_title('Flux Range Changes in o2 vs WT:\nConstrained vs Newly Flexible Reactions', fontsize=11, fontweight='bold')
ax.legend(); ax.grid(axis='y', alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(FIG, 'FigFVA3_constrained_flexible.png'), dpi=200, bbox_inches='tight')
plt.close()
print("  -> FigFVA3_constrained_flexible.png")

# ── Fig D: Lysine pathway FVA range heatmap ──────────────────────────────────
key_rxns_short = {
    'R00451[K,p]': 'R00451\n(DAP decarboxylase)',
    'R00480[K,p]': 'R00480\n(Aspartate kinase)',
    'R00715[K,c]': 'R00715c\n(LKR cytosol)',
    'R00715[K,p]': 'R00715p\n(LKR plastid)',
    'R00716[K,c]': 'R00716\n(Alt-LKR)',
    'cmTransport_C00047[K]': 'cmTransport\n(Lys)',
    'Sink1[K]': 'Sink1\n(Lys export)',
}

fig, axes = plt.subplots(1, 2, figsize=(14, 6))
for ax_idx, (geno, data_dict) in enumerate([('WT', wt_data), ('O2', o2_data)]):
    mat = np.zeros((len(key_rxns_short), len(DAP_LIST)))
    for j, dap in enumerate(DAP_LIST):
        df_g = data_dict[dap]
        for i, rxn in enumerate(key_rxns_short.keys()):
            if rxn in df_g.index:
                mat[i, j] = df_g.loc[rxn, 'max_rate']
    im = axes[ax_idx].imshow(mat, aspect='auto', cmap='YlOrRd')
    axes[ax_idx].set_xticks(range(len(DAP_LIST)))
    axes[ax_idx].set_xticklabels(['DAP %d' % d for d in DAP_LIST], rotation=45, ha='right', fontsize=9)
    axes[ax_idx].set_yticks(range(len(key_rxns_short)))
    axes[ax_idx].set_yticklabels(list(key_rxns_short.values()), fontsize=9)
    axes[ax_idx].set_title('max_rate: %s' % geno, fontsize=11, fontweight='bold')
    plt.colorbar(im, ax=axes[ax_idx], label='Max flux (mmol g DW⁻¹ h⁻¹)')
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat[i, j]
            axes[ax_idx].text(j, i, '%.3f' % v if v < 1 else '%.1f' % v,
                              ha='center', va='center', fontsize=7,
                              color='white' if v > mat.max() * 0.6 else 'black')

plt.suptitle('Flux Variability Analysis: Lysine Pathway\nMaximum Achievable Flux (100% Optimality)', fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join(FIG, 'FigFVA4_lysine_pathway_heatmap.png'), dpi=200, bbox_inches='tight')
plt.close()
print("  -> FigFVA4_lysine_pathway_heatmap.png")

# ── Fig E: LKR catabolism capacity comparison ────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

for ax_idx, rxn in enumerate(['R00715[K,c]', 'R00715[K,p]']):
    ax = axes[ax_idx]
    wt_v = [wt_data[d].loc[rxn, 'max_rate'] if rxn in wt_data[d].index else 0 for d in DAP_LIST]
    o2_v = [o2_data[d].loc[rxn, 'max_rate'] if rxn in o2_data[d].index else 0 for d in DAP_LIST]
    ax.plot(DAP_LIST, wt_v, 'o-', color=COLORS['WT'], linewidth=2, markersize=7, label='Wild Type')
    ax.plot(DAP_LIST, o2_v, 's--', color=COLORS['O2'], linewidth=2, markersize=7, label='o2 mutant')
    ax.set_xlabel('Days After Pollination'); ax.set_ylabel('Max LKR Flux (mmol g DW⁻¹ h⁻¹)')
    compartment = 'Cytosolic' if '[K,c]' in rxn else 'Plastidic'
    ax.set_title('LKR Catabolic Capacity (%s)\n%s' % (compartment, rxn), fontsize=10, fontweight='bold')
    ax.legend(); ax.grid(alpha=0.3)
    ax.set_xticks(DAP_LIST)

plt.suptitle('Reduced Lysine Catabolism Capacity in o2 Mutant\n(FVA Maximum Rates at 100% Optimality)',
             fontsize=11, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join(FIG, 'FigFVA5_LKR_capacity.png'), dpi=200, bbox_inches='tight')
plt.close()
print("  -> FigFVA5_LKR_capacity.png")

# ── Print key numbers for manuscript ─────────────────────────────────────────
print()
print("=== KEY NUMBERS FOR MANUSCRIPT ===")
print("\nR00451 (DAP decarboxylase) max achievable flux:")
for d in DAP_LIST:
    wt_v = wt_data[d].loc['R00451[K,p]', 'max_rate'] if 'R00451[K,p]' in wt_data[d].index else 0
    o2_v = o2_data[d].loc['R00451[K,p]', 'max_rate'] if 'R00451[K,p]' in o2_data[d].index else 0
    ratio = o2_v / wt_v if wt_v > 1e-10 else float('inf')
    print("  DAP %2d: WT_max=%.5f  O2_max=%.5f  ratio=%.2f" % (d, wt_v, o2_v, ratio))

print("\ncmTransport_C00047 (lysine transport) max flux:")
for d in DAP_LIST:
    wt_v = wt_data[d].loc['cmTransport_C00047[K]', 'max_rate'] if 'cmTransport_C00047[K]' in wt_data[d].index else 0
    o2_v = o2_data[d].loc['cmTransport_C00047[K]', 'max_rate'] if 'cmTransport_C00047[K]' in o2_data[d].index else 0
    print("  DAP %2d: WT_max=%.5f  O2_max=%.5f" % (d, wt_v, o2_v))

print("\nNewly constrained reactions in O2 (reactions that lose flexibility):")
print(df_ess[['DAP','Newly_fixed_in_O2','Newly_flexible_in_O2']].to_string(index=False))

print("\nDone.")
