#!/usr/bin/env python3
"""Compute manuscript metrics: Rewiring Index, Carbon Allocation Ratio, Flux Fraction to Lysine."""
import pandas as pd
import numpy as np
import csv, io, re

BASE = "/mnt/nrdstor/ssbio/nampimah/ENDOSPERM/2026/Endosperm_Model_MARY"
DAP_VALUES = [6, 8, 10, 12, 15, 18, 22, 30]

def read_rxn_csv(path):
    """Read a CSV where reaction names contain unquoted commas like R00307[K,c].
    Strategy: merge the first two columns (reaction prefix + compartment suffix)
    back into the reaction name, then treat remaining columns as numeric data.
    """
    with open(path) as f:
        reader = csv.reader(f)
        raw_header = next(reader)
        # The header row has no comma in field 0 ('Reaction'), so it parses normally
        # raw_header[0] = 'Reaction', raw_header[1] = 'WT_DAP6', etc.
        # But for data rows, field 0 = 'R00307[K', field 1 = 'c]', field 2 = WT_DAP6 value
        # Detect if the header looks misaligned by checking if raw_header[1] is numeric
        records = []
        header_merged = [raw_header[0]] + raw_header[1:]  # 'Reaction', 'WT_DAP6', ...
        for row in reader:
            if len(row) == 0:
                continue
            # Try to parse row[1] as a float — if it fails, the reaction name has a comma
            try:
                float(row[1])
                # row[1] is numeric → reaction name is just row[0], data starts at row[1]
                rxn_name = row[0]
                data = row[1:]
            except (ValueError, IndexError):
                # row[1] is 'c]' or 'p]' etc — merge first two fields as reaction name
                rxn_name = row[0] + ',' + row[1] if len(row) > 1 else row[0]
                data = row[2:]
            records.append([rxn_name] + data)
    # Build DataFrame with merged header (trim to min length)
    min_len = min(len(r) for r in records)
    records = [r[:min_len] for r in records]
    # Column count may differ from header; use header up to min_len
    cols = header_merged[:min_len]
    df = pd.DataFrame(records, columns=cols)
    # Convert numeric columns
    for c in df.columns[1:]:
        df[c] = pd.to_numeric(df[c], errors='coerce')
    return df

# ── 1. REWIRING INDEX ─────────────────────────────────────────────────────────
df = read_rxn_csv(f"{BASE}/analysis_results/06_differential_flux_all.csv")
ri = {}
for d in DAP_VALUES:
    ri[d] = df[f"Diff_DAP{d}"].abs().sum()

print("=== REWIRING INDEX  RI = sum(|v_O2 - v_WT|) ===")
for d, v in ri.items():
    print(f"  DAP {d:2d}: RI = {v:.4f}")

ri_df = pd.DataFrame({'DAP': list(ri.keys()), 'Rewiring_Index': list(ri.values())})
ri_df.to_csv(f"{BASE}/analysis_results/13_rewiring_index.csv", index=False)
print(f"  -> Saved to 13_rewiring_index.csv\n")

# ── 2. CARBON ALLOCATION TO LYSINE ────────────────────────────────────────────
lys = read_rxn_csv(f"{BASE}/analysis_results/02_lysine_pathway_fluxes.csv")
bm  = pd.read_csv(f"{BASE}/analysis_results/01_biomass_comparison.csv")

r00451 = lys[lys['Reaction'].str.contains('R00451', na=False)]
glc    = df[df['Reaction'].str.contains('Exchange_C00031', na=False)]

print(f"R00451 rows: {r00451['Reaction'].values}")
print(f"Glucose exchange rows: {glc['Reaction'].values}")
print()

rows = []
print("=== CARBON ALLOCATION TO LYSINE (% glucose carbon directed to Lys) ===")
print(f"{'DAP':>5}  {'WT_R00451':>12}  {'O2_R00451':>12}  {'WT_glc_flux':>12}  {'O2_glc_flux':>12}  {'WT_%C_Lys':>10}  {'O2_%C_Lys':>10}")
for d in DAP_VALUES:
    wt_lys = float(r00451[f'WT_DAP{d}'].values[0]) if not r00451.empty else np.nan
    o2_lys = float(r00451[f'O2_DAP{d}'].values[0]) if not r00451.empty else np.nan
    if not glc.empty and f'WT_DAP{d}' in glc.columns:
        wt_glc = abs(float(glc[f'WT_DAP{d}'].values[0]))
        o2_glc = abs(float(glc[f'O2_DAP{d}'].values[0]))
    else:
        wt_glc = np.nan; o2_glc = np.nan
    wt_pct = (wt_lys / wt_glc) * 100 if (not np.isnan(wt_glc) and wt_glc > 0) else np.nan
    o2_pct = (o2_lys / o2_glc) * 100 if (not np.isnan(o2_glc) and o2_glc > 0) else np.nan
    print(f"  {d:3d}  {wt_lys:12.6f}  {o2_lys:12.6f}  {wt_glc:12.4f}  {o2_glc:12.4f}  {wt_pct:10.4f}  {o2_pct:10.4f}")
    rows.append({'DAP': d, 'WT_R00451': wt_lys, 'O2_R00451': o2_lys,
                 'WT_glc': wt_glc, 'O2_glc': o2_glc, 'WT_pct_C_to_Lys': wt_pct,
                 'O2_pct_C_to_Lys': o2_pct})
print()

pd.DataFrame(rows).to_csv(f"{BASE}/analysis_results/14_carbon_allocation_to_lysine.csv", index=False)
print("  -> Saved to 14_carbon_allocation_to_lysine.csv\n")

# ── 3. FLUX FRACTION TO LYSINE from aspartate branch ──────────────────────────
r00480 = lys[lys['Reaction'].str.contains('R00480', na=False)]
print(f"R00480 (aspartate kinase) rows: {r00480['Reaction'].values}")
print()
print("=== FLUX FRACTION TO LYSINE (R00451 / R00480, aspartate-branch lysine fraction) ===")
print(f"{'DAP':>5}  {'WT_R00480':>12}  {'WT_R00451':>12}  {'WT_frac':>8}  {'O2_R00480':>12}  {'O2_R00451':>12}  {'O2_frac':>8}")
frac_rows = []
for d in DAP_VALUES:
    wt_kin = float(r00480[f'WT_DAP{d}'].values[0]) if not r00480.empty else np.nan
    o2_kin = float(r00480[f'O2_DAP{d}'].values[0]) if not r00480.empty else np.nan
    wt_lys = float(r00451[f'WT_DAP{d}'].values[0]) if not r00451.empty else np.nan
    o2_lys = float(r00451[f'O2_DAP{d}'].values[0]) if not r00451.empty else np.nan
    wt_frac = wt_lys / wt_kin if (not np.isnan(wt_kin) and wt_kin > 1e-12) else np.nan
    o2_frac = o2_lys / o2_kin if (not np.isnan(o2_kin) and o2_kin > 1e-12) else np.nan
    print(f"  {d:3d}  {wt_kin:12.6f}  {wt_lys:12.6f}  {wt_frac:8.4f}  {o2_kin:12.6f}  {o2_lys:12.6f}  {o2_frac:8.4f}")
    frac_rows.append({'DAP': d, 'WT_AK_flux': wt_kin, 'WT_DAPDC_flux': wt_lys,
                      'WT_lys_fraction': wt_frac, 'O2_AK_flux': o2_kin,
                      'O2_DAPDC_flux': o2_lys, 'O2_lys_fraction': o2_frac})
print()
pd.DataFrame(frac_rows).to_csv(f"{BASE}/analysis_results/15_aspartate_branch_lysine_fraction.csv", index=False)
print("  -> Saved to 15_aspartate_branch_lysine_fraction.csv\n")

print("=== ALL METRICS COMPUTED SUCCESSFULLY ===")
