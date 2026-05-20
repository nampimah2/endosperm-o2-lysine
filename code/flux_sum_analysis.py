#!/usr/bin/env python3
"""
==============================================================================
Flux Sum Analysis: O2 Mutant vs Wild Type Maize Endosperm
==============================================================================

Flux Sum (Φ_i) for metabolite i is defined as:
    Φ_i = 0.5 * Σ_j |S(i,j) * v(j)|

This equals the total production (= total consumption at steady state) flux
through each metabolite, quantifying its turnover rate.

By comparing flux sums between O2 mutant and Wild Type across developmental
stages, we identify metabolites with significantly altered turnover that may
represent targets for increasing lysine accumulation in wild type.

Analyses performed:
  1. Compute flux sums for all metabolites (WT and O2, all DAPs)
  2. Differential flux sum analysis (O2 − WT)
  3. Rank metabolites by flux sum change
  4. Pathway-level flux sum analysis
  5. Identify metabolite targets for lysine enhancement
  6. Generate publication-quality figures

Reference:
  Kim & Reed (2012) "Refining metabolic models and accounting for regulatory
  effects" - Flux sum as a measure of metabolite turnover.
==============================================================================
"""

import os
import re
import csv
import math
from collections import defaultdict

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

# ============================================================================
# CONFIGURATION
# ============================================================================
BASE_DIR = "/mnt/nrdstor/ssbio/nampimah/ENDOSPERM/2026/Endosperm_Model_MARY"
OUTPUT_DIR = os.path.join(BASE_DIR, "flux_sum_analysis")
FIGURE_DIR = os.path.join(OUTPUT_DIR, "figures")
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(FIGURE_DIR, exist_ok=True)

DAP_VALUES = [6, 8, 10, 12, 15, 18, 22, 30]
WT_DIRS = {dap: os.path.join(BASE_DIR, "Wild_Type", f"B{dap}") for dap in DAP_VALUES}
O2_DIRS = {dap: os.path.join(BASE_DIR, "O2_mutant", f"O{dap}") for dap in DAP_VALUES}

# Metabolite names (KEGG IDs)
METABOLITE_NAMES = {
    "C00047": "L-Lysine",
    "C00049": "L-Aspartate",
    "C00441": "L-Aspartate-4-semialdehyde",
    "C03340": "Dihydrodipicolinate",
    "C00680": "meso-Diaminopimelate",
    "C00449": "L-Saccharopine",
    "C04076": "L-2-Aminoadipate-6-semialdehyde",
    "C00956": "L-2-Aminoadipate",
    "C00022": "Pyruvate",
    "C00024": "Acetyl-CoA",
    "C00026": "2-Oxoglutarate",
    "C00158": "Citrate",
    "C00149": "Malate",
    "C00042": "Succinate",
    "C00036": "Oxaloacetate",
    "C00311": "Isocitrate",
    "C00122": "Fumarate",
    "C00074": "Phosphoenolpyruvate",
    "C00031": "D-Glucose",
    "C00089": "Sucrose",
    "C00267": "alpha-D-Glucose",
    "C00668": "alpha-D-Glucose-6P",
    "C00085": "D-Fructose-6P",
    "C00354": "D-Fructose-1,6-bisP",
    "C00118": "D-Glyceraldehyde-3P",
    "C00065": "L-Serine",
    "C00037": "Glycine",
    "C00041": "L-Alanine",
    "C00025": "L-Glutamate",
    "C00064": "L-Glutamine",
    "C00123": "L-Leucine",
    "C00183": "L-Valine",
    "C00407": "L-Isoleucine",
    "C00062": "L-Arginine",
    "C00078": "L-Tryptophan",
    "C00079": "L-Phenylalanine",
    "C00082": "L-Tyrosine",
    "C00073": "L-Methionine",
    "C00097": "L-Cysteine",
    "C00135": "L-Histidine",
    "C00148": "L-Proline",
    "C00188": "L-Threonine",
    "C00152": "L-Asparagine",
    "C00263": "L-Homoserine",
    "C00001": "H2O",
    "C00002": "ATP",
    "C00003": "NAD+",
    "C00004": "NADH",
    "C00005": "NADPH",
    "C00006": "NADP+",
    "C00008": "ADP",
    "C00009": "Orthophosphate",
    "C00013": "Diphosphate",
    "C00020": "AMP",
    "C00035": "GDP",
    "C00044": "GTP",
    "C00010": "CoA",
    "C00011": "CO2",
    "C00007": "Oxygen",
    "C00080": "H+",
    "C00095": "D-Fructose",
    "C00103": "D-Glucose-1P",
    "C00111": "Glycerone-phosphate",
    "C00197": "3-Phospho-D-glycerate",
    "C00236": "3-Phospho-D-glyceroyl-phosphate",
    "C01159": "2,3-Diphospho-D-glycerate",
    "C00631": "2-Phospho-D-glycerate",
    "C05345": "Beta-D-Fructose-6P",
    "C00229": "Acyl-carrier protein",
    "C00416": "Phosphatidate",
    "C00173": "Acyl-[acp]",
    "C00681": "1-Acyl-sn-glycerol-3P",
    "C00033": "Acetate",
    "C00048": "Glyoxylate",
    "C00050": "Glutathione-SH",
    "C00051": "Glutathione-SS",
    "C00127": "Glutathione-disulfide",
    "C00077": "L-Ornithine",
    "C00134": "Putrescine",
    "C00315": "Spermidine",
    "C00750": "Spermine",
}

# Pathway groupings for metabolites
METABOLITE_PATHWAYS = {
    "Lysine_Biosynthesis": ["C00049", "C00441", "C03340", "C00680", "C00047"],
    "Lysine_Degradation": ["C00047", "C00449", "C04076", "C00956"],
    "Aspartate_Family": ["C00049", "C00152", "C00263", "C00188", "C00073", "C00047"],
    "TCA_Cycle": ["C00022", "C00024", "C00158", "C00311", "C00026", "C00042", "C00122", "C00149", "C00036"],
    "Glycolysis": ["C00031", "C00267", "C00668", "C00085", "C00354", "C00118", "C00111",
                    "C00197", "C00236", "C00631", "C00074", "C00022"],
    "Energy_Cofactors": ["C00002", "C00008", "C00020", "C00003", "C00004", "C00005", "C00006"],
    "Amino_Acids": ["C00041", "C00025", "C00064", "C00065", "C00037", "C00123",
                     "C00183", "C00407", "C00062", "C00078", "C00079", "C00082",
                     "C00097", "C00135", "C00148"],
    "Polyamines": ["C00077", "C00134", "C00315", "C00750"],
}

# Cofactors / currency metabolites to optionally filter
CURRENCY_METABOLITES = {"C00001", "C00080", "C00007", "C00011"}  # H2O, H+, O2, CO2


# ============================================================================
# PARSING UTILITY FUNCTIONS
# ============================================================================
def parse_results(filepath):
    """Parse GAMS FBA results file → (flux_dict, objective_value)."""
    fluxes = {}
    obj_value = None
    if not os.path.exists(filepath):
        return None, None
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith("The max Biomass value is"):
                try:
                    obj_value = float(line.split(":")[-1].strip())
                except ValueError:
                    obj_value = None
            elif line and not line.startswith("The"):
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        fluxes[parts[0]] = float(parts[1])
                    except ValueError:
                        pass
    return fluxes, obj_value


def parse_sij(filepath):
    """Parse stoichiometric matrix → dict: (metabolite, reaction) → coefficient."""
    sij = {}
    if not os.path.exists(filepath):
        return sij
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line == '/' or line == '':
                continue
            match = re.match(r"'([^']+)'\.'([^']+)'\s+([-\d.eE+]+)", line)
            if match:
                met, rxn, coeff = match.group(1), match.group(2), float(match.group(3))
                sij[(met, rxn)] = coeff
    return sij


def get_metabolite_name(met_id):
    """Get human-readable name from full metabolite ID like 'C00047[K,c]'."""
    # Extract base KEGG ID
    match = re.match(r"(M?C\d+)", met_id)
    if match:
        base_id = match.group(1)
        name = METABOLITE_NAMES.get(base_id, base_id)
        # Extract compartment
        comp_match = re.search(r"\[K,(\w+)\]", met_id)
        compartment = comp_match.group(1) if comp_match else ""
        comp_map = {"c": "cyt", "p": "pla", "m": "mit", "v": "vac", "L": "ext"}
        comp_label = comp_map.get(compartment, compartment)
        return f"{name} [{comp_label}]"
    return met_id


def get_base_kegg_id(met_id):
    """Extract base KEGG ID (e.g. 'C00047') from full ID 'C00047[K,c]'."""
    match = re.match(r"(M?C\d+)", met_id)
    return match.group(1) if match else met_id


def get_compartment(met_id):
    """Extract compartment from metabolite ID."""
    match = re.search(r"\[K,(\w+)\]", met_id)
    return match.group(1) if match else ""


# ============================================================================
# FLUX SUM COMPUTATION
# ============================================================================
def compute_flux_sums(sij, fluxes):
    """
    Compute flux sum for every metabolite:
        Φ_i = 0.5 * Σ_j |S(i,j) * v(j)|

    Returns dict: metabolite_id → flux_sum_value
    """
    # Build metabolite → list of (reaction, stoich_coeff)
    met_reactions = defaultdict(list)
    for (met, rxn), coeff in sij.items():
        met_reactions[met].append((rxn, coeff))

    flux_sums = {}
    for met, rxn_list in met_reactions.items():
        total = 0.0
        for rxn, coeff in rxn_list:
            v = fluxes.get(rxn, 0.0)
            total += abs(coeff * v)
        flux_sums[met] = 0.5 * total

    return flux_sums


def compute_flux_sum_contributions(sij, fluxes, met_id):
    """
    Compute individual reaction contributions to a metabolite's flux sum.
    Returns list of (reaction, S(i,j), v(j), |S*v|, direction) sorted by |S*v|.
    """
    contributions = []
    for (met, rxn), coeff in sij.items():
        if met == met_id:
            v = fluxes.get(rxn, 0.0)
            sv = coeff * v
            direction = "producing" if sv > 0 else ("consuming" if sv < 0 else "inactive")
            contributions.append((rxn, coeff, v, abs(sv), direction))
    contributions.sort(key=lambda x: x[3], reverse=True)
    return contributions


# ============================================================================
# MAIN ANALYSIS
# ============================================================================
print("=" * 80)
print("FLUX SUM ANALYSIS: O2 MUTANT vs WILD TYPE MAIZE ENDOSPERM")
print("=" * 80)

# --- Step 1: Load all FBA results ---
print("\n--- Loading FBA results ---")
wt_fluxes = {}
o2_fluxes = {}
wt_biomass = {}
o2_biomass = {}

for dap in DAP_VALUES:
    wt_path = os.path.join(WT_DIRS[dap], "results_FBA.txt")
    o2_path = os.path.join(O2_DIRS[dap], "results_FBA.txt")

    f, obj = parse_results(wt_path)
    if f is not None:
        wt_fluxes[dap] = f
        wt_biomass[dap] = obj
        print(f"  WT DAP {dap:2d}: {len(f)} reactions, Biomass = {obj:.8f}")

    f, obj = parse_results(o2_path)
    if f is not None:
        o2_fluxes[dap] = f
        o2_biomass[dap] = obj
        print(f"  O2 DAP {dap:2d}: {len(f)} reactions, Biomass = {obj:.8f}")

# --- Step 2: Load stoichiometric matrices ---
print("\n--- Loading stoichiometric matrices ---")
wt_sij = {}
o2_sij = {}
for dap in DAP_VALUES:
    wt_sij[dap] = parse_sij(os.path.join(WT_DIRS[dap], "sij.txt"))
    o2_sij[dap] = parse_sij(os.path.join(O2_DIRS[dap], "sij.txt"))
    print(f"  DAP {dap:2d}: WT S(i,j) = {len(wt_sij[dap]):5d} entries, "
          f"O2 S(i,j) = {len(o2_sij[dap]):5d} entries")

# --- Step 3: Compute flux sums for all metabolites at all DAPs ---
print("\n--- Computing flux sums ---")
wt_flux_sums = {}
o2_flux_sums = {}

for dap in DAP_VALUES:
    wt_flux_sums[dap] = compute_flux_sums(wt_sij[dap], wt_fluxes[dap])
    o2_flux_sums[dap] = compute_flux_sums(o2_sij[dap], o2_fluxes[dap])
    print(f"  DAP {dap:2d}: WT = {len(wt_flux_sums[dap])} metabolites, "
          f"O2 = {len(o2_flux_sums[dap])} metabolites")

# Collect the union of all metabolites across all DAPs
all_mets = set()
for dap in DAP_VALUES:
    all_mets.update(wt_flux_sums[dap].keys())
    all_mets.update(o2_flux_sums[dap].keys())
print(f"\n  Total unique metabolites: {len(all_mets)}")


# ============================================================================
# OUTPUT 1: Complete Flux Sum Table (all metabolites × all DAPs)
# ============================================================================
print("\n" + "=" * 80)
print("OUTPUT 1: COMPLETE FLUX SUM TABLE")
print("=" * 80)

rows = []
for met in sorted(all_mets):
    row = {
        "Metabolite_ID": met,
        "Metabolite_Name": get_metabolite_name(met),
        "Base_KEGG_ID": get_base_kegg_id(met),
        "Compartment": get_compartment(met),
    }
    for dap in DAP_VALUES:
        wt_fs = wt_flux_sums[dap].get(met, 0.0)
        o2_fs = o2_flux_sums[dap].get(met, 0.0)
        diff = o2_fs - wt_fs
        row[f"WT_DAP{dap}"] = wt_fs
        row[f"O2_DAP{dap}"] = o2_fs
        row[f"Diff_DAP{dap}"] = diff
    rows.append(row)

df_all = pd.DataFrame(rows)

# Add summary statistics
for prefix in ["WT", "O2", "Diff"]:
    cols = [f"{prefix}_DAP{d}" for d in DAP_VALUES]
    df_all[f"{prefix}_Mean"] = df_all[cols].mean(axis=1)
    df_all[f"{prefix}_Max"] = df_all[cols].max(axis=1)
    if prefix == "Diff":
        df_all["Diff_MaxAbs"] = df_all[cols].abs().max(axis=1)
        df_all["Diff_SumAbs"] = df_all[cols].abs().sum(axis=1)

outfile1 = os.path.join(OUTPUT_DIR, "01_flux_sums_complete.csv")
df_all.to_csv(outfile1, index=False, float_format="%.10f")
print(f"  Saved: {outfile1}  ({len(df_all)} metabolites)")


# ============================================================================
# OUTPUT 2: Differential Flux Sum Analysis (excluding currency metabolites)
# ============================================================================
print("\n" + "=" * 80)
print("OUTPUT 2: DIFFERENTIAL FLUX SUM ANALYSIS")
print("=" * 80)

# Filter out currency metabolites for ranking
df_nocurrency = df_all[~df_all["Base_KEGG_ID"].isin(CURRENCY_METABOLITES)].copy()
df_nocurrency = df_nocurrency.sort_values("Diff_MaxAbs", ascending=False)

outfile2 = os.path.join(OUTPUT_DIR, "02_differential_flux_sums_ranked.csv")
df_nocurrency.to_csv(outfile2, index=False, float_format="%.10f")
print(f"  Saved: {outfile2}  ({len(df_nocurrency)} metabolites, currency excluded)")

# Print top 30
print(f"\n  Top 30 metabolites by max |ΔΦ| (excluding H2O, H+, O2, CO2):")
print(f"  {'Rank':>4s}  {'Metabolite':<40s}  {'MaxAbs_Diff':>12s}  {'Mean_Diff':>12s}  {'Direction':>10s}")
print(f"  {'-'*4}  {'-'*40}  {'-'*12}  {'-'*12}  {'-'*10}")
for i, (_, row) in enumerate(df_nocurrency.head(30).iterrows()):
    name = row["Metabolite_Name"][:40]
    direction = "INCREASED" if row["Diff_Mean"] > 0 else "DECREASED"
    print(f"  {i+1:4d}  {name:<40s}  {row['Diff_MaxAbs']:12.6f}  {row['Diff_Mean']:+12.6f}  {direction:>10s}")


# ============================================================================
# OUTPUT 3: Lysine-Related Metabolite Flux Sums (detailed)
# ============================================================================
print("\n" + "=" * 80)
print("OUTPUT 3: LYSINE-RELATED METABOLITE FLUX SUMS")
print("=" * 80)

lysine_kegg_ids = ["C00049", "C00441", "C03340", "C00680", "C00047",
                   "C00449", "C04076", "C00956", "C00263", "C00188",
                   "C00073", "C00152"]

df_lysine = df_all[df_all["Base_KEGG_ID"].isin(lysine_kegg_ids)].copy()
df_lysine = df_lysine.sort_values("Diff_MaxAbs", ascending=False)

outfile3 = os.path.join(OUTPUT_DIR, "03_lysine_pathway_flux_sums.csv")
df_lysine.to_csv(outfile3, index=False, float_format="%.10f")
print(f"  Saved: {outfile3}  ({len(df_lysine)} metabolites)")

for _, row in df_lysine.iterrows():
    print(f"\n  {row['Metabolite_Name']} ({row['Metabolite_ID']})")
    for dap in DAP_VALUES:
        wt = row[f"WT_DAP{dap}"]
        o2 = row[f"O2_DAP{dap}"]
        diff = row[f"Diff_DAP{dap}"]
        marker = " ***" if abs(diff) > 0.001 else ""
        print(f"    DAP {dap:2d}: WT={wt:10.6f}  O2={o2:10.6f}  Diff={diff:+10.6f}{marker}")


# ============================================================================
# OUTPUT 4: Reaction Contributions to Key Metabolite Flux Sums
# ============================================================================
print("\n" + "=" * 80)
print("OUTPUT 4: REACTION CONTRIBUTIONS TO LYSINE FLUX SUM")
print("=" * 80)

# Focus on lysine (C00047) in all compartments
lysine_mets = [m for m in all_mets if "C00047" in m]
print(f"  Lysine metabolite instances: {lysine_mets}")

contribution_rows = []
for dap in DAP_VALUES:
    for met in lysine_mets:
        wt_contribs = compute_flux_sum_contributions(wt_sij[dap], wt_fluxes[dap], met)
        o2_contribs = compute_flux_sum_contributions(o2_sij[dap], o2_fluxes[dap], met)

        # Merge WT and O2 contributions
        all_rxns = set()
        wt_map = {}
        o2_map = {}
        for rxn, coeff, v, sv, direction in wt_contribs:
            wt_map[rxn] = (coeff, v, sv, direction)
            all_rxns.add(rxn)
        for rxn, coeff, v, sv, direction in o2_contribs:
            o2_map[rxn] = (coeff, v, sv, direction)
            all_rxns.add(rxn)

        for rxn in sorted(all_rxns):
            wt_info = wt_map.get(rxn, (0, 0, 0, "absent"))
            o2_info = o2_map.get(rxn, (0, 0, 0, "absent"))
            contribution_rows.append({
                "DAP": dap,
                "Metabolite": met,
                "Reaction": rxn,
                "Stoich_Coeff": wt_info[0] if wt_info[0] != 0 else o2_info[0],
                "WT_Flux": wt_info[1],
                "WT_abs_Sv": wt_info[2],
                "WT_Direction": wt_info[3],
                "O2_Flux": o2_info[1],
                "O2_abs_Sv": o2_info[2],
                "O2_Direction": o2_info[3],
                "Diff_abs_Sv": o2_info[2] - wt_info[2],
            })

df_contributions = pd.DataFrame(contribution_rows)
outfile4 = os.path.join(OUTPUT_DIR, "04_lysine_flux_sum_contributions.csv")
df_contributions.to_csv(outfile4, index=False, float_format="%.10f")
print(f"  Saved: {outfile4}  ({len(df_contributions)} entries)")


# ============================================================================
# OUTPUT 5: Pathway-Level Flux Sum Comparison
# ============================================================================
print("\n" + "=" * 80)
print("OUTPUT 5: PATHWAY-LEVEL FLUX SUM COMPARISON")
print("=" * 80)

pathway_rows = []
for pathway_name, kegg_ids in METABOLITE_PATHWAYS.items():
    for dap in DAP_VALUES:
        wt_pathway_total = 0.0
        o2_pathway_total = 0.0
        met_count = 0
        for met in all_mets:
            base = get_base_kegg_id(met)
            if base in kegg_ids:
                wt_fs = wt_flux_sums[dap].get(met, 0.0)
                o2_fs = o2_flux_sums[dap].get(met, 0.0)
                wt_pathway_total += wt_fs
                o2_pathway_total += o2_fs
                met_count += 1

        diff = o2_pathway_total - wt_pathway_total
        fold_change = o2_pathway_total / wt_pathway_total if wt_pathway_total > 1e-12 else float('nan')
        pathway_rows.append({
            "Pathway": pathway_name,
            "DAP": dap,
            "Num_Metabolites": met_count,
            "WT_FluxSum_Total": wt_pathway_total,
            "O2_FluxSum_Total": o2_pathway_total,
            "Difference": diff,
            "Fold_Change": fold_change,
        })

df_pathways = pd.DataFrame(pathway_rows)
outfile5 = os.path.join(OUTPUT_DIR, "05_pathway_flux_sum_comparison.csv")
df_pathways.to_csv(outfile5, index=False, float_format="%.8f")
print(f"  Saved: {outfile5}")

# Print pathway summary
for pathway in METABOLITE_PATHWAYS.keys():
    sub = df_pathways[df_pathways["Pathway"] == pathway]
    avg_diff = sub["Difference"].mean()
    direction = "HIGHER in O2" if avg_diff > 0 else "LOWER in O2"
    print(f"  {pathway:<25s}: Avg ΔΦ = {avg_diff:+10.6f} ({direction})")


# ============================================================================
# OUTPUT 6: Metabolite Target Identification
# ============================================================================
print("\n" + "=" * 80)
print("OUTPUT 6: METABOLITE TARGET IDENTIFICATION FOR LYSINE ENHANCEMENT")
print("=" * 80)

# Strategy: Identify metabolites where O2 has consistently different flux sums
# and that are biologically connected to lysine metabolism.

# Score each metabolite by:
# 1. Consistency of direction across DAPs
# 2. Magnitude of change
# 3. Biological relevance (closer to lysine pathway = higher weight)

target_rows = []
for _, row in df_nocurrency.iterrows():
    met_id = row["Metabolite_ID"]
    base_kegg = row["Base_KEGG_ID"]
    name = row["Metabolite_Name"]
    compartment = row["Compartment"]

    diffs = [row[f"Diff_DAP{d}"] for d in DAP_VALUES]
    non_zero_diffs = [d for d in diffs if abs(d) > 1e-10]

    if len(non_zero_diffs) == 0:
        continue

    # Direction consistency
    n_positive = sum(1 for d in non_zero_diffs if d > 0)
    n_negative = sum(1 for d in non_zero_diffs if d < 0)
    consistency = max(n_positive, n_negative) / len(non_zero_diffs) if non_zero_diffs else 0

    # Predominant direction
    predominant_dir = "INCREASED" if n_positive > n_negative else "DECREASED"

    # Mean absolute difference
    mean_abs_diff = sum(abs(d) for d in diffs) / len(diffs)
    max_abs_diff = max(abs(d) for d in diffs)

    # Biological relevance score (heuristic)
    relevance = 1.0
    lysine_related = ["C00049", "C00441", "C03340", "C00680", "C00047",
                      "C00449", "C04076", "C00956", "C00263", "C00188", "C00073"]
    tca_related = ["C00022", "C00024", "C00158", "C00311", "C00026",
                   "C00042", "C00122", "C00149", "C00036"]
    aa_biosynthesis = ["C00025", "C00064", "C00152", "C00062", "C00065",
                       "C00037", "C00041", "C00123", "C00183", "C00407",
                       "C00078", "C00079", "C00082", "C00097", "C00135", "C00148"]
    energy = ["C00002", "C00008", "C00020", "C00003", "C00004", "C00005", "C00006"]
    carbon = ["C00031", "C00089", "C00267", "C00668", "C00085", "C00354",
              "C00118", "C00074"]

    if base_kegg in lysine_related:
        relevance = 5.0
    elif base_kegg in tca_related:
        relevance = 3.0
    elif base_kegg in aa_biosynthesis:
        relevance = 2.5
    elif base_kegg in energy:
        relevance = 2.0
    elif base_kegg in carbon:
        relevance = 2.0

    # Composite score
    composite_score = mean_abs_diff * consistency * relevance

    # Determine target strategy
    if predominant_dir == "INCREASED":
        strategy = "Enhance turnover (increase supply/demand)"
    else:
        strategy = "Reduce turnover (limit competing usage)"

    # Pathway membership
    pathways_found = []
    for pw_name, pw_keggs in METABOLITE_PATHWAYS.items():
        if base_kegg in pw_keggs:
            pathways_found.append(pw_name)

    target_rows.append({
        "Rank": 0,
        "Metabolite_ID": met_id,
        "Metabolite_Name": name,
        "Base_KEGG_ID": base_kegg,
        "Compartment": compartment,
        "Direction_in_O2": predominant_dir,
        "Consistency": consistency,
        "DAPs_with_change": len(non_zero_diffs),
        "Mean_Abs_Diff": mean_abs_diff,
        "Max_Abs_Diff": max_abs_diff,
        "Relevance_Weight": relevance,
        "Composite_Score": composite_score,
        "Suggested_Strategy": strategy,
        "Pathways": "; ".join(pathways_found) if pathways_found else "Other",
    })
    # Add individual DAP values
    for dap in DAP_VALUES:
        target_rows[-1][f"Diff_DAP{dap}"] = row[f"Diff_DAP{dap}"]

df_targets = pd.DataFrame(target_rows)
df_targets = df_targets.sort_values("Composite_Score", ascending=False).reset_index(drop=True)
df_targets["Rank"] = df_targets.index + 1

outfile6 = os.path.join(OUTPUT_DIR, "06_metabolite_targets_ranked.csv")
df_targets.to_csv(outfile6, index=False, float_format="%.8f")
print(f"  Saved: {outfile6}  ({len(df_targets)} metabolite targets)")

# Print top 30 targets
print(f"\n  Top 30 Metabolite Targets for Lysine Enhancement in WT:")
print(f"  {'Rank':>4s}  {'Metabolite':<35s}  {'Comp':>4s}  {'Direction':>10s}  {'Score':>10s}  {'Strategy':<40s}")
print(f"  {'-'*4}  {'-'*35}  {'-'*4}  {'-'*10}  {'-'*10}  {'-'*40}")
for _, row in df_targets.head(30).iterrows():
    name = (row["Metabolite_Name"])[:35]
    print(f"  {int(row['Rank']):4d}  {name:<35s}  {row['Compartment']:>4s}  "
          f"{row['Direction_in_O2']:>10s}  {row['Composite_Score']:10.6f}  "
          f"{row['Suggested_Strategy']:<40s}")


# ============================================================================
# OUTPUT 7: Summary Statistics
# ============================================================================
print("\n" + "=" * 80)
print("OUTPUT 7: SUMMARY STATISTICS")
print("=" * 80)

summary_file = os.path.join(OUTPUT_DIR, "07_summary_statistics.txt")
with open(summary_file, 'w') as f:
    f.write("FLUX SUM ANALYSIS: SUMMARY STATISTICS\n")
    f.write("=" * 70 + "\n\n")

    for dap in DAP_VALUES:
        n_mets = len(wt_flux_sums[dap])
        wt_total = sum(wt_flux_sums[dap].values())
        o2_total = sum(o2_flux_sums[dap].values())

        # Count metabolites with increased/decreased flux sum in O2
        n_increased = 0
        n_decreased = 0
        n_unchanged = 0
        for met in all_mets:
            wt_fs = wt_flux_sums[dap].get(met, 0.0)
            o2_fs = o2_flux_sums[dap].get(met, 0.0)
            if o2_fs - wt_fs > 1e-10:
                n_increased += 1
            elif wt_fs - o2_fs > 1e-10:
                n_decreased += 1
            else:
                n_unchanged += 1

        line = (f"DAP {dap:2d}:  Metabolites={n_mets}  "
                f"WT_total_Φ={wt_total:.4f}  O2_total_Φ={o2_total:.4f}  "
                f"Increased={n_increased}  Decreased={n_decreased}  "
                f"Unchanged={n_unchanged}\n")
        f.write(line)
        print(f"  {line.strip()}")

    f.write("\n\nTOP 20 METABOLITE TARGETS:\n")
    f.write("-" * 70 + "\n")
    for _, row in df_targets.head(20).iterrows():
        f.write(f"  {int(row['Rank']):3d}. {row['Metabolite_Name']:<35s} "
                f"[{row['Compartment']}]  "
                f"{row['Direction_in_O2']:>10s}  Score={row['Composite_Score']:.6f}\n")
        f.write(f"       Strategy: {row['Suggested_Strategy']}\n")
        f.write(f"       Pathways: {row['Pathways']}\n")

print(f"  Saved: {summary_file}")


# ============================================================================
# FIGURES
# ============================================================================
print("\n" + "=" * 80)
print("GENERATING FIGURES")
print("=" * 80)

# Color scheme
WT_COLOR = "#2166ac"
O2_COLOR = "#b2182b"

# --- Figure 1: Top 20 Metabolites by Max Absolute Flux Sum Difference ---
fig, ax = plt.subplots(figsize=(12, 8))
top20 = df_nocurrency.head(20).copy()
top20 = top20.iloc[::-1]  # Reverse for horizontal bar
y_labels = [f"{row['Metabolite_Name']}" for _, row in top20.iterrows()]
y_pos = range(len(y_labels))
bars = ax.barh(y_pos, top20["Diff_Mean"].values, color=[O2_COLOR if v > 0 else WT_COLOR for v in top20["Diff_Mean"].values],
               edgecolor='black', linewidth=0.5)
ax.set_yticks(y_pos)
ax.set_yticklabels(y_labels, fontsize=9)
ax.set_xlabel("Mean Flux Sum Difference (O2 − WT)", fontsize=11)
ax.set_title("Top 20 Metabolites by Flux Sum Change (O2 vs WT)", fontsize=13, fontweight='bold')
ax.axvline(x=0, color='black', linestyle='-', linewidth=0.8)
ax.legend(handles=[plt.Rectangle((0,0),1,1, fc=O2_COLOR, label="Higher in O2"),
                    plt.Rectangle((0,0),1,1, fc=WT_COLOR, label="Higher in WT")],
          loc='lower right', fontsize=9)
plt.tight_layout()
fig.savefig(os.path.join(FIGURE_DIR, "fig1_top20_flux_sum_diff.png"), dpi=300, bbox_inches='tight')
plt.close(fig)
print("  Saved: fig1_top20_flux_sum_diff.png")

# --- Figure 2: Temporal Flux Sum Profiles for Lysine Pathway Metabolites ---
lysine_plot_keggs = ["C00049", "C00441", "C00047", "C00449", "C00680"]
lysine_plot_mets = {}
for kegg_id in lysine_plot_keggs:
    for met in sorted(all_mets):
        if kegg_id in met and get_compartment(met) in ("c", "p"):
            if kegg_id not in lysine_plot_mets:
                lysine_plot_mets[kegg_id] = []
            lysine_plot_mets[kegg_id].append(met)

n_plots = sum(len(v) for v in lysine_plot_mets.values())
fig, axes = plt.subplots(n_plots, 1, figsize=(10, 3.5 * n_plots), sharex=True)
if n_plots == 1:
    axes = [axes]

idx = 0
for kegg_id in lysine_plot_keggs:
    for met in lysine_plot_mets.get(kegg_id, []):
        ax = axes[idx]
        wt_vals = [wt_flux_sums[d].get(met, 0) for d in DAP_VALUES]
        o2_vals = [o2_flux_sums[d].get(met, 0) for d in DAP_VALUES]
        ax.plot(DAP_VALUES, wt_vals, 'o-', color=WT_COLOR, linewidth=2, markersize=6, label="Wild Type")
        ax.plot(DAP_VALUES, o2_vals, 's-', color=O2_COLOR, linewidth=2, markersize=6, label="O2 Mutant")
        ax.set_ylabel("Flux Sum (Φ)", fontsize=10)
        ax.set_title(f"{get_metabolite_name(met)}", fontsize=11, fontweight='bold')
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        idx += 1

axes[-1].set_xlabel("Days After Pollination (DAP)", fontsize=11)
axes[-1].set_xticks(DAP_VALUES)
plt.tight_layout()
fig.savefig(os.path.join(FIGURE_DIR, "fig2_lysine_pathway_flux_sums_temporal.png"), dpi=300, bbox_inches='tight')
plt.close(fig)
print("  Saved: fig2_lysine_pathway_flux_sums_temporal.png")

# --- Figure 3: Pathway-Level Flux Sum Comparison (heatmap-style) ---
pathway_names = list(METABOLITE_PATHWAYS.keys())
fig, ax = plt.subplots(figsize=(12, 6))
bar_width = 0.35
x_positions = np.arange(len(DAP_VALUES))

# Pick a few key pathways
key_pathways = ["Lysine_Biosynthesis", "Lysine_Degradation", "TCA_Cycle", "Aspartate_Family", "Glycolysis"]
fig, axes = plt.subplots(len(key_pathways), 1, figsize=(10, 3 * len(key_pathways)), sharex=True)
for pidx, pw_name in enumerate(key_pathways):
    ax = axes[pidx]
    pw_data = df_pathways[df_pathways["Pathway"] == pw_name]
    wt_totals = [pw_data[pw_data["DAP"] == d]["WT_FluxSum_Total"].values[0] for d in DAP_VALUES]
    o2_totals = [pw_data[pw_data["DAP"] == d]["O2_FluxSum_Total"].values[0] for d in DAP_VALUES]

    ax.bar(x_positions - bar_width/2, wt_totals, bar_width, color=WT_COLOR, label="Wild Type", edgecolor='black', linewidth=0.5)
    ax.bar(x_positions + bar_width/2, o2_totals, bar_width, color=O2_COLOR, label="O2 Mutant", edgecolor='black', linewidth=0.5)
    ax.set_ylabel("Σ Flux Sum", fontsize=10)
    ax.set_title(pw_name.replace("_", " "), fontsize=11, fontweight='bold')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, axis='y')

axes[-1].set_xticks(x_positions)
axes[-1].set_xticklabels([f"DAP {d}" for d in DAP_VALUES], fontsize=10)
axes[-1].set_xlabel("Developmental Stage", fontsize=11)
plt.tight_layout()
fig.savefig(os.path.join(FIGURE_DIR, "fig3_pathway_flux_sum_comparison.png"), dpi=300, bbox_inches='tight')
plt.close(fig)
print("  Saved: fig3_pathway_flux_sum_comparison.png")

# --- Figure 4: Heatmap of Flux Sum Differences ---
# Top 30 non-currency metabolites
top30_for_heatmap = df_nocurrency.head(30)
diff_cols = [f"Diff_DAP{d}" for d in DAP_VALUES]
heatmap_data = top30_for_heatmap[diff_cols].values
y_labels_hm = [f"{row['Metabolite_Name']}" for _, row in top30_for_heatmap.iterrows()]

fig, ax = plt.subplots(figsize=(10, 10))
vmax = np.max(np.abs(heatmap_data))
im = ax.imshow(heatmap_data, aspect='auto', cmap='RdBu_r', vmin=-vmax, vmax=vmax)
ax.set_xticks(range(len(DAP_VALUES)))
ax.set_xticklabels([f"DAP {d}" for d in DAP_VALUES], fontsize=10)
ax.set_yticks(range(len(y_labels_hm)))
ax.set_yticklabels(y_labels_hm, fontsize=8)
ax.set_title("Flux Sum Difference (O2 − WT) Across Development", fontsize=13, fontweight='bold')
ax.set_xlabel("Developmental Stage", fontsize=11)
cbar = plt.colorbar(im, ax=ax, shrink=0.8)
cbar.set_label("ΔΦ (O2 − WT)", fontsize=10)
plt.tight_layout()
fig.savefig(os.path.join(FIGURE_DIR, "fig4_flux_sum_heatmap.png"), dpi=300, bbox_inches='tight')
plt.close(fig)
print("  Saved: fig4_flux_sum_heatmap.png")

# --- Figure 5: Metabolite Target Scores ---
top20_targets = df_targets.head(20).copy()
top20_targets = top20_targets.iloc[::-1]
fig, ax = plt.subplots(figsize=(10, 8))
colors = [O2_COLOR if row["Direction_in_O2"] == "INCREASED" else WT_COLOR for _, row in top20_targets.iterrows()]
bars = ax.barh(range(len(top20_targets)),
               top20_targets["Composite_Score"].values,
               color=colors, edgecolor='black', linewidth=0.5)
ax.set_yticks(range(len(top20_targets)))
ax.set_yticklabels([f"{row['Metabolite_Name']}" for _, row in top20_targets.iterrows()], fontsize=9)
ax.set_xlabel("Composite Target Score", fontsize=11)
ax.set_title("Top 20 Metabolite Targets for Lysine Enhancement", fontsize=13, fontweight='bold')
ax.legend(handles=[plt.Rectangle((0,0),1,1, fc=O2_COLOR, label="Turnover ↑ in O2"),
                    plt.Rectangle((0,0),1,1, fc=WT_COLOR, label="Turnover ↓ in O2")],
          loc='lower right', fontsize=9)
plt.tight_layout()
fig.savefig(os.path.join(FIGURE_DIR, "fig5_metabolite_target_scores.png"), dpi=300, bbox_inches='tight')
plt.close(fig)
print("  Saved: fig5_metabolite_target_scores.png")

# --- Figure 6: Lysine Flux Sum vs Biomass Scatter ---
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# WT
lys_mets_cyt = [m for m in all_mets if "C00047" in m and get_compartment(m) == "c"]
if lys_mets_cyt:
    lys_met = lys_mets_cyt[0]
    wt_lys_fs = [wt_flux_sums[d].get(lys_met, 0) for d in DAP_VALUES]
    o2_lys_fs = [o2_flux_sums[d].get(lys_met, 0) for d in DAP_VALUES]
    wt_bm = [wt_biomass.get(d, 0) for d in DAP_VALUES]
    o2_bm = [o2_biomass.get(d, 0) for d in DAP_VALUES]

    axes[0].scatter(wt_bm, wt_lys_fs, c=WT_COLOR, s=80, edgecolors='black', zorder=5)
    for i, dap in enumerate(DAP_VALUES):
        axes[0].annotate(f"DAP{dap}", (wt_bm[i], wt_lys_fs[i]), fontsize=8, ha='left', va='bottom')
    axes[0].set_xlabel("Biomass Flux", fontsize=11)
    axes[0].set_ylabel("Lysine Flux Sum (Φ)", fontsize=11)
    axes[0].set_title("Wild Type", fontsize=12, fontweight='bold')
    axes[0].grid(True, alpha=0.3)

    axes[1].scatter(o2_bm, o2_lys_fs, c=O2_COLOR, s=80, edgecolors='black', zorder=5)
    for i, dap in enumerate(DAP_VALUES):
        axes[1].annotate(f"DAP{dap}", (o2_bm[i], o2_lys_fs[i]), fontsize=8, ha='left', va='bottom')
    axes[1].set_xlabel("Biomass Flux", fontsize=11)
    axes[1].set_ylabel("Lysine Flux Sum (Φ)", fontsize=11)
    axes[1].set_title("O2 Mutant", fontsize=12, fontweight='bold')
    axes[1].grid(True, alpha=0.3)

plt.suptitle("Lysine Turnover vs Biomass Production", fontsize=13, fontweight='bold')
plt.tight_layout()
fig.savefig(os.path.join(FIGURE_DIR, "fig6_lysine_fluxsum_vs_biomass.png"), dpi=300, bbox_inches='tight')
plt.close(fig)
print("  Saved: fig6_lysine_fluxsum_vs_biomass.png")


# ============================================================================
# FINAL INTERPRETATION
# ============================================================================
print("\n" + "=" * 80)
print("INTERPRETATION & RECOMMENDATIONS")
print("=" * 80)

# Identify consistent patterns
increased_in_o2 = df_targets[
    (df_targets["Direction_in_O2"] == "INCREASED") &
    (df_targets["Consistency"] >= 0.6) &
    (df_targets["Mean_Abs_Diff"] > 1e-6)
].head(15)

decreased_in_o2 = df_targets[
    (df_targets["Direction_in_O2"] == "DECREASED") &
    (df_targets["Consistency"] >= 0.6) &
    (df_targets["Mean_Abs_Diff"] > 1e-6)
].head(15)

print("\n  Metabolites with CONSISTENTLY HIGHER turnover in O2 mutant:")
print("  (These suggest increased metabolic activity that may support lysine)")
for _, row in increased_in_o2.iterrows():
    print(f"    - {row['Metabolite_Name']:<35s} [{row['Compartment']}]  "
          f"Mean ΔΦ = {row['Mean_Abs_Diff']:+.6f}  Pathways: {row['Pathways']}")

print("\n  Metabolites with CONSISTENTLY LOWER turnover in O2 mutant:")
print("  (These suggest reduced competing pathways in O2)")
for _, row in decreased_in_o2.iterrows():
    print(f"    - {row['Metabolite_Name']:<35s} [{row['Compartment']}]  "
          f"Mean ΔΦ = {row['Mean_Abs_Diff']:+.6f}  Pathways: {row['Pathways']}")

print("\n  POTENTIAL STRATEGIES to increase lysine in Wild Type:")
print("  Based on metabolic reprogramming observed in O2 mutant:")
print("  1. Increase turnover of metabolites that are higher in O2")
print("     (enhance supply of lysine precursors)")
print("  2. Reduce turnover of metabolites that are lower in O2")
print("     (limit flux through competing pathways)")
print("  3. Focus on metabolites with high composite scores")
print("     (high magnitude, consistent direction, pathway relevance)")

print(f"\n{'='*80}")
print(f"ANALYSIS COMPLETE. All outputs saved to: {OUTPUT_DIR}")
print(f"{'='*80}")
