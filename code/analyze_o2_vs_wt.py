#!/usr/bin/env python3
"""
==============================================================================
Comparative Metabolic Analysis: O2 Mutant vs Wild Type Maize Endosperm
Focus: Metabolic Reprogramming for Enhanced Lysine (C00047) Production
==============================================================================

This script performs a series of analyses on FBA results from context-specific
GAMS models across 8 developmental stages (DAP 6,8,10,12,15,18,22,30).

Analyses:
1. Biomass comparison across DAPs
2. Lysine pathway flux comparison (biosynthesis, degradation, transport)
3. Key amino acid pathway flux analysis
4. Flux Variability Analysis (FVA) for lysine reactions
5. Differential flux analysis (all reactions)
6. Pathway enrichment of differentially active reactions
7. Carbon flux redistribution analysis
8. Temporal dynamics visualization
==============================================================================
"""

import os
import re
import sys
import csv
import math
from collections import defaultdict

BASE_DIR = "/mnt/nrdstor/ssbio/nampimah/ENDOSPERM/2026/Endosperm_Model_MARY"
OUTPUT_DIR = os.path.join(BASE_DIR, "analysis_results")
os.makedirs(OUTPUT_DIR, exist_ok=True)

DAP_VALUES = [6, 8, 10, 12, 15, 18, 22, 30]
WT_DIRS = {dap: os.path.join(BASE_DIR, "Wild_Type", f"B{dap}") for dap in DAP_VALUES}
O2_DIRS = {dap: os.path.join(BASE_DIR, "O2_mutant", f"O{dap}") for dap in DAP_VALUES}

# ============================================================================
# Key metabolite KEGG IDs and names
# ============================================================================
METABOLITE_NAMES = {
    "C00047": "L-Lysine",
    "C00049": "L-Aspartate",
    "C00441": "L-Aspartate-4-semialdehyde",
    "C03340": "Dihydrodipicolinate",
    "C00680": "meso-Diaminopimelate",
    "C00449": "L-Saccharopine",
    "C04076": "L-2-Aminoadipate-6-semialdehyde",
    "C00956": "L-2-Aminoadipate",
    "C00322": "2-Oxoglutarate",  # also alpha-ketoglutarate
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
}

# ============================================================================
# Lysine pathway reactions with biological annotation
# ============================================================================
LYSINE_PATHWAY_REACTIONS = {
    # ---- LYSINE BIOSYNTHESIS (DAP pathway in plastid) ----
    "R00480[K,p]":   "Aspartate kinase (Asp → Asp-4-P) [Biosynthesis]",
    "R00355[K,p]":   "Aspartate aminotransferase [Biosynthesis]",
    "R01954[K,p]":   "Aspartate-semialdehyde dehydrogenase [Biosynthesis]",
    "R02291[K,p]":   "Asp-semialdehyde production [Biosynthesis]",
    "R02292[K,p]":   "Dihydrodipicolinate synthase (DHDPS) [Biosynthesis]",
    "R04198[K,p]":   "Dihydrodipicolinate reductase (DHDPR) [Biosynthesis]",
    "R02735[K,p]":   "meso-DAP production [Biosynthesis]",
    "R00451[K,p]":   "DAP decarboxylase (meso-DAP → Lysine) [Biosynthesis]",
    
    # ---- LYSINE DEGRADATION (Saccharopine pathway) ----
    "R00715[K,c]":   "Lysine-ketoglutarate reductase/LKR (Lys → Saccharopine) [Degradation-cytosol]",
    "R00715[K,p]":   "Lysine-ketoglutarate reductase/LKR (Lys → Saccharopine) [Degradation-plastid]",
    "R00716[K,c]":   "LKR alternative (Lys → Saccharopine) [Degradation-cytosol]",
    "R02313[K,c]":   "Saccharopine dehydrogenase/SDH (Saccharopine → AAS) [Degradation]",
    "R03102[K,c]":   "Aminoadipate-semialdehyde dehydrogenase [Degradation]",
    "R03658[K,m]":   "Lysine degradation (mitochondrial) [Degradation]",
    "R03658[K,p]":   "Lysine degradation (plastid) [Degradation]",
    
    # ---- TRANSPORT ----
    "cpTransport_C00047[K]":  "Lysine cytosol↔plastid transport",
    "ExB_C00047[K]":          "Lysine cytosol→boundary export",
    "Exchange_C00047[K,L]":   "Lysine exchange reaction",
    "cpTransport_C00449[K]":  "Saccharopine transport",
    
    # ---- BIOMASS ----
    "Seed_Biomass[K]":        "Seed Biomass reaction (lysine drain)",
}

# Broader set: amino acid biosynthesis, central carbon, TCA, etc.
PATHWAY_CATEGORIES = {
    "Lysine_Biosynthesis": ["R00480", "R00355", "R01954", "R02291", "R02292", "R04198", "R02735", "R00451"],
    "Lysine_Degradation": ["R00715", "R00716", "R02313", "R03102", "R03103", "R03658", "R02317"],
    "Aspartate_Pathway": ["R00355", "R00480", "R01954", "R00486", "MR00850", "R05577", "R00483"],
    "TCA_Cycle": ["R00351", "R01325", "R00342", "R00267", "R00405", "R00432", "R01082", "R00621"],
    "Glycolysis": ["R00299", "R01786", "R04779", "R01070", "R01015", "R01061", "R00200", "R00658"],
    "Amino_Acid_Transport": ["cpTransport_C00047", "ExB_C00047", "Exchange_C00047",
                              "ExB_C00049", "cpTransport_C00049"],
    "Zein_Storage_Protein": ["MR01119", "MR01138", "MR01151", "MR01188", "MR01196", "MR01279"],
    "Starch_Metabolism": ["R00292", "R00959", "R02421", "R00948", "R00300", "R06050", "R00010"],
    "Pentose_Phosphate": ["R01528", "R01529", "R01641", "R01056", "R01827", "R00835"],
}


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================
def parse_results(filepath):
    """Parse GAMS FBA results file into dict of reaction fluxes and objective value."""
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
    """Parse stoichiometric matrix file into dict: (metabolite, reaction) -> coefficient."""
    sij = {}
    if not os.path.exists(filepath):
        return sij
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line == '/' or line == '':
                continue
            # Format: 'metabolite'.'reaction'  coefficient
            match = re.match(r"'([^']+)'\.'([^']+)'\s+([-\d.eE+]+)", line)
            if match:
                met, rxn, coeff = match.group(1), match.group(2), float(match.group(3))
                sij[(met, rxn)] = coeff
    return sij


def parse_bounds(filepath):
    """Parse bounds file into dict: reaction -> value."""
    bounds = {}
    if not os.path.exists(filepath):
        return bounds
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line == '/' or line == '':
                continue
            match = re.match(r"'([^']+)'\s+([-\d.eE+]+)", line)
            if match:
                bounds[match.group(1)] = float(match.group(2))
    return bounds


def get_reactions_for_metabolite(sij, met_id_pattern):
    """Find all reactions involving a metabolite (any compartment)."""
    rxn_roles = {}  # rxn -> list of (met_full, coeff)
    for (met, rxn), coeff in sij.items():
        if met_id_pattern in met:
            if rxn not in rxn_roles:
                rxn_roles[rxn] = []
            rxn_roles[rxn].append((met, coeff))
    return rxn_roles


# ============================================================================
# LOAD ALL RESULTS
# ============================================================================
print("=" * 80)
print("LOADING FBA RESULTS")
print("=" * 80)

wt_results = {}
o2_results = {}
wt_biomass = {}
o2_biomass = {}

for dap in DAP_VALUES:
    wt_path = os.path.join(WT_DIRS[dap], "results_FBA.txt")
    o2_path = os.path.join(O2_DIRS[dap], "results_FBA.txt")
    
    wt_fluxes, wt_obj = parse_results(wt_path)
    o2_fluxes, o2_obj = parse_results(o2_path)
    
    if wt_fluxes is not None:
        wt_results[dap] = wt_fluxes
        wt_biomass[dap] = wt_obj
        print(f"  WT DAP {dap:2d}: Loaded {len(wt_fluxes)} reactions, Biomass = {wt_obj:.8f}")
    else:
        print(f"  WT DAP {dap:2d}: MISSING")
    
    if o2_fluxes is not None:
        o2_results[dap] = o2_fluxes
        o2_biomass[dap] = o2_obj
        print(f"  O2 DAP {dap:2d}: Loaded {len(o2_fluxes)} reactions, Biomass = {o2_obj:.8f}")
    else:
        print(f"  O2 DAP {dap:2d}: MISSING")


# ============================================================================
# ANALYSIS 1: BIOMASS COMPARISON
# ============================================================================
print("\n" + "=" * 80)
print("ANALYSIS 1: BIOMASS (Seed_Biomass) COMPARISON ACROSS DAPs")
print("=" * 80)

biomass_file = os.path.join(OUTPUT_DIR, "01_biomass_comparison.csv")
with open(biomass_file, 'w') as f:
    f.write("DAP,WT_Biomass,O2_Biomass,Difference,Fold_Change,Percent_Change\n")
    for dap in DAP_VALUES:
        wt_b = wt_biomass.get(dap, 0)
        o2_b = o2_biomass.get(dap, 0)
        diff = o2_b - wt_b
        fc = o2_b / wt_b if wt_b != 0 else float('inf')
        pct = (diff / wt_b * 100) if wt_b != 0 else float('inf')
        print(f"  DAP {dap:2d}: WT = {wt_b:.8f}  O2 = {o2_b:.8f}  Diff = {diff:+.8f}  FC = {fc:.4f}  ({pct:+.2f}%)")
        f.write(f"{dap},{wt_b:.8f},{o2_b:.8f},{diff:.8f},{fc:.6f},{pct:.4f}\n")
print(f"  Saved to: {biomass_file}")


# ============================================================================
# ANALYSIS 2: LYSINE PATHWAY FLUX ANALYSIS
# ============================================================================
print("\n" + "=" * 80)
print("ANALYSIS 2: LYSINE PATHWAY FLUX ANALYSIS")
print("=" * 80)

lysine_file = os.path.join(OUTPUT_DIR, "02_lysine_pathway_fluxes.csv")
with open(lysine_file, 'w') as f:
    header = "Reaction,Annotation"
    for dap in DAP_VALUES:
        header += f",WT_DAP{dap},O2_DAP{dap},Diff_DAP{dap}"
    f.write(header + "\n")
    
    for rxn, annotation in LYSINE_PATHWAY_REACTIONS.items():
        line = f"{rxn},{annotation}"
        for dap in DAP_VALUES:
            wt_v = wt_results.get(dap, {}).get(rxn, 0)
            o2_v = o2_results.get(dap, {}).get(rxn, 0)
            diff = o2_v - wt_v
            line += f",{wt_v:.8f},{o2_v:.8f},{diff:.8f}"
        f.write(line + "\n")

print(f"  Saved to: {lysine_file}")

# Print key lysine reactions summary
print("\n  --- Key Lysine Reactions at Each DAP ---")
key_rxns = [
    ("R00451[K,p]", "DAP decarboxylase (Lys biosynthesis)"),
    ("R00715[K,c]", "LKR (Lys degradation, cytosol)"),
    ("R00715[K,p]", "LKR (Lys degradation, plastid)"),
    ("R00716[K,c]", "LKR alt (Lys degradation)"),
    ("R02292[K,p]", "DHDPS (pathway entry)"),
    ("ExB_C00047[K]", "Lys export to boundary"),
    ("Exchange_C00047[K,L]", "Lys exchange"),
    ("Seed_Biomass[K]", "Biomass"),
]

for rxn, name in key_rxns:
    print(f"\n  {rxn} ({name}):")
    for dap in DAP_VALUES:
        wt_v = wt_results.get(dap, {}).get(rxn, 0)
        o2_v = o2_results.get(dap, {}).get(rxn, 0)
        diff = o2_v - wt_v
        arrow = "↑" if diff > 1e-10 else ("↓" if diff < -1e-10 else "=")
        print(f"    DAP {dap:2d}: WT={wt_v:12.8f}  O2={o2_v:12.8f}  Diff={diff:+12.8f} {arrow}")


# ============================================================================
# ANALYSIS 3: LYSINE NET PRODUCTION, CONSUMPTION & BALANCE
# ============================================================================
print("\n" + "=" * 80)
print("ANALYSIS 3: LYSINE NET FLUX BALANCE (production vs consumption)")
print("=" * 80)

# Load one sij to identify all lysine reactions
sij_wt = parse_sij(os.path.join(WT_DIRS[10], "sij.txt"))
lys_reactions = get_reactions_for_metabolite(sij_wt, "C00047")

balance_file = os.path.join(OUTPUT_DIR, "03_lysine_flux_balance.csv")
with open(balance_file, 'w') as f:
    f.write("DAP,WT_Lys_Production,WT_Lys_Consumption,WT_Net,O2_Lys_Production,O2_Lys_Consumption,O2_Net\n")
    
    for dap in DAP_VALUES:
        wt_prod, wt_cons = 0, 0
        o2_prod, o2_cons = 0, 0
        for rxn, met_list in lys_reactions.items():
            for met, coeff in met_list:
                wt_flux = wt_results.get(dap, {}).get(rxn, 0)
                o2_flux = o2_results.get(dap, {}).get(rxn, 0)
                wt_contribution = coeff * wt_flux
                o2_contribution = coeff * o2_flux
                if wt_contribution > 0:
                    wt_prod += wt_contribution
                else:
                    wt_cons += wt_contribution
                if o2_contribution > 0:
                    o2_prod += o2_contribution
                else:
                    o2_cons += o2_contribution
        
        wt_net = wt_prod + wt_cons
        o2_net = o2_prod + o2_cons
        print(f"  DAP {dap:2d}: WT prod={wt_prod:10.6f} cons={wt_cons:10.6f} net={wt_net:+10.6f}")
        print(f"          O2 prod={o2_prod:10.6f} cons={o2_cons:10.6f} net={o2_net:+10.6f}")
        f.write(f"{dap},{wt_prod:.8f},{wt_cons:.8f},{wt_net:.8f},{o2_prod:.8f},{o2_cons:.8f},{o2_net:.8f}\n")

print(f"  Saved to: {balance_file}")


# ============================================================================
# ANALYSIS 4: LYSINE BIOMASS COEFFICIENT COMPARISON
# ============================================================================
print("\n" + "=" * 80)
print("ANALYSIS 4: LYSINE STOICHIOMETRIC COEFFICIENT IN BIOMASS ACROSS DAPs")
print("=" * 80)

coeff_file = os.path.join(OUTPUT_DIR, "04_lysine_biomass_coefficient.csv")
with open(coeff_file, 'w') as f:
    f.write("DAP,WT_Lys_Coeff,O2_Lys_Coeff,Difference,WT_Asp_Coeff,O2_Asp_Coeff\n")
    for dap in DAP_VALUES:
        wt_sij = parse_sij(os.path.join(WT_DIRS[dap], "sij.txt"))
        o2_sij = parse_sij(os.path.join(O2_DIRS[dap], "sij.txt"))
        
        wt_lys_coeff = wt_sij.get(("C00047[K,B]", "Seed_Biomass[K]"), 0)
        o2_lys_coeff = o2_sij.get(("C00047[K,B]", "Seed_Biomass[K]"), 0)
        wt_asp_coeff = wt_sij.get(("C00049[K,B]", "Seed_Biomass[K]"), 0)
        o2_asp_coeff = o2_sij.get(("C00049[K,B]", "Seed_Biomass[K]"), 0)
        diff = o2_lys_coeff - wt_lys_coeff
        
        print(f"  DAP {dap:2d}: WT Lys coeff = {wt_lys_coeff:.6f}  O2 Lys coeff = {o2_lys_coeff:.6f}  "
              f"Diff = {diff:+.6f}  | WT Asp = {wt_asp_coeff:.6f}  O2 Asp = {o2_asp_coeff:.6f}")
        f.write(f"{dap},{wt_lys_coeff:.8f},{o2_lys_coeff:.8f},{diff:.8f},{wt_asp_coeff:.8f},{o2_asp_coeff:.8f}\n")

print(f"  Saved to: {coeff_file}")


# ============================================================================
# ANALYSIS 5: COMPREHENSIVE AMINO ACID FLUX IN BIOMASS
# ============================================================================
print("\n" + "=" * 80)
print("ANALYSIS 5: AMINO ACID BIOMASS COEFFICIENTS (WT vs O2)")
print("=" * 80)

amino_acids = {
    "C00041": "Alanine", "C00049": "Aspartate", "C00152": "Asparagine",
    "C00025": "Glutamate", "C00064": "Glutamine", "C00037": "Glycine",
    "C00135": "Histidine", "C00407": "Isoleucine", "C00123": "Leucine",
    "C00047": "Lysine", "C00073": "Methionine", "C00079": "Phenylalanine",
    "C00148": "Proline", "C00065": "Serine", "C00188": "Threonine",
    "C00078": "Tryptophan", "C00082": "Tyrosine", "C00183": "Valine",
    "C00062": "Arginine", "C00097": "Cysteine",
}

aa_file = os.path.join(OUTPUT_DIR, "05_amino_acid_biomass_coefficients.csv")
with open(aa_file, 'w') as f:
    header = "KEGG_ID,Amino_Acid"
    for dap in DAP_VALUES:
        header += f",WT_DAP{dap},O2_DAP{dap},Diff_DAP{dap}"
    f.write(header + "\n")
    
    for aa_id, aa_name in sorted(amino_acids.items(), key=lambda x: x[1]):
        line = f"{aa_id},{aa_name}"
        for dap in DAP_VALUES:
            wt_sij = parse_sij(os.path.join(WT_DIRS[dap], "sij.txt"))
            o2_sij = parse_sij(os.path.join(O2_DIRS[dap], "sij.txt"))
            wt_c = wt_sij.get((f"{aa_id}[K,B]", "Seed_Biomass[K]"), 0)
            o2_c = o2_sij.get((f"{aa_id}[K,B]", "Seed_Biomass[K]"), 0)
            diff = o2_c - wt_c
            line += f",{wt_c:.8f},{o2_c:.8f},{diff:.8f}"
        f.write(line + "\n")
        
        # Print lysine specifically
        if aa_id == "C00047":
            print(f"  ** {aa_name} (C00047) - THE KEY METABOLITE **")
            for dap in DAP_VALUES:
                wt_sij = parse_sij(os.path.join(WT_DIRS[dap], "sij.txt"))
                o2_sij = parse_sij(os.path.join(O2_DIRS[dap], "sij.txt"))
                wt_c = wt_sij.get((f"{aa_id}[K,B]", "Seed_Biomass[K]"), 0)
                o2_c = o2_sij.get((f"{aa_id}[K,B]", "Seed_Biomass[K]"), 0)
                print(f"    DAP {dap:2d}: WT = {wt_c:.6f}  O2 = {o2_c:.6f}  Diff = {o2_c - wt_c:+.6f}")

print(f"  Saved to: {aa_file}")


# ============================================================================
# ANALYSIS 6: DIFFERENTIAL FLUX ANALYSIS (all reactions)
# ============================================================================
print("\n" + "=" * 80)
print("ANALYSIS 6: DIFFERENTIAL FLUX ANALYSIS (Top Changed Reactions)")
print("=" * 80)

diff_file = os.path.join(OUTPUT_DIR, "06_differential_flux_all.csv")

# Collect all reactions
all_rxns = set()
for dap in DAP_VALUES:
    all_rxns.update(wt_results.get(dap, {}).keys())
    all_rxns.update(o2_results.get(dap, {}).keys())

with open(diff_file, 'w') as f:
    header = "Reaction"
    for dap in DAP_VALUES:
        header += f",WT_DAP{dap},O2_DAP{dap},Diff_DAP{dap}"
    header += ",MaxAbsDiff,AvgAbsDiff,ConsistentDirection"
    f.write(header + "\n")
    
    rxn_stats = []
    for rxn in sorted(all_rxns):
        diffs = []
        vals = []
        for dap in DAP_VALUES:
            wt_v = wt_results.get(dap, {}).get(rxn, 0)
            o2_v = o2_results.get(dap, {}).get(rxn, 0)
            d = o2_v - wt_v
            diffs.append(d)
            vals.append((wt_v, o2_v, d))
        
        max_abs = max(abs(d) for d in diffs)
        avg_abs = sum(abs(d) for d in diffs) / len(diffs)
        
        # Check if direction is consistent (all up or all down, ignoring zeros)
        non_zero_diffs = [d for d in diffs if abs(d) > 1e-12]
        if non_zero_diffs:
            all_pos = all(d > 0 for d in non_zero_diffs)
            all_neg = all(d < 0 for d in non_zero_diffs)
            consistent = "UP" if all_pos else ("DOWN" if all_neg else "MIXED")
        else:
            consistent = "NONE"
        
        rxn_stats.append((rxn, vals, max_abs, avg_abs, consistent))
    
    # Sort by max absolute difference
    rxn_stats.sort(key=lambda x: -x[2])
    
    for rxn, vals, max_abs, avg_abs, consistent in rxn_stats:
        line = rxn
        for wt_v, o2_v, d in vals:
            line += f",{wt_v:.8f},{o2_v:.8f},{d:.8f}"
        line += f",{max_abs:.8f},{avg_abs:.8f},{consistent}"
        f.write(line + "\n")

print(f"  Saved to: {diff_file}")

# Print top 30 most changed reactions
print("\n  --- Top 30 Most Differentially Active Reactions ---")
print(f"  {'Reaction':<35s} {'MaxAbsDiff':>12s} {'AvgAbsDiff':>12s} {'Direction':>10s}")
print(f"  {'-'*35} {'-'*12} {'-'*12} {'-'*10}")
for rxn, vals, max_abs, avg_abs, consistent in rxn_stats[:30]:
    print(f"  {rxn:<35s} {max_abs:12.6f} {avg_abs:12.6f} {consistent:>10s}")


# ============================================================================
# ANALYSIS 7: PATHWAY-LEVEL FLUX COMPARISON
# ============================================================================
print("\n" + "=" * 80)
print("ANALYSIS 7: PATHWAY-LEVEL FLUX COMPARISON")
print("=" * 80)

pathway_file = os.path.join(OUTPUT_DIR, "07_pathway_flux_summary.csv")
with open(pathway_file, 'w') as f:
    f.write("Pathway,DAP,WT_TotalAbsFlux,O2_TotalAbsFlux,Difference,Percent_Change\n")
    
    for pathway_name, rxn_ids in PATHWAY_CATEGORIES.items():
        print(f"\n  {pathway_name}:")
        for dap in DAP_VALUES:
            wt_total = 0
            o2_total = 0
            for rxn_base in rxn_ids:
                # Match all compartment variants
                for rxn_full in all_rxns:
                    if rxn_base in rxn_full:
                        wt_total += abs(wt_results.get(dap, {}).get(rxn_full, 0))
                        o2_total += abs(o2_results.get(dap, {}).get(rxn_full, 0))
            
            diff = o2_total - wt_total
            pct = (diff / wt_total * 100) if wt_total > 1e-12 else 0
            arrow = "↑" if diff > 1e-10 else ("↓" if diff < -1e-10 else "=")
            print(f"    DAP {dap:2d}: WT={wt_total:10.6f}  O2={o2_total:10.6f}  Diff={diff:+10.6f} ({pct:+.1f}%) {arrow}")
            f.write(f"{pathway_name},{dap},{wt_total:.8f},{o2_total:.8f},{diff:.8f},{pct:.4f}\n")

print(f"\n  Saved to: {pathway_file}")


# ============================================================================
# ANALYSIS 8: FLUX BOUND DIFFERENCES (constraint analysis)
# ============================================================================
print("\n" + "=" * 80)
print("ANALYSIS 8: CONSTRAINT DIFFERENCES (v_min/v_max bounds)")
print("=" * 80)

bounds_file = os.path.join(OUTPUT_DIR, "08_bound_differences.csv")

# Focus on lysine-related reactions
with open(bounds_file, 'w') as f:
    f.write("Reaction,DAP,WT_vmin,O2_vmin,WT_vmax,O2_vmax,vmin_diff,vmax_diff\n")
    
    print("  Key lysine pathway reaction bounds:")
    for rxn in LYSINE_PATHWAY_REACTIONS.keys():
        printed_header = False
        for dap in DAP_VALUES:
            wt_vmin = parse_bounds(os.path.join(WT_DIRS[dap], "v_min.txt"))
            wt_vmax = parse_bounds(os.path.join(WT_DIRS[dap], "v_max.txt"))
            o2_vmin = parse_bounds(os.path.join(O2_DIRS[dap], "v_min.txt"))
            o2_vmax = parse_bounds(os.path.join(O2_DIRS[dap], "v_max.txt"))
            
            wt_lo = wt_vmin.get(rxn, 0)
            wt_hi = wt_vmax.get(rxn, 0)
            o2_lo = o2_vmin.get(rxn, 0)
            o2_hi = o2_vmax.get(rxn, 0)
            
            if abs(wt_lo - o2_lo) > 1e-10 or abs(wt_hi - o2_hi) > 1e-10:
                if not printed_header:
                    print(f"\n  {rxn} ({LYSINE_PATHWAY_REACTIONS[rxn]}):")
                    printed_header = True
                print(f"    DAP {dap:2d}: WT bounds=[{wt_lo:.6f}, {wt_hi:.6f}]  O2 bounds=[{o2_lo:.6f}, {o2_hi:.6f}]")
            
            f.write(f"{rxn},{dap},{wt_lo:.8f},{o2_lo:.8f},{wt_hi:.8f},{o2_hi:.8f},"
                    f"{o2_lo-wt_lo:.8f},{o2_hi-wt_hi:.8f}\n")

print(f"\n  Saved to: {bounds_file}")


# ============================================================================
# ANALYSIS 9: ASPARTATE FAMILY AMINO ACID COMPETITION
# ============================================================================
print("\n" + "=" * 80)
print("ANALYSIS 9: ASPARTATE-FAMILY AMINO ACID COMPETITION ANALYSIS")
print("=" * 80)
print("  (Lysine, Threonine, Methionine, Isoleucine all branch from aspartate)")

asp_family = {
    "C00047": "Lysine",
    "C00188": "Threonine",
    "C00073": "Methionine",
    "C00407": "Isoleucine",
    "C00049": "Aspartate",
    "C00263": "Homoserine",
}

asp_file = os.path.join(OUTPUT_DIR, "09_aspartate_family_analysis.csv")
with open(asp_file, 'w') as f:
    header = "MetaboliteID,Name"
    for dap in DAP_VALUES:
        header += f",WT_BiomassCoeff_DAP{dap},O2_BiomassCoeff_DAP{dap}"
    f.write(header + "\n")
    
    for met_id, met_name in asp_family.items():
        line = f"{met_id},{met_name}"
        print(f"\n  {met_name} ({met_id}):")
        for dap in DAP_VALUES:
            wt_sij = parse_sij(os.path.join(WT_DIRS[dap], "sij.txt"))
            o2_sij = parse_sij(os.path.join(O2_DIRS[dap], "sij.txt"))
            wt_c = wt_sij.get((f"{met_id}[K,B]", "Seed_Biomass[K]"), 0)
            o2_c = o2_sij.get((f"{met_id}[K,B]", "Seed_Biomass[K]"), 0)
            diff = o2_c - wt_c
            arrow = "↑" if diff > 1e-8 else ("↓" if diff < -1e-8 else "=")
            print(f"    DAP {dap:2d}: WT coeff = {wt_c:.6f}  O2 coeff = {o2_c:.6f}  {arrow}")
            line += f",{wt_c:.8f},{o2_c:.8f}"
        f.write(line + "\n")

    # Also check fluxes through key branch point reactions
    print("\n  --- Branch Point Reaction Fluxes ---")
    branch_rxns = {
        "R00480[K,p]": "Aspartate kinase → Lysine branch",
        "R01773[K,p]": "Aspartate-semialdehyde → Lysine branch (to DHDP)",
        "R00355[K,p]": "Aspartate aminotransferase",
        "R01954[K,p]": "Towards Thr/Met/Ile (homoserine branch)",
    }
    for rxn, name in branch_rxns.items():
        print(f"\n  {rxn} ({name}):")
        for dap in DAP_VALUES:
            wt_v = wt_results.get(dap, {}).get(rxn, 0)
            o2_v = o2_results.get(dap, {}).get(rxn, 0)
            diff = o2_v - wt_v
            arrow = "↑" if diff > 1e-10 else ("↓" if diff < -1e-10 else "=")
            print(f"    DAP {dap:2d}: WT={wt_v:10.6f}  O2={o2_v:10.6f}  Diff={diff:+10.6f} {arrow}")

print(f"\n  Saved to: {asp_file}")


# ============================================================================
# ANALYSIS 10: REACTIONS UNIQUE TO ONE GENOTYPE (active in one, zero in other)
# ============================================================================
print("\n" + "=" * 80)
print("ANALYSIS 10: UNIQUELY ACTIVE REACTIONS (nonzero in one genotype only)")
print("=" * 80)

unique_file = os.path.join(OUTPUT_DIR, "10_uniquely_active_reactions.csv")
with open(unique_file, 'w') as f:
    f.write("DAP,Reaction,Active_In,Flux_Value\n")
    
    for dap in DAP_VALUES:
        wt_active = {r for r, v in wt_results.get(dap, {}).items() if abs(v) > 1e-10}
        o2_active = {r for r, v in o2_results.get(dap, {}).items() if abs(v) > 1e-10}
        
        only_wt = wt_active - o2_active
        only_o2 = o2_active - wt_active
        
        print(f"\n  DAP {dap:2d}: {len(wt_active)} active in WT, {len(o2_active)} active in O2")
        print(f"          {len(only_wt)} unique to WT, {len(only_o2)} unique to O2")
        
        # Check if any lysine-related reactions are unique
        lys_rxns_set = set(LYSINE_PATHWAY_REACTIONS.keys())
        lys_only_wt = only_wt & lys_rxns_set
        lys_only_o2 = only_o2 & lys_rxns_set
        if lys_only_wt:
            print(f"          Lysine rxns only in WT: {lys_only_wt}")
        if lys_only_o2:
            print(f"          Lysine rxns only in O2: {lys_only_o2}")
        
        for rxn in sorted(only_wt):
            f.write(f"{dap},{rxn},WT_only,{wt_results[dap][rxn]:.8f}\n")
        for rxn in sorted(only_o2):
            f.write(f"{dap},{rxn},O2_only,{o2_results[dap][rxn]:.8f}\n")

print(f"\n  Saved to: {unique_file}")


# ============================================================================
# ANALYSIS 11: TEMPORAL FLUX PROFILES FOR KEY REACTIONS
# ============================================================================
print("\n" + "=" * 80)
print("ANALYSIS 11: TEMPORAL FLUX PROFILES FOR KEY REACTIONS")
print("=" * 80)

temporal_file = os.path.join(OUTPUT_DIR, "11_temporal_flux_profiles.csv")

# Focus on reactions important for lysine metabolism
key_temporal_rxns = list(LYSINE_PATHWAY_REACTIONS.keys())
# Add some central carbon metabolism reactions
extra_rxns = ["R00200[K,c]", "R00200[K,p]", "R00206[K,c]", "R00206[K,p]",
              "R00224[K,c]", "R00258[K,c]"]  
# Add phloem imports
for rxn in sorted(all_rxns):
    if "PhloemImport" in rxn or "PhloemTransport" in rxn:
        extra_rxns.append(rxn)
key_temporal_rxns.extend(extra_rxns)

with open(temporal_file, 'w') as f:
    header = "Reaction"
    for dap in DAP_VALUES:
        header += f",WT_DAP{dap},O2_DAP{dap}"
    f.write(header + "\n")
    
    for rxn in key_temporal_rxns:
        line = rxn
        for dap in DAP_VALUES:
            wt_v = wt_results.get(dap, {}).get(rxn, 0)
            o2_v = o2_results.get(dap, {}).get(rxn, 0)
            line += f",{wt_v:.8f},{o2_v:.8f}"
        f.write(line + "\n")

print(f"  Saved to: {temporal_file}")


# ============================================================================
# ANALYSIS 12: SUMMARY STATISTICS
# ============================================================================
print("\n" + "=" * 80)
print("ANALYSIS 12: SUMMARY STATISTICS")
print("=" * 80)

summary_file = os.path.join(OUTPUT_DIR, "12_summary_statistics.txt")
with open(summary_file, 'w') as f:
    def write_both(text):
        print(text)
        f.write(text + "\n")
    
    write_both("\n  === OVERALL SUMMARY ===")
    write_both(f"  Total DAP time points: {len(DAP_VALUES)}")
    write_both(f"  Total reactions in model: ~{len(all_rxns)}")
    
    # Count significantly different reactions per DAP
    write_both("\n  Significantly different reactions (|diff| > 1e-6) per DAP:")
    for dap in DAP_VALUES:
        n_diff = 0
        for rxn in all_rxns:
            wt_v = wt_results.get(dap, {}).get(rxn, 0)
            o2_v = o2_results.get(dap, {}).get(rxn, 0)
            if abs(o2_v - wt_v) > 1e-6:
                n_diff += 1
        write_both(f"    DAP {dap:2d}: {n_diff} reactions")
    
    # Lysine flux through biomass
    write_both("\n  === LYSINE DEMAND THROUGH BIOMASS ===")
    write_both("  (Biomass_flux × |Lys_coefficient| = Lysine consumed for biomass)")
    for dap in DAP_VALUES:
        wt_sij = parse_sij(os.path.join(WT_DIRS[dap], "sij.txt"))
        o2_sij = parse_sij(os.path.join(O2_DIRS[dap], "sij.txt"))
        wt_lys_c = abs(wt_sij.get(("C00047[K,B]", "Seed_Biomass[K]"), 0))
        o2_lys_c = abs(o2_sij.get(("C00047[K,B]", "Seed_Biomass[K]"), 0))
        wt_bm = wt_biomass.get(dap, 0)
        o2_bm = o2_biomass.get(dap, 0)
        wt_lys_demand = wt_lys_c * wt_bm
        o2_lys_demand = o2_lys_c * o2_bm
        fc = o2_lys_demand / wt_lys_demand if wt_lys_demand > 0 else float('inf')
        write_both(f"    DAP {dap:2d}: WT Lys demand = {wt_lys_demand:.8f}  O2 = {o2_lys_demand:.8f}  "
                   f"FC = {fc:.4f}")
    
    # LKR degradation flux
    write_both("\n  === LKR (LYSINE DEGRADATION) FLUX ===")
    write_both("  R00715 is Lysine-ketoglutarate reductase (main degradation)")
    for dap in DAP_VALUES:
        wt_lkr_c = wt_results.get(dap, {}).get("R00715[K,c]", 0)
        o2_lkr_c = o2_results.get(dap, {}).get("R00715[K,c]", 0)
        wt_lkr_p = wt_results.get(dap, {}).get("R00715[K,p]", 0)
        o2_lkr_p = o2_results.get(dap, {}).get("R00715[K,p]", 0)
        write_both(f"    DAP {dap:2d}: WT(c)={wt_lkr_c:10.6f} O2(c)={o2_lkr_c:10.6f}  |  "
                   f"WT(p)={wt_lkr_p:10.6f} O2(p)={o2_lkr_p:10.6f}")
    
    # Biosynthesis flux
    write_both("\n  === DAP DECARBOXYLASE (R00451 - Final Lysine Biosynthesis Step) ===")
    for dap in DAP_VALUES:
        wt_v = wt_results.get(dap, {}).get("R00451[K,p]", 0)
        o2_v = o2_results.get(dap, {}).get("R00451[K,p]", 0)
        diff = o2_v - wt_v
        fc = o2_v / wt_v if wt_v > 1e-12 else float('inf')
        write_both(f"    DAP {dap:2d}: WT={wt_v:10.6f}  O2={o2_v:10.6f}  FC={fc:.4f}")

print(f"  Saved to: {summary_file}")


# ============================================================================
# GENERATE GNUPLOT-COMPATIBLE DATA FOR VISUALIZATION
# ============================================================================
print("\n" + "=" * 80)
print("GENERATING PLOT DATA FILES")
print("=" * 80)

# Plot 1: Biomass over time
plot1_file = os.path.join(OUTPUT_DIR, "plot_biomass.dat")
with open(plot1_file, 'w') as f:
    f.write("# DAP  WT_Biomass  O2_Biomass\n")
    for dap in DAP_VALUES:
        f.write(f"{dap}\t{wt_biomass.get(dap, 0):.8f}\t{o2_biomass.get(dap, 0):.8f}\n")

# Plot 2: Lysine biosynthesis (R00451)
plot2_file = os.path.join(OUTPUT_DIR, "plot_lysine_biosynthesis.dat")
with open(plot2_file, 'w') as f:
    f.write("# DAP  WT_R00451  O2_R00451\n")
    for dap in DAP_VALUES:
        wt_v = wt_results.get(dap, {}).get("R00451[K,p]", 0)
        o2_v = o2_results.get(dap, {}).get("R00451[K,p]", 0)
        f.write(f"{dap}\t{wt_v:.8f}\t{o2_v:.8f}\n")

# Plot 3: LKR degradation
plot3_file = os.path.join(OUTPUT_DIR, "plot_lysine_degradation.dat")
with open(plot3_file, 'w') as f:
    f.write("# DAP  WT_LKR_cyt  O2_LKR_cyt  WT_LKR_plast  O2_LKR_plast\n")
    for dap in DAP_VALUES:
        f.write(f"{dap}\t"
                f"{wt_results.get(dap, {}).get('R00715[K,c]', 0):.8f}\t"
                f"{o2_results.get(dap, {}).get('R00715[K,c]', 0):.8f}\t"
                f"{wt_results.get(dap, {}).get('R00715[K,p]', 0):.8f}\t"
                f"{o2_results.get(dap, {}).get('R00715[K,p]', 0):.8f}\n")

# Plot 4: DHDPS flux
plot4_file = os.path.join(OUTPUT_DIR, "plot_DHDPS.dat")
with open(plot4_file, 'w') as f:
    f.write("# DAP  WT_DHDPS  O2_DHDPS\n")
    for dap in DAP_VALUES:
        wt_v = wt_results.get(dap, {}).get("R02292[K,p]", 0)
        o2_v = o2_results.get(dap, {}).get("R02292[K,p]", 0)
        f.write(f"{dap}\t{wt_v:.8f}\t{o2_v:.8f}\n")

# Plot 5: Exchange/export flux
plot5_file = os.path.join(OUTPUT_DIR, "plot_lysine_exchange.dat")
with open(plot5_file, 'w') as f:
    f.write("# DAP  WT_Exchange  O2_Exchange  WT_ExB  O2_ExB\n")
    for dap in DAP_VALUES:
        f.write(f"{dap}\t"
                f"{wt_results.get(dap, {}).get('Exchange_C00047[K,L]', 0):.8f}\t"
                f"{o2_results.get(dap, {}).get('Exchange_C00047[K,L]', 0):.8f}\t"
                f"{wt_results.get(dap, {}).get('ExB_C00047[K]', 0):.8f}\t"
                f"{o2_results.get(dap, {}).get('ExB_C00047[K]', 0):.8f}\n")

print("  Plot data files saved in analysis_results/")


print("\n" + "=" * 80)
print("ALL ANALYSES COMPLETE!")
print(f"Results saved in: {OUTPUT_DIR}")
print("=" * 80)
print("\nGenerated files:")
for fname in sorted(os.listdir(OUTPUT_DIR)):
    fpath = os.path.join(OUTPUT_DIR, fname)
    size = os.path.getsize(fpath)
    print(f"  {fname:<50s} ({size:,} bytes)")
print("")
