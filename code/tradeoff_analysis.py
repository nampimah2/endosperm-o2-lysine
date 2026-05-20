#!/usr/bin/env python3
"""
==============================================================================
Trade-off & Resource Reallocation Analysis: O2 DAP 22 Focus
==============================================================================

Answers the key question: WHERE does the o2 mutant at DAP 22 pull its carbon
and nitrogen resources from to achieve massively higher lysine flux, and what
pathways are sacrificed (trade-offs) compared to:
  - All WT DAPs
  - All other O2 DAPs

Analyses:
  1. Flux sum matrix for ALL 16 conditions (8 DAPs × 2 genotypes)
  2. Carbon budget: glucose import → glycolysis → TCA → amino acids
  3. Nitrogen budget: glutamate/glutamine → transamination → lysine
  4. Trade-off identification: what DECREASES when lysine INCREASES at O2 DAP22
  5. Bound difference analysis: which regulatory constraints change
  6. Cross-condition comparison (O2 DAP22 vs all others)
  7. Zein/storage protein vs lysine flux competition
  8. Publication-quality figures

Reference: Kim & Reed (2012), Flux Sum Analysis
==============================================================================
"""

import os
import re
import csv
from collections import defaultdict
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# ============================================================================
# CONFIGURATION
# ============================================================================
BASE_DIR = "/mnt/nrdstor/ssbio/nampimah/ENDOSPERM/2026/Endosperm_Model_MARY"
OUTPUT_DIR = os.path.join(BASE_DIR, "tradeoff_analysis")
FIGURE_DIR = os.path.join(OUTPUT_DIR, "figures")
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(FIGURE_DIR, exist_ok=True)

DAP_VALUES = [6, 8, 10, 12, 15, 18, 22, 30]
WT_DIRS = {dap: os.path.join(BASE_DIR, "Wild_Type", f"B{dap}") for dap in DAP_VALUES}
O2_DIRS = {dap: os.path.join(BASE_DIR, "O2_mutant", f"O{dap}") for dap in DAP_VALUES}

# ============================================================================
# METABOLITE ANNOTATIONS
# ============================================================================
METABOLITE_NAMES = {
    "C00047": "L-Lysine", "C00049": "L-Aspartate", "C00441": "Aspartate-4-semialdehyde",
    "C03340": "Dihydrodipicolinate", "C00680": "meso-DAP", "C00449": "Saccharopine",
    "C04076": "Aminoadipate-semialdehyde", "C00956": "Aminoadipate",
    "C00022": "Pyruvate", "C00024": "Acetyl-CoA", "C00026": "2-Oxoglutarate",
    "C00158": "Citrate", "C00149": "Malate", "C00042": "Succinate",
    "C00036": "Oxaloacetate", "C00311": "Isocitrate", "C00122": "Fumarate",
    "C00074": "PEP", "C00031": "Glucose", "C00089": "Sucrose",
    "C00267": "alpha-Glucose", "C00668": "Glucose-6P", "C00085": "Fructose-6P",
    "C00354": "Fructose-1,6-bisP", "C00118": "GAP", "C00065": "Serine",
    "C00037": "Glycine", "C00041": "Alanine", "C00025": "Glutamate",
    "C00064": "Glutamine", "C00123": "Leucine", "C00183": "Valine",
    "C00407": "Isoleucine", "C00062": "Arginine", "C00078": "Tryptophan",
    "C00079": "Phenylalanine", "C00082": "Tyrosine", "C00073": "Methionine",
    "C00097": "Cysteine", "C00135": "Histidine", "C00148": "Proline",
    "C00188": "Threonine", "C00152": "Asparagine", "C00263": "Homoserine",
    "C00001": "H2O", "C00002": "ATP", "C00003": "NAD+", "C00004": "NADH",
    "C00005": "NADPH", "C00006": "NADP+", "C00008": "ADP",
    "C00009": "Pi", "C00013": "PPi", "C00020": "AMP",
    "C00010": "CoA", "C00011": "CO2", "C00007": "O2", "C00080": "H+",
    "C00103": "Glucose-1P", "C00111": "DHAP", "C00197": "3PGA",
    "C00236": "1,3-bisP-glycerate", "C00631": "2PGA",
    "C00077": "Ornithine", "C00134": "Putrescine",
    "C00315": "Spermidine", "C00750": "Spermine",
    "C00033": "Acetate", "C00048": "Glyoxylate",
    "C00095": "Fructose",
    "C00169": "Carbamoyl-P", "C00327": "Citrulline",
}

# Carbon source / sink metabolites
CARBON_METABOLITES = {
    "Glucose_import": ["C00031", "C00267"],
    "Sucrose_import": ["C00089"],
    "Glycolysis": ["C00668", "C00085", "C00354", "C00118", "C00111",
                    "C00197", "C00236", "C00631", "C00074", "C00022"],
    "TCA_cycle": ["C00024", "C00158", "C00311", "C00026", "C00042",
                  "C00122", "C00149", "C00036"],
    "Pyruvate_node": ["C00022", "C00024", "C00036"],
}

# Nitrogen source / sink metabolites
NITROGEN_METABOLITES = {
    "Nitrogen_donors": ["C00025", "C00064"],  # Glu, Gln
    "Aspartate_family": ["C00049", "C00152", "C00263", "C00188", "C00073", "C00047"],
    "Lysine_pathway": ["C00049", "C00441", "C03340", "C00680", "C00047"],
    "Lysine_degradation": ["C00047", "C00449", "C04076", "C00956"],
    "BCAA": ["C00123", "C00183", "C00407"],
    "Aromatic_AA": ["C00078", "C00079", "C00082"],
    "Other_AA": ["C00041", "C00065", "C00037", "C00062", "C00097",
                 "C00135", "C00148"],
    "Polyamines": ["C00077", "C00134", "C00315", "C00750"],
}

CURRENCY_METABOLITES = {"C00001", "C00080", "C00007", "C00011"}

# Important reactions for trade-off
LYSINE_REACTIONS = {
    "R00480": "Aspartate kinase",
    "R02292": "DHDPS (dihydrodipicolinate synthase)",
    "R00451": "DAP decarboxylase",
    "R01213": "DAP aminotransferase",
    "R00715": "Aspartate-semialdehyde dehydrogenase",
    "R04475": "Dihydrodipicolinate reductase",
}

ZEIN_KEYWORDS = ["Zein", "zein", "ZEIN", "Prot_", "protein", "ProtSyn"]


# ============================================================================
# PARSING FUNCTIONS
# ============================================================================
def parse_results(filepath):
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


def parse_bounds(filepath):
    d = {}
    if not os.path.exists(filepath):
        return d
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line in ('/', ''):
                continue
            m = re.match(r"'([^']+)'\s+([-\d.eE+]+)", line)
            if m:
                d[m.group(1)] = float(m.group(2))
    return d


def parse_metabolites(filepath):
    mets = []
    if not os.path.exists(filepath):
        return mets
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line in ('/', ''):
                continue
            m = re.match(r"'([^']+)'", line)
            if m:
                mets.append(m.group(1))
    return mets


def parse_reactions(filepath):
    rxns = []
    if not os.path.exists(filepath):
        return rxns
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line in ('/', ''):
                continue
            m = re.match(r"'([^']+)'", line)
            if m:
                rxns.append(m.group(1))
    return rxns


def get_name(met_id):
    match = re.match(r"(M?C\d+)", met_id)
    if match:
        base_id = match.group(1)
        name = METABOLITE_NAMES.get(base_id, base_id)
        comp_match = re.search(r"\[K,(\w+)\]", met_id)
        comp = comp_match.group(1) if comp_match else ""
        comp_map = {"c": "cyt", "p": "pla", "m": "mit", "v": "vac", "L": "ext", "B": "bnd"}
        return f"{name} [{comp_map.get(comp, comp)}]"
    return met_id


def get_base_kegg(met_id):
    match = re.match(r"(M?C\d+)", met_id)
    return match.group(1) if match else met_id


def get_compartment(met_id):
    match = re.search(r"\[K,(\w+)\]", met_id)
    return match.group(1) if match else ""


def compute_flux_sums(sij, fluxes):
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


def compute_net_flux(sij, fluxes, met_id):
    """Net flux through a metabolite: Σ S(i,j) * v(j) — should be ~0 at steady state."""
    net = 0.0
    for (met, rxn), coeff in sij.items():
        if met == met_id:
            v = fluxes.get(rxn, 0.0)
            net += coeff * v
    return net


def compute_production_consumption(sij, fluxes, met_id):
    """Separate production and consumption fluxes for a metabolite."""
    production = 0.0
    consumption = 0.0
    details = []
    for (met, rxn), coeff in sij.items():
        if met == met_id:
            v = fluxes.get(rxn, 0.0)
            sv = coeff * v
            if sv > 0:
                production += sv
                details.append((rxn, coeff, v, sv, "producing"))
            elif sv < 0:
                consumption += abs(sv)
                details.append((rxn, coeff, v, sv, "consuming"))
    details.sort(key=lambda x: abs(x[3]), reverse=True)
    return production, consumption, details


# ============================================================================
# MAIN ANALYSIS
# ============================================================================
print("=" * 80)
print("TRADE-OFF & RESOURCE REALLOCATION ANALYSIS")
print("Focus: O2 DAP 22 — Where do resources come from for lysine?")
print("=" * 80)

# --- Step 1: Load ALL data ---
print("\n[1] Loading all FBA results and model data...")
wt_fluxes = {}
o2_fluxes = {}
wt_biomass = {}
o2_biomass = {}
wt_sij = {}
o2_sij = {}
wt_vmin = {}
wt_vmax = {}
o2_vmin = {}
o2_vmax = {}

for dap in DAP_VALUES:
    # WT
    f, obj = parse_results(os.path.join(WT_DIRS[dap], "results_FBA.txt"))
    if f:
        wt_fluxes[dap] = f
        wt_biomass[dap] = obj
        print(f"  WT DAP {dap:2d}: {len(f):5d} rxns, Biomass={obj:.8f}")
    wt_sij[dap] = parse_sij(os.path.join(WT_DIRS[dap], "sij.txt"))
    wt_vmin[dap] = parse_bounds(os.path.join(WT_DIRS[dap], "v_min.txt"))
    wt_vmax[dap] = parse_bounds(os.path.join(WT_DIRS[dap], "v_max.txt"))

    # O2
    f, obj = parse_results(os.path.join(O2_DIRS[dap], "results_FBA.txt"))
    if f:
        o2_fluxes[dap] = f
        o2_biomass[dap] = obj
        print(f"  O2 DAP {dap:2d}: {len(f):5d} rxns, Biomass={obj:.8f}")
    o2_sij[dap] = parse_sij(os.path.join(O2_DIRS[dap], "sij.txt"))
    o2_vmin[dap] = parse_bounds(os.path.join(O2_DIRS[dap], "v_min.txt"))
    o2_vmax[dap] = parse_bounds(os.path.join(O2_DIRS[dap], "v_max.txt"))


# --- Step 2: Compute flux sums for all conditions ---
print("\n[2] Computing flux sums for all 16 conditions...")
wt_fs = {}
o2_fs = {}
for dap in DAP_VALUES:
    if dap in wt_fluxes:
        wt_fs[dap] = compute_flux_sums(wt_sij[dap], wt_fluxes[dap])
    if dap in o2_fluxes:
        o2_fs[dap] = compute_flux_sums(o2_sij[dap], o2_fluxes[dap])

# Get all metabolites present in DAP 22
all_mets_22 = set(wt_fs[22].keys()) | set(o2_fs[22].keys())
print(f"  Total metabolites at DAP 22: {len(all_mets_22)}")


# --- Step 3: Carbon budget analysis ---
print("\n[3] Carbon budget analysis...")

# Track import/export reactions (boundary fluxes — reactions with [K,B] metabolites)
carbon_budget_rows = []
for dap in DAP_VALUES:
    for genotype, fluxes, sij in [("WT", wt_fluxes.get(dap, {}), wt_sij.get(dap, {})),
                                    ("O2", o2_fluxes.get(dap, {}), o2_sij.get(dap, {}))]:
        if not fluxes:
            continue

        row = {"DAP": dap, "Genotype": genotype}

        # Compute flux sums for key carbon metabolites in cytosol & plastid
        for category, kegg_ids in CARBON_METABOLITES.items():
            cat_total = 0.0
            for kid in kegg_ids:
                for met_id in all_mets_22:
                    if get_base_kegg(met_id) == kid and get_compartment(met_id) in ('c', 'p'):
                        fs_dict = wt_fs[dap] if genotype == "WT" else o2_fs[dap]
                        cat_total += fs_dict.get(met_id, 0.0)
            row[category] = cat_total

        carbon_budget_rows.append(row)

carbon_df = pd.DataFrame(carbon_budget_rows)
carbon_df.to_csv(os.path.join(OUTPUT_DIR, "01_carbon_budget.csv"), index=False)
print("  Saved: 01_carbon_budget.csv")


# --- Step 4: Nitrogen budget analysis ---
print("\n[4] Nitrogen budget analysis...")
nitrogen_budget_rows = []
for dap in DAP_VALUES:
    for genotype, fluxes, sij in [("WT", wt_fluxes.get(dap, {}), wt_sij.get(dap, {})),
                                    ("O2", o2_fluxes.get(dap, {}), o2_sij.get(dap, {}))]:
        if not fluxes:
            continue

        row = {"DAP": dap, "Genotype": genotype}

        for category, kegg_ids in NITROGEN_METABOLITES.items():
            cat_total = 0.0
            for kid in kegg_ids:
                for met_id in all_mets_22:
                    if get_base_kegg(met_id) == kid and get_compartment(met_id) in ('c', 'p'):
                        fs_dict = wt_fs[dap] if genotype == "WT" else o2_fs[dap]
                        cat_total += fs_dict.get(met_id, 0.0)
            row[category] = cat_total

        nitrogen_budget_rows.append(row)

nitrogen_df = pd.DataFrame(nitrogen_budget_rows)
nitrogen_df.to_csv(os.path.join(OUTPUT_DIR, "02_nitrogen_budget.csv"), index=False)
print("  Saved: 02_nitrogen_budget.csv")


# --- Step 5: Identify trade-offs at DAP 22 ---
print("\n[5] Identifying trade-offs: What DECREASES in O2 DAP22 when lysine INCREASES...")

tradeoff_rows = []
for met_id in sorted(all_mets_22):
    base_kegg = get_base_kegg(met_id)
    if base_kegg in CURRENCY_METABOLITES:
        continue
    comp = get_compartment(met_id)
    if comp == 'B':  # skip boundary
        continue

    wt_val = wt_fs[22].get(met_id, 0.0)
    o2_val = o2_fs[22].get(met_id, 0.0)
    diff = o2_val - wt_val

    # Compute fold change safely
    if wt_val > 1e-10:
        fold = o2_val / wt_val
    elif o2_val > 1e-10:
        fold = float('inf')
    else:
        fold = 1.0

    # Also get flux sums at all other DAPs for context
    wt_mean = np.mean([wt_fs[d].get(met_id, 0.0) for d in DAP_VALUES if d in wt_fs])
    o2_mean = np.mean([o2_fs[d].get(met_id, 0.0) for d in DAP_VALUES if d in o2_fs])

    tradeoff_rows.append({
        "Metabolite_ID": met_id,
        "Name": get_name(met_id),
        "Compartment": comp,
        "WT_DAP22_FluxSum": wt_val,
        "O2_DAP22_FluxSum": o2_val,
        "Diff_O2_minus_WT": diff,
        "Fold_Change": fold,
        "WT_AllDAP_Mean": wt_mean,
        "O2_AllDAP_Mean": o2_mean,
        "Direction": "INCREASE" if diff > 1e-8 else ("DECREASE" if diff < -1e-8 else "UNCHANGED"),
    })

tradeoff_df = pd.DataFrame(tradeoff_rows)

# Separate increases and decreases
increases = tradeoff_df[tradeoff_df["Direction"] == "INCREASE"].sort_values("Diff_O2_minus_WT", ascending=False)
decreases = tradeoff_df[tradeoff_df["Direction"] == "DECREASE"].sort_values("Diff_O2_minus_WT", ascending=True)

print(f"  Metabolites that INCREASE in O2 DAP22: {len(increases)}")
print(f"  Metabolites that DECREASE in O2 DAP22: {len(decreases)}")
print(f"  Unchanged: {len(tradeoff_df[tradeoff_df['Direction'] == 'UNCHANGED'])}")

tradeoff_df.sort_values("Diff_O2_minus_WT", key=abs, ascending=False).to_csv(
    os.path.join(OUTPUT_DIR, "03_tradeoff_all_metabolites_dap22.csv"), index=False)
print("  Saved: 03_tradeoff_all_metabolites_dap22.csv")

# Top increases (what O2 gains)
print("\n  === TOP 30 INCREASES in O2 DAP22 (resources gained) ===")
for _, r in increases.head(30).iterrows():
    print(f"    {r['Name']:<40s}  WT={r['WT_DAP22_FluxSum']:.6f}  O2={r['O2_DAP22_FluxSum']:.6f}  Diff={r['Diff_O2_minus_WT']:+.6f}  Fold={r['Fold_Change']:.1f}x")

# Top decreases (what O2 sacrifices)
print("\n  === TOP 30 DECREASES in O2 DAP22 (resources sacrificed — TRADE-OFFS) ===")
for _, r in decreases.head(30).iterrows():
    print(f"    {r['Name']:<40s}  WT={r['WT_DAP22_FluxSum']:.6f}  O2={r['O2_DAP22_FluxSum']:.6f}  Diff={r['Diff_O2_minus_WT']:+.6f}")


# --- Step 6: Amino acid-specific trade-offs ---
print("\n[6] Amino acid trade-offs at DAP 22...")
all_aa_kegg = set()
for cat, ids in NITROGEN_METABOLITES.items():
    all_aa_kegg.update(ids)

aa_tradeoff_rows = []
for met_id in sorted(all_mets_22):
    base_kegg = get_base_kegg(met_id)
    if base_kegg not in all_aa_kegg:
        continue
    comp = get_compartment(met_id)
    if comp not in ('c', 'p'):
        continue

    row = {"Metabolite_ID": met_id, "Name": get_name(met_id), "Compartment": comp}
    for dap in DAP_VALUES:
        row[f"WT_DAP{dap}"] = wt_fs[dap].get(met_id, 0.0) if dap in wt_fs else 0.0
        row[f"O2_DAP{dap}"] = o2_fs[dap].get(met_id, 0.0) if dap in o2_fs else 0.0
        row[f"Diff_DAP{dap}"] = row[f"O2_DAP{dap}"] - row[f"WT_DAP{dap}"]

    aa_tradeoff_rows.append(row)

aa_df = pd.DataFrame(aa_tradeoff_rows)
aa_df.to_csv(os.path.join(OUTPUT_DIR, "04_amino_acid_tradeoffs.csv"), index=False)
print("  Saved: 04_amino_acid_tradeoffs.csv")

print("\n  Amino acid flux sum differences at DAP 22 (cytosol + plastid, O2 - WT):")
aa_22 = []
for _, r in aa_df.iterrows():
    aa_22.append((r["Name"], r["Diff_DAP22"], r["WT_DAP22"], r["O2_DAP22"]))
aa_22.sort(key=lambda x: x[1], reverse=True)
for name, diff, wt, o2 in aa_22:
    marker = "▲" if diff > 1e-8 else ("▼" if diff < -1e-8 else "=")
    print(f"    {marker} {name:<35s}  WT={wt:.6f}  O2={o2:.6f}  Diff={diff:+.6f}")


# --- Step 7: Bound differences at DAP 22 ---
print("\n[7] Flux bound differences between WT B22 and O2 O22...")

wt_max_22 = wt_vmax[22]
o2_max_22 = o2_vmax[22]
wt_min_22 = wt_vmin[22]
o2_min_22 = o2_vmin[22]

bound_diff_rows = []
all_rxns = set(wt_max_22.keys()) | set(o2_max_22.keys())
for rxn in sorted(all_rxns):
    wmax = wt_max_22.get(rxn, 0.0)
    omax = o2_max_22.get(rxn, 0.0)
    wmin = wt_min_22.get(rxn, 0.0)
    omin = o2_min_22.get(rxn, 0.0)
    if abs(wmax - omax) > 1e-10 or abs(wmin - omin) > 1e-10:
        # Also get actual flux values
        wt_v = wt_fluxes[22].get(rxn, 0.0)
        o2_v = o2_fluxes[22].get(rxn, 0.0)
        bound_diff_rows.append({
            "Reaction": rxn,
            "WT_v_min": wmin, "WT_v_max": wmax,
            "O2_v_min": omin, "O2_v_max": omax,
            "Delta_v_max": omax - wmax,
            "Delta_v_min": omin - wmin,
            "WT_flux": wt_v,
            "O2_flux": o2_v,
            "Flux_diff": o2_v - wt_v,
        })

bound_df = pd.DataFrame(bound_diff_rows)
bound_df.to_csv(os.path.join(OUTPUT_DIR, "05_bound_differences_dap22.csv"), index=False)
print(f"  Reactions with different bounds: {len(bound_df)}")
print("  Saved: 05_bound_differences_dap22.csv")

# Identify which bound changes actually AFFECT flux (bounds that are binding)
binding_changes = bound_df[abs(bound_df["Flux_diff"]) > 1e-6].sort_values("Flux_diff", key=abs, ascending=False)
print(f"  Of these, {len(binding_changes)} have actual flux differences")
print("\n  Top 30 bound changes WITH flux impact:")
for _, r in binding_changes.head(30).iterrows():
    print(f"    {r['Reaction']:<40s}  WT_v=[{r['WT_v_min']:.4f}, {r['WT_v_max']:.4f}]  O2_v=[{r['O2_v_min']:.4f}, {r['O2_v_max']:.4f}]  Flux: WT={r['WT_flux']:.6f} O2={r['O2_flux']:.6f}  Δ={r['Flux_diff']:+.6f}")


# --- Step 8: Cross-condition comparison: O2 DAP22 uniqueness ---
print("\n[8] Cross-condition comparison: What makes O2 DAP22 unique?")

# Build a matrix of flux sums for KEY metabolites across ALL 16 conditions
key_kegg_ids = {
    "C00047": "Lysine", "C00049": "Aspartate", "C00441": "Asp-semialdehyde",
    "C03340": "DHDP", "C00680": "meso-DAP",
    "C00025": "Glutamate", "C00064": "Glutamine",
    "C00022": "Pyruvate", "C00024": "Acetyl-CoA", "C00036": "OAA",
    "C00026": "2-Oxoglutarate", "C00158": "Citrate",
    "C00149": "Malate", "C00042": "Succinate",
    "C00668": "Glucose-6P", "C00074": "PEP",
    "C00188": "Threonine", "C00073": "Methionine", "C00263": "Homoserine",
    "C00123": "Leucine", "C00183": "Valine", "C00407": "Isoleucine",
    "C00079": "Phenylalanine", "C00082": "Tyrosine", "C00078": "Tryptophan",
    "C00041": "Alanine", "C00065": "Serine", "C00037": "Glycine",
    "C00148": "Proline", "C00062": "Arginine",
}

cross_rows = []
for kegg_id, name in sorted(key_kegg_ids.items(), key=lambda x: x[1]):
    # Find in cytosol and plastid
    for comp_code, comp_name in [("c", "cyt"), ("p", "pla")]:
        target_met = f"{kegg_id}[K,{comp_code}]"
        row = {"Metabolite": f"{name} [{comp_name}]", "KEGG_ID": kegg_id, "Compartment": comp_code}

        for dap in DAP_VALUES:
            row[f"WT_DAP{dap}"] = wt_fs[dap].get(target_met, 0.0) if dap in wt_fs else 0.0
            row[f"O2_DAP{dap}"] = o2_fs[dap].get(target_met, 0.0) if dap in o2_fs else 0.0

        # O2 DAP22-specific metrics
        o2_22 = row["O2_DAP22"]
        wt_22 = row["WT_DAP22"]
        o2_others = [row[f"O2_DAP{d}"] for d in DAP_VALUES if d != 22]
        wt_all = [row[f"WT_DAP{d}"] for d in DAP_VALUES]

        row["O2_DAP22_vs_WT_DAP22"] = o2_22 - wt_22
        row["O2_DAP22_vs_O2_mean_others"] = o2_22 - np.mean(o2_others)
        row["O2_DAP22_vs_WT_mean_all"] = o2_22 - np.mean(wt_all)

        cross_rows.append(row)

cross_df = pd.DataFrame(cross_rows)
cross_df.to_csv(os.path.join(OUTPUT_DIR, "06_cross_condition_key_metabolites.csv"), index=False)
print("  Saved: 06_cross_condition_key_metabolites.csv")


# --- Step 9: Detailed lysine reaction flux comparison ---
print("\n[9] Detailed lysine pathway reaction fluxes across all conditions...")

lys_rxn_rows = []
# Find all lysine-related reactions (those involving lysine pathway metabolites)
lys_pathway_mets = {"C00049", "C00441", "C03340", "C00680", "C00047", "C00449", "C04076", "C00956"}
lys_rxn_ids = set()

# From DAP 22 sij, find all reactions involving lysine pathway metabolites in plastid
for (met, rxn), coeff in wt_sij[22].items():
    base = get_base_kegg(met)
    if base in lys_pathway_mets and get_compartment(met) in ('c', 'p'):
        lys_rxn_ids.add(rxn)

print(f"  Found {len(lys_rxn_ids)} reactions involving lysine pathway metabolites")

for rxn in sorted(lys_rxn_ids):
    row = {"Reaction": rxn}
    # Get reaction name if we have it
    base_rxn = rxn.split("[")[0] if "[" in rxn else rxn
    row["Annotation"] = LYSINE_REACTIONS.get(base_rxn, "")

    for dap in DAP_VALUES:
        row[f"WT_DAP{dap}"] = wt_fluxes[dap].get(rxn, 0.0) if dap in wt_fluxes else 0.0
        row[f"O2_DAP{dap}"] = o2_fluxes[dap].get(rxn, 0.0) if dap in o2_fluxes else 0.0

    row["WT_DAP22"] = wt_fluxes[22].get(rxn, 0.0)
    row["O2_DAP22"] = o2_fluxes[22].get(rxn, 0.0)
    row["Diff_DAP22"] = row["O2_DAP22"] - row["WT_DAP22"]

    lys_rxn_rows.append(row)

lys_rxn_df = pd.DataFrame(lys_rxn_rows)
lys_rxn_df.sort_values("Diff_DAP22", key=abs, ascending=False).to_csv(
    os.path.join(OUTPUT_DIR, "07_lysine_reaction_fluxes.csv"), index=False)
print("  Saved: 07_lysine_reaction_fluxes.csv")


# --- Step 10: Storage protein / zein analysis ---
print("\n[10] Storage protein (zein) flux analysis...")

# Find reactions possibly related to storage proteins / zein
zein_rxns = {}
for rxn_id in wt_fluxes[22]:
    # Look for protein synthesis reactions
    if any(kw.lower() in rxn_id.lower() for kw in ZEIN_KEYWORDS):
        zein_rxns[rxn_id] = {"WT_22": wt_fluxes[22].get(rxn_id, 0.0),
                             "O2_22": o2_fluxes[22].get(rxn_id, 0.0)}

# Also scan for import/export reactions involving amino acids at boundary
amino_import_rxns = {}
for (met, rxn), coeff in wt_sij[22].items():
    if get_compartment(met) == 'B':
        base = get_base_kegg(met)
        if base in all_aa_kegg:
            wt_v = wt_fluxes[22].get(rxn, 0.0)
            o2_v = o2_fluxes[22].get(rxn, 0.0)
            if abs(wt_v) > 1e-10 or abs(o2_v) > 1e-10:
                amino_import_rxns[rxn] = {
                    "Metabolite": get_name(met),
                    "WT_flux": wt_v,
                    "O2_flux": o2_v,
                    "Diff": o2_v - wt_v,
                }

print(f"  Storage protein-related reactions found: {len(zein_rxns)}")
print(f"  Amino acid boundary reactions with flux: {len(amino_import_rxns)}")

# Save boundary amino acid fluxes
aa_boundary_rows = []
for rxn, info in sorted(amino_import_rxns.items(), key=lambda x: abs(x[1]["Diff"]), reverse=True):
    aa_boundary_rows.append({"Reaction": rxn, **info})
aa_boundary_df = pd.DataFrame(aa_boundary_rows)
if not aa_boundary_df.empty:
    aa_boundary_df.to_csv(os.path.join(OUTPUT_DIR, "08_amino_acid_boundary_fluxes.csv"), index=False)
    print("  Saved: 08_amino_acid_boundary_fluxes.csv")
    print("\n  Top boundary amino acid flux changes at DAP 22:")
    for _, r in aa_boundary_df.head(20).iterrows():
        print(f"    {r['Reaction']:<40s}  {r['Metabolite']:<25s}  WT={r['WT_flux']:.6f}  O2={r['O2_flux']:.6f}  Diff={r['Diff']:+.6f}")


# --- Step 11: Comprehensive summary report ---
print("\n[11] Generating comprehensive summary report...")

report_lines = []
report_lines.append("=" * 80)
report_lines.append("TRADE-OFF ANALYSIS: WHERE DOES O2 DAP22 GET RESOURCES FOR LYSINE?")
report_lines.append("=" * 80)

# Biomass comparison
report_lines.append("\n## 1. BIOMASS COMPARISON")
report_lines.append(f"{'DAP':<6} {'WT Biomass':>15} {'O2 Biomass':>15} {'Diff':>15} {'O2/WT Ratio':>15}")
report_lines.append("-" * 68)
for dap in DAP_VALUES:
    wt_b = wt_biomass.get(dap, 0.0)
    o2_b = o2_biomass.get(dap, 0.0)
    ratio = o2_b / wt_b if wt_b > 0 else float('inf')
    report_lines.append(f"{dap:<6} {wt_b:>15.8f} {o2_b:>15.8f} {o2_b - wt_b:>+15.8f} {ratio:>15.4f}")

# Carbon budget comparison
report_lines.append("\n## 2. CARBON BUDGET (Flux Sum Totals)")
report_lines.append("Comparing WT DAP22 vs O2 DAP22:")
for _, r in carbon_df[(carbon_df["DAP"] == 22)].iterrows():
    report_lines.append(f"  {r['Genotype']}:")
    for cat in CARBON_METABOLITES:
        report_lines.append(f"    {cat:<25s}: {r[cat]:.6f}")

# Nitrogen budget comparison
report_lines.append("\n## 3. NITROGEN BUDGET (Flux Sum Totals)")
report_lines.append("Comparing WT DAP22 vs O2 DAP22:")
for _, r in nitrogen_df[(nitrogen_df["DAP"] == 22)].iterrows():
    report_lines.append(f"  {r['Genotype']}:")
    for cat in NITROGEN_METABOLITES:
        report_lines.append(f"    {cat:<25s}: {r[cat]:.6f}")

# KEY TRADE-OFFS
report_lines.append("\n## 4. KEY TRADE-OFFS AT DAP 22")
report_lines.append("What DECREASES in O2 when lysine INCREASES:")
report_lines.append(f"{'Metabolite':<45s} {'WT':>12} {'O2':>12} {'O2-WT':>12}")
report_lines.append("-" * 82)
for _, r in decreases.head(40).iterrows():
    report_lines.append(f"{r['Name']:<45s} {r['WT_DAP22_FluxSum']:>12.6f} {r['O2_DAP22_FluxSum']:>12.6f} {r['Diff_O2_minus_WT']:>+12.6f}")

report_lines.append("\nWhat INCREASES in O2 (resources channeled toward):")
report_lines.append(f"{'Metabolite':<45s} {'WT':>12} {'O2':>12} {'O2-WT':>12}")
report_lines.append("-" * 82)
for _, r in increases.head(40).iterrows():
    report_lines.append(f"{r['Name']:<45s} {r['WT_DAP22_FluxSum']:>12.6f} {r['O2_DAP22_FluxSum']:>12.6f} {r['Diff_O2_minus_WT']:>+12.6f}")

# DAP 22 UNIQUENESS
report_lines.append("\n## 5. O2 DAP22 UNIQUENESS")
report_lines.append("Key metabolites: O2 DAP22 vs O2 mean (other DAPs) vs WT mean (all):")

# Filter to cytosol only for cleaner display
cross_cyt = cross_df[cross_df["Compartment"] == "c"].copy()
cross_cyt_sorted = cross_cyt.reindex(cross_cyt["O2_DAP22_vs_WT_DAP22"].abs().sort_values(ascending=False).index)

report_lines.append(f"{'Metabolite':<30s} {'O2 DAP22':>10} {'WT DAP22':>10} {'O2-WT':>10} {'O2_22 vs O2_mean':>16} {'O2_22 vs WT_mean':>16}")
report_lines.append("-" * 95)
for _, r in cross_cyt_sorted.head(30).iterrows():
    report_lines.append(f"{r['Metabolite']:<30s} {r['O2_DAP22']:>10.6f} {r['WT_DAP22']:>10.6f} {r['O2_DAP22_vs_WT_DAP22']:>+10.6f} {r['O2_DAP22_vs_O2_mean_others']:>+16.6f} {r['O2_DAP22_vs_WT_mean_all']:>+16.6f}")

# BOUND CHANGES
report_lines.append("\n## 6. KEY REGULATORY (BOUND) CHANGES AT DAP 22")
report_lines.append(f"Total reactions with different bounds: {len(bound_df)}")
report_lines.append(f"Reactions with both bound change AND flux change: {len(binding_changes)}")
report_lines.append(f"\n{'Reaction':<40s} {'WT bounds':>25} {'O2 bounds':>25} {'WT flux':>10} {'O2 flux':>10} {'Δflux':>10}")
report_lines.append("-" * 122)
for _, r in binding_changes.head(40).iterrows():
    report_lines.append(f"{r['Reaction']:<40s} [{r['WT_v_min']:.4f},{r['WT_v_max']:.4f}]{' ':>5} [{r['O2_v_min']:.4f},{r['O2_v_max']:.4f}]{' ':>5} {r['WT_flux']:>10.6f} {r['O2_flux']:>10.6f} {r['Flux_diff']:>+10.6f}")

report_text = "\n".join(report_lines)
with open(os.path.join(OUTPUT_DIR, "09_tradeoff_summary_report.txt"), 'w') as f:
    f.write(report_text)
print("  Saved: 09_tradeoff_summary_report.txt")
print("\n" + report_text)


# =============================================================================
# FIGURES
# =============================================================================
print("\n[12] Generating publication-quality figures...")

# Figure style
plt.rcParams.update({
    'font.size': 10,
    'font.family': 'sans-serif',
    'axes.linewidth': 1.2,
    'xtick.major.width': 1.0,
    'ytick.major.width': 1.0,
    'pdf.fonttype': 42,
})
WT_COLOR = '#2166ac'
O2_COLOR = '#b2182b'
COLORS_WT = plt.cm.Blues(np.linspace(0.3, 0.9, len(DAP_VALUES)))
COLORS_O2 = plt.cm.Reds(np.linspace(0.3, 0.9, len(DAP_VALUES)))


# --- Figure 1: Carbon & Nitrogen Budget Comparison at DAP 22 ---
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Carbon budget
ax = axes[0]
categories_c = list(CARBON_METABOLITES.keys())
wt_vals_c = [carbon_df[(carbon_df["DAP"] == 22) & (carbon_df["Genotype"] == "WT")][cat].values[0] for cat in categories_c]
o2_vals_c = [carbon_df[(carbon_df["DAP"] == 22) & (carbon_df["Genotype"] == "O2")][cat].values[0] for cat in categories_c]

x = np.arange(len(categories_c))
w = 0.35
ax.bar(x - w/2, wt_vals_c, w, color=WT_COLOR, label='Wild Type', edgecolor='black', linewidth=0.5)
ax.bar(x + w/2, o2_vals_c, w, color=O2_COLOR, label='O2 Mutant', edgecolor='black', linewidth=0.5)
ax.set_ylabel("Flux Sum (mmol/gDW/h)")
ax.set_title("A) Carbon Budget — DAP 22", fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels([c.replace("_", "\n") for c in categories_c], fontsize=8, rotation=45, ha='right')
ax.legend(fontsize=9)
ax.grid(axis='y', alpha=0.3)

# Nitrogen budget
ax = axes[1]
categories_n = list(NITROGEN_METABOLITES.keys())
wt_vals_n = [nitrogen_df[(nitrogen_df["DAP"] == 22) & (nitrogen_df["Genotype"] == "WT")][cat].values[0] for cat in categories_n]
o2_vals_n = [nitrogen_df[(nitrogen_df["DAP"] == 22) & (nitrogen_df["Genotype"] == "O2")][cat].values[0] for cat in categories_n]

x = np.arange(len(categories_n))
ax.bar(x - w/2, wt_vals_n, w, color=WT_COLOR, label='Wild Type', edgecolor='black', linewidth=0.5)
ax.bar(x + w/2, o2_vals_n, w, color=O2_COLOR, label='O2 Mutant', edgecolor='black', linewidth=0.5)
ax.set_ylabel("Flux Sum (mmol/gDW/h)")
ax.set_title("B) Nitrogen Budget — DAP 22", fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels([c.replace("_", "\n") for c in categories_n], fontsize=8, rotation=45, ha='right')
ax.legend(fontsize=9)
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
for ext in ['png', 'pdf']:
    fig.savefig(os.path.join(FIGURE_DIR, f"Fig1_carbon_nitrogen_budget.{ext}"), dpi=600, bbox_inches='tight')
plt.close()
print("  Fig1: Carbon & Nitrogen budget")


# --- Figure 2: Trade-off Waterfall — Top gains and losses at DAP 22 ---
fig, ax = plt.subplots(figsize=(14, 8))

n_show = 20
top_inc = increases.head(n_show)
top_dec = decreases.head(n_show)

# Combine: decreases on left (negative), increases on right (positive)
names_combined = list(top_dec["Name"].values) + list(top_inc["Name"].values)
diffs_combined = list(top_dec["Diff_O2_minus_WT"].values) + list(top_inc["Diff_O2_minus_WT"].values)

colors = [O2_COLOR if d < 0 else WT_COLOR for d in diffs_combined]  # Red = decrease, Blue = increase
# Actually swap: red for O2 loss (decrease), blue for O2 gain doesn't make sense
# Use: RED = decrease in O2 (trade-off), GREEN = increase in O2 (gain)
colors = ['#d73027' if d < 0 else '#1a9850' for d in diffs_combined]

y_pos = np.arange(len(names_combined))
ax.barh(y_pos, diffs_combined, color=colors, edgecolor='black', linewidth=0.3, height=0.7)
ax.set_yticks(y_pos)
ax.set_yticklabels(names_combined, fontsize=8)
ax.axvline(x=0, color='black', linewidth=1)
ax.set_xlabel("Flux Sum Difference (O2 − WT) at DAP 22")
ax.set_title("Trade-off Map: Metabolite Flux Sum Changes in O2 DAP22", fontweight='bold')
ax.annotate("← Resources LOST in O2\n(Trade-offs)", xy=(0.15, 0.02), xycoords='axes fraction',
            fontsize=9, color='#d73027', fontweight='bold', ha='center')
ax.annotate("Resources GAINED in O2 →\n(Channeled to lysine)", xy=(0.85, 0.02), xycoords='axes fraction',
            fontsize=9, color='#1a9850', fontweight='bold', ha='center')
ax.grid(axis='x', alpha=0.3)
ax.invert_yaxis()

plt.tight_layout()
for ext in ['png', 'pdf']:
    fig.savefig(os.path.join(FIGURE_DIR, f"Fig2_tradeoff_waterfall.{ext}"), dpi=600, bbox_inches='tight')
plt.close()
print("  Fig2: Trade-off waterfall")


# --- Figure 3: Amino acid flux sums - Heatmap across all conditions ---
fig, ax = plt.subplots(figsize=(16, 10))

# Build heatmap matrix: amino acids (cytosol only) × conditions
aa_cyt = aa_df[aa_df["Compartment"] == "c"].copy()
if aa_cyt.empty:
    aa_cyt = aa_df[aa_df["Compartment"] == "p"].copy()

conditions = [f"WT_{d}" for d in DAP_VALUES] + [f"O2_{d}" for d in DAP_VALUES]
heatmap_data = []
row_labels = []
for _, r in aa_cyt.iterrows():
    row_vals = []
    for dap in DAP_VALUES:
        row_vals.append(r.get(f"WT_DAP{dap}", 0.0))
    for dap in DAP_VALUES:
        row_vals.append(r.get(f"O2_DAP{dap}", 0.0))
    heatmap_data.append(row_vals)
    row_labels.append(r["Name"])

if heatmap_data:
    hm_array = np.array(heatmap_data)

    # Log-scale for better visualization (add small epsilon)
    hm_log = np.log10(hm_array + 1e-10)

    im = ax.imshow(hm_log, aspect='auto', cmap='RdYlBu_r', interpolation='nearest')
    ax.set_xticks(range(len(conditions)))
    ax.set_xticklabels(conditions, rotation=90, fontsize=8)
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=8)
    ax.set_title("Amino Acid Flux Sums Across All Conditions (log10 scale)", fontweight='bold')

    # Add vertical separator between WT and O2
    ax.axvline(x=len(DAP_VALUES) - 0.5, color='white', linewidth=2)

    # Highlight DAP 22 columns
    wt_22_idx = DAP_VALUES.index(22)
    o2_22_idx = len(DAP_VALUES) + DAP_VALUES.index(22)
    for idx in [wt_22_idx, o2_22_idx]:
        ax.axvline(x=idx - 0.5, color='gold', linewidth=1.5, linestyle='--')
        ax.axvline(x=idx + 0.5, color='gold', linewidth=1.5, linestyle='--')

    plt.colorbar(im, ax=ax, label='log10(Flux Sum)', shrink=0.8)
    ax.text(wt_22_idx, -1.5, "DAP22", ha='center', fontsize=8, fontweight='bold', color='gold')
    ax.text(o2_22_idx, -1.5, "DAP22", ha='center', fontsize=8, fontweight='bold', color='gold')

plt.tight_layout()
for ext in ['png', 'pdf']:
    fig.savefig(os.path.join(FIGURE_DIR, f"Fig3_amino_acid_heatmap.{ext}"), dpi=600, bbox_inches='tight')
plt.close()
print("  Fig3: Amino acid heatmap")


# --- Figure 4: Temporal profiles of key metabolites (O2 vs WT) ---
key_targets = [
    ("C00047", "c", "L-Lysine [cyt]"),
    ("C00049", "c", "L-Aspartate [cyt]"),
    ("C00025", "c", "L-Glutamate [cyt]"),
    ("C00036", "c", "OAA [cyt]"),
    ("C00022", "c", "Pyruvate [cyt]"),
    ("C00024", "c", "Acetyl-CoA [cyt]"),
    ("C00026", "c", "2-Oxoglutarate [cyt]"),
    ("C00188", "c", "L-Threonine [cyt]"),
    ("C00073", "c", "L-Methionine [cyt]"),
    ("C00047", "p", "L-Lysine [pla]"),
    ("C00049", "p", "L-Aspartate [pla]"),
    ("C00680", "p", "meso-DAP [pla]"),
]

fig, axes = plt.subplots(4, 3, figsize=(16, 14))
axes = axes.flatten()

for idx, (kegg, comp, label) in enumerate(key_targets):
    if idx >= len(axes):
        break
    ax = axes[idx]
    met_id = f"{kegg}[K,{comp}]"

    wt_vals = [wt_fs[d].get(met_id, 0.0) for d in DAP_VALUES]
    o2_vals = [o2_fs[d].get(met_id, 0.0) for d in DAP_VALUES]

    ax.plot(DAP_VALUES, wt_vals, 'o-', color=WT_COLOR, label='WT', linewidth=2, markersize=6)
    ax.plot(DAP_VALUES, o2_vals, 's-', color=O2_COLOR, label='O2', linewidth=2, markersize=6)

    # Highlight DAP 22
    ax.axvline(x=22, color='gray', linewidth=0.8, linestyle='--', alpha=0.5)
    dap22_idx = DAP_VALUES.index(22)
    ax.plot(22, o2_vals[dap22_idx], 'D', color='gold', markersize=10, zorder=5, markeredgecolor='black')

    ax.set_title(label, fontweight='bold', fontsize=10)
    ax.set_xlabel("DAP")
    ax.set_ylabel("Flux Sum")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

plt.suptitle("Temporal Flux Sum Profiles: Key Metabolites (O2 vs WT)", fontweight='bold', fontsize=13, y=1.01)
plt.tight_layout()
for ext in ['png', 'pdf']:
    fig.savefig(os.path.join(FIGURE_DIR, f"Fig4_temporal_profiles.{ext}"), dpi=600, bbox_inches='tight')
plt.close()
print("  Fig4: Temporal profiles")


# --- Figure 5: Carbon & Nitrogen Budget Temporal Profiles ---
fig, axes = plt.subplots(2, 3, figsize=(16, 10))

# Select key budget categories to plot
budget_cats = [
    ("carbon", "Glycolysis", "Glycolysis Carbon Flux"),
    ("carbon", "TCA_cycle", "TCA Cycle Carbon Flux"),
    ("carbon", "Pyruvate_node", "Pyruvate Node Flux"),
    ("nitrogen", "Nitrogen_donors", "N Donors (Glu + Gln)"),
    ("nitrogen", "Lysine_pathway", "Lysine Pathway"),
    ("nitrogen", "Aspartate_family", "Aspartate Family"),
]

for idx, (budget_type, cat, title) in enumerate(budget_cats):
    ax = axes.flatten()[idx]
    df_budget = carbon_df if budget_type == "carbon" else nitrogen_df

    wt_vals = [df_budget[(df_budget["DAP"] == d) & (df_budget["Genotype"] == "WT")][cat].values[0] for d in DAP_VALUES]
    o2_vals = [df_budget[(df_budget["DAP"] == d) & (df_budget["Genotype"] == "O2")][cat].values[0] for d in DAP_VALUES]

    ax.plot(DAP_VALUES, wt_vals, 'o-', color=WT_COLOR, label='WT', linewidth=2, markersize=6)
    ax.plot(DAP_VALUES, o2_vals, 's-', color=O2_COLOR, label='O2', linewidth=2, markersize=6)
    ax.axvline(x=22, color='gray', linewidth=0.8, linestyle='--', alpha=0.5)

    dap22_idx = DAP_VALUES.index(22)
    ax.plot(22, o2_vals[dap22_idx], 'D', color='gold', markersize=10, zorder=5, markeredgecolor='black')

    ax.set_title(title, fontweight='bold', fontsize=10)
    ax.set_xlabel("DAP")
    ax.set_ylabel("Flux Sum")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

plt.suptitle("Temporal Budget Profiles: Carbon & Nitrogen Pathways", fontweight='bold', fontsize=13, y=1.01)
plt.tight_layout()
for ext in ['png', 'pdf']:
    fig.savefig(os.path.join(FIGURE_DIR, f"Fig5_budget_temporal.{ext}"), dpi=600, bbox_inches='tight')
plt.close()
print("  Fig5: Budget temporal profiles")


# --- Figure 6: O2 DAP22 Uniqueness Radar / Bar Chart ---
fig, ax = plt.subplots(figsize=(14, 8))

# For key metabolites, show: O2_DAP22 vs WT_DAP22 vs O2_mean_others
radar_mets = [
    "Lysine [cyt]", "Aspartate [cyt]", "Glutamate [cyt]", "OAA [cyt]",
    "Pyruvate [cyt]", "Acetyl-CoA [cyt]", "2-Oxoglutarate [cyt]",
    "Threonine [cyt]", "Methionine [cyt]", "Homoserine [cyt]",
    "Leucine [cyt]", "Valine [cyt]", "Isoleucine [cyt]",
    "Phenylalanine [cyt]", "Alanine [cyt]", "Proline [cyt]",
]

radar_data = []
radar_labels = []
for mname in radar_mets:
    row = cross_df[cross_df["Metabolite"] == mname]
    if not row.empty:
        r = row.iloc[0]
        radar_data.append({
            "name": mname,
            "O2_DAP22": r["O2_DAP22"],
            "WT_DAP22": r["WT_DAP22"],
            "O2_others_mean": r["O2_DAP22"] - r["O2_DAP22_vs_O2_mean_others"],  # = O2 mean others
        })
        radar_labels.append(mname.replace(" [cyt]", ""))

if radar_data:
    x = np.arange(len(radar_labels))
    w = 0.25
    o2_22_vals = [d["O2_DAP22"] for d in radar_data]
    wt_22_vals = [d["WT_DAP22"] for d in radar_data]
    o2_other_vals = [d["O2_others_mean"] for d in radar_data]

    ax.bar(x - w, wt_22_vals, w, color=WT_COLOR, label='WT DAP22', edgecolor='black', linewidth=0.3)
    ax.bar(x, o2_other_vals, w, color='#fdae61', label='O2 Mean (other DAPs)', edgecolor='black', linewidth=0.3)
    ax.bar(x + w, o2_22_vals, w, color=O2_COLOR, label='O2 DAP22', edgecolor='black', linewidth=0.3)

    ax.set_xticks(x)
    ax.set_xticklabels(radar_labels, rotation=45, ha='right', fontsize=9)
    ax.set_ylabel("Flux Sum (mmol/gDW/h)")
    ax.set_title("O2 DAP22 Uniqueness: Key Amino Acid Flux Sums", fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
for ext in ['png', 'pdf']:
    fig.savefig(os.path.join(FIGURE_DIR, f"Fig6_o2_dap22_uniqueness.{ext}"), dpi=600, bbox_inches='tight')
plt.close()
print("  Fig6: O2 DAP22 uniqueness")


# --- Figure 7: Lysine pathway reaction flux profiles ---
fig, axes = plt.subplots(2, 3, figsize=(16, 10))

# Key lysine reactions to plot
key_lys_rxns = [
    ("R00480[K,p]", "Aspartate kinase [pla]"),
    ("R02292[K,p]", "DHDPS [pla]"),
    ("R00451[K,p]", "DAP decarboxylase [pla]"),
    ("R00480[K,c]", "Aspartate kinase [cyt]"),
    ("R02292[K,c]", "DHDPS [cyt]"),
    ("R00451[K,c]", "DAP decarboxylase [cyt]"),
]

for idx, (rxn_id, title) in enumerate(key_lys_rxns):
    if idx >= len(axes.flatten()):
        break
    ax = axes.flatten()[idx]

    wt_vals = [wt_fluxes[d].get(rxn_id, 0.0) for d in DAP_VALUES]
    o2_vals = [o2_fluxes[d].get(rxn_id, 0.0) for d in DAP_VALUES]

    ax.plot(DAP_VALUES, wt_vals, 'o-', color=WT_COLOR, label='WT', linewidth=2, markersize=6)
    ax.plot(DAP_VALUES, o2_vals, 's-', color=O2_COLOR, label='O2', linewidth=2, markersize=6)
    ax.axvline(x=22, color='gray', linewidth=0.8, linestyle='--', alpha=0.5)

    dap22_idx = DAP_VALUES.index(22)
    ax.plot(22, o2_vals[dap22_idx], 'D', color='gold', markersize=10, zorder=5, markeredgecolor='black')

    ax.set_title(title, fontweight='bold', fontsize=10)
    ax.set_xlabel("DAP")
    ax.set_ylabel("Flux (mmol/gDW/h)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

plt.suptitle("Lysine Pathway Reaction Fluxes Across Development", fontweight='bold', fontsize=13, y=1.01)
plt.tight_layout()
for ext in ['png', 'pdf']:
    fig.savefig(os.path.join(FIGURE_DIR, f"Fig7_lysine_reaction_profiles.{ext}"), dpi=600, bbox_inches='tight')
plt.close()
print("  Fig7: Lysine reaction profiles")


# --- Figure 8: Stacked area - Pathway trade-offs over development ---
fig, axes = plt.subplots(1, 2, figsize=(16, 7))

pathway_categories = {
    "Lysine": ["C00047", "C00441", "C03340", "C00680"],
    "Aspartate_family": ["C00049", "C00152", "C00263", "C00188", "C00073"],
    "BCAA": ["C00123", "C00183", "C00407"],
    "Aromatic_AA": ["C00078", "C00079", "C00082"],
    "Glutamate_family": ["C00025", "C00064", "C00062", "C00148", "C00077"],
    "Other_AA": ["C00041", "C00065", "C00037", "C00097", "C00135"],
}

for ax_idx, (genotype_label, fs_dict) in enumerate([("Wild Type", wt_fs), ("O2 Mutant", o2_fs)]):
    ax = axes[ax_idx]

    pathway_vals = {pname: [] for pname in pathway_categories}
    for dap in DAP_VALUES:
        if dap not in fs_dict:
            for pname in pathway_categories:
                pathway_vals[pname].append(0.0)
            continue
        for pname, kegg_ids in pathway_categories.items():
            total = 0.0
            for kid in kegg_ids:
                for comp in ['c', 'p']:
                    met_id = f"{kid}[K,{comp}]"
                    total += fs_dict[dap].get(met_id, 0.0)
            pathway_vals[pname].append(total)

    # Stacked area
    colors_stack = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#a65628']
    bottom = np.zeros(len(DAP_VALUES))
    for pidx, (pname, vals) in enumerate(pathway_vals.items()):
        vals_arr = np.array(vals)
        ax.fill_between(DAP_VALUES, bottom, bottom + vals_arr,
                        alpha=0.7, color=colors_stack[pidx % len(colors_stack)], label=pname.replace("_", " "))
        bottom += vals_arr

    ax.axvline(x=22, color='black', linewidth=1.5, linestyle='--', alpha=0.7)
    ax.set_xlabel("Days After Pollination")
    ax.set_ylabel("Total Flux Sum")
    ax.set_title(f"{genotype_label}", fontweight='bold', fontsize=12)
    ax.legend(fontsize=8, loc='upper left')
    ax.grid(alpha=0.3)

plt.suptitle("Amino Acid Pathway Flux Distribution Over Development\n(Stacked: How resources are allocated between pathways)",
             fontweight='bold', fontsize=13, y=1.03)
plt.tight_layout()
for ext in ['png', 'pdf']:
    fig.savefig(os.path.join(FIGURE_DIR, f"Fig8_pathway_stacked_area.{ext}"), dpi=600, bbox_inches='tight')
plt.close()
print("  Fig8: Pathway stacked area")


# --- Figure 9: Difference heatmap — O2 vs WT for key metabolites across DAPs ---
fig, ax = plt.subplots(figsize=(12, 10))

diff_data = []
diff_labels = []
for kegg_id, name in sorted(key_kegg_ids.items(), key=lambda x: x[1]):
    met_id = f"{kegg_id}[K,c]"
    row_vals = []
    has_data = False
    for dap in DAP_VALUES:
        wt_v = wt_fs[dap].get(met_id, 0.0) if dap in wt_fs else 0.0
        o2_v = o2_fs[dap].get(met_id, 0.0) if dap in o2_fs else 0.0
        diff_val = o2_v - wt_v
        row_vals.append(diff_val)
        if abs(diff_val) > 1e-10:
            has_data = True
    if has_data:
        diff_data.append(row_vals)
        diff_labels.append(name)

if diff_data:
    diff_array = np.array(diff_data)

    # Use diverging colormap
    vmax = np.max(np.abs(diff_array))
    im = ax.imshow(diff_array, aspect='auto', cmap='RdBu_r', vmin=-vmax, vmax=vmax, interpolation='nearest')
    ax.set_xticks(range(len(DAP_VALUES)))
    ax.set_xticklabels([f"DAP {d}" for d in DAP_VALUES], fontsize=10)
    ax.set_yticks(range(len(diff_labels)))
    ax.set_yticklabels(diff_labels, fontsize=9)
    ax.set_title("Flux Sum Difference (O2 − WT) Across Development\n(Red = Higher in O2, Blue = Higher in WT)",
                 fontweight='bold')

    # Highlight DAP 22 column
    dap22_idx = DAP_VALUES.index(22)
    ax.axvline(x=dap22_idx - 0.5, color='gold', linewidth=2, linestyle='--')
    ax.axvline(x=dap22_idx + 0.5, color='gold', linewidth=2, linestyle='--')

    plt.colorbar(im, ax=ax, label='Flux Sum Diff (O2 − WT)', shrink=0.8)

plt.tight_layout()
for ext in ['png', 'pdf']:
    fig.savefig(os.path.join(FIGURE_DIR, f"Fig9_difference_heatmap.{ext}"), dpi=600, bbox_inches='tight')
plt.close()
print("  Fig9: Difference heatmap")


# --- Figure 10: Aspartate family competition diagram ---
fig, axes = plt.subplots(2, 3, figsize=(16, 10))

asp_family_mets = [
    ("C00049", "c", "L-Aspartate [cyt]"),
    ("C00263", "c", "L-Homoserine [cyt]"),
    ("C00188", "c", "L-Threonine [cyt]"),
    ("C00073", "c", "L-Methionine [cyt]"),
    ("C00047", "c", "L-Lysine [cyt]"),
    ("C00047", "p", "L-Lysine [pla]"),
]

for idx, (kegg, comp, label) in enumerate(asp_family_mets):
    ax = axes.flatten()[idx]
    met_id = f"{kegg}[K,{comp}]"

    wt_vals = [wt_fs[d].get(met_id, 0.0) for d in DAP_VALUES]
    o2_vals = [o2_fs[d].get(met_id, 0.0) for d in DAP_VALUES]

    ax.fill_between(DAP_VALUES, wt_vals, alpha=0.3, color=WT_COLOR)
    ax.fill_between(DAP_VALUES, o2_vals, alpha=0.3, color=O2_COLOR)
    ax.plot(DAP_VALUES, wt_vals, 'o-', color=WT_COLOR, label='WT', linewidth=2, markersize=5)
    ax.plot(DAP_VALUES, o2_vals, 's-', color=O2_COLOR, label='O2', linewidth=2, markersize=5)
    ax.axvline(x=22, color='gray', linewidth=0.8, linestyle='--', alpha=0.5)

    ax.set_title(label, fontweight='bold', fontsize=10)
    ax.set_xlabel("DAP")
    ax.set_ylabel("Flux Sum")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

plt.suptitle("Aspartate Family Competition:\nLysine vs Threonine/Methionine/Homoserine",
             fontweight='bold', fontsize=13, y=1.03)
plt.tight_layout()
for ext in ['png', 'pdf']:
    fig.savefig(os.path.join(FIGURE_DIR, f"Fig10_aspartate_family_competition.{ext}"), dpi=600, bbox_inches='tight')
plt.close()
print("  Fig10: Aspartate family competition")


print("\n" + "=" * 80)
print("ANALYSIS COMPLETE")
print(f"Output directory: {OUTPUT_DIR}")
print(f"  CSV files: 9")
print(f"  Figures: 10 (PNG + PDF)")
print("=" * 80)
