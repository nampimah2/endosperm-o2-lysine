#!/usr/bin/env python3
"""
Gene/Reaction Target Identification for Lysine (R00451) Upregulation in Wild-Type Maize Endosperm
=================================================================================================

This script identifies specific reactions in the wild-type model that, when modified,
lead to increased flux through R00451 (DAP decarboxylase, the final step of lysine biosynthesis).

Analyses:
1. Baseline FBA - establish reference fluxes
2. FSEOF (Flux Scanning based on Enforced Objective Flux) - identify reactions co-varying with R00451
3. Single-reaction overexpression screen - increase v_max for each reaction
4. Single-reaction knockout screen - set flux bounds to 0
5. O2-bound mimicry - apply O2 mutant bounds to WT model
6. Combination target analysis - test top hits together
7. Verification - confirm targets produce R00451 upregulation

Author: Computational analysis pipeline
Date: March 2026
"""

import os
import sys
import csv
import re
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

BASE_DIR = '/mnt/nrdstor/ssbio/nampimah/ENDOSPERM/2026/Endosperm_Model_MARY'
OUTPUT_DIR = os.path.join(BASE_DIR, 'gene_target_analysis')
os.makedirs(OUTPUT_DIR, exist_ok=True)

WT_DIRS = {
    6: 'Wild_Type/B6', 8: 'Wild_Type/B8', 10: 'Wild_Type/B10', 12: 'Wild_Type/B12',
    15: 'Wild_Type/B15', 18: 'Wild_Type/B18', 22: 'Wild_Type/B22', 30: 'Wild_Type/B30'
}
O2_DIRS = {
    6: 'O2_mutant/O6', 8: 'O2_mutant/O8', 10: 'O2_mutant/O10', 12: 'O2_mutant/O12',
    15: 'O2_mutant/O15', 18: 'O2_mutant/O18', 22: 'O2_mutant/O22', 30: 'O2_mutant/O30'
}

LYSINE_TARGET = "R00451[K,p]"  # DAP decarboxylase - final lysine biosynthesis step
BIOMASS_RXN = "Seed_Biomass[K]"

# ============================================================
# MODEL PARSING
# ============================================================

def parse_model(model_dir):
    """Parse GAMS model files into numpy arrays for LP solving."""
    full_dir = os.path.join(BASE_DIR, model_dir)
    
    # Parse metabolites
    metabolites = []
    with open(os.path.join(full_dir, 'metabolites.txt')) as f:
        for line in f:
            line = line.strip()
            if line and line != '/':
                met = line.strip("'").strip(",").strip("'")
                if met:
                    metabolites.append(met)
    
    # Parse reactions
    reactions = []
    with open(os.path.join(full_dir, 'reactions.txt')) as f:
        for line in f:
            line = line.strip()
            if line and line != '/':
                rxn = line.strip("'").strip(",").strip("'")
                if rxn:
                    reactions.append(rxn)
    
    met_idx = {m: i for i, m in enumerate(metabolites)}
    rxn_idx = {r: i for i, r in enumerate(reactions)}
    
    n_met = len(metabolites)
    n_rxn = len(reactions)
    
    # Parse stoichiometric matrix (sparse)
    rows, cols, vals = [], [], []
    with open(os.path.join(full_dir, 'sij.txt')) as f:
        for line in f:
            line = line.strip()
            if line and line != '/':
                # Format: 'met'.'rxn'  value
                match = re.match(r"'([^']+)'\s*\.\s*'([^']+)'\s+([-\d.eE+]+)", line)
                if match:
                    met_name = match.group(1)
                    rxn_name = match.group(2)
                    val = float(match.group(3))
                    if met_name in met_idx and rxn_name in rxn_idx:
                        rows.append(met_idx[met_name])
                        cols.append(rxn_idx[rxn_name])
                        vals.append(val)
    
    S = csr_matrix((vals, (rows, cols)), shape=(n_met, n_rxn))
    
    # Parse bounds
    v_max = np.zeros(n_rxn)
    with open(os.path.join(full_dir, 'v_max.txt')) as f:
        for line in f:
            line = line.strip()
            if line and line != '/':
                match = re.match(r"'([^']+)'\s+([-\d.eE+]+)", line)
                if match:
                    rxn_name = match.group(1)
                    val = float(match.group(2))
                    if rxn_name in rxn_idx:
                        v_max[rxn_idx[rxn_name]] = val
    
    v_min = np.zeros(n_rxn)
    with open(os.path.join(full_dir, 'v_min.txt')) as f:
        for line in f:
            line = line.strip()
            if line and line != '/':
                match = re.match(r"'([^']+)'\s+([-\d.eE+]+)", line)
                if match:
                    rxn_name = match.group(1)
                    val = float(match.group(2))
                    if rxn_name in rxn_idx:
                        v_min[rxn_idx[rxn_name]] = val
    
    return {
        'S': S,
        'v_min': v_min,
        'v_max': v_max,
        'metabolites': metabolites,
        'reactions': reactions,
        'met_idx': met_idx,
        'rxn_idx': rxn_idx,
        'n_met': n_met,
        'n_rxn': n_rxn
    }


def solve_fba(model, objective_rxn=None, objective_sense='max',
              extra_lb=None, extra_ub=None):
    """
    Solve FBA using scipy linprog (minimization).
    For maximization, we negate the objective.
    
    extra_lb, extra_ub: dict of {rxn_index: bound_value} for modified bounds
    """
    n_rxn = model['n_rxn']
    
    # Objective: minimize c^T * v
    c = np.zeros(n_rxn)
    if objective_rxn is not None:
        if isinstance(objective_rxn, str):
            obj_idx = model['rxn_idx'].get(objective_rxn)
            if obj_idx is None:
                return None, None, 'infeasible'
        else:
            obj_idx = objective_rxn
        
        if objective_sense == 'max':
            c[obj_idx] = -1  # minimize negative = maximize
        else:
            c[obj_idx] = 1
    
    # Bounds
    lb = model['v_min'].copy()
    ub = model['v_max'].copy()
    
    if extra_lb:
        for idx, val in extra_lb.items():
            lb[idx] = max(lb[idx], val) if val > lb[idx] else val
    if extra_ub:
        for idx, val in extra_ub.items():
            ub[idx] = val
    
    bounds = list(zip(lb, ub))
    
    # Equality constraint: S * v = 0
    A_eq = model['S'].toarray()
    b_eq = np.zeros(model['n_met'])
    
    try:
        result = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method='highs')
        if result.success:
            obj_val = -result.fun if objective_sense == 'max' else result.fun
            return obj_val, result.x, 'optimal'
        else:
            return None, None, result.message
    except Exception as e:
        return None, None, str(e)


# ============================================================
# ANALYSIS 1: BASELINE FBA
# ============================================================

def analysis_baseline(models_wt):
    """Compute baseline FBA for all WT models."""
    print("\n" + "="*70)
    print("ANALYSIS 1: BASELINE FBA FOR ALL WT MODELS")
    print("="*70)
    
    results = {}
    for dap in sorted(models_wt.keys()):
        model = models_wt[dap]
        biomass_idx = model['rxn_idx'].get(BIOMASS_RXN)
        target_idx = model['rxn_idx'].get(LYSINE_TARGET)
        
        obj_val, fluxes, status = solve_fba(model, BIOMASS_RXN, 'max')
        
        if status == 'optimal':
            r00451_flux = fluxes[target_idx] if target_idx is not None else 0
            biomass_flux = fluxes[biomass_idx] if biomass_idx is not None else 0
            results[dap] = {
                'biomass': biomass_flux,
                'R00451': r00451_flux,
                'fluxes': fluxes,
                'status': status
            }
            print(f"  DAP {dap:2d}: Biomass = {biomass_flux:.8e}, R00451 = {r00451_flux:.8e}")
        else:
            results[dap] = {'biomass': 0, 'R00451': 0, 'fluxes': None, 'status': status}
            print(f"  DAP {dap:2d}: INFEASIBLE ({status})")
    
    # Save baseline
    with open(os.path.join(OUTPUT_DIR, '01_baseline_fba.csv'), 'w') as f:
        f.write("DAP,Biomass_Flux,R00451_Flux,Status\n")
        for dap in sorted(results.keys()):
            r = results[dap]
            f.write(f"{dap},{r['biomass']:.10e},{r['R00451']:.10e},{r['status']}\n")
    
    return results


# ============================================================
# ANALYSIS 2: FSEOF (Flux Scanning based on Enforced Objective Flux)
# ============================================================

def analysis_fseof(models_wt, baseline_results, n_steps=10):
    """
    FSEOF: Fix biomass at suboptimal level, then force increasing R00451 flux.
    Identify reactions whose flux must change to accommodate increased R00451.
    """
    print("\n" + "="*70)
    print("ANALYSIS 2: FSEOF - FLUX SCANNING FOR R00451 TARGETS")
    print("="*70)
    
    all_fseof = {}
    
    for dap in sorted(models_wt.keys()):
        model = models_wt[dap]
        base = baseline_results[dap]
        if base['status'] != 'optimal':
            continue
        
        biomass_idx = model['rxn_idx'][BIOMASS_RXN]
        target_idx = model['rxn_idx'].get(LYSINE_TARGET)
        if target_idx is None:
            continue
        
        base_biomass = base['biomass']
        base_r00451 = base['R00451']
        
        print(f"\n  DAP {dap}: Baseline R00451 = {base_r00451:.6e}, Biomass = {base_biomass:.6e}")
        
        # Step 1: Find max R00451 when biomass is fixed at 90% of optimal
        biomass_floor = 0.90 * base_biomass
        
        # Fix biomass at floor
        extra_lb_bio = {biomass_idx: biomass_floor}
        
        max_r00451_val, max_fluxes, status = solve_fba(
            model, LYSINE_TARGET, 'max',
            extra_lb=extra_lb_bio
        )
        
        if status != 'optimal' or max_r00451_val is None:
            print(f"    Could not find max R00451 at 90% biomass: {status}")
            # Try with 80% biomass
            biomass_floor = 0.80 * base_biomass
            extra_lb_bio = {biomass_idx: biomass_floor}
            max_r00451_val, max_fluxes, status = solve_fba(
                model, LYSINE_TARGET, 'max',
                extra_lb=extra_lb_bio
            )
            if status != 'optimal' or max_r00451_val is None:
                print(f"    Still infeasible at 80% biomass: {status}")
                continue
            print(f"    Using 80% biomass floor: {biomass_floor:.6e}")
        
        print(f"    Max R00451 at constrained biomass: {max_r00451_val:.6e}")
        
        if max_r00451_val <= base_r00451 * 1.01:
            print(f"    No room to increase R00451")
            continue
        
        # Step 2: Scan increasing R00451 levels
        r00451_levels = np.linspace(base_r00451, max_r00451_val, n_steps + 1)
        flux_profiles = {}
        
        for step_i, r00451_target in enumerate(r00451_levels):
            extra_lb_target = {biomass_idx: biomass_floor, target_idx: r00451_target}
            
            obj_val, fluxes, status = solve_fba(
                model, BIOMASS_RXN, 'max',
                extra_lb=extra_lb_target
            )
            
            if status == 'optimal' and fluxes is not None:
                flux_profiles[step_i] = fluxes
        
        if len(flux_profiles) < 3:
            print(f"    Not enough feasible steps ({len(flux_profiles)})")
            continue
        
        # Step 3: Identify co-varying reactions
        # A reaction is a target if its flux monotonically increases or decreases with R00451
        n_rxn = model['n_rxn']
        step_keys = sorted(flux_profiles.keys())
        
        covarying = {}
        for j in range(n_rxn):
            fluxes_j = [flux_profiles[s][j] for s in step_keys]
            
            # Skip near-zero reactions
            max_abs = max(abs(f) for f in fluxes_j)
            if max_abs < 1e-10:
                continue
            
            # Compute correlation with R00451 level
            r00451_vals = [r00451_levels[s] if s < len(r00451_levels) else r00451_levels[-1] for s in step_keys]
            
            if np.std(fluxes_j) < 1e-12:
                continue
            
            corr = np.corrcoef(r00451_vals[:len(fluxes_j)], fluxes_j)[0, 1]
            
            # Flux change from baseline to max
            delta = fluxes_j[-1] - fluxes_j[0]
            fold = fluxes_j[-1] / fluxes_j[0] if abs(fluxes_j[0]) > 1e-12 else float('inf')
            
            if abs(corr) > 0.7 and abs(delta) > 1e-8:
                rxn_name = model['reactions'][j]
                direction = "UP" if corr > 0 else "DOWN"
                covarying[rxn_name] = {
                    'correlation': corr,
                    'direction': direction,
                    'baseline_flux': fluxes_j[0],
                    'max_flux': fluxes_j[-1],
                    'delta': delta,
                    'fold_change': fold
                }
        
        all_fseof[dap] = covarying
        print(f"    Found {len(covarying)} co-varying reactions")
        
        # Print top hits
        sorted_targets = sorted(covarying.items(), key=lambda x: abs(x[1]['correlation']), reverse=True)
        for rxn, info in sorted_targets[:10]:
            print(f"      {rxn:40s} corr={info['correlation']:+.3f} dir={info['direction']:4s} delta={info['delta']:+.4e}")
    
    # Save FSEOF results
    with open(os.path.join(OUTPUT_DIR, '02_fseof_targets.csv'), 'w') as f:
        f.write("DAP,Reaction,Correlation,Direction,Baseline_Flux,Max_Flux,Delta,Fold_Change\n")
        for dap in sorted(all_fseof.keys()):
            for rxn, info in sorted(all_fseof[dap].items(), key=lambda x: abs(x[1]['correlation']), reverse=True):
                f.write(f"{dap},{rxn},{info['correlation']:.6f},{info['direction']},"
                       f"{info['baseline_flux']:.8e},{info['max_flux']:.8e},"
                       f"{info['delta']:.8e},{info['fold_change']:.6f}\n")
    
    return all_fseof


# ============================================================
# ANALYSIS 3: SINGLE-REACTION OVEREXPRESSION SCREEN
# ============================================================

def analysis_overexpression_screen(models_wt, baseline_results, factor=2.0):
    """
    For each reaction, increase v_max by a factor and check if R00451 flux increases.
    """
    print("\n" + "="*70)
    print(f"ANALYSIS 3: SINGLE-REACTION OVEREXPRESSION SCREEN ({factor}x v_max)")
    print("="*70)
    
    all_results = {}
    
    for dap in sorted(models_wt.keys()):
        model = models_wt[dap]
        base = baseline_results[dap]
        if base['status'] != 'optimal':
            continue
        
        target_idx = model['rxn_idx'].get(LYSINE_TARGET)
        biomass_idx = model['rxn_idx'][BIOMASS_RXN]
        if target_idx is None:
            continue
        
        base_r00451 = base['R00451']
        base_biomass = base['biomass']
        n_rxn = model['n_rxn']
        
        print(f"\n  DAP {dap}: Screening {n_rxn} reactions...")
        
        hits = {}
        tested = 0
        
        for j in range(n_rxn):
            rxn_name = model['reactions'][j]
            
            # Skip biomass and exchange reactions for overexpression
            if rxn_name == BIOMASS_RXN:
                continue
            
            orig_vmax = model['v_max'][j]
            orig_vmin = model['v_min'][j]
            
            # Skip reactions with zero bounds (not expressed)
            if abs(orig_vmax) < 1e-12 and abs(orig_vmin) < 1e-12:
                continue
            
            # Increase upper bound
            new_vmax = orig_vmax * factor if orig_vmax > 0 else orig_vmax
            # Also increase lower bound magnitude for reversible reactions
            new_vmin = orig_vmin * factor if orig_vmin < 0 else orig_vmin
            
            extra_ub = {j: new_vmax}
            extra_lb = {}
            if new_vmin != orig_vmin:
                extra_lb[j] = new_vmin
            
            obj_val, fluxes, status = solve_fba(
                model, BIOMASS_RXN, 'max',
                extra_lb=extra_lb if extra_lb else None,
                extra_ub=extra_ub
            )
            
            tested += 1
            
            if status == 'optimal' and fluxes is not None:
                new_r00451 = fluxes[target_idx]
                new_biomass = fluxes[biomass_idx]
                delta_r00451 = new_r00451 - base_r00451
                
                if delta_r00451 > 1e-10:  # Meaningful increase
                    fold = new_r00451 / base_r00451 if abs(base_r00451) > 1e-15 else float('inf')
                    hits[rxn_name] = {
                        'new_R00451': new_r00451,
                        'delta_R00451': delta_r00451,
                        'fold_change': fold,
                        'new_biomass': new_biomass,
                        'biomass_change': new_biomass - base_biomass,
                        'orig_vmax': orig_vmax,
                        'new_vmax': new_vmax
                    }
        
        all_results[dap] = hits
        print(f"    Tested {tested} reactions, found {len(hits)} hits")
        
        # Print top hits
        sorted_hits = sorted(hits.items(), key=lambda x: x[1]['delta_R00451'], reverse=True)
        for rxn, info in sorted_hits[:15]:
            print(f"      {rxn:40s} R00451: {info['new_R00451']:.4e} (FC={info['fold_change']:.2f}, "
                  f"dBiomass={info['biomass_change']:+.4e})")
    
    # Save results
    with open(os.path.join(OUTPUT_DIR, '03_overexpression_screen.csv'), 'w') as f:
        f.write("DAP,Reaction,New_R00451,Delta_R00451,Fold_Change,New_Biomass,Biomass_Change,Orig_Vmax,New_Vmax\n")
        for dap in sorted(all_results.keys()):
            for rxn, info in sorted(all_results[dap].items(), key=lambda x: x[1]['delta_R00451'], reverse=True):
                f.write(f"{dap},{rxn},{info['new_R00451']:.10e},{info['delta_R00451']:.10e},"
                       f"{info['fold_change']:.6f},{info['new_biomass']:.10e},"
                       f"{info['biomass_change']:.10e},{info['orig_vmax']:.10e},{info['new_vmax']:.10e}\n")
    
    return all_results


# ============================================================
# ANALYSIS 4: SINGLE-REACTION KNOCKOUT SCREEN
# ============================================================

def analysis_knockout_screen(models_wt, baseline_results):
    """
    For each reaction, set bounds to zero (knockout) and check if R00451 flux changes.
    Some knockouts may redirect flux toward lysine biosynthesis.
    """
    print("\n" + "="*70)
    print("ANALYSIS 4: SINGLE-REACTION KNOCKOUT SCREEN")
    print("="*70)
    
    all_results = {}
    
    for dap in sorted(models_wt.keys()):
        model = models_wt[dap]
        base = baseline_results[dap]
        if base['status'] != 'optimal':
            continue
        
        target_idx = model['rxn_idx'].get(LYSINE_TARGET)
        biomass_idx = model['rxn_idx'][BIOMASS_RXN]
        if target_idx is None:
            continue
        
        base_r00451 = base['R00451']
        base_biomass = base['biomass']
        n_rxn = model['n_rxn']
        
        print(f"\n  DAP {dap}: Knockout screening {n_rxn} reactions...")
        
        hits = {}
        lethal = []
        tested = 0
        
        for j in range(n_rxn):
            rxn_name = model['reactions'][j]
            
            if rxn_name == BIOMASS_RXN or rxn_name == LYSINE_TARGET:
                continue
            
            # Skip already inactive reactions
            if base['fluxes'] is not None and abs(base['fluxes'][j]) < 1e-12:
                continue
            
            extra_lb = {j: 0.0}
            extra_ub = {j: 0.0}
            
            obj_val, fluxes, status = solve_fba(
                model, BIOMASS_RXN, 'max',
                extra_lb=extra_lb,
                extra_ub=extra_ub
            )
            
            tested += 1
            
            if status == 'optimal' and fluxes is not None:
                new_r00451 = fluxes[target_idx]
                new_biomass = fluxes[biomass_idx]
                delta_r00451 = new_r00451 - base_r00451
                
                if abs(new_biomass) < 1e-10:
                    lethal.append(rxn_name)
                    continue
                
                if delta_r00451 > 1e-10:
                    fold = new_r00451 / base_r00451 if abs(base_r00451) > 1e-15 else float('inf')
                    hits[rxn_name] = {
                        'new_R00451': new_r00451,
                        'delta_R00451': delta_r00451,
                        'fold_change': fold,
                        'new_biomass': new_biomass,
                        'biomass_pct': (new_biomass / base_biomass * 100) if base_biomass > 0 else 0,
                        'original_flux': base['fluxes'][j]
                    }
            else:
                lethal.append(rxn_name)
        
        all_results[dap] = {'hits': hits, 'lethal': lethal}
        print(f"    Tested {tested} active reactions, found {len(hits)} positive hits, {len(lethal)} lethal KOs")
        
        sorted_hits = sorted(hits.items(), key=lambda x: x[1]['delta_R00451'], reverse=True)
        for rxn, info in sorted_hits[:15]:
            print(f"      {rxn:40s} R00451: {info['new_R00451']:.4e} (FC={info['fold_change']:.2f}, "
                  f"Biomass={info['biomass_pct']:.1f}%)")
    
    # Save results
    with open(os.path.join(OUTPUT_DIR, '04_knockout_screen.csv'), 'w') as f:
        f.write("DAP,Reaction,New_R00451,Delta_R00451,Fold_Change,New_Biomass,Biomass_Pct_of_WT,Original_Flux\n")
        for dap in sorted(all_results.keys()):
            for rxn, info in sorted(all_results[dap]['hits'].items(), key=lambda x: x[1]['delta_R00451'], reverse=True):
                f.write(f"{dap},{rxn},{info['new_R00451']:.10e},{info['delta_R00451']:.10e},"
                       f"{info['fold_change']:.6f},{info['new_biomass']:.10e},"
                       f"{info['biomass_pct']:.2f},{info['original_flux']:.10e}\n")
    
    with open(os.path.join(OUTPUT_DIR, '04b_lethal_knockouts.csv'), 'w') as f:
        f.write("DAP,Reaction\n")
        for dap in sorted(all_results.keys()):
            for rxn in all_results[dap]['lethal']:
                f.write(f"{dap},{rxn}\n")
    
    return all_results


# ============================================================
# ANALYSIS 5: O2-BOUND MIMICRY
# ============================================================

def analysis_o2_mimicry(models_wt, models_o2, baseline_results):
    """
    For each reaction where O2 bounds differ from WT, apply O2 bounds
    to the WT model one at a time and check R00451 impact.
    """
    print("\n" + "="*70)
    print("ANALYSIS 5: O2-BOUND MIMICRY (APPLY O2 BOUNDS TO WT)")
    print("="*70)
    
    all_results = {}
    
    for dap in sorted(models_wt.keys()):
        wt = models_wt[dap]
        o2 = models_o2[dap]
        base = baseline_results[dap]
        if base['status'] != 'optimal':
            continue
        
        target_idx = wt['rxn_idx'].get(LYSINE_TARGET)
        biomass_idx = wt['rxn_idx'][BIOMASS_RXN]
        if target_idx is None:
            continue
        
        base_r00451 = base['R00451']
        base_biomass = base['biomass']
        
        print(f"\n  DAP {dap}: Testing O2 bound application...")
        
        hits = {}
        tested = 0
        
        for j in range(wt['n_rxn']):
            rxn_name = wt['reactions'][j]
            
            if rxn_name == BIOMASS_RXN:
                continue
            
            # Check if bounds differ
            wt_vmax = wt['v_max'][j]
            wt_vmin = wt['v_min'][j]
            
            # Find corresponding reaction in O2
            o2_j = o2['rxn_idx'].get(rxn_name)
            if o2_j is None:
                continue
            
            o2_vmax = o2['v_max'][o2_j]
            o2_vmin = o2['v_min'][o2_j]
            
            # Skip if bounds are identical
            if abs(wt_vmax - o2_vmax) < 1e-10 and abs(wt_vmin - o2_vmin) < 1e-10:
                continue
            
            # Apply O2 bounds to WT model
            extra_ub = {j: o2_vmax}
            extra_lb = {j: o2_vmin}
            
            obj_val, fluxes, status = solve_fba(
                wt, BIOMASS_RXN, 'max',
                extra_lb=extra_lb,
                extra_ub=extra_ub
            )
            
            tested += 1
            
            if status == 'optimal' and fluxes is not None:
                new_r00451 = fluxes[target_idx]
                new_biomass = fluxes[biomass_idx]
                delta_r00451 = new_r00451 - base_r00451
                
                if abs(delta_r00451) > 1e-10:
                    fold = new_r00451 / base_r00451 if abs(base_r00451) > 1e-15 else float('inf')
                    hits[rxn_name] = {
                        'new_R00451': new_r00451,
                        'delta_R00451': delta_r00451,
                        'fold_change': fold,
                        'new_biomass': new_biomass,
                        'biomass_change': new_biomass - base_biomass,
                        'wt_vmax': wt_vmax, 'o2_vmax': o2_vmax,
                        'wt_vmin': wt_vmin, 'o2_vmin': o2_vmin,
                        'vmax_change_pct': ((o2_vmax - wt_vmax) / wt_vmax * 100) if abs(wt_vmax) > 1e-10 else float('inf')
                    }
        
        all_results[dap] = hits
        print(f"    Tested {tested} reactions with different bounds, found {len(hits)} hits")
        
        # Show top positive hits
        pos_hits = {k: v for k, v in hits.items() if v['delta_R00451'] > 0}
        neg_hits = {k: v for k, v in hits.items() if v['delta_R00451'] < 0}
        
        print(f"    {len(pos_hits)} increase R00451, {len(neg_hits)} decrease R00451")
        sorted_hits = sorted(pos_hits.items(), key=lambda x: x[1]['delta_R00451'], reverse=True)
        for rxn, info in sorted_hits[:10]:
            print(f"      {rxn:40s} R00451: {info['new_R00451']:.4e} (FC={info['fold_change']:.2f}, "
                  f"vmax: {info['wt_vmax']:.3f}->{info['o2_vmax']:.3f})")
    
    # Save results
    with open(os.path.join(OUTPUT_DIR, '05_o2_bound_mimicry.csv'), 'w') as f:
        f.write("DAP,Reaction,New_R00451,Delta_R00451,Fold_Change,New_Biomass,Biomass_Change,"
               "WT_Vmax,O2_Vmax,WT_Vmin,O2_Vmin,Vmax_Change_Pct\n")
        for dap in sorted(all_results.keys()):
            for rxn, info in sorted(all_results[dap].items(), key=lambda x: abs(x[1]['delta_R00451']), reverse=True):
                f.write(f"{dap},{rxn},{info['new_R00451']:.10e},{info['delta_R00451']:.10e},"
                       f"{info['fold_change']:.6f},{info['new_biomass']:.10e},{info['biomass_change']:.10e},"
                       f"{info['wt_vmax']:.10e},{info['o2_vmax']:.10e},"
                       f"{info['wt_vmin']:.10e},{info['o2_vmin']:.10e},{info['vmax_change_pct']:.4f}\n")
    
    return all_results


# ============================================================
# ANALYSIS 6: COMBINATION TARGET ANALYSIS
# ============================================================

def analysis_combinations(models_wt, baseline_results, overexp_results, ko_results, mimicry_results):
    """
    Test top targets in combination (pairs and triples).
    """
    print("\n" + "="*70)
    print("ANALYSIS 6: COMBINATION TARGET ANALYSIS")
    print("="*70)
    
    all_combo_results = {}
    
    for dap in sorted(models_wt.keys()):
        model = models_wt[dap]
        base = baseline_results[dap]
        if base['status'] != 'optimal':
            continue
        
        target_idx = model['rxn_idx'].get(LYSINE_TARGET)
        biomass_idx = model['rxn_idx'][BIOMASS_RXN]
        if target_idx is None:
            continue
        
        base_r00451 = base['R00451']
        base_biomass = base['biomass']
        
        # Collect unique top targets from all screens
        top_targets = set()
        
        # From overexpression
        if dap in overexp_results:
            for rxn in sorted(overexp_results[dap].keys(),
                            key=lambda x: overexp_results[dap][x]['delta_R00451'], reverse=True)[:5]:
                top_targets.add(('OE', rxn))
        
        # From knockout
        if dap in ko_results and 'hits' in ko_results[dap]:
            for rxn in sorted(ko_results[dap]['hits'].keys(),
                            key=lambda x: ko_results[dap]['hits'][x]['delta_R00451'], reverse=True)[:5]:
                top_targets.add(('KO', rxn))
        
        # From O2 mimicry
        if dap in mimicry_results:
            for rxn in sorted(mimicry_results[dap].keys(),
                            key=lambda x: mimicry_results[dap][x]['delta_R00451'], reverse=True)[:5]:
                top_targets.add(('MIM', rxn))
        
        if len(top_targets) < 2:
            continue
        
        print(f"\n  DAP {dap}: Testing {len(top_targets)} unique targets in combinations...")
        
        target_list = list(top_targets)
        combos = {}
        
        # Test all pairs
        for i in range(len(target_list)):
            for j_idx in range(i + 1, min(len(target_list), i + 10)):
                t1_type, t1_rxn = target_list[i]
                t2_type, t2_rxn = target_list[j_idx]
                
                extra_lb = {}
                extra_ub = {}
                
                for t_type, t_rxn in [(t1_type, t1_rxn), (t2_type, t2_rxn)]:
                    rxn_j = model['rxn_idx'].get(t_rxn)
                    if rxn_j is None:
                        continue
                    
                    if t_type == 'OE':
                        extra_ub[rxn_j] = model['v_max'][rxn_j] * 2.0
                        if model['v_min'][rxn_j] < 0:
                            extra_lb[rxn_j] = model['v_min'][rxn_j] * 2.0
                    elif t_type == 'KO':
                        extra_lb[rxn_j] = 0.0
                        extra_ub[rxn_j] = 0.0
                    elif t_type == 'MIM':
                        o2 = models_o2.get(dap)
                        if o2:
                            o2_j = o2['rxn_idx'].get(t_rxn)
                            if o2_j is not None:
                                extra_ub[rxn_j] = o2['v_max'][o2_j]
                                extra_lb[rxn_j] = o2['v_min'][o2_j]
                
                obj_val, fluxes, status = solve_fba(
                    model, BIOMASS_RXN, 'max',
                    extra_lb=extra_lb if extra_lb else None,
                    extra_ub=extra_ub
                )
                
                if status == 'optimal' and fluxes is not None:
                    new_r00451 = fluxes[target_idx]
                    new_biomass = fluxes[biomass_idx]
                    delta_r00451 = new_r00451 - base_r00451
                    
                    combo_name = f"{t1_type}:{t1_rxn} + {t2_type}:{t2_rxn}"
                    if delta_r00451 > 1e-10:
                        fold = new_r00451 / base_r00451 if abs(base_r00451) > 1e-15 else float('inf')
                        combos[combo_name] = {
                            'new_R00451': new_r00451,
                            'delta_R00451': delta_r00451,
                            'fold_change': fold,
                            'new_biomass': new_biomass,
                            'biomass_pct': (new_biomass / base_biomass * 100) if base_biomass > 0 else 0
                        }
        
        all_combo_results[dap] = combos
        print(f"    Found {len(combos)} positive combinations")
        
        sorted_combos = sorted(combos.items(), key=lambda x: x[1]['delta_R00451'], reverse=True)
        for combo, info in sorted_combos[:10]:
            print(f"      {combo[:60]:60s} R00451 FC={info['fold_change']:.2f}, Biomass={info['biomass_pct']:.1f}%")
    
    # Save
    with open(os.path.join(OUTPUT_DIR, '06_combination_targets.csv'), 'w') as f:
        f.write("DAP,Combination,New_R00451,Delta_R00451,Fold_Change,New_Biomass,Biomass_Pct\n")
        for dap in sorted(all_combo_results.keys()):
            for combo, info in sorted(all_combo_results[dap].items(), key=lambda x: x[1]['delta_R00451'], reverse=True):
                f.write(f"{dap},\"{combo}\",{info['new_R00451']:.10e},{info['delta_R00451']:.10e},"
                       f"{info['fold_change']:.6f},{info['new_biomass']:.10e},{info['biomass_pct']:.2f}\n")
    
    return all_combo_results


# ============================================================
# ANALYSIS 7: COMPREHENSIVE TARGET SUMMARY & VERIFICATION
# ============================================================

def analysis_summary_and_verification(models_wt, baseline_results,
                                       overexp_results, ko_results,
                                       mimicry_results, fseof_results):
    """
    Create a master ranked list of targets, cross-referencing all screens.
    Verify top targets by re-running FBA with modifications applied.
    """
    print("\n" + "="*70)
    print("ANALYSIS 7: MASTER TARGET RANKING & VERIFICATION")
    print("="*70)
    
    # Aggregate targets across all screens and DAPs
    target_scores = defaultdict(lambda: {
        'overexp_hits': 0, 'ko_hits': 0, 'mimicry_hits': 0, 'fseof_hits': 0,
        'max_fold_overexp': 0, 'max_fold_ko': 0, 'max_fold_mimicry': 0,
        'best_dap_overexp': None, 'best_dap_ko': None, 'best_dap_mimicry': None,
        'daps_hit': set(), 'screen_types': set()
    })
    
    # Score overexpression hits
    for dap, hits in overexp_results.items():
        for rxn, info in hits.items():
            target_scores[rxn]['overexp_hits'] += 1
            target_scores[rxn]['daps_hit'].add(dap)
            target_scores[rxn]['screen_types'].add('OE')
            if info['fold_change'] > target_scores[rxn]['max_fold_overexp']:
                target_scores[rxn]['max_fold_overexp'] = info['fold_change']
                target_scores[rxn]['best_dap_overexp'] = dap
    
    # Score knockout hits
    for dap, data in ko_results.items():
        if 'hits' not in data:
            continue
        for rxn, info in data['hits'].items():
            target_scores[rxn]['ko_hits'] += 1
            target_scores[rxn]['daps_hit'].add(dap)
            target_scores[rxn]['screen_types'].add('KO')
            if info['fold_change'] > target_scores[rxn]['max_fold_ko']:
                target_scores[rxn]['max_fold_ko'] = info['fold_change']
                target_scores[rxn]['best_dap_ko'] = dap
    
    # Score mimicry hits  
    for dap, hits in mimicry_results.items():
        for rxn, info in hits.items():
            if info['delta_R00451'] > 0:
                target_scores[rxn]['mimicry_hits'] += 1
                target_scores[rxn]['daps_hit'].add(dap)
                target_scores[rxn]['screen_types'].add('MIM')
                if info['fold_change'] > target_scores[rxn]['max_fold_mimicry']:
                    target_scores[rxn]['max_fold_mimicry'] = info['fold_change']
                    target_scores[rxn]['best_dap_mimicry'] = dap
    
    # Score FSEOF hits
    for dap, hits in fseof_results.items():
        for rxn, info in hits.items():
            if info['direction'] == 'UP':
                target_scores[rxn]['fseof_hits'] += 1
                target_scores[rxn]['daps_hit'].add(dap)
                target_scores[rxn]['screen_types'].add('FSEOF')
    
    # Compute composite score
    for rxn, scores in target_scores.items():
        n_screens = len(scores['screen_types'])
        n_daps = len(scores['daps_hit'])
        max_fold = max(scores['max_fold_overexp'], scores['max_fold_ko'], scores['max_fold_mimicry'])
        
        # Composite score: number of screens * number of DAPs * log(max fold change + 1)
        scores['composite'] = n_screens * n_daps * np.log10(max_fold + 1) if max_fold > 0 else 0
        scores['n_screens'] = n_screens
        scores['n_daps'] = n_daps
        scores['max_fold_any'] = max_fold
    
    # Rank by composite score
    ranked = sorted(target_scores.items(), key=lambda x: x[1]['composite'], reverse=True)
    
    print(f"\n  Total unique targets identified: {len(ranked)}")
    print(f"\n  Top 30 targets by composite score:")
    print(f"  {'Rank':>4s} {'Reaction':40s} {'Score':>6s} {'Screens':>7s} {'DAPs':>4s} {'MaxFC':>8s} {'Types'}")
    print(f"  {'-'*4} {'-'*40} {'-'*6} {'-'*7} {'-'*4} {'-'*8} {'-'*20}")
    
    for rank, (rxn, scores) in enumerate(ranked[:30], 1):
        types = ','.join(sorted(scores['screen_types']))
        print(f"  {rank:4d} {rxn:40s} {scores['composite']:6.2f} {scores['n_screens']:7d} "
              f"{scores['n_daps']:4d} {scores['max_fold_any']:8.2f} {types}")
    
    # VERIFICATION: Re-run top targets and confirm R00451 upregulation
    print(f"\n  VERIFICATION: Confirming top 20 targets across all DAPs...")
    
    verification = {}
    top_20 = [rxn for rxn, _ in ranked[:20]]
    
    for rxn in top_20:
        verification[rxn] = {}
        
        for dap in sorted(models_wt.keys()):
            model = models_wt[dap]
            base = baseline_results[dap]
            if base['status'] != 'optimal':
                continue
            
            target_idx = model['rxn_idx'].get(LYSINE_TARGET)
            biomass_idx = model['rxn_idx'][BIOMASS_RXN]
            rxn_j = model['rxn_idx'].get(rxn)
            
            if target_idx is None or rxn_j is None:
                continue
            
            base_r00451 = base['R00451']
            base_biomass = base['biomass']
            
            # Test overexpression (2x)
            extra_ub = {rxn_j: model['v_max'][rxn_j] * 2.0}
            extra_lb = {}
            if model['v_min'][rxn_j] < 0:
                extra_lb[rxn_j] = model['v_min'][rxn_j] * 2.0
            
            obj_val, fluxes, status = solve_fba(
                model, BIOMASS_RXN, 'max',
                extra_lb=extra_lb if extra_lb else None,
                extra_ub=extra_ub
            )
            
            if status == 'optimal' and fluxes is not None:
                new_r00451 = fluxes[target_idx]
                new_biomass = fluxes[biomass_idx]
                verification[rxn][dap] = {
                    'base_R00451': base_r00451,
                    'new_R00451': new_r00451,
                    'fold_change': new_r00451 / base_r00451 if abs(base_r00451) > 1e-15 else float('inf'),
                    'base_biomass': base_biomass,
                    'new_biomass': new_biomass,
                    'biomass_retained_pct': (new_biomass / base_biomass * 100) if base_biomass > 0 else 0,
                    'verified': new_r00451 > base_r00451 * 1.001  # >0.1% increase
                }
    
    # Save master ranking
    with open(os.path.join(OUTPUT_DIR, '07_master_target_ranking.csv'), 'w') as f:
        f.write("Rank,Reaction,Composite_Score,N_Screens,N_DAPs,Max_Fold_Change,"
               "OE_Hits,KO_Hits,MIM_Hits,FSEOF_Hits,Screen_Types,DAPs_Hit\n")
        for rank, (rxn, scores) in enumerate(ranked, 1):
            types = '+'.join(sorted(scores['screen_types']))
            daps = '+'.join(str(d) for d in sorted(scores['daps_hit']))
            f.write(f"{rank},{rxn},{scores['composite']:.6f},{scores['n_screens']},"
                   f"{scores['n_daps']},{scores['max_fold_any']:.6f},"
                   f"{scores['overexp_hits']},{scores['ko_hits']},"
                   f"{scores['mimicry_hits']},{scores['fseof_hits']},"
                   f"{types},{daps}\n")
    
    # Save verification
    with open(os.path.join(OUTPUT_DIR, '08_target_verification.csv'), 'w') as f:
        f.write("Reaction,DAP,Base_R00451,New_R00451,Fold_Change,Base_Biomass,New_Biomass,Biomass_Retained_Pct,Verified\n")
        for rxn in top_20:
            for dap in sorted(verification.get(rxn, {}).keys()):
                v = verification[rxn][dap]
                f.write(f"{rxn},{dap},{v['base_R00451']:.10e},{v['new_R00451']:.10e},"
                       f"{v['fold_change']:.6f},{v['base_biomass']:.10e},{v['new_biomass']:.10e},"
                       f"{v['biomass_retained_pct']:.2f},{v['verified']}\n")
    
    # Print verification summary
    print(f"\n  Verification Summary (2x overexpression):")
    print(f"  {'Reaction':40s} {'DAPs Verified':>13s} {'Avg FC':>8s} {'Avg Biomass%':>12s}")
    print(f"  {'-'*40} {'-'*13} {'-'*8} {'-'*12}")
    
    for rxn in top_20:
        vals = verification.get(rxn, {})
        if not vals:
            continue
        verified_daps = sum(1 for v in vals.values() if v['verified'])
        total_daps = len(vals)
        avg_fc = np.mean([v['fold_change'] for v in vals.values() if v['fold_change'] < 1e10])
        avg_bio = np.mean([v['biomass_retained_pct'] for v in vals.values()])
        print(f"  {rxn:40s} {verified_daps:4d}/{total_daps:<4d}     {avg_fc:8.2f} {avg_bio:11.1f}%")
    
    return ranked, verification


# ============================================================
# ANALYSIS 8: MECHANISTIC PATHWAY ANALYSIS OF TOP TARGETS
# ============================================================

def analysis_pathway_mechanism(models_wt, baseline_results, top_targets, verification):
    """
    For top verified targets, trace flux changes through the metabolic network
    to explain HOW they increase R00451.
    """
    print("\n" + "="*70)
    print("ANALYSIS 8: MECHANISTIC ANALYSIS OF TOP TARGETS")
    print("="*70)
    
    LYSINE_RELATED = [
        'R00480[K,p]', 'R00355[K,p]', 'R02291[K,p]', 'R02292[K,p]',
        'R04198[K,p]', 'R02735[K,p]', 'R00451[K,p]',  # biosynthesis
        'R00715[K,c]', 'R00716[K,c]', 'R02313[K,c]', 'R03102[K,c]', 'R03658[K,m]',  # degradation
        'cpTransport_C00047[K]', 'ExB_C00047[K]', 'Exchange_C00047[K,L]',  # transport/exchange
        'Seed_Biomass[K]'
    ]
    
    # Central metabolism reactions to track
    CENTRAL_METABOLISM = [
        'R00351[K,c]',  # Pyruvate kinase
        'R00200[K,m]',  # Citrate synthase
        'R01082[K,c]',  # Aspartate aminotransferase cytosol
        'R00355[K,p]',  # Aspartate kinase
        'R00258[K,c]',  # Alanine aminotransferase
    ]
    
    ALL_TRACKED = list(set(LYSINE_RELATED + CENTRAL_METABOLISM))
    
    results = {}
    
    # Focus on DAP 15 and 22 (most relevant) plus one early (DAP 10)
    focus_daps = [10, 15, 22]
    
    for target_rxn, _ in top_targets[:10]:
        results[target_rxn] = {}
        print(f"\n  Target: {target_rxn}")
        
        for dap in focus_daps:
            if dap not in models_wt:
                continue
            model = models_wt[dap]
            base = baseline_results[dap]
            if base['status'] != 'optimal' or base['fluxes'] is None:
                continue
            
            rxn_j = model['rxn_idx'].get(target_rxn)
            if rxn_j is None:
                continue
            
            # Apply 2x overexpression
            extra_ub = {rxn_j: model['v_max'][rxn_j] * 2.0}
            extra_lb = {}
            if model['v_min'][rxn_j] < 0:
                extra_lb[rxn_j] = model['v_min'][rxn_j] * 2.0
            
            obj_val, fluxes, status = solve_fba(
                model, BIOMASS_RXN, 'max',
                extra_lb=extra_lb if extra_lb else None,
                extra_ub=extra_ub
            )
            
            if status != 'optimal' or fluxes is None:
                continue
            
            flux_changes = {}
            for tracked_rxn in ALL_TRACKED:
                t_idx = model['rxn_idx'].get(tracked_rxn)
                if t_idx is not None:
                    orig = base['fluxes'][t_idx]
                    new = fluxes[t_idx]
                    delta = new - orig
                    if abs(orig) > 1e-15:
                        fc = new / orig
                    elif abs(new) > 1e-15:
                        fc = float('inf')
                    else:
                        fc = 1.0
                    flux_changes[tracked_rxn] = {
                        'original': orig, 'new': new, 'delta': delta, 'fold_change': fc
                    }
            
            results[target_rxn][dap] = flux_changes
            
            # Print notable changes
            notable = {k: v for k, v in flux_changes.items() if abs(v['delta']) > 1e-10}
            if notable:
                print(f"    DAP {dap}: {len(notable)} pathway reactions changed")
                for r, v in sorted(notable.items(), key=lambda x: abs(x[1]['delta']), reverse=True)[:5]:
                    print(f"      {r:40s} {v['original']:+.4e} -> {v['new']:+.4e} (FC={v['fold_change']:.2f})")
    
    # Save
    with open(os.path.join(OUTPUT_DIR, '09_mechanistic_analysis.csv'), 'w') as f:
        f.write("Target_Reaction,DAP,Tracked_Reaction,Original_Flux,New_Flux,Delta,Fold_Change\n")
        for target, dap_data in results.items():
            for dap, flux_data in sorted(dap_data.items()):
                for tracked, info in sorted(flux_data.items()):
                    f.write(f"{target},{dap},{tracked},{info['original']:.10e},"
                           f"{info['new']:.10e},{info['delta']:.10e},{info['fold_change']:.6f}\n")
    
    return results


# ============================================================
# MAIN EXECUTION
# ============================================================

def main():
    print("=" * 70)
    print("GENE TARGET IDENTIFICATION FOR LYSINE (R00451) UPREGULATION")
    print("IN WILD-TYPE MAIZE ENDOSPERM")
    print("=" * 70)
    
    # Load all models
    print("\nLoading WT models...")
    models_wt = {}
    for dap, path in sorted(WT_DIRS.items()):
        print(f"  Loading B{dap}...", end=" ")
        models_wt[dap] = parse_model(path)
        print(f"({models_wt[dap]['n_rxn']} rxns, {models_wt[dap]['n_met']} mets)")
    
    print("\nLoading O2 models...")
    global models_o2
    models_o2 = {}
    for dap, path in sorted(O2_DIRS.items()):
        print(f"  Loading O{dap}...", end=" ")
        models_o2[dap] = parse_model(path)
        print(f"({models_o2[dap]['n_rxn']} rxns, {models_o2[dap]['n_met']} mets)")
    
    # Verify target reaction exists
    test_model = models_wt[22]
    if LYSINE_TARGET not in test_model['rxn_idx']:
        print(f"\nERROR: Target reaction {LYSINE_TARGET} not found!")
        print("Available R00451 reactions:")
        for r in test_model['reactions']:
            if 'R00451' in r:
                print(f"  {r}")
        return
    
    print(f"\nTarget reaction: {LYSINE_TARGET} (index {test_model['rxn_idx'][LYSINE_TARGET]})")
    print(f"Biomass reaction: {BIOMASS_RXN} (index {test_model['rxn_idx'][BIOMASS_RXN]})")
    
    # Run analyses
    baseline = analysis_baseline(models_wt)
    
    fseof = analysis_fseof(models_wt, baseline, n_steps=10)
    
    overexp = analysis_overexpression_screen(models_wt, baseline, factor=2.0)
    
    ko = analysis_knockout_screen(models_wt, baseline)
    
    mimicry = analysis_o2_mimicry(models_wt, models_o2, baseline)
    
    combos = analysis_combinations(models_wt, baseline, overexp, ko, mimicry)
    
    ranked, verification = analysis_summary_and_verification(
        models_wt, baseline, overexp, ko, mimicry, fseof)
    
    mechanism = analysis_pathway_mechanism(models_wt, baseline, ranked, verification)
    
    # Final summary
    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)
    print(f"\nOutput files in: {OUTPUT_DIR}/")
    print("  01_baseline_fba.csv              - Baseline WT FBA results")
    print("  02_fseof_targets.csv             - FSEOF co-varying reactions")
    print("  03_overexpression_screen.csv     - Overexpression screen hits")
    print("  04_knockout_screen.csv           - Knockout screen hits")
    print("  04b_lethal_knockouts.csv         - Lethal knockouts")
    print("  05_o2_bound_mimicry.csv          - O2-bound mimicry results")
    print("  06_combination_targets.csv       - Combination target results")
    print("  07_master_target_ranking.csv     - Master ranked target list")
    print("  08_target_verification.csv       - Verification of top targets")
    print("  09_mechanistic_analysis.csv      - Mechanistic pathway analysis")


if __name__ == '__main__':
    main()
