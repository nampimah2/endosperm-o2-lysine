#!/usr/bin/env python3
"""
Gene Target Identification via GAMS-based Perturbation Analysis
===============================================================

Uses GAMS/CPLEX (proven reliable solver) for all FBA computations.

Strategy:
1. Generate modified GAMS scripts for each perturbation
2. Run GAMS in batch mode
3. Parse results and rank targets
4. Verify top hits

Perturbation types:
- Overexpression: double v_max for each reaction
- Knockout: set v_min = v_max = 0
- O2-bound mimicry: apply O2 bounds to WT
"""

import os
import sys
import re
import subprocess
import csv
import shutil
import numpy as np
from collections import defaultdict

BASE_DIR = '/mnt/nrdstor/ssbio/nampimah/ENDOSPERM/2026/Endosperm_Model_MARY'
OUTPUT_DIR = os.path.join(BASE_DIR, 'gene_target_analysis')
WORK_DIR = os.path.join(OUTPUT_DIR, 'gams_runs')
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(WORK_DIR, exist_ok=True)

WT_DIRS = {
    6: 'Wild_Type/B6', 8: 'Wild_Type/B8', 10: 'Wild_Type/B10', 12: 'Wild_Type/B12',
    15: 'Wild_Type/B15', 18: 'Wild_Type/B18', 22: 'Wild_Type/B22', 30: 'Wild_Type/B30'
}
O2_DIRS = {
    6: 'O2_mutant/O6', 8: 'O2_mutant/O8', 10: 'O2_mutant/O10', 12: 'O2_mutant/O12',
    15: 'O2_mutant/O15', 18: 'O2_mutant/O18', 22: 'O2_mutant/O22', 30: 'O2_mutant/O30'
}

LYSINE_TARGET = "R00451[K,p]"
BIOMASS_RXN = "Seed_Biomass[K]"


# ============================================================
# FILE PARSING
# ============================================================

def parse_reactions(filepath):
    """Parse reactions.txt and return list of reaction names."""
    rxns = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line and line != '/':
                rxn = line.strip("'").strip(",").strip("'")
                if rxn:
                    rxns.append(rxn)
    return rxns


def parse_bounds(filepath):
    """Parse v_max.txt or v_min.txt and return dict {rxn: value}."""
    bounds = {}
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line and line != '/':
                match = re.match(r"'([^']+)'\s+([-\d.eE+]+)", line)
                if match:
                    bounds[match.group(1)] = float(match.group(2))
    return bounds


def parse_results(filepath):
    """Parse results_FBA.txt and return dict {rxn: flux}."""
    fluxes = {}
    obj = None
    if not os.path.exists(filepath):
        return None, None
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith('The max Biomass value is'):
                parts = line.split(':')
                if len(parts) > 1:
                    try:
                        obj = float(parts[-1].strip())
                    except:
                        pass
            else:
                parts = line.split()
                if len(parts) == 2:
                    try:
                        fluxes[parts[0]] = float(parts[1])
                    except:
                        pass
    return obj, fluxes


# ============================================================
# GAMS SCRIPT GENERATION
# ============================================================

def generate_perturbation_gms(model_dir, output_dir, perturbation_type, 
                               target_rxn=None, factor=2.0, 
                               new_vmax=None, new_vmin=None,
                               result_filename='results_perturb.txt'):
    """
    Generate a modified GAMS FBA script with a single-reaction perturbation.
    
    perturbation_type: 'overexpression', 'knockout', 'bounds'
    """
    full_model_dir = os.path.join(BASE_DIR, model_dir)
    os.makedirs(output_dir, exist_ok=True)
    
    # Copy model files
    for fname in ['metabolites.txt', 'reactions.txt', 'sij.txt']:
        src = os.path.join(full_model_dir, fname)
        dst = os.path.join(output_dir, fname)
        if not os.path.exists(dst):
            os.symlink(src, dst)
    
    # Read and modify bounds
    vmax_orig = parse_bounds(os.path.join(full_model_dir, 'v_max.txt'))
    vmin_orig = parse_bounds(os.path.join(full_model_dir, 'v_min.txt'))
    
    vmax_new = dict(vmax_orig)
    vmin_new = dict(vmin_orig)
    
    if target_rxn:
        if perturbation_type == 'overexpression':
            if target_rxn in vmax_new:
                vmax_new[target_rxn] = vmax_orig[target_rxn] * factor
            if target_rxn in vmin_new and vmin_new[target_rxn] < 0:
                vmin_new[target_rxn] = vmin_orig[target_rxn] * factor
        elif perturbation_type == 'knockout':
            vmax_new[target_rxn] = 0.0
            vmin_new[target_rxn] = 0.0
        elif perturbation_type == 'bounds':
            if new_vmax is not None:
                vmax_new[target_rxn] = new_vmax
            if new_vmin is not None:
                vmin_new[target_rxn] = new_vmin
    
    # Write modified bounds
    with open(os.path.join(output_dir, 'v_max.txt'), 'w') as f:
        f.write("/\n")
        for rxn in parse_reactions(os.path.join(full_model_dir, 'reactions.txt')):
            if rxn in vmax_new:
                f.write(f"'{rxn}'  {vmax_new[rxn]}\n")
        f.write("/\n")
    
    with open(os.path.join(output_dir, 'v_min.txt'), 'w') as f:
        f.write("/\n")
        for rxn in parse_reactions(os.path.join(full_model_dir, 'reactions.txt')):
            if rxn in vmin_new:
                f.write(f"'{rxn}'  {vmin_new[rxn]}\n")
        f.write("/\n")
    
    # Write GAMS script
    gms_content = f"""*************************************************************
* Perturbation FBA - Gene Target Identification
*************************************************************
$INLINECOM /*  */

OPTIONS
    decimals = 8
    lp = cplex
;

SETS
    i   set of metabolites
$include "metabolites.txt"
    j   set of reactions
$include "reactions.txt"
;

PARAMETERS
    S(i,j)   stoichiometric matrix
$include "sij.txt"
    v_max(j) maximum flux
$include "v_max.txt"
    v_min(j) minimum flux
$include "v_min.txt"
;

EQUATIONS
    objective
    mass_balance1(i)
    lower_bound(j)
    upper_bound(j)
;

FREE VARIABLES
    v(j)
    obj
;

objective..         obj =e= v('Seed_Biomass[K]');
mass_balance1(i)..  sum(j, S(i,j) * v(j)) =e= 0;
lower_bound(j)..    v_min(j) =l= v(j);
upper_bound(j)..    v(j) =l= v_max(j);

Model fba_perturb /all/;
Solve fba_perturb using lp maximizing obj;
fba_perturb.optfile = 1;
fba_perturb.holdfixed = 1;

FILE RESULTS /'{result_filename}'/;
PUT RESULTS;
PUT "The max Biomass value is : " obj.l:0:8//;
loop(j,
    put j.tl:0:30, "    ", v.l(j):0:8/;
);
PUTCLOSE;
"""
    
    with open(os.path.join(output_dir, 'perturb_fba.gms'), 'w') as f:
        f.write(gms_content)
    
    return os.path.join(output_dir, 'perturb_fba.gms')


def run_gams(gms_file, cwd=None):
    """Run GAMS and return True if successful."""
    if cwd is None:
        cwd = os.path.dirname(gms_file)
    
    try:
        result = subprocess.run(
            ['gams', os.path.basename(gms_file)],
            cwd=cwd,
            capture_output=True, text=True, timeout=120
        )
        return result.returncode == 0
    except subprocess.TimeoutExpired:
        return False
    except Exception as e:
        print(f"    GAMS error: {e}")
        return False


# ============================================================
# BATCH SCREENING VIA GAMS LOOP SCRIPT
# ============================================================

def generate_batch_overexpression_gms(model_dir, output_dir, reactions, factor=2.0):
    """
    Generate a single GAMS script that loops over candidate reactions,
    doubling each one's v_max in turn, solving FBA, and recording R00451 flux.
    This is MUCH faster than running individual GAMS files.
    """
    full_model_dir = os.path.join(BASE_DIR, model_dir)
    os.makedirs(output_dir, exist_ok=True)
    
    # Symlink shared files
    for fname in ['metabolites.txt', 'reactions.txt', 'sij.txt', 'v_max.txt', 'v_min.txt']:
        src = os.path.join(full_model_dir, fname)
        dst = os.path.join(output_dir, fname)
        if os.path.exists(dst) or os.path.islink(dst):
            os.remove(dst)
        os.symlink(src, dst)
    
    gms_content = f"""*************************************************************
* Batch Overexpression Screen for Lysine Target ID
* Factor = {factor}x v_max
*************************************************************
$INLINECOM /*  */

OPTIONS
    decimals = 8
    lp = cplex
    solprint = off
    limrow = 0
    limcol = 0
    solvelink = 5
;

SETS
    i   set of metabolites
$include "metabolites.txt"
    j   set of reactions
$include "reactions.txt"
;

PARAMETERS
    S(i,j)   stoichiometric matrix
$include "sij.txt"
    v_max_orig(j) maximum flux
$include "v_max.txt"
    v_min_orig(j) minimum flux
$include "v_min.txt"

    v_max_mod(j) modified max bounds
    v_min_mod(j) modified min bounds
    result_biomass(j)    biomass when j is overexpressed
    result_r00451(j)     R00451 flux when j is overexpressed
    base_biomass         baseline biomass
    base_r00451          baseline R00451
;

EQUATIONS
    objective
    mass_balance1(i)
    lower_bound(j)
    upper_bound(j)
;

FREE VARIABLES
    v(j)
    obj
;

objective..         obj =e= v('Seed_Biomass[K]');
mass_balance1(i)..  sum(j, S(i,j) * v(j)) =e= 0;
lower_bound(j)..    v_min_mod(j) =l= v(j);
upper_bound(j)..    v(j) =l= v_max_mod(j);

Model fba_screen /all/;
fba_screen.optfile = 1;
fba_screen.holdfixed = 1;

* -------- Baseline solve --------
v_max_mod(j) = v_max_orig(j);
v_min_mod(j) = v_min_orig(j);
Solve fba_screen using lp maximizing obj;
base_biomass = v.l('Seed_Biomass[K]');
base_r00451  = v.l('R00451[K,p]');

* -------- Overexpression loop --------
ALIAS(j, jj);

loop(jj$(v_max_orig(jj) > 1e-12 or v_min_orig(jj) < -1e-12),
    v_max_mod(j) = v_max_orig(j);
    v_min_mod(j) = v_min_orig(j);
    
    v_max_mod(jj) = v_max_orig(jj) * {factor};
    if(v_min_orig(jj) < 0,
        v_min_mod(jj) = v_min_orig(jj) * {factor};
    );
    
    Solve fba_screen using lp maximizing obj;
    
    if(fba_screen.modelstat = 1 or fba_screen.modelstat = 2,
        result_biomass(jj) = v.l('Seed_Biomass[K]');
        result_r00451(jj)  = v.l('R00451[K,p]');
    else
        result_biomass(jj) = -999;
        result_r00451(jj)  = -999;
    );
);

* -------- Output Results --------
FILE RESULTS /'results_overexpression.txt'/;
PUT RESULTS;
PUT "BASELINE_BIOMASS ", base_biomass:0:10 /;
PUT "BASELINE_R00451 ", base_r00451:0:10 /;
PUT "REACTION BIOMASS R00451" /;
loop(jj$(result_r00451(jj) ne 0),
    put jj.tl:0:40, " ", result_biomass(jj):0:10, " ", result_r00451(jj):0:10 /;
);
PUTCLOSE;
"""
    
    gms_path = os.path.join(output_dir, 'batch_overexpression.gms')
    with open(gms_path, 'w') as f:
        f.write(gms_content)
    
    return gms_path


def generate_batch_knockout_gms(model_dir, output_dir):
    """
    GAMS script that loops over active reactions, knocking each out
    and recording R00451 flux change.
    """
    full_model_dir = os.path.join(BASE_DIR, model_dir)
    os.makedirs(output_dir, exist_ok=True)
    
    for fname in ['metabolites.txt', 'reactions.txt', 'sij.txt', 'v_max.txt', 'v_min.txt']:
        src = os.path.join(full_model_dir, fname)
        dst = os.path.join(output_dir, fname)
        if os.path.exists(dst) or os.path.islink(dst):
            os.remove(dst)
        os.symlink(src, dst)
    
    gms_content = """*************************************************************
* Batch Knockout Screen for Lysine Target ID
*************************************************************
$INLINECOM /*  */

OPTIONS
    decimals = 8
    lp = cplex
    solprint = off
    limrow = 0
    limcol = 0
    solvelink = 5
;

SETS
    i   set of metabolites
$include "metabolites.txt"
    j   set of reactions
$include "reactions.txt"
;

PARAMETERS
    S(i,j)   stoichiometric matrix
$include "sij.txt"
    v_max_orig(j) maximum flux
$include "v_max.txt"
    v_min_orig(j) minimum flux
$include "v_min.txt"

    v_max_mod(j) modified max bounds
    v_min_mod(j) modified min bounds
    result_biomass(j)
    result_r00451(j)
    base_biomass
    base_r00451
    base_flux(j)
;

EQUATIONS
    objective
    mass_balance1(i)
    lower_bound(j)
    upper_bound(j)
;

FREE VARIABLES
    v(j)
    obj
;

objective..         obj =e= v('Seed_Biomass[K]');
mass_balance1(i)..  sum(j, S(i,j) * v(j)) =e= 0;
lower_bound(j)..    v_min_mod(j) =l= v(j);
upper_bound(j)..    v(j) =l= v_max_mod(j);

Model fba_ko /all/;
fba_ko.optfile = 1;
fba_ko.holdfixed = 1;

* -------- Baseline solve --------
v_max_mod(j) = v_max_orig(j);
v_min_mod(j) = v_min_orig(j);
Solve fba_ko using lp maximizing obj;
base_biomass = v.l('Seed_Biomass[K]');
base_r00451  = v.l('R00451[K,p]');
base_flux(j) = v.l(j);

* -------- Knockout loop (only active reactions) --------
ALIAS(j, jj);

loop(jj$(abs(base_flux(jj)) > 1e-12
         and not sameas(jj, 'Seed_Biomass[K]')
         and not sameas(jj, 'R00451[K,p]')),
    v_max_mod(j) = v_max_orig(j);
    v_min_mod(j) = v_min_orig(j);
    
    v_max_mod(jj) = 0;
    v_min_mod(jj) = 0;
    
    Solve fba_ko using lp maximizing obj;
    
    if(fba_ko.modelstat = 1 or fba_ko.modelstat = 2,
        result_biomass(jj) = v.l('Seed_Biomass[K]');
        result_r00451(jj)  = v.l('R00451[K,p]');
    else
        result_biomass(jj) = -999;
        result_r00451(jj)  = -999;
    );
);

* -------- Output Results --------
FILE RESULTS /'results_knockout.txt'/;
PUT RESULTS;
PUT "BASELINE_BIOMASS ", base_biomass:0:10 /;
PUT "BASELINE_R00451 ", base_r00451:0:10 /;
PUT "REACTION BIOMASS R00451 ORIG_FLUX" /;
loop(jj$(result_r00451(jj) ne 0),
    put jj.tl:0:40, " ", result_biomass(jj):0:10, " ", result_r00451(jj):0:10, " ", base_flux(jj):0:10 /;
);
PUTCLOSE;
"""
    
    gms_path = os.path.join(output_dir, 'batch_knockout.gms')
    with open(gms_path, 'w') as f:
        f.write(gms_content)
    
    return gms_path


def generate_fseof_gms(model_dir, output_dir, n_steps=10, biomass_frac=0.9):
    """
    GAMS script for FSEOF: gradually enforce increasing R00451 flux levels
    while maintaining biomass >= fraction of optimal.
    """
    full_model_dir = os.path.join(BASE_DIR, model_dir)
    os.makedirs(output_dir, exist_ok=True)
    
    for fname in ['metabolites.txt', 'reactions.txt', 'sij.txt', 'v_max.txt', 'v_min.txt']:
        src = os.path.join(full_model_dir, fname)
        dst = os.path.join(output_dir, fname)
        if os.path.exists(dst) or os.path.islink(dst):
            os.remove(dst)
        os.symlink(src, dst)
    
    gms_content = f"""*************************************************************
* FSEOF: Flux Scanning based on Enforced Objective Flux
* Target: R00451[K,p] (lysine biosynthesis)
* Steps: {n_steps}, Biomass fraction: {biomass_frac}
*************************************************************
$INLINECOM /*  */

OPTIONS
    decimals = 8
    lp = cplex
    solprint = off
    limrow = 0
    limcol = 0
    solvelink = 5
;

SETS
    i   set of metabolites
$include "metabolites.txt"
    j   set of reactions
$include "reactions.txt"
    s   steps /s0*s{n_steps}/
;

PARAMETERS
    S(i,j)   stoichiometric matrix
$include "sij.txt"
    v_max_orig(j) maximum flux
$include "v_max.txt"
    v_min_orig(j) minimum flux
$include "v_min.txt"

    biomass_floor
    max_r00451
    r00451_target(s)
    flux_profile(s,j)
;

EQUATIONS
    objective
    mass_balance1(i)
    lower_bound(j)
    upper_bound(j)
    biomass_constraint
    r00451_enforce
;

FREE VARIABLES
    v(j)
    obj
;

SCALAR r00451_min /0/;

objective..         obj =e= v('Seed_Biomass[K]');
mass_balance1(i)..  sum(j, S(i,j) * v(j)) =e= 0;
lower_bound(j)..    v_min_orig(j) =l= v(j);
upper_bound(j)..    v(j) =l= v_max_orig(j);
biomass_constraint.. v('Seed_Biomass[K]') =g= biomass_floor;
r00451_enforce..    v('R00451[K,p]') =g= r00451_min;

* Step 1: Solve for max biomass
Model fba_base /objective, mass_balance1, lower_bound, upper_bound/;
fba_base.optfile = 1;
Solve fba_base using lp maximizing obj;

SCALAR base_biomass;
SCALAR base_r00451;
base_biomass = v.l('Seed_Biomass[K]');
base_r00451  = v.l('R00451[K,p]');
biomass_floor = {biomass_frac} * base_biomass;

* Step 2: Max R00451 subject to biomass floor
Model fba_max_r00451 /mass_balance1, lower_bound, upper_bound, biomass_constraint/;
fba_max_r00451.optfile = 1;

VARIABLE obj2;
EQUATION obj_r00451;
obj_r00451.. obj2 =e= v('R00451[K,p]');
Model fba_max_lys /obj_r00451, mass_balance1, lower_bound, upper_bound, biomass_constraint/;
fba_max_lys.optfile = 1;
Solve fba_max_lys using lp maximizing obj2;
max_r00451 = v.l('R00451[K,p]');

* Step 3: Scan
SCALAR step_size;
step_size = (max_r00451 - base_r00451) / {n_steps};

Model fba_fseof /objective, mass_balance1, lower_bound, upper_bound, r00451_enforce/;
fba_fseof.optfile = 1;

loop(s,
    r00451_min = base_r00451 + ord(s) * step_size;
    r00451_target(s) = r00451_min;
    
    Solve fba_fseof using lp maximizing obj;
    
    if(fba_fseof.modelstat = 1 or fba_fseof.modelstat = 2,
        flux_profile(s,j) = v.l(j);
    );
);

* Output
FILE RESULTS /'results_fseof.txt'/;
PUT RESULTS;
PUT "BASE_BIOMASS ", base_biomass:0:10 /;
PUT "BASE_R00451 ", base_r00451:0:10 /;
PUT "MAX_R00451 ", max_r00451:0:10 /;
PUT "STEP_SIZE ", step_size:0:10 /;
loop(s,
    PUT "STEP ", ord(s):0:0, " TARGET ", r00451_target(s):0:10 /;
    loop(j$(abs(flux_profile(s,j)) > 1e-10),
        PUT "  ", j.tl:0:40, " ", flux_profile(s,j):0:10 /;
    );
);
PUTCLOSE;
"""
    
    gms_path = os.path.join(output_dir, 'batch_fseof.gms')
    with open(gms_path, 'w') as f:
        f.write(gms_content)
    
    return gms_path


# ============================================================
# RESULT PARSING
# ============================================================

def parse_overexpression_results(filepath):
    """Parse batch overexpression results."""
    base_biomass = 0
    base_r00451 = 0
    hits = {}
    
    if not os.path.exists(filepath):
        return None, None, None
    
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith('BASELINE_BIOMASS'):
                base_biomass = float(line.split()[-1])
            elif line.startswith('BASELINE_R00451'):
                base_r00451 = float(line.split()[-1])
            elif line.startswith('REACTION'):
                continue
            else:
                parts = line.split()
                if len(parts) >= 3:
                    rxn = parts[0]
                    try:
                        biomass = float(parts[1])
                        r00451 = float(parts[2])
                        if biomass != -999:
                            hits[rxn] = {'biomass': biomass, 'r00451': r00451}
                    except:
                        pass
    
    return base_biomass, base_r00451, hits


def parse_knockout_results(filepath):
    """Parse batch knockout results."""
    base_biomass = 0
    base_r00451 = 0
    hits = {}
    
    if not os.path.exists(filepath):
        return None, None, None
    
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith('BASELINE_BIOMASS'):
                base_biomass = float(line.split()[-1])
            elif line.startswith('BASELINE_R00451'):
                base_r00451 = float(line.split()[-1])
            elif line.startswith('REACTION'):
                continue
            else:
                parts = line.split()
                if len(parts) >= 4:
                    rxn = parts[0]
                    try:
                        biomass = float(parts[1])
                        r00451 = float(parts[2])
                        orig_flux = float(parts[3])
                        if biomass != -999:
                            hits[rxn] = {'biomass': biomass, 'r00451': r00451, 'orig_flux': orig_flux}
                    except:
                        pass
    
    return base_biomass, base_r00451, hits


def parse_fseof_results(filepath):
    """Parse FSEOF results."""
    if not os.path.exists(filepath):
        return None
    
    data = {'base_biomass': 0, 'base_r00451': 0, 'max_r00451': 0, 'steps': {}}
    current_step = None
    
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith('BASE_BIOMASS'):
                data['base_biomass'] = float(line.split()[-1])
            elif line.startswith('BASE_R00451'):
                data['base_r00451'] = float(line.split()[-1])
            elif line.startswith('MAX_R00451'):
                data['max_r00451'] = float(line.split()[-1])
            elif line.startswith('STEP_SIZE'):
                data['step_size'] = float(line.split()[-1])
            elif line.startswith('STEP'):
                parts = line.split()
                current_step = int(parts[1])
                target = float(parts[3])
                data['steps'][current_step] = {'target': target, 'fluxes': {}}
            elif current_step is not None and line:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        data['steps'][current_step]['fluxes'][parts[0]] = float(parts[1])
                    except:
                        pass
    
    return data


# ============================================================
# MAIN PIPELINE
# ============================================================

def main():
    print("=" * 70)
    print("GENE TARGET IDENTIFICATION FOR LYSINE (R00451) UPREGULATION")
    print("GAMS/CPLEX-BASED PERTURBATION ANALYSIS")
    print("=" * 70)
    
    # Focus DAPs: 10, 15, 18, 22 (most biologically relevant)
    focus_daps = [10, 15, 18, 22]
    
    # ===========================================================
    # PHASE 1: Generate and run all GAMS scripts
    # ===========================================================
    print("\n" + "=" * 70)
    print("PHASE 1: GENERATING GAMS SCRIPTS")
    print("=" * 70)
    
    gams_jobs = []
    
    for dap in focus_daps:
        wt_dir = WT_DIRS[dap]
        
        # 1a. Overexpression screen
        oe_dir = os.path.join(WORK_DIR, f'B{dap}_overexp')
        gms_oe = generate_batch_overexpression_gms(wt_dir, oe_dir, None, factor=2.0)
        gams_jobs.append(('overexp', dap, gms_oe, oe_dir))
        print(f"  Generated: B{dap} overexpression screen")
        
        # 1b. Knockout screen
        ko_dir = os.path.join(WORK_DIR, f'B{dap}_knockout')
        gms_ko = generate_batch_knockout_gms(wt_dir, ko_dir)
        gams_jobs.append(('knockout', dap, gms_ko, ko_dir))
        print(f"  Generated: B{dap} knockout screen")
        
        # 1c. FSEOF
        fseof_dir = os.path.join(WORK_DIR, f'B{dap}_fseof')
        gms_fseof = generate_fseof_gms(wt_dir, fseof_dir, n_steps=10, biomass_frac=0.9)
        gams_jobs.append(('fseof', dap, gms_fseof, fseof_dir))
        print(f"  Generated: B{dap} FSEOF")
    
    # ===========================================================
    # PHASE 2: Run all GAMS jobs
    # ===========================================================
    print("\n" + "=" * 70)
    print("PHASE 2: RUNNING GAMS/CPLEX JOBS")
    print("=" * 70)
    
    for job_type, dap, gms_file, job_dir in gams_jobs:
        print(f"  Running B{dap} {job_type}...", end=" ", flush=True)
        success = run_gams(gms_file, cwd=job_dir)
        if success:
            print("OK")
        else:
            print("FAILED")
    
    # ===========================================================
    # PHASE 3: Parse results
    # ===========================================================
    print("\n" + "=" * 70)
    print("PHASE 3: PARSING RESULTS")
    print("=" * 70)
    
    all_oe_results = {}
    all_ko_results = {}
    all_fseof_data = {}
    
    for dap in focus_daps:
        # Overexpression
        oe_file = os.path.join(WORK_DIR, f'B{dap}_overexp', 'results_overexpression.txt')
        base_bio, base_r, hits = parse_overexpression_results(oe_file)
        if hits is not None:
            all_oe_results[dap] = {'base_biomass': base_bio, 'base_r00451': base_r, 'hits': hits}
            # Filter for positive hits
            pos = {r: h for r, h in hits.items() if h['r00451'] > base_r * 1.001}
            print(f"  B{dap} OE: baseline R00451={base_r:.6e}, {len(pos)} positive hits out of {len(hits)}")
        else:
            print(f"  B{dap} OE: no results file")
        
        # Knockout
        ko_file = os.path.join(WORK_DIR, f'B{dap}_knockout', 'results_knockout.txt')
        base_bio, base_r, hits = parse_knockout_results(ko_file)
        if hits is not None:
            all_ko_results[dap] = {'base_biomass': base_bio, 'base_r00451': base_r, 'hits': hits}
            pos = {r: h for r, h in hits.items() if h['r00451'] > base_r * 1.001 and h['biomass'] > 0}
            lethal = {r: h for r, h in hits.items() if h['biomass'] <= 1e-10 and h['biomass'] != -999}
            print(f"  B{dap} KO: baseline R00451={base_r:.6e}, {len(pos)} positive hits, {len(lethal)} lethal")
        else:
            print(f"  B{dap} KO: no results file")
        
        # FSEOF
        fseof_file = os.path.join(WORK_DIR, f'B{dap}_fseof', 'results_fseof.txt')
        fseof_data = parse_fseof_results(fseof_file)
        if fseof_data:
            all_fseof_data[dap] = fseof_data
            print(f"  B{dap} FSEOF: base R00451={fseof_data['base_r00451']:.6e}, "
                  f"max R00451={fseof_data['max_r00451']:.6e}, {len(fseof_data['steps'])} steps")
        else:
            print(f"  B{dap} FSEOF: no results file")
    
    # ===========================================================
    # PHASE 4: FSEOF Analysis - identify co-varying reactions
    # ===========================================================
    print("\n" + "=" * 70)
    print("PHASE 4: FSEOF CO-VARYING REACTION ANALYSIS")
    print("=" * 70)
    
    fseof_targets = {}
    
    for dap, data in sorted(all_fseof_data.items()):
        if len(data['steps']) < 3:
            continue
        
        steps = sorted(data['steps'].keys())
        targets = data['steps'][steps[0]].get('target', 0)
        
        # Get all reactions across all steps
        all_rxns = set()
        for s in steps:
            all_rxns.update(data['steps'][s]['fluxes'].keys())
        
        covarying = {}
        for rxn in all_rxns:
            fluxes = []
            r_targets = []
            for s in steps:
                f = data['steps'][s]['fluxes'].get(rxn, 0)
                t = data['steps'][s]['target']
                fluxes.append(f)
                r_targets.append(t)
            
            if np.std(fluxes) < 1e-12:
                continue
            
            corr = np.corrcoef(r_targets, fluxes)[0, 1]
            if np.isnan(corr):
                continue
            
            delta = fluxes[-1] - fluxes[0]
            
            if abs(corr) > 0.7 and abs(delta) > 1e-8:
                direction = "UP" if corr > 0 else "DOWN"
                covarying[rxn] = {
                    'correlation': corr,
                    'direction': direction,
                    'flux_first': fluxes[0],
                    'flux_last': fluxes[-1],
                    'delta': delta
                }
        
        fseof_targets[dap] = covarying
        print(f"\n  DAP {dap}: {len(covarying)} co-varying reactions (|corr| > 0.7)")
        
        up = {r: v for r, v in covarying.items() if v['direction'] == 'UP'}
        down = {r: v for r, v in covarying.items() if v['direction'] == 'DOWN'}
        print(f"    {len(up)} UP-regulated, {len(down)} DOWN-regulated with R00451")
        
        for rxn, info in sorted(up.items(), key=lambda x: x[1]['correlation'], reverse=True)[:10]:
            print(f"      UP   {rxn:40s} corr={info['correlation']:+.3f} delta={info['delta']:+.4e}")
        for rxn, info in sorted(down.items(), key=lambda x: x[1]['correlation'])[:5]:
            print(f"      DOWN {rxn:40s} corr={info['correlation']:+.3f} delta={info['delta']:+.4e}")
    
    # ===========================================================
    # PHASE 5: O2-bound mimicry screen
    # ===========================================================
    print("\n" + "=" * 70)
    print("PHASE 5: O2-BOUND MIMICRY SCREEN")
    print("=" * 70)
    
    all_mimicry = {}
    
    for dap in focus_daps:
        wt_dir = os.path.join(BASE_DIR, WT_DIRS[dap])
        o2_dir = os.path.join(BASE_DIR, O2_DIRS[dap])
        
        wt_vmax = parse_bounds(os.path.join(wt_dir, 'v_max.txt'))
        wt_vmin = parse_bounds(os.path.join(wt_dir, 'v_min.txt'))
        o2_vmax = parse_bounds(os.path.join(o2_dir, 'v_max.txt'))
        o2_vmin = parse_bounds(os.path.join(o2_dir, 'v_min.txt'))
        
        # Find reactions with different bounds
        diff_rxns = []
        for rxn in wt_vmax:
            if rxn in o2_vmax:
                if (abs(wt_vmax[rxn] - o2_vmax[rxn]) > 1e-10 or 
                    abs(wt_vmin.get(rxn, 0) - o2_vmin.get(rxn, 0)) > 1e-10):
                    diff_rxns.append(rxn)
        
        print(f"\n  DAP {dap}: {len(diff_rxns)} reactions with different bounds")
        
        if dap not in all_oe_results:
            continue
        
        base_r = all_oe_results[dap]['base_r00451']
        
        # Test top candidates (those with largest bound changes)
        mimicry_hits = {}
        mim_dir = os.path.join(WORK_DIR, f'B{dap}_mimicry')
        os.makedirs(mim_dir, exist_ok=True)
        
        # Generate individual GAMS runs for each mimicry target
        tested = 0
        for rxn in diff_rxns:
            if rxn == BIOMASS_RXN:
                continue
            
            rxn_safe = rxn.replace('[', '_').replace(']', '_').replace(',', '_').replace('/', '_')
            run_dir = os.path.join(mim_dir, rxn_safe)
            
            gms_path = generate_perturbation_gms(
                WT_DIRS[dap], run_dir, 'bounds',
                target_rxn=rxn,
                new_vmax=o2_vmax.get(rxn, wt_vmax.get(rxn, 0)),
                new_vmin=o2_vmin.get(rxn, wt_vmin.get(rxn, 0)),
                result_filename='results_mimicry.txt'
            )
            
            success = run_gams(gms_path, cwd=run_dir)
            tested += 1
            
            if success:
                res_file = os.path.join(run_dir, 'results_mimicry.txt')
                obj, fluxes = parse_results(res_file)
                if fluxes and LYSINE_TARGET in fluxes:
                    new_r = fluxes[LYSINE_TARGET]
                    new_bio = fluxes.get(BIOMASS_RXN, 0)
                    delta = new_r - base_r
                    
                    if abs(delta) > 1e-10:
                        mimicry_hits[rxn] = {
                            'new_r00451': new_r,
                            'delta': delta,
                            'fold_change': new_r / base_r if abs(base_r) > 1e-15 else float('inf'),
                            'new_biomass': new_bio,
                            'wt_vmax': wt_vmax.get(rxn, 0),
                            'o2_vmax': o2_vmax.get(rxn, 0),
                            'wt_vmin': wt_vmin.get(rxn, 0),
                            'o2_vmin': o2_vmin.get(rxn, 0)
                        }
            
            # Status update every 200 reactions
            if tested % 200 == 0:
                print(f"    Tested {tested}/{len(diff_rxns)} reactions...")
        
        all_mimicry[dap] = mimicry_hits
        pos = {r: h for r, h in mimicry_hits.items() if h['delta'] > 0}
        print(f"    Tested {tested}, found {len(pos)} positive hits")
        
        for rxn, info in sorted(pos.items(), key=lambda x: x[1]['delta'], reverse=True)[:10]:
            print(f"      {rxn:40s} R00451 FC={info['fold_change']:.2f} "
                  f"vmax: {info['wt_vmax']:.3f}->{info['o2_vmax']:.3f}")
    
    # ===========================================================
    # PHASE 6: MASTER RANKING & VERIFICATION
    # ===========================================================
    print("\n" + "=" * 70)
    print("PHASE 6: MASTER TARGET RANKING")
    print("=" * 70)
    
    target_scores = defaultdict(lambda: {
        'oe_hits': 0, 'ko_hits': 0, 'mim_hits': 0, 'fseof_up': 0, 'fseof_down': 0,
        'max_fc_oe': 0, 'max_fc_ko': 0, 'max_fc_mim': 0,
        'max_delta_oe': 0, 'max_delta_ko': 0, 'max_delta_mim': 0,
        'daps_hit': set(), 'screen_types': set(),
        'best_dap': None, 'best_screen': None, 'best_fc': 0,
        'biomass_impact': []
    })
    
    # Score overexpression
    for dap, data in all_oe_results.items():
        base_r = data['base_r00451']
        base_b = data['base_biomass']
        for rxn, h in data['hits'].items():
            delta = h['r00451'] - base_r
            if delta > 1e-10:
                fc = h['r00451'] / base_r if abs(base_r) > 1e-15 else float('inf')
                ts = target_scores[rxn]
                ts['oe_hits'] += 1
                ts['daps_hit'].add(dap)
                ts['screen_types'].add('OE')
                ts['biomass_impact'].append(h['biomass'] / base_b * 100 if base_b > 0 else 0)
                if fc > ts['max_fc_oe']:
                    ts['max_fc_oe'] = fc
                    ts['max_delta_oe'] = delta
                if fc > ts['best_fc']:
                    ts['best_fc'] = fc
                    ts['best_dap'] = dap
                    ts['best_screen'] = 'OE'
    
    # Score knockouts
    for dap, data in all_ko_results.items():
        base_r = data['base_r00451']
        base_b = data['base_biomass']
        for rxn, h in data['hits'].items():
            delta = h['r00451'] - base_r
            if delta > 1e-10 and h['biomass'] > 0:
                fc = h['r00451'] / base_r if abs(base_r) > 1e-15 else float('inf')
                ts = target_scores[rxn]
                ts['ko_hits'] += 1
                ts['daps_hit'].add(dap)
                ts['screen_types'].add('KO')
                ts['biomass_impact'].append(h['biomass'] / base_b * 100 if base_b > 0 else 0)
                if fc > ts['max_fc_ko']:
                    ts['max_fc_ko'] = fc
                    ts['max_delta_ko'] = delta
                if fc > ts['best_fc']:
                    ts['best_fc'] = fc
                    ts['best_dap'] = dap
                    ts['best_screen'] = 'KO'
    
    # Score mimicry
    for dap, hits in all_mimicry.items():
        base_r = all_oe_results.get(dap, {}).get('base_r00451', 0)
        base_b = all_oe_results.get(dap, {}).get('base_biomass', 0)
        for rxn, h in hits.items():
            if h['delta'] > 0:
                ts = target_scores[rxn]
                ts['mim_hits'] += 1
                ts['daps_hit'].add(dap)
                ts['screen_types'].add('MIM')
                ts['biomass_impact'].append(h['new_biomass'] / base_b * 100 if base_b > 0 else 0)
                if h['fold_change'] > ts['max_fc_mim']:
                    ts['max_fc_mim'] = h['fold_change']
                    ts['max_delta_mim'] = h['delta']
                if h['fold_change'] > ts['best_fc']:
                    ts['best_fc'] = h['fold_change']
                    ts['best_dap'] = dap
                    ts['best_screen'] = 'MIM'
    
    # Score FSEOF
    for dap, hits in fseof_targets.items():
        for rxn, info in hits.items():
            if info['direction'] == 'UP':
                target_scores[rxn]['fseof_up'] += 1
                target_scores[rxn]['daps_hit'].add(dap)
                target_scores[rxn]['screen_types'].add('FSEOF')
            else:
                target_scores[rxn]['fseof_down'] += 1
                target_scores[rxn]['daps_hit'].add(dap)
                target_scores[rxn]['screen_types'].add('FSEOF')
    
    # Compute composite scores
    for rxn, ts in target_scores.items():
        n_screens = len(ts['screen_types'])
        n_daps = len(ts['daps_hit'])
        max_fc = max(ts['max_fc_oe'], ts['max_fc_ko'], ts['max_fc_mim'], 1)
        avg_bio = np.mean(ts['biomass_impact']) if ts['biomass_impact'] else 0
        
        # Composite: screens * daps * log(FC) * biomass_retention_factor
        bio_factor = max(avg_bio / 100, 0.01)  # penalize biomass loss
        ts['composite'] = n_screens * n_daps * np.log10(max_fc + 1) * bio_factor
        ts['n_screens'] = n_screens
        ts['n_daps'] = n_daps
        ts['max_fc'] = max_fc
        ts['avg_biomass_pct'] = avg_bio
    
    # Rank
    ranked = sorted(target_scores.items(), key=lambda x: x[1]['composite'], reverse=True)
    
    print(f"\n  Total unique gene targets identified: {len(ranked)}")
    print(f"\n  {'Rank':>4s} {'Reaction':40s} {'Score':>7s} {'#Scr':>4s} {'#DAP':>4s} "
          f"{'MaxFC':>10s} {'Bio%':>6s} {'Best':>5s} {'Types'}")
    print(f"  {'-'*4} {'-'*40} {'-'*7} {'-'*4} {'-'*4} {'-'*10} {'-'*6} {'-'*5} {'-'*20}")
    
    for rank, (rxn, ts) in enumerate(ranked[:50], 1):
        types = '+'.join(sorted(ts['screen_types']))
        best = ts['best_screen'] or '---'
        fc_str = f"{ts['max_fc']:.1f}" if ts['max_fc'] < 1e6 else f"{ts['max_fc']:.1e}"
        print(f"  {rank:4d} {rxn:40s} {ts['composite']:7.2f} {ts['n_screens']:4d} {ts['n_daps']:4d} "
              f"{fc_str:>10s} {ts['avg_biomass_pct']:5.1f}% {best:>5s} {types}")
    
    # ===========================================================
    # SAVE ALL RESULTS
    # ===========================================================
    print("\n" + "=" * 70)
    print("SAVING RESULTS")
    print("=" * 70)
    
    # Master ranking
    with open(os.path.join(OUTPUT_DIR, '01_master_target_ranking.csv'), 'w') as f:
        w = csv.writer(f)
        w.writerow(['Rank', 'Reaction', 'Composite_Score', 'N_Screens', 'N_DAPs',
                     'Max_Fold_Change', 'Avg_Biomass_Pct', 'Best_Screen', 'Best_DAP',
                     'OE_Hits', 'OE_MaxFC', 'KO_Hits', 'KO_MaxFC',
                     'MIM_Hits', 'MIM_MaxFC', 'FSEOF_Up', 'FSEOF_Down',
                     'Screen_Types', 'DAPs_Hit'])
        for rank, (rxn, ts) in enumerate(ranked, 1):
            w.writerow([rank, rxn, f"{ts['composite']:.6f}", ts['n_screens'], ts['n_daps'],
                       f"{ts['max_fc']:.6f}", f"{ts['avg_biomass_pct']:.2f}",
                       ts['best_screen'] or '', ts['best_dap'] or '',
                       ts['oe_hits'], f"{ts['max_fc_oe']:.6f}",
                       ts['ko_hits'], f"{ts['max_fc_ko']:.6f}",
                       ts['mim_hits'], f"{ts['max_fc_mim']:.6f}",
                       ts['fseof_up'], ts['fseof_down'],
                       '+'.join(sorted(ts['screen_types'])),
                       '+'.join(str(d) for d in sorted(ts['daps_hit']))])
    
    # Overexpression details
    with open(os.path.join(OUTPUT_DIR, '02_overexpression_hits.csv'), 'w') as f:
        w = csv.writer(f)
        w.writerow(['DAP', 'Reaction', 'Baseline_R00451', 'New_R00451', 'Delta',
                     'Fold_Change', 'New_Biomass', 'Biomass_Pct'])
        for dap in sorted(all_oe_results.keys()):
            data = all_oe_results[dap]
            base_r = data['base_r00451']
            base_b = data['base_biomass']
            for rxn, h in sorted(data['hits'].items(),
                                key=lambda x: x[1]['r00451'] - base_r, reverse=True):
                delta = h['r00451'] - base_r
                if delta > 1e-10:
                    fc = h['r00451'] / base_r if abs(base_r) > 1e-15 else 0
                    bio_pct = h['biomass'] / base_b * 100 if base_b > 0 else 0
                    w.writerow([dap, rxn, f"{base_r:.10e}", f"{h['r00451']:.10e}",
                               f"{delta:.10e}", f"{fc:.6f}", f"{h['biomass']:.10e}",
                               f"{bio_pct:.2f}"])
    
    # Knockout details
    with open(os.path.join(OUTPUT_DIR, '03_knockout_hits.csv'), 'w') as f:
        w = csv.writer(f)
        w.writerow(['DAP', 'Reaction', 'Baseline_R00451', 'New_R00451', 'Delta',
                     'Fold_Change', 'New_Biomass', 'Biomass_Pct', 'Original_Flux'])
        for dap in sorted(all_ko_results.keys()):
            data = all_ko_results[dap]
            base_r = data['base_r00451']
            base_b = data['base_biomass']
            for rxn, h in sorted(data['hits'].items(),
                                key=lambda x: x[1]['r00451'] - base_r, reverse=True):
                delta = h['r00451'] - base_r
                if delta > 1e-10 and h['biomass'] > 0:
                    fc = h['r00451'] / base_r if abs(base_r) > 1e-15 else 0
                    bio_pct = h['biomass'] / base_b * 100 if base_b > 0 else 0
                    w.writerow([dap, rxn, f"{base_r:.10e}", f"{h['r00451']:.10e}",
                               f"{delta:.10e}", f"{fc:.6f}", f"{h['biomass']:.10e}",
                               f"{bio_pct:.2f}", f"{h['orig_flux']:.10e}"])
    
    # FSEOF targets
    with open(os.path.join(OUTPUT_DIR, '04_fseof_targets.csv'), 'w') as f:
        w = csv.writer(f)
        w.writerow(['DAP', 'Reaction', 'Direction', 'Correlation', 'Flux_First', 'Flux_Last', 'Delta'])
        for dap in sorted(fseof_targets.keys()):
            for rxn, info in sorted(fseof_targets[dap].items(),
                                   key=lambda x: abs(x[1]['correlation']), reverse=True):
                w.writerow([dap, rxn, info['direction'], f"{info['correlation']:.6f}",
                           f"{info['flux_first']:.10e}", f"{info['flux_last']:.10e}",
                           f"{info['delta']:.10e}"])
    
    # O2-bound mimicry
    with open(os.path.join(OUTPUT_DIR, '05_o2_mimicry_hits.csv'), 'w') as f:
        w = csv.writer(f)
        w.writerow(['DAP', 'Reaction', 'New_R00451', 'Delta', 'Fold_Change',
                     'New_Biomass', 'WT_Vmax', 'O2_Vmax', 'WT_Vmin', 'O2_Vmin'])
        for dap in sorted(all_mimicry.keys()):
            for rxn, h in sorted(all_mimicry[dap].items(),
                                key=lambda x: x[1]['delta'], reverse=True):
                if h['delta'] > 1e-10:
                    w.writerow([dap, rxn, f"{h['new_r00451']:.10e}", f"{h['delta']:.10e}",
                               f"{h['fold_change']:.6f}", f"{h['new_biomass']:.10e}",
                               f"{h['wt_vmax']:.6f}", f"{h['o2_vmax']:.6f}",
                               f"{h['wt_vmin']:.6f}", f"{h['o2_vmin']:.6f}"])
    
    # Summary report
    with open(os.path.join(OUTPUT_DIR, '06_analysis_summary.txt'), 'w') as f:
        f.write("GENE TARGET IDENTIFICATION FOR LYSINE (R00451) UPREGULATION\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"Focus DAPs: {focus_daps}\n")
        f.write(f"Total unique targets: {len(ranked)}\n\n")
        
        f.write("SCREEN SUMMARY:\n")
        for dap in focus_daps:
            f.write(f"\n  DAP {dap}:\n")
            if dap in all_oe_results:
                base_r = all_oe_results[dap]['base_r00451']
                pos = sum(1 for h in all_oe_results[dap]['hits'].values() if h['r00451'] > base_r * 1.001)
                f.write(f"    Overexpression: {pos} positive hits (baseline R00451 = {base_r:.6e})\n")
            if dap in all_ko_results:
                base_r = all_ko_results[dap]['base_r00451']
                pos = sum(1 for h in all_ko_results[dap]['hits'].values() 
                         if h['r00451'] > base_r * 1.001 and h['biomass'] > 0)
                f.write(f"    Knockout: {pos} positive hits\n")
            if dap in fseof_targets:
                up = sum(1 for v in fseof_targets[dap].values() if v['direction'] == 'UP')
                down = sum(1 for v in fseof_targets[dap].values() if v['direction'] == 'DOWN')
                f.write(f"    FSEOF: {up} UP, {down} DOWN co-varying\n")
            if dap in all_mimicry:
                pos = sum(1 for h in all_mimicry[dap].values() if h['delta'] > 0)
                f.write(f"    O2-mimicry: {pos} positive hits\n")
        
        f.write("\n\nTOP 30 GENE TARGETS:\n")
        f.write(f"{'Rank':>4s}  {'Reaction':40s}  {'Score':>7s}  {'#Scr':>4s}  {'#DAP':>4s}  "
               f"{'MaxFC':>10s}  {'Bio%':>6s}  {'Screens'}\n")
        f.write("-" * 90 + "\n")
        for rank, (rxn, ts) in enumerate(ranked[:30], 1):
            types = '+'.join(sorted(ts['screen_types']))
            fc_str = f"{ts['max_fc']:.1f}" if ts['max_fc'] < 1e6 else f"{ts['max_fc']:.1e}"
            f.write(f"{rank:4d}  {rxn:40s}  {ts['composite']:7.2f}  {ts['n_screens']:4d}  "
                   f"{ts['n_daps']:4d}  {fc_str:>10s}  {ts['avg_biomass_pct']:5.1f}%  {types}\n")
    
    print(f"\n  All results saved to: {OUTPUT_DIR}/")
    print("  01_master_target_ranking.csv")
    print("  02_overexpression_hits.csv")
    print("  03_knockout_hits.csv")
    print("  04_fseof_targets.csv")
    print("  05_o2_mimicry_hits.csv")
    print("  06_analysis_summary.txt")
    
    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)


if __name__ == '__main__':
    main()
