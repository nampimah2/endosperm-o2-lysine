#!/usr/bin/env python3
"""
Preprocessing script for OptKnock analysis.
Reads model reactions and classifies them as protected (non-knockout) or
candidate (knockout-eligible) reactions.

Generates:
  - protected_reactions.txt  (GAMS set include file)
  - candidates_reactions.txt (GAMS set include file)

Protected reactions (NOT available for knockout):
  - Seed_Biomass[K]        : biomass objective
  - R00451[K,p]            : lysine biosynthesis target (DAP decarboxylase)
  - Exchange_*             : model exchange/boundary reactions
  - ExB_*                  : boundary export reactions
  - Exe*                   : extracellular exchange reactions
  - Trans*                 : tissue-level transport reactions
  - PhloemTransport_*      : phloem transport
  - PhloemImport_*         : phloem import

Candidate reactions (available for knockout):
  - MR* / R*               : internal metabolic reactions
  - cpTransport_*          : cytosol-plastid transport
  - cmTransport_*          : cytosol-mitochondria transport
  - cvTransport_*          : cytosol-vacuole transport
  - cxTransport_*          : cytosol-peroxisome transport
  - pmTransport_*          : plasma membrane transport
  - Other intracellular metabolic/transport reactions

Usage:
    python preprocess_optknock.py <data_dir>
    e.g., python preprocess_optknock.py ../Wild_Type/B18
"""

import sys
import os
import re


def read_reactions(data_dir):
    """Read reaction names from reactions.txt"""
    reactions = []
    filepath = os.path.join(data_dir, 'reactions.txt')
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('/') or line == '':
                continue
            # Remove trailing /
            if line.endswith('/'):
                line = line[:-1].strip()
            if line == '':
                continue
            # Remove quotes
            rxn = line.strip("'").strip('"')
            if rxn:
                reactions.append(rxn)
    return reactions


def read_bounds(data_dir, filename):
    """Read flux bounds from v_min.txt or v_max.txt"""
    bounds = {}
    filepath = os.path.join(data_dir, filename)
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('/') or line == '' or line == '/':
                continue
            if line.endswith('/'):
                line = line[:-1].strip()
            parts = line.split('\t')
            if len(parts) >= 2:
                rxn = parts[0].strip("'").strip('"')
                val = float(parts[1])
                bounds[rxn] = val
    return bounds


def classify_reactions(reactions, v_min, v_max):
    """Classify reactions as protected or candidate"""
    # Protected prefixes/patterns
    protected_prefixes = [
        'Exchange_',
        'ExB_',
        'Exe',
        'Trans',       # Trans1[K], Trans2[K], etc.
        'PhloemTransport_',
        'PhloemImport_',
    ]
    
    # Explicitly protected reactions
    protected_exact = {
        'Seed_Biomass[K]',
        'R00451[K,p]',       # DAP decarboxylase (lysine target)
    }
    
    protected = []
    candidates = []
    inactive = []  # Reactions with zero bounds (v_min=0, v_max=0)
    
    for rxn in reactions:
        # Check explicit protection
        if rxn in protected_exact:
            protected.append(rxn)
            continue
        
        # Check prefix-based protection
        is_protected = False
        for prefix in protected_prefixes:
            if rxn.startswith(prefix):
                is_protected = True
                break
        
        # Also protect Trans reactions (Trans6[K], Trans13[K], etc.)
        if re.match(r'^Trans\d+\[', rxn):
            is_protected = True
        
        if is_protected:
            protected.append(rxn)
        else:
            # Check if reaction is inactive (zero bounds)
            vmin = v_min.get(rxn, 0)
            vmax = v_max.get(rxn, 0)
            if abs(vmin) < 1e-10 and abs(vmax) < 1e-10:
                inactive.append(rxn)
                # Inactive reactions are technically protected (no point knocking out)
                protected.append(rxn)
            else:
                candidates.append(rxn)
    
    return protected, candidates, inactive


def write_gams_set(filename, reactions, set_name=""):
    """Write a GAMS set include file"""
    with open(filename, 'w') as f:
        f.write('/\n')
        for rxn in reactions:
            f.write(f"'{rxn}'\n")
        f.write('/\n')


def main():
    if len(sys.argv) < 2:
        print("Usage: python preprocess_optknock.py <data_dir>")
        print("  e.g., python preprocess_optknock.py ../Wild_Type/B18")
        sys.exit(1)
    
    data_dir = sys.argv[1]
    
    print(f"Reading model data from: {data_dir}")
    reactions = read_reactions(data_dir)
    v_min = read_bounds(data_dir, 'v_min.txt')
    v_max = read_bounds(data_dir, 'v_max.txt')
    
    print(f"  Total reactions: {len(reactions)}")
    
    protected, candidates, inactive = classify_reactions(reactions, v_min, v_max)
    
    print(f"  Protected reactions: {len(protected)}")
    print(f"    (of which inactive/zero-bounds: {len(inactive)})")
    print(f"  Candidate reactions: {len(candidates)}")
    
    # Write GAMS set files
    write_gams_set('protected_reactions.txt', protected)
    write_gams_set('candidate_reactions.txt', candidates)
    
    print(f"\nGenerated:")
    print(f"  protected_reactions.txt  ({len(protected)} reactions)")
    print(f"  candidate_reactions.txt  ({len(candidates)} reactions)")
    
    # Print summary by category
    cat_counts = {}
    for rxn in protected:
        if rxn in {'Seed_Biomass[K]', 'R00451[K,p]'}:
            cat = 'Explicit'
        elif rxn.startswith('Exchange_'):
            cat = 'Exchange'
        elif rxn.startswith('ExB_'):
            cat = 'ExB (boundary export)'
        elif rxn.startswith('Exe'):
            cat = 'Exe (extracellular)'
        elif rxn.startswith('Trans'):
            cat = 'Trans (tissue transport)'
        elif rxn.startswith('Phloem'):
            cat = 'Phloem'
        elif rxn in inactive:
            cat = 'Inactive (zero bounds)'
        else:
            cat = 'Other'
        cat_counts[cat] = cat_counts.get(cat, 0) + 1
    
    print("\nProtected reaction categories:")
    for cat, count in sorted(cat_counts.items()):
        print(f"  {cat}: {count}")
    
    cand_counts = {}
    for rxn in candidates:
        if rxn.startswith('MR'):
            cat = 'MR (metabolic)'
        elif rxn.startswith('R'):
            cat = 'R (metabolic)'
        elif 'Transport' in rxn:
            cat = 'Intracellular transport'
        else:
            cat = 'Other'
        cand_counts[cat] = cand_counts.get(cat, 0) + 1
    
    print("\nCandidate reaction categories:")
    for cat, count in sorted(cand_counts.items()):
        print(f"  {cat}: {count}")


if __name__ == '__main__':
    main()
