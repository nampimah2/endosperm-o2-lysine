#!/usr/bin/env python3
"""
Map gene target reactions to their specific genes using GPR rules
from endosperm_GSM.xlsx.

Reads:
  - endosperm_GSM.xlsx (Reactions, GPR, Pathway columns)
  - gene_target_analysis/07_master_target_ranking.csv
  - gene_target_analysis/03_overexpression_screen.csv
  - gene_target_analysis/04_knockout_screen.csv
  - gene_target_analysis/05_o2_bound_mimicry.csv
  - gene_target_analysis/08_target_verification.csv

Outputs:
  - gene_target_analysis/gene_mapping_all_targets.csv
  - gene_target_analysis/gene_mapping_summary.csv
  - gene_target_analysis/gene_mapping_report.txt
"""

import pandas as pd
import re
import os

BASEDIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTDIR  = os.path.join(BASEDIR, 'gene_target_analysis')

# ──────────────────────────────────────────────────────────────
# 1. Load GPR data
# ──────────────────────────────────────────────────────────────
print("Loading endosperm_GSM.xlsx ...")
gsm = pd.read_excel(os.path.join(BASEDIR, 'endosperm_GSM.xlsx'), sheet_name='Sheet1')
gsm.columns = gsm.columns.str.strip()
print(f"  Total reactions in GSM: {len(gsm)}")
print(f"  Reactions with GPR: {gsm['GPR'].notna().sum()}")

# Build reaction -> GPR dictionary
rxn_gpr = {}
rxn_pathway = {}
for _, row in gsm.iterrows():
    rxn = str(row['Reactions']).strip()
    gpr = row['GPR']
    pathway = row['Pathway']
    rxn_gpr[rxn] = str(gpr).strip() if pd.notna(gpr) else None
    rxn_pathway[rxn] = str(pathway).strip() if pd.notna(pathway) else None

# ──────────────────────────────────────────────────────────────
# 2. Helper: parse genes from a GPR rule string
# ──────────────────────────────────────────────────────────────
GENE_PATTERN = re.compile(r'((?:GRMZM\d+G\d+|AC\d+\.\d+_FG\w+|AC\d+\.\d+_FGP\d+))')

def parse_genes(gpr_string):
    """Extract all unique gene IDs from a GPR rule string."""
    if gpr_string is None or gpr_string == 'nan' or gpr_string == '':
        return []
    genes = GENE_PATTERN.findall(gpr_string)
    return sorted(set(genes))

def classify_gpr(gpr_string):
    """Classify GPR rule as 'single_gene', 'isozymes (or)', 'complex (and)', 'mixed', or 'none'."""
    if gpr_string is None or gpr_string == 'nan' or gpr_string == '':
        return 'none'
    genes = parse_genes(gpr_string)
    if len(genes) == 0:
        return 'none'
    if len(genes) == 1:
        return 'single_gene'
    has_and = ' and ' in gpr_string
    has_or  = ' or '  in gpr_string
    if has_and and has_or:
        return 'complex_and_or'
    elif has_and:
        return 'enzyme_complex'
    elif has_or:
        return 'isozymes'
    return 'single_gene'

# ──────────────────────────────────────────────────────────────
# 3. Fuzzy reaction matching
# ──────────────────────────────────────────────────────────────
def find_reaction_in_gsm(rxn_name):
    """
    Try to find the reaction in the GSM dictionary.
    The reaction names in the target analysis may differ slightly
    from those in the GSM (e.g., compartment annotations).
    """
    rxn_name = rxn_name.strip()
    # Direct match
    if rxn_name in rxn_gpr:
        return rxn_name
    # Try without trailing whitespace variants
    for key in rxn_gpr:
        if key.strip() == rxn_name:
            return key
    # Try matching the base reaction ID (before [)
    base = rxn_name.split('[')[0]
    matches = [k for k in rxn_gpr if k.startswith(base + '[')]
    if len(matches) == 1:
        return matches[0]
    elif len(matches) > 1:
        # Prefer exact compartment match
        for m in matches:
            if m == rxn_name:
                return m
        return matches[0]  # return first match
    # Try substring match
    matches = [k for k in rxn_gpr if base in k]
    if matches:
        return matches[0]
    return None

# ──────────────────────────────────────────────────────────────
# 4. Load robust CSV (handles commas in reaction names)
# ──────────────────────────────────────────────────────────────
def load_csv_robust(filepath):
    """Load CSV handling reaction names that contain commas."""
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
        if not lines:
            return pd.DataFrame()
        header = lines[0].strip().split(',')
        n_cols = len(header)
        rows = []
        for line in lines[1:]:
            line = line.strip()
            if not line:
                continue
            parts = line.split(',')
            if len(parts) > n_cols:
                excess = len(parts) - n_cols
                first_col = ','.join(parts[:excess+1])
                rest = parts[excess+1:]
                rows.append([first_col] + rest)
            else:
                rows.append(parts)
        df = pd.DataFrame(rows, columns=header)
        return df
    except Exception:
        return pd.read_csv(filepath)

# ──────────────────────────────────────────────────────────────
# 5. Collect all target reactions from analysis files
# ──────────────────────────────────────────────────────────────
print("\nLoading gene target analysis files ...")

all_targets = {}  # reaction_name -> {screen_type, fold_change, ...}

# Master ranking
try:
    master = load_csv_robust(os.path.join(OUTDIR, '07_master_target_ranking.csv'))
    print(f"  Master ranking: {len(master)} targets")
    for _, row in master.iterrows():
        rxn = str(row.iloc[0]).strip()
        if rxn not in all_targets:
            all_targets[rxn] = {
                'source': 'master_ranking',
                'rank': row.iloc[1] if master.shape[1] > 1 else '',
                'max_fc': row.iloc[3] if master.shape[1] > 3 else '',
                'screen_type': row.iloc[2] if master.shape[1] > 2 else '',
            }
except Exception as e:
    print(f"  Warning loading master ranking: {e}")

# Knockout screen
try:
    ko = load_csv_robust(os.path.join(OUTDIR, '04_knockout_screen.csv'))
    print(f"  Knockout screen: {len(ko)} targets")
    for _, row in ko.iterrows():
        rxn = str(row.iloc[0]).strip()
        if rxn not in all_targets:
            all_targets[rxn] = {'source': 'knockout', 'screen_type': 'KO'}
except Exception as e:
    print(f"  Warning loading knockout screen: {e}")

# Overexpression screen
try:
    oe = load_csv_robust(os.path.join(OUTDIR, '03_overexpression_screen.csv'))
    print(f"  Overexpression screen: {len(oe)} targets")
    for _, row in oe.iterrows():
        rxn = str(row.iloc[0]).strip()
        if rxn not in all_targets:
            all_targets[rxn] = {'source': 'overexpression', 'screen_type': 'OE'}
except Exception as e:
    print(f"  Warning loading OE screen: {e}")

# Mimicry screen
try:
    mim = load_csv_robust(os.path.join(OUTDIR, '05_o2_bound_mimicry.csv'))
    print(f"  Mimicry screen: {len(mim)} targets")
    for _, row in mim.iterrows():
        rxn = str(row.iloc[0]).strip()
        if rxn not in all_targets:
            all_targets[rxn] = {'source': 'mimicry', 'screen_type': 'MIM'}
except Exception as e:
    print(f"  Warning loading mimicry screen: {e}")

# Verification
try:
    verif = load_csv_robust(os.path.join(OUTDIR, '08_target_verification.csv'))
    print(f"  Verified targets: {len(verif)} targets")
except Exception as e:
    print(f"  Warning loading verification: {e}")
    verif = pd.DataFrame()

print(f"\nTotal unique target reactions: {len(all_targets)}")

# ──────────────────────────────────────────────────────────────
# 6. Key target reactions (from manuscript analysis)
# ──────────────────────────────────────────────────────────────
key_targets = [
    # Tier 1: High-impact knockouts
    ('Exe2[K]',               'KO',  'Lipid export',                      '46,434x'),
    ('Trans2[K]',             'KO',  'Lipid transport',                   '46,434x'),
    ('cpTransport_C00007[K]', 'KO',  'O2 transport (plastid-cytosol)',    '23,470x'),
    ('R00216[K,p]',           'KO',  'Pyruvate dehydrogenase (plastid)',  '16,234x'),
    ('R00282[K,c]',           'KO',  'Acetyl-CoA carboxylase (cytosol)',  '28x'),
    # Tier 2: Overexpression targets
    ('MR00558[K,p]',          'OE',  'Lipid biosynthesis',                '1.096x'),
    ('R01308[K,p]',           'OE',  'Shikimate pathway',                 '1.096x'),
    ('R01280[K,p]',           'OE',  'Chorismate synthesis',              '1.096x'),
    ('R00883[K,c]',           'OE',  'Phosphofructokinase',               '1.096x'),
    ('R00512[K,c]',           'OE',  'Pyruvate kinase (cytosol)',         '1.096x'),
    ('R01015[K,c]',           'OE',  'Glycolysis',                        '1.096x'),
    ('R00667[K,m]',           'OE',  'TCA cycle (mitochondrial)',         '1.096x'),
    # Tier 3: Mimicry targets
    ('MR00422[K,p]',          'MIM', 'Lipid metabolism (plastid)',         '1.096x'),
    ('MR00732[K,m]',          'MIM', 'Mitochondrial metabolism',           '1.096x'),
    ('R00342[K,p]',           'MIM', 'Transketolase (plastid)',            '1.096x'),
    # Verified targets
    ('cpTransport_C00712[K]', 'VER', 'Nicotinamide transport',             '1.096x'),
    # Lysine pathway (for reference)
    ('R00451[K,p]',           'REF', 'DAP decarboxylase (lysine)',         'N/A'),
    ('R02292[K,p]',           'REF', 'DHDPS (lysine committed step)',      'N/A'),
]

# ──────────────────────────────────────────────────────────────
# 7. Map reactions to genes
# ──────────────────────────────────────────────────────────────
print("\n" + "="*80)
print("MAPPING TARGET REACTIONS TO GENES")
print("="*80)

results = []
report_lines = []

report_lines.append("="*80)
report_lines.append("GENE TARGET MAPPING REPORT")
report_lines.append("Mapping metabolic model reactions to maize (Zea mays) gene IDs")
report_lines.append("Source: endosperm_GSM.xlsx GPR rules")
report_lines.append("="*80)
report_lines.append("")

# --- A. Key targets first ---
report_lines.append("="*80)
report_lines.append("SECTION A: KEY ENGINEERING TARGETS")
report_lines.append("="*80)

for rxn_name, screen_type, function, fold_change in key_targets:
    gsm_key = find_reaction_in_gsm(rxn_name)
    gpr = rxn_gpr.get(gsm_key, None) if gsm_key else None
    genes = parse_genes(gpr)
    gpr_type = classify_gpr(gpr)
    pathway = rxn_pathway.get(gsm_key, None) if gsm_key else None
    
    found = gsm_key is not None
    n_genes = len(genes)
    
    results.append({
        'Reaction': rxn_name,
        'GSM_Match': gsm_key if gsm_key else 'NOT_FOUND',
        'Screen_Type': screen_type,
        'Function': function,
        'Fold_Change': fold_change,
        'Has_GPR': 'Yes' if gpr and gpr != 'nan' else 'No',
        'GPR_Type': gpr_type,
        'N_Genes': n_genes,
        'Gene_IDs': '; '.join(genes) if genes else 'No GPR / Transport reaction',
        'Pathway': pathway if pathway else 'N/A',
        'GPR_Rule': gpr if gpr and gpr != 'nan' else 'N/A',
        'Category': 'Key_Target',
    })
    
    report_lines.append("")
    report_lines.append(f"─── {rxn_name} ({screen_type}) ───")
    report_lines.append(f"  Function:    {function}")
    report_lines.append(f"  Fold Change: {fold_change}")
    report_lines.append(f"  GSM Match:   {gsm_key if gsm_key else 'NOT FOUND in GSM'}")
    report_lines.append(f"  Pathway:     {pathway if pathway else 'N/A'}")
    if genes:
        report_lines.append(f"  GPR Type:    {gpr_type}")
        report_lines.append(f"  # Genes:     {n_genes}")
        report_lines.append(f"  Genes:")
        for g in genes:
            report_lines.append(f"    - {g}")
    else:
        report_lines.append(f"  GPR:         No gene association (transport/exchange reaction)")
    
    status = f"  {rxn_name:30s} -> {n_genes:3d} genes" if genes else f"  {rxn_name:30s} -> No GPR"
    print(status)

# --- B. All targets from analysis ---
report_lines.append("")
report_lines.append("="*80)
report_lines.append("SECTION B: ALL TARGETS FROM SCREENING ANALYSIS")
report_lines.append("="*80)

mapped_count = 0
no_gpr_count = 0
not_found_count = 0

for rxn_name, info in sorted(all_targets.items()):
    if any(rxn_name == kt[0] for kt in key_targets):
        continue  # already processed above
    
    gsm_key = find_reaction_in_gsm(rxn_name)
    gpr = rxn_gpr.get(gsm_key, None) if gsm_key else None
    genes = parse_genes(gpr)
    gpr_type = classify_gpr(gpr)
    pathway = rxn_pathway.get(gsm_key, None) if gsm_key else None
    
    if genes:
        mapped_count += 1
    elif gsm_key:
        no_gpr_count += 1
    else:
        not_found_count += 1
    
    results.append({
        'Reaction': rxn_name,
        'GSM_Match': gsm_key if gsm_key else 'NOT_FOUND',
        'Screen_Type': info.get('screen_type', 'unknown'),
        'Function': '',
        'Fold_Change': str(info.get('max_fc', '')),
        'Has_GPR': 'Yes' if gpr and gpr != 'nan' else 'No',
        'GPR_Type': gpr_type,
        'N_Genes': len(genes),
        'Gene_IDs': '; '.join(genes) if genes else 'No GPR',
        'Pathway': pathway if pathway else 'N/A',
        'GPR_Rule': gpr if gpr and gpr != 'nan' else 'N/A',
        'Category': 'Screening_Target',
    })

report_lines.append(f"\n  Targets with gene associations: {mapped_count}")
report_lines.append(f"  Targets without GPR (transport/exchange): {no_gpr_count}")
report_lines.append(f"  Targets not found in GSM: {not_found_count}")

# ──────────────────────────────────────────────────────────────
# 8. Save outputs
# ──────────────────────────────────────────────────────────────
# Full results CSV
df_results = pd.DataFrame(results)
outfile1 = os.path.join(OUTDIR, 'gene_mapping_all_targets.csv')
df_results.to_csv(outfile1, index=False)
print(f"\nSaved: {outfile1}")

# Summary CSV (key targets only, with genes expanded)
summary_rows = []
for rxn_name, screen_type, function, fold_change in key_targets:
    gsm_key = find_reaction_in_gsm(rxn_name)
    gpr = rxn_gpr.get(gsm_key, None) if gsm_key else None
    genes = parse_genes(gpr)
    pathway = rxn_pathway.get(gsm_key, None) if gsm_key else None
    gpr_type = classify_gpr(gpr)
    
    if genes:
        for gene in genes:
            summary_rows.append({
                'Reaction': rxn_name,
                'Screen_Type': screen_type,
                'Function': function,
                'Fold_Change': fold_change,
                'Gene_ID': gene,
                'GPR_Type': gpr_type,
                'Pathway': pathway if pathway else 'N/A',
            })
    else:
        summary_rows.append({
            'Reaction': rxn_name,
            'Screen_Type': screen_type,
            'Function': function,
            'Fold_Change': fold_change,
            'Gene_ID': 'NO_GPR',
            'GPR_Type': 'none',
            'Pathway': pathway if pathway else 'N/A',
        })

df_summary = pd.DataFrame(summary_rows)
outfile2 = os.path.join(OUTDIR, 'gene_mapping_summary.csv')
df_summary.to_csv(outfile2, index=False)
print(f"Saved: {outfile2}")

# Report
report_lines.append("")
report_lines.append("="*80)
report_lines.append("SUMMARY STATISTICS")
report_lines.append("="*80)

key_with_genes = sum(1 for r in results if r['Category'] == 'Key_Target' and r['N_Genes'] > 0)
key_total = sum(1 for r in results if r['Category'] == 'Key_Target')
all_genes = set()
for r in results:
    if r['Gene_IDs'] and r['Gene_IDs'] not in ('No GPR', 'No GPR / Transport reaction'):
        for g in r['Gene_IDs'].split('; '):
            all_genes.add(g.strip())

report_lines.append(f"  Key targets with gene associations: {key_with_genes}/{key_total}")
report_lines.append(f"  Total unique genes across all targets: {len(all_genes)}")
report_lines.append(f"  Total target reactions: {len(results)}")

outfile3 = os.path.join(OUTDIR, 'gene_mapping_report.txt')
with open(outfile3, 'w') as f:
    f.write('\n'.join(report_lines))
print(f"Saved: {outfile3}")

# ──────────────────────────────────────────────────────────────
# 9. Print key target summary to console
# ──────────────────────────────────────────────────────────────
print("\n" + "="*80)
print("KEY TARGET GENES SUMMARY")
print("="*80)
for rxn_name, screen_type, function, fold_change in key_targets:
    gsm_key = find_reaction_in_gsm(rxn_name)
    gpr = rxn_gpr.get(gsm_key, None) if gsm_key else None
    genes = parse_genes(gpr)
    n = len(genes)
    if genes:
        gene_str = ', '.join(genes[:5])
        if n > 5:
            gene_str += f' ... (+{n-5} more)'
        print(f"\n  {rxn_name} ({screen_type}, {fold_change})")
        print(f"    {function}")
        print(f"    {n} gene(s): {gene_str}")
    else:
        print(f"\n  {rxn_name} ({screen_type}, {fold_change})")
        print(f"    {function}")
        print(f"    ** No GPR association (transport/exchange reaction) **")

print(f"\nTotal unique genes across all key targets: ", end='')
key_genes = set()
for rxn_name, _, _, _ in key_targets:
    gsm_key = find_reaction_in_gsm(rxn_name)
    gpr = rxn_gpr.get(gsm_key, None) if gsm_key else None
    genes = parse_genes(gpr)
    key_genes.update(genes)
print(f"{len(key_genes)}")
print("\nDone!")
