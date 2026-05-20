# OptKnock: Bilevel Optimization for Metabolite Overproduction

Implementation of the OptKnock algorithm (Burgard et al., 2003) in GAMS.
Identifies gene knockout strategies that couple growth to overproduction
of a target metabolite.

## Reference

Burgard AP, Pharkya P, Maranas CD (2003). "OptKnock: A bilevel programming
framework for identifying gene knockout strategies for microbial strain
optimization." *Biotechnol Bioeng* 84(6):647-657.

---

## Quick Start

### 1. Prepare your input files

Place these **5 files** in a single directory (e.g., `./my_model/`):

| File | Format | Description |
|------|--------|-------------|
| `reactions.txt` | GAMS set | List of reaction identifiers |
| `metabolites.txt` | GAMS set | List of metabolite identifiers |
| `sij.txt` | GAMS parameter | Stoichiometric matrix S(i,j) |
| `upper_bound.txt` | GAMS parameter | Upper flux bounds per reaction |
| `lower_bound.txt` | GAMS parameter | Lower flux bounds per reaction |

#### File format examples

**reactions.txt** (GAMS set — begins and ends with `/`):
```
/
'R00001'
'R00002'
'Biomass'
'EX_glucose_e'
/
```

**metabolites.txt** (same format):
```
/
'glucose_c'
'pyruvate_c'
/
```

**sij.txt** (GAMS parameter — `'metabolite'.'reaction'  coefficient`):
```
/
'glucose_c'.'R00001'  -1
'pyruvate_c'.'R00001'   2
/
```

**upper_bound.txt** / **lower_bound.txt** (GAMS parameter — `'reaction'  value`):
```
/
'R00001'	1000
'R00002'	0
'Biomass'	1000
/
```

### 2. Run OptKnock

#### Option A: Using the shell wrapper (recommended)

```bash
bash run_optknock_generic.sh \
    --datadir  ./my_model \
    --biomass  "Biomass" \
    --target   "EX_lysine_e" \
    --label    WT_run1 \
    --max-ko   3 \
    --min-bio  0.10
```

#### Option B: Directly with GAMS

```bash
# Step 1: Generate protected_reactions.txt (optional but recommended)
python3 preprocess_optknock_generic.py ./my_model \
    --biomass "Biomass" --target "EX_lysine_e"

# Step 2: Create the config file with reaction names
#   (GAMS CLI cannot handle brackets/commas in values,
#    so reaction names go in this include file)
cat > optknock_config.inc <<'EOF'
$setglobal BIOMASS Biomass
$setglobal TARGET  EX_lysine_e
EOF

# Step 3: Run GAMS
gams OptKnock.gms \
    --DATADIR="./my_model" \
    --MAX_KO=3 \
    --MIN_BIO=0.10 \
    --DAPNAME="WT" \
    lo=3
```

### 3. View results

- `results_optknock.txt` — Human-readable report with knockout details
- `results_optknock_summary.csv` — Machine-readable summary per K value
- `results_optknock_knockouts.csv` — Knocked-out reactions per K value

---

## Command-Line Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--DATADIR` | `.` | Path to folder with the 5 input files |
| `--MAX_KO` | `5` | Maximum simultaneous knockouts (K=1..MAX_KO) |
| `--MIN_BIO` | `0.01` | Minimum biomass as fraction of wild-type |
| `--BIGM` | `1000` | Big-M constant for McCormick linearisation |
| `--TIMELIM` | `3600` | Solver time limit in seconds |
| `--DAPNAME` | `run` | Label for this run (appears in CSV output) |

**Reaction names** (`BIOMASS`, `TARGET`) and **bound filenames** (`UBFILE`,
`LBFILE`) are set in `optknock_config.inc` rather than on the command line,
because GAMS CLI cannot handle brackets or commas.

```
$setglobal BIOMASS  Seed_Biomass[K]
$setglobal TARGET   R00451[K,p]
$setglobal UBFILE   upper_bound.txt
$setglobal LBFILE   lower_bound.txt
```

---

## Protected Reactions

Reactions that should **not** be knocked out (exchange reactions, transport
boundaries, biomass, target) are listed in `protected_reactions.txt`.

- **If you provide this file** in the working directory, OptKnock uses it.
- **If you omit it**, only the biomass and target reactions are protected
  (every other reaction becomes a knockout candidate).

The preprocessing script `preprocess_optknock_generic.py` can auto-generate
this file. It protects reactions based on name prefixes (configurable via
`--protect-prefixes`). Default prefixes: `Exchange_, ExB_, Exe, Trans,
PhloemTransport_, PhloemImport_`.

---

## Output Interpretation

For each K value (1 to MAX_KO), OptKnock reports three verification checks:

| Check | What it means |
|-------|--------------|
| **A. Max biomass FBA** | Standard FBA with knockouts applied → confirms cell can grow |
| **B. Pessimistic** | Minimise target at maximum biomass → the *guaranteed minimum* target production |
| **C. Optimistic** | Maximise target at maximum biomass → the *best possible* target production |

**True growth coupling** occurs only when the **pessimistic** target flux
exceeds the wild-type target flux (i.e., the cell *must* produce more target
to achieve maximum growth).

---

## Files

| File | Description |
|------|-------------|
| `OptKnock.gms` | Main GAMS model (generic, no hardcoded reaction names) |
| `preprocess_optknock_generic.py` | Classifies reactions as protected vs candidate |
| `run_optknock_generic.sh` | Shell wrapper for preprocessing + GAMS execution |
| `cplex.opt` | CPLEX solver options (threads, time limit, gap tolerance) |
| `analyze_optknock_results.py` | Post-processing: figures, LaTeX tables, rankings |

---

## Requirements

- **GAMS** (tested with 24.7+) with **CPLEX** solver
- **Python 3** (for preprocessing and analysis scripts)
- **matplotlib** (optional, for figure generation in the analysis script)

---

## Maize Endosperm Example

For the maize endosperm model in this project:

```bash
bash run_optknock_generic.sh \
    --datadir  ../Wild_Type/B18 \
    --biomass  "Seed_Biomass[K]" \
    --target   "R00451[K,p]" \
    --label    B18 \
    --max-ko   5 \
    --min-bio  0.01 \
    --ub-file  v_max.txt \
    --lb-file  v_min.txt
```

Note: the original model uses `v_max.txt` / `v_min.txt` instead of
`upper_bound.txt` / `lower_bound.txt`. Use `--ub-file` and `--lb-file`
to override the default filenames.
