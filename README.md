# Endosperm O2-Lysine Metabolic Model

**Title:** Metabolic Reprogramming in the *opaque2* Maize Endosperm: A Genome-Scale Model Reveals Coordinated Flux Rewiring for Lysine Overaccumulation

**Authors:** Mary Nampimah, Rajdeep Khangura, Luis Avila-Ospina, Bertrand Hirel, Alexander Johnson, Sue Rhee, Surinder Chopra, Siela Maximova, Cathie Martin, Oliver Fiehn

---

## Repository Structure

```
├── manuscript/          # LaTeX manuscript and compiled PDF
├── supplementary/       # Supplementary tables (LaTeX + PDF) and CSV data
│   └── data/            # Tabular results from all analyses
├── models/              # GAMS genome-scale metabolic models
│   ├── Wild_Type/       # WT endosperm (8 developmental stages: DAP 6–30)
│   └── O2_mutant/       # o2 mutant endosperm (8 developmental stages)
└── code/                # Python/shell analysis scripts
    ├── gene_target_analysis/
    └── optknock/
```

---

## Models

Genome-scale metabolic model (GSMM) of maize (*Zea mays*) endosperm with:
- **3,024 reactions**, ~3,460 metabolites, **6 compartments**
- **16 conditions**: 8 developmental stages (DAP 6, 8, 10, 12, 15, 18, 22, 30) × 2 genotypes (WT, *o2*)
- Solver: **CPLEX** via GAMS (`.gms` files)

Each model folder (e.g., `models/Wild_Type/B6/`) contains:
| File | Description |
|------|-------------|
| `Endosperm_Model_FBA.gms` | FBA model (primal) |
| `FCA_FVA.gms` | Flux Coupling Analysis + FVA |
| `Scenario_FBA.gms` | Scenario testing |
| `reactions.txt` | Reaction list with GPR associations |
| `metabolites.txt` | Metabolite list |
| `v_max.txt` / `v_min.txt` | FVA flux ranges |
| `results_FBA.txt` | FBA optimal flux solution |

---

## Code

| Script | Description |
|--------|-------------|
| `analyze_o2_vs_wt.py` | Main comparative flux analysis (WT vs *o2*) |
| `analyze_fva.py` | Flux Variability Analysis processing |
| `flux_sum_analysis.py` | Metabolite-level flux sum calculations |
| `tradeoff_analysis.py` | Trade-off and rewiring index computation |
| `compute_manuscript_metrics.py` | Rewiring Index, Carbon Allocation, Flux Fraction |
| `generate_figures.py` | Publication-quality figure generation |
| `run_all_models.sh` | Shell script to run all GAMS models |
| `gene_target_analysis/identify_lysine_targets.py` | FSEOF + KO/OE screening |
| `gene_target_analysis/map_reactions_to_genes.py` | GPR gene mapping |
| `optknock/OptKnock.gms` | OptKnock formulation (generic) |
| `optknock/OptKnock_Endosperm.gms` | OptKnock for endosperm model |

### Requirements
```
conda create -n endosperm python=3.10
conda activate endosperm
pip install cobra pandas numpy matplotlib seaborn openpyxl scipy
```
> Note: GAMS models require a licensed GAMS installation with CPLEX solver.

---

## Supplementary Data

All supplementary tables (S1–S22) are compiled in `supplementary/supplementary_tables.tex` / `.pdf`.

Raw CSV data in `supplementary/data/`:
| Folder | Contents |
|--------|----------|
| `analysis_results/` | Biomass flux, lysine pathway, differential flux tables |
| `flux_sum_analysis/` | Metabolite flux sums across conditions |
| `gene_target_analysis/` | KO/OE screening results, gene mappings |
| `flux_coupling_analysis/` | Reaction coupling coefficients |

---

## Citation

> Nampimah M, *et al.* (2026). Metabolic Reprogramming in the *opaque2* Maize Endosperm: A Genome-Scale Model Reveals Coordinated Flux Rewiring for Lysine Overaccumulation. *Manuscript in preparation.*

---

## License

This repository is made available for academic use. Model files and analysis code are provided under the MIT License.
