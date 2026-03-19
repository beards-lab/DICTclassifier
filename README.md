# Predicting Metabolism-Mediated Drug-Induced Cardiotoxicity Using Genome-Scale Metabolic Flux Analysis

This repository contains the MATLAB code for reproducing all computational results and figures in:

> Schenk, N.A., Sturgess, A., Mahoney, C., & Beard, D.A. "Predicting Metabolism-Mediated Mechanisms of Drug-Induced Cardiotoxicity Using Genome-Scale Metabolic Flux Analysis." *University of Michigan, Department of Molecular & Integrative Physiology.*

## Overview

We developed a computational pipeline that integrates genome-scale metabolic modeling with machine learning to predict drug-induced cardiotoxicity (DICT). The approach:

1. Builds a thermodynamically constrained cardiomyocyte metabolic model (iCardio) using the RECON3D reconstruction
2. Integrates RNA-seq data from drug-treated iPSC-derived cardiomyocytes (6 cell lines, 54 drugs) to compute metabolic flux perturbations
3. Trains ensemble machine learning classifiers on reaction-level flux features, validated against FDA Adverse Event Reporting System (FAERS) cardiotoxicity labels
4. Identifies metabolic subsystems driving toxicity predictions and evaluates drug combination (anthracycline + olmesartan) effects on cardiotoxicity risk

## Requirements

- **MATLAB R2024b** (or later)
- **Required Toolboxes:**
  - [Optimization Toolbox](https://www.mathworks.com/products/optimization.html) — `solve`, `linprog`, `optimvar`, `optimproblem`
  - [Statistics and Machine Learning Toolbox](https://www.mathworks.com/products/statistics.html) — `fitcensemble`, `rocmetrics`, `tsne`, `cvpartition`, `bayesopt`, `predict`, `swarmchart`, `boxplot`
  - [Bioinformatics Toolbox](https://www.mathworks.com/products/bioinfo.html) — `mattest`, `mafdr`
  - [Parallel Computing Toolbox](https://www.mathworks.com/products/parallel-computing.html) — `parfor`, `parpool` (optional, for accelerated computation)
- **GPU** (optional) — Driver 10 uses `gpuArray` for permutation testing; falls back to CPU if unavailable

## Input Data

Place the following in the project directory before running:

| Directory/File | Description |
|---|---|
| `data/drugMap.xlsx` | Drug name ↔ 3-letter code mapping table |
| `data/GEGcsv/` | RNA-seq differential expression files (`*_mapped.csv`) for 54 drugs × 6 cell lines (266 files) |
| `data/FAERS/` | FDA Adverse Event Reporting System data for ROR calculation |

> **Note:** RNA-seq data paths are configured at the top of each script. Update `pathRNA` to point to your local data directory.

## Pipeline

Run scripts sequentially in the order below. Each script is self-contained — all helper functions are defined locally within the file.

| Step | Script | Description | Output | Figures |
|---|---|---|---|---|
| 1 | `step01_calculate_ROR.m` | Compute Reporting Odds Ratios from FAERS data | `out/ROR_results.mat` | — |
| 2 | `step02_optimize_GEM.m` | Thermodynamic optimization of iCardio GEM | `out/iCardio_optimized.mat` | Fig 1B |
| 3 | `step03_calculate_delta_flux.m` | Delta-flux for 54 drugs vs DMSO (all cell lines) | `out/deltaFlux_stats.mat` | Fig 2, 3 (data) |
| 4 | `step04_train_ML_models.m` | Train ensemble ML classifiers (100×4-fold CV) | `out/DICT_L0X/c1_10_6_30hr/` | — |
| 5 | `step05_evaluate_model_metrics.m` | Calculate comprehensive metrics (22 metrics, ROC, PRC) | `out/trainedModels_comprehensive.mat` | — |
| 6 | `step06_generate_publication_figures.m` | Generate Figures 4–6 (ROR, t-SNE, correlograms, model performance) | `figures/` | Fig 4, 5, 6 |
| 7 | `step07_subsystem_importance.m` | Predictor importance & subsystem sensitivity sweep | `out/driver7_results.mat` | Fig 7A |
| 8 | `step08_subsystem_sweep_mass_balanced.m` | Mass-balanced subsystem sweep (null-space projection) | `out/driver8_results.mat` | Fig 7B |
| 9 | `step09_drug_combination_flux.m` | Delta-flux for anthracycline + OLM combinations | `out/RxnStats_*.mat` | — |
| 10 | `step10_combination_ML_inference.m` | ML inference + beeswarm/boxplot for combo scoring | `out/driver10_*.mat` | Fig 8A |
| 11 | `step11_dose_response_grid_flux.m` | Dose-response grid flux (21×21 concentration grid) | `outPf_GL/` | — |
| 12 | `step12_dose_response_heatmaps.m` | Grid search heatmaps with ML risk scoring | `figures/` | Fig 8B |

### Quick Start

```matlab
% Step 1: Configure paths
% Open each script and set pathRNA to your data directory

% Step 2: Run the pipeline sequentially
step01_calculate_ROR
step02_optimize_GEM
step03_calculate_delta_flux
step04_train_ML_models
step05_evaluate_model_metrics
step06_generate_publication_figures
step07_subsystem_importance
step08_subsystem_sweep_mass_balanced
step09_drug_combination_flux
step10_combination_ML_inference
step11_dose_response_grid_flux
step12_dose_response_heatmaps
```

> **Runtime note:** Steps 3, 4, 9, and 11 are computationally intensive and benefit significantly from the Parallel Computing Toolbox. On a 34-core workstation, the full pipeline takes approximately 24–48 hours. Steps 6, 7, 8, 10, and 12 complete in minutes.

## Suggested File Renaming

For clarity, we recommend renaming the original driver scripts as follows:

| Original Name | Suggested Name |
|---|---|
| `driver1_calc_ROR-11.m` | `step01_calculate_ROR.m` |
| `driver2_optimize_iCardio-5.m` | `step02_optimize_GEM.m` |
| `driver3_deltaFlux_SJ0-6.m` | `step03_calculate_delta_flux.m` |
| `driver4_train_ML_models-2.m` | `step04_train_ML_models.m` |
| `driver5_calc_comprehensive_metrics-13.m` | `step05_evaluate_model_metrics.m` |
| `driver6_make_publication_figs-3.m` | `step06_generate_publication_figures.m` |
| `driver7_parseModels_final-7.m` | `step07_subsystem_importance.m` |
| `driver8_parseSubSystems_MassBalanced-14.m` | `step08_subsystem_sweep_mass_balanced.m` |
| `driver9_drugConvComp_ver2-12.m` | `step09_drug_combination_flux.m` |
| `driver10_v3-4.m` | `step10_combination_ML_inference.m` |
| `driver11_drugCombComp_NoahsFinalAttempt-8.m` | `step11_dose_response_grid_flux.m` |
| `driver12_vP-9.m` | `step12_dose_response_heatmaps.m` |

## Project Structure

```
├── README.md
├── step01_calculate_ROR.m
├── step02_optimize_GEM.m
├── step03_calculate_delta_flux.m
├── step04_train_ML_models.m
├── step05_evaluate_model_metrics.m
├── step06_generate_publication_figures.m
├── step07_subsystem_importance.m
├── step08_subsystem_sweep_mass_balanced.m
├── step09_drug_combination_flux.m
├── step10_combination_ML_inference.m
├── step11_dose_response_grid_flux.m
├── step12_dose_response_heatmaps.m
├── data/
│   ├── drugMap.xlsx
│   ├── FAERS/
│   └── GEGcsv/
├── out/                    % Generated intermediate results
└── figures/                % Generated figures
```

## Citation

If you use this code, please cite:

```
Schenk, N.A., Sturgess, A., Mahoney, C., & Beard, D.A. (2026). Predicting Metabolism-Mediated 
Mechanisms of Drug-Induced Cardiotoxicity Using Genome-Scale Metabolic Flux Analysis. 
University of Michigan, Department of Molecular & Integrative Physiology.
```

## License

[Choose your license — MIT, GPL, BSD, etc.]

## Contact

Noah A. Schenk — naschenk@umich.edu  
University of Michigan, Department of Molecular & Integrative Physiology
