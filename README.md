# Decoding the GPRC5A Paradox in Pancreatic Ductal Adenocarcinoma

**A subtype-stratified, treatment-deconfounded, multi-omic investigation**

> *Mark Barsoum Markarian — Department of Anatomy, Cell Biology and Physiological Sciences, Faculty of Medicine, American University of Beirut*  
> *Correspondence: mkm25@aub.edu.lb*

---

## Background

A prior batch-harmonized machine learning framework ([Markarian 2025, bioRxiv](https://doi.org/10.1101/2025.11.14.688421)) identified GPRC5A as prognostically relevant in PDAC but found **reduced expression in deceased patients** — the opposite of its established oncogenic role. This repository contains the full analytical pipeline to resolve that paradox through five mechanistic aims.

## Study Design

| Aim | Question | Data | Method |
|-----|----------|------|--------|
| 1 | Does molecular subtype mixing explain the paradox? | TCGA-PAAD (n=177) | Moffitt 2015 classifier, subtype-stratified Cox |
| 2 | Does gemcitabine confound the survival signal? | TCGA-PAAD (n=177) | Treatment-stratified Cox, multivariable adjustment |
| 3 | Is post-transcriptional regulation a driver? | CPTAC-PAAD (n=140) | Spearman RNA–protein correlation, genome-wide benchmarking |
| 4 | Are somatic mutations structurally organized? | TCGA-PAAD + AlphaFold2 | GDC mutation extraction, pLDDT domain mapping |
| 5 | Can ML predict GPRC5A functional role state? | TCGA-PAAD | Leakage-free RF / XGBoost / Logistic LOOCV pipeline |

## Key Results

- **Aim 1:** GPRC5A is oncogenic in classical subtype (HR=1.53, p=0.0017) with opposing KM directionality in basal-like (log-rank p=0.022) — subtype mixing explains the bulk paradox (Simpson's paradox)
- **Aim 2:** Gemcitabine attenuates the HR to non-significance in treated patients; fully adjusted model retains significance (HR=1.44, p=3.89×10⁻⁶) — gemcitabine is a secondary, not primary, confound
- **Aim 3:** RNA–protein Spearman r=0.571 (84.6th genome-wide percentile) — post-transcriptional repression is not a major driver
- **Aim 4:** Zero somatic GPRC5A mutations across 177 TCGA-PAAD samples — dysregulation is regulatory, not structural
- **Aim 5:** Random Forest role-state classifier — held-out test AUC=0.833, LOOCV AUC=0.758; classical co-expression features dominate over GPRC5A expression itself (GPRC5A ranks 23rd in importance)

## Repository Structure

```
gprc5a-paradox-pdac/
├── R/
│   ├── utils_clinical.R                  # Shared utilities (TCGAbiolinks version-agnostic)
│   ├── aim1_subtype_stratification.R     # Moffitt classifier, subtype-stratified Cox
│   ├── aim2.R                            # Treatment deconfounding, multivariable Cox
│   ├── aim3.R                            # CPTAC RNA-protein correlation
│   ├── aim4.R                            # AlphaFold2 structural mapping, somatic mutations
│   └── aim5.R                            # Leakage-free ML role-state classifier
├── results/
│   ├── figures/                          # PDF figures (all aims)
│   └── tables/                           # CSV result tables (all aims)
├── data/                                 # Not tracked — see Data Access below
└── docs/
```

## Reproducing the Analysis

### Requirements

```r
# Core
library(TCGAbiolinks)
library(DESeq2)
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(here)

# ML (Aim 5)
library(caret)
library(randomForest)
library(xgboost)
library(glmnet)
library(pROC)
```

### Execution Order

```r
source("R/utils_clinical.R")   # Always source first
source("R/aim1_subtype_stratification.R")
source("R/aim2.R")
source("R/aim3.R")
source("R/aim4.R")
source("R/aim5.R")
```

Each script is self-contained and writes outputs to `results/figures/` and `results/tables/`. Run from project root (`here::here()` handles paths).

### Data Access

Raw data are not tracked due to data use agreements.

| Dataset | Access |
|---------|--------|
| TCGA-PAAD RNA-seq, clinical, mutations | [GDC Portal](https://portal.gdc.cancer.gov) via `TCGAbiolinks` (auto-downloaded on first run) |
| CPTAC-PAAD proteomics | [CPTAC Data Portal](https://cptac-data-portal.georgetown.edu) |
| AlphaFold2 GPRC5A structure (Q8NFJ5) | [AlphaFold DB](https://alphafold.ebi.ac.uk/entry/Q8NFJ5) |

## Citation

If you use this code, please cite the companion preprint:

> Markarian MB. Batch-harmonized machine learning framework for cross-cohort RNA biomarker discovery in pancreatic adenocarcinoma. *bioRxiv*. 2025. https://doi.org/10.1101/2025.11.14.688421

The GPRC5A paradox manuscript is in preparation.

## License

MIT

## Keywords

`GPRC5A` · `PDAC` · `pancreatic cancer` · `molecular subtypes` · `Moffitt classifier` · `gemcitabine` · `CPTAC` · `AlphaFold2` · `machine learning` · `Simpson's paradox` · `oncogenic switching` · `R` · `bioinformatics`
