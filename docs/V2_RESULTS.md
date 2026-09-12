# Version 2 manuscript-ready results and interpretation

## Proposed title

**GPRC5A is cancer-cell enriched but not a transportable stand-alone prognostic biomarker in pancreatic ductal adenocarcinoma: a five-cohort reanalysis**

## Proposed abstract

### Background

GPRC5A is highly expressed in pancreatic ductal adenocarcinoma (PDAC) and has been proposed as both a therapeutic target and prognostic biomarker. An earlier analysis from this repository incorrectly interpreted subtype-stratified results as a direction-reversing Simpson's paradox. We reassessed the distinction between target expression and prognostic transportability.

### Methods

We analysed overall survival in five independent primary-tumour cohorts spanning RNA sequencing and three microarray platforms (TCGA-PAAD, CPTAC-PAAD, GSE85916, GSE57495, and GSE62452). Expression was standardized within each cohort and was never batch-corrected across platforms. Cox proportional-hazards models used observed follow-up time and censoring. Cohort log hazard ratios were combined by random-effects meta-analysis with Hartung-Knapp inference, a prediction interval, and leave-one-cohort-out sensitivity analysis. We formally tested a GPRC5A-by-Moffitt-subtype interaction in TCGA. Target-expression evidence was evaluated using paired tumour-adjacent tissue, paired laser-captured cancer-stroma, and matched CPTAC RNA-protein data.

### Results

The survival analysis included 481 patients and 304 deaths after restricting TCGA to explicitly ductal primary tumours. Cohort hazard ratios per within-cohort standard deviation ranged from 0.92 in CPTAC-PAAD to 1.39 in GSE57495. The pooled association was inconclusive (HR 1.16, 95% CI 0.94-1.43, P=0.118), with moderate heterogeneity (I²=28.2%) and a 95% prediction interval crossing the null (0.80-1.69). Leave-one-cohort-out estimates were also inconclusive. GSE62452 showed non-proportional hazards (P=0.023). Both TCGA subtype-specific point estimates were adverse, and the adjusted expression-by-subtype interaction was not significant (P=0.146), refuting a direction-reversing subtype paradox. In contrast, GPRC5A expression was higher in tumour than adjacent tissue in 45 matched pairs (mean 2.20-fold, paired Wilcoxon P=6.03×10⁻⁸) and higher in laser-captured cancer cells than paired stroma in 19 untreated tumours (mean 5.67-fold, P=2.31×10⁻⁴). RNA and protein abundance were correlated in 135 matched CPTAC PDAC tumours (Spearman rho=0.567, P=7.40×10⁻¹³).

### Conclusions

GPRC5A shows a reproducible PDAC target-expression phenotype but lacks a stable, transportable bulk-tumour survival association. Target candidacy should therefore be separated from stand-alone prognostic use. Direct validation of surface localization, normal-tissue safety, and therapeutic response is required.

## Exact primary survival results

| Cohort | N | Deaths | HR per SD | 95% CI | P | C-index | PH P |
|---|---:|---:|---:|---:|---:|---:|---:|
| TCGA-PAAD | 145 | 84 | 1.270 | 0.990-1.629 | 0.060 | 0.583 | 0.278 |
| CPTAC-PAAD | 129 | 72 | 0.918 | 0.715-1.177 | 0.498 | 0.516 | 0.512 |
| GSE85916 | 79 | 57 | 1.132 | 0.851-1.507 | 0.394 | 0.524 | 0.522 |
| GSE57495 | 63 | 42 | 1.393 | 1.003-1.934 | 0.048 | 0.560 | 0.795 |
| GSE62452 | 65 | 49 | 1.242 | 0.911-1.695 | 0.171 | 0.516 | 0.023 |

Random-effects pooled HR 1.162 (95% CI 0.942-1.434), P=0.118; tau²=0.0081; I²=28.2%; prediction interval 0.800-1.689.

## Robustness and limitations

- Age/stage/grade-adjusted sensitivity models leave the key pattern unchanged: all cohort estimates remain weak or imprecise, and the TCGA estimate attenuates after adjustment. GSE85916 lacks suitable deposited clinical covariates and is unadjusted.
- Both annotated GSE85916 probes independently give similar null estimates, reducing concern that probe averaging created its result.
- No cross-platform ComBat step is used. Such correction would entangle cohort, platform, case mix, and outcome and could create rather than solve bias.
- The GSE62452 12-month piecewise analysis is a post-hoc proportional-hazards diagnostic, not a confirmatory subgroup finding.
- Cohorts differ in platform, sampling, stage distribution, treatment era, and follow-up. Those differences are part of the transportability question but prevent a causal explanation for heterogeneity.
- Adjacent tissue is not healthy pancreas. Laser-capture results localize the transcript more directly but include only 19 patients.
- RNA-protein correlation does not prove membrane accessibility, target safety, or response to inhibition.
- The work is retrospective, single-gene, and computational. Prospective or experimental validation remains necessary.

## Claims that must not reappear

Do not describe the results as Simpson's paradox, subtype switching, tumour-suppressive GPRC5A in basal-like PDAC, causal gemcitabine deconfounding, independent machine-learning prediction, or clinical/therapeutic validation.
