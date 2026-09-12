# Scientific correction notice

The original repository narrative described a subtype-dependent reversal of the
GPRC5A-survival association in pancreatic ductal adenocarcinoma (PDAC) and
called it Simpson's paradox. That interpretation is not supported by the saved
continuous Cox models.

- Overall: HR 1.362, 95% CI 1.173-1.582, p = 5.26e-05.
- Classical: HR 1.530, 95% CI 1.174-1.995, p = 0.00165.
- Basal-like: HR 1.261, 95% CI 1.063-1.496, p = 0.00772.
- GPRC5A-by-subtype interaction: HR 1.192, 95% CI 0.871-1.632,
  p = 0.272.

Both subtype-specific estimates point toward higher hazard with higher GPRC5A,
and the interaction test does not support different effects by subtype. A
separate check of the repository's high-versus-low split likewise shows worse,
not better, survival for the GPRC5A-high Basal-like group.

The previous "role-state" classifier is also not evidence that GPRC5A switches
between oncogenic and suppressive states. Its target label was constructed from
subtype, GPRC5A expression and vital status, while closely related variables
were then used as predictors. This makes the task circular and clinically
non-independent. The held-out test contained only approximately 12 patients.

The gemcitabine analysis is retained only as a historical exploratory analysis.
There was one treatment-naive patient, treatment timing was not modelled, and
baseline Cox adjustment for a post-diagnosis treatment indicator cannot support
a causal deconfounding claim.

The project is being rebuilt as a multi-cohort assessment of GPRC5A prognostic
heterogeneity and cellular compartment expression. The scripts in `R/aim1_*`
through `R/aim5.R` are preserved for auditability but are superseded by the
versioned workflow under `R/v2/`.

## What remains potentially useful

- GPRC5A is a biologically plausible and experimentally supported PDAC target.
- TCGA contains an adverse expression-survival association worth testing for
  transportability.
- Independent survival cohorts can determine whether the association is
  reproducible or cohort-dependent.
- Laser-capture, single-cell or spatial data can determine whether bulk GPRC5A
  reflects malignant epithelium, stroma or sample composition.
- Matched RNA-protein data can assess whether transcript abundance represents
  protein abundance, but correlation alone is not prognostic validation.

No result in this repository should currently be interpreted as a clinically
validated prognostic biomarker or as evidence of a GPRC5A functional switch.
