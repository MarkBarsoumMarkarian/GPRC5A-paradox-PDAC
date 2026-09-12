# GPRC5A PDAC version 2 workflow

Run from the repository root:

```bash
Rscript R/v2/run_v2.R
```

The script downloads or reuses public inputs, verifies cohort eligibility,
performs native-platform survival analyses, runs a random-effects meta-analysis,
tests the GPRC5A-by-subtype interaction in TCGA, evaluates paired tumour and
cellular-compartment expression, and calculates matched CPTAC RNA-protein
correlation.

The only non-base recommended package required is `survival`.

Large inputs are cached under `data/v2/cache/` and are ignored by Git. Small
auditable tables and publication figures are written under `results/v2/`.

See `docs/V2_PROTOCOL.md` and `CORRECTION_NOTICE.md` before interpreting any
result.
