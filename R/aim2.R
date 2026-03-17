# =============================================================================
# GPRC5A Paradox in PDAC — Aim 2: Gemcitabine Deconfounding
# Project: Decoding the GPRC5A Paradox in Pancreatic Ductal Adenocarcinoma
# Author:  Mark Barsoum Markarian
# Date:    2026
#
# Description:
#   Tests whether gemcitabine treatment artificially elevates GPRC5A expression
#   in surviving patients, confounding the survival-expression association found
#   in the original batch-harmonized ML study (Markarian 2025).
#
#   Logic: Zhou et al. (2016) showed gemcitabine upregulates GPRC5A in PDAC.
#   If more surviving patients received gemcitabine (i.e., they lived long
#   enough to complete therapy), their GPRC5A would be artificially elevated —
#   making "Alive" patients appear to have higher GPRC5A, which inverts the
#   true oncogenic signal.
#
#   This script:
#     1. Extracts treatment history from TCGA-PAAD clinical annotations
#     2. Stratifies patients: treatment-naive vs. gemcitabine-treated
#     3. Re-tests GPRC5A vs. survival within each treatment stratum
#     4. Builds multivariable Cox models to isolate GPRC5A's independent effect
#     5. Tests whether treatment status mediates/moderates the GPRC5A signal
#
# Prerequisite: aim1_subtype_stratification.R must have been run first
#               (uses cached TCGA data and subtype assignments)
#
# Outputs (saved to results/):
#   - figures/aim2_treatment_boxplot.pdf
#   - figures/aim2_km_by_treatment_stratum.pdf
#   - figures/aim2_forest_multivariable.pdf
#   - figures/aim2_mediation_diagram.pdf
#   - tables/aim2_treatment_stratified_cox.csv
#   - tables/aim2_multivariable_cox.csv
#   - tables/aim2_treatment_summary.csv
# =============================================================================

# ── 0. Libraries ──────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(broom)
  library(ggpubr)
  library(patchwork)
  library(stringr)
  library(forcats)
  library(RColorBrewer)
  library(mediation)   # for mediation analysis
})

# ── 1. Configuration ──────────────────────────────────────────────────────────
set.seed(42)
RESULTS_DIR  <- here::here("results")
FIGURES_DIR  <- file.path(RESULTS_DIR, "figures")
TABLES_DIR   <- file.path(RESULTS_DIR, "tables")
DATA_DIR     <- here::here("data")

dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES_DIR,  showWarnings = FALSE)

message("=== Aim 2: Gemcitabine Deconfounding Analysis ===")

# ── 2. Load Cached Data from Aim 1 ────────────────────────────────────────────
vst_cache      <- file.path(DATA_DIR, "tcga_paad_vst.rds")
clinical_cache <- file.path(DATA_DIR, "tcga_paad_clinical.rds")
subtype_cache  <- file.path(TABLES_DIR, "aim1_subtype_assignments.csv")

if (!file.exists(vst_cache)) {
  stop("Run aim1_subtype_stratification.R first to generate cached data.")
}

vst_mat      <- readRDS(vst_cache)
clinical_raw <- readRDS(clinical_cache)
subtype_df   <- read.csv(subtype_cache)

message("Loaded cached TCGA-PAAD data.")

# ── 3. Build Clinical + Treatment Dataframe ────────────────────────────────────
# Uses shared utility — handles TCGAbiolinks column name differences
# and extracts treatment flags from all available clinical fields
source(here::here("R", "utils_clinical.R"))

message("Building clinical + treatment dataframe...")
clin_df <- build_clinical_df(clinical_raw)

# Report treatment distribution
tx_summary <- clin_df %>%
  count(treatment_group) %>%
  mutate(pct = round(100 * n / sum(n), 1))
message("Treatment distribution:")
print(tx_summary)

# ── 4. Build Master Analysis Dataframe ────────────────────────────────────────

# ── 5. Merge with Subtype Labels and GPRC5A Expression ────────────────────────

# Build GPRC5A expression lookup from the VST matrix
gprc5a_vec <- as.numeric(vst_mat["GPRC5A", ])
gprc5a_lookup <- data.frame(
  barcode     = colnames(vst_mat),
  gprc5a_expr = gprc5a_vec,
  stringsAsFactors = FALSE
)

# Subset subtype columns explicitly (avoids dplyr::select namespace conflicts)
subtype_cols <- c("barcode", "subtype_moffitt", "classical_score",
                  "basal_score", "subtype_score")
subtype_sub  <- subtype_df[, intersect(subtype_cols, colnames(subtype_df))]

analysis_df <- clin_df %>%
  dplyr::inner_join(subtype_sub,    by = "barcode") %>%
  dplyr::inner_join(gprc5a_lookup,  by = "barcode") %>%
  dplyr::mutate(
    gprc5a_high = ifelse(gprc5a_expr >= median(gprc5a_expr, na.rm = TRUE),
                         "High", "Low"),
    gem_treated = ifelse(received_gemcitabine, "Gemcitabine", "No Gemcitabine")
  ) %>%
  dplyr::filter(!is.na(gprc5a_expr))

message(sprintf("\nAnalysis cohort: n=%d", nrow(analysis_df)))
message(sprintf("  Gemcitabine-treated:  n=%d (%.1f%%)",
                sum(analysis_df$received_gemcitabine, na.rm=TRUE),
                100 * mean(analysis_df$received_gemcitabine, na.rm=TRUE)))
message(sprintf("  Treatment-naive:      n=%d (%.1f%%)",
                sum(analysis_df$treatment_naive, na.rm=TRUE),
                100 * mean(analysis_df$treatment_naive, na.rm=TRUE)))

# Save treatment summary
write.csv(
  analysis_df %>%
    group_by(treatment_group, survival_group) %>%
    summarise(n = n(),
              median_gprc5a = round(median(gprc5a_expr), 3),
              median_OS_mo  = round(median(OS_months), 1),
              .groups = "drop"),
  file.path(TABLES_DIR, "aim2_treatment_summary.csv"), row.names = FALSE
)

# ── 5. Core Test: GPRC5A Expression by Treatment × Survival ───────────────────
# The key question: does treatment status explain the paradoxical direction?
# If gemcitabine-treated survivors have artificially elevated GPRC5A,
# we expect: Alive+Gem > Alive+NoGem, and Dead+Gem > Dead+NoGem

message("Generating treatment × survival expression analysis...")

palette_treat   <- c("Gemcitabine"    = "#E67E22",
                     "No Gemcitabine" = "#8E44AD")
palette_survive <- c("Alive" = "#27AE60", "Deceased" = "#E74C3C")

# Figure 1: 2×2 expression grid (treatment × survival)
p1 <- ggplot(analysis_df,
             aes(x = gem_treated, y = gprc5a_expr, fill = survival_group)) +
  geom_boxplot(alpha = 0.8, outlier.shape = 21, outlier.size = 1.5,
               width = 0.6, position = position_dodge(0.7)) +
  geom_jitter(aes(color = survival_group),
              position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7),
              alpha = 0.35, size = 1.2) +
  scale_fill_manual(values  = palette_survive) +
  scale_color_manual(values = palette_survive) +
  stat_compare_means(aes(group = survival_group), method = "wilcox.test",
                     label = "p.format", label.y.npc = 0.92, size = 3.5) +
  labs(
    title    = "GPRC5A Expression: Treatment × Survival Status",
    subtitle = "Key test: is GPRC5A elevated in alive+gemcitabine vs. alive+no-gemcitabine?\nIf yes, treatment confounds the survival signal.",
    x        = "Gemcitabine Treatment",
    y        = "GPRC5A Expression (VST)",
    fill     = "Survival", color = "Survival"
  ) +
  theme_bw(base_size = 13) +
  theme(plot.title    = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(color = "grey40", size = 10),
        legend.position = "bottom")

# Figure 2: All 4 treatment groups (None / Gem / Other Chemo / Radiation)
p2 <- ggplot(analysis_df,
             aes(x = fct_reorder(treatment_group, gprc5a_expr, median),
                 y = gprc5a_expr, fill = treatment_group)) +
  geom_boxplot(alpha = 0.8, outlier.shape = 21, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1) +
  scale_fill_brewer(palette = "Set2") +
  coord_flip() +
  facet_wrap(~survival_group) +
  stat_compare_means(method = "kruskal.test", label.y.npc = 0.05, size = 3.5) +
  labs(
    title = "GPRC5A Expression Across All Treatment Groups",
    x     = NULL, y = "GPRC5A Expression (VST)",
    fill  = "Treatment"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "none")

pdf(file.path(FIGURES_DIR, "aim2_treatment_boxplot.pdf"), width = 14, height = 7)
  print(p1 + p2 + plot_layout(widths = c(1, 1.4)))
dev.off()

# ── 6. Kaplan-Meier Within Treatment Strata ────────────────────────────────────
message("Generating KM curves within treatment strata...")

km_plots <- list()

strata <- list(
  "All patients"        = analysis_df,
  "Treatment-naive"     = analysis_df[analysis_df$treatment_naive == TRUE, ],
  "Gemcitabine-treated" = analysis_df[analysis_df$received_gemcitabine == TRUE, ]
)

for (stratum_name in names(strata)) {
  sub <- strata[[stratum_name]]
  n   <- nrow(sub)
  
  # Skip strata too small to fit a KM curve
  if (n < 10) {
    message(sprintf("  Skipping '%s' — only n=%d samples", stratum_name, n))
    next
  }
  
  # Re-compute gprc5a_high within this stratum so the median split
  # reflects the local distribution, and both groups are guaranteed to exist
  sub <- sub %>%
    dplyr::mutate(
      gprc5a_high = ifelse(gprc5a_expr >= median(gprc5a_expr, na.rm = TRUE),
                           "GPRC5A High", "GPRC5A Low")
    )
  
  # Check both groups exist after the split
  groups_present <- sort(unique(sub$gprc5a_high))
  n_groups       <- length(groups_present)
  
  if (n_groups < 2) {
    message(sprintf("  Skipping '%s' — only one GPRC5A group present", stratum_name))
    next
  }
  
  # Build palette and legend labels from the actual groups present
  full_palette <- c("GPRC5A High" = "#E74C3C", "GPRC5A Low" = "#2E75B6")
  this_palette <- full_palette[groups_present]
  this_labs    <- groups_present   # labels match strata names exactly
  
  fit <- survfit(Surv(OS_months, OS_event) ~ gprc5a_high, data = sub)
  
  km_plots[[stratum_name]] <- tryCatch(
    ggsurvplot(
      fit,
      data              = sub,
      pval              = TRUE,
      pval.method       = TRUE,
      conf.int          = TRUE,
      risk.table        = TRUE,
      risk.table.height = 0.28,
      palette           = unname(this_palette),
      title             = sprintf("GPRC5A Prognosis — %s (n=%d)", stratum_name, n),
      subtitle          = "Paradox resolves here if treatment-naive HR flips direction",
      xlab              = "Time (months)",
      ylab              = "Overall Survival",
      legend.labs       = this_labs,
      legend.title      = "GPRC5A",
      ggtheme           = theme_bw(base_size = 12),
      font.title        = c(13, "bold"),
      surv.median.line  = "hv"
    ),
    error = function(e) {
      message(sprintf("  KM plot failed for '%s': %s", stratum_name, e$message))
      NULL
    }
  )
}

pdf(file.path(FIGURES_DIR, "aim2_km_by_treatment_stratum.pdf"), width = 11, height = 8)
  for (nm in names(km_plots)) {
    if (!is.null(km_plots[[nm]])) print(km_plots[[nm]])
  }
dev.off()

# ── 7. Multivariable Cox Models ────────────────────────────────────────────────
# Test GPRC5A's independent prognostic value after adjusting for:
# (a) gemcitabine treatment
# (b) subtype
# (c) both together

message("Running multivariable Cox models...")

cox_models <- list(
  "Univariate_GPRC5A" = coxph(
    Surv(OS_months, OS_event) ~ gprc5a_expr,
    data = analysis_df
  ),
  "Adjusted_Treatment" = coxph(
    Surv(OS_months, OS_event) ~ gprc5a_expr + received_gemcitabine,
    data = analysis_df
  ),
  "Adjusted_Subtype" = coxph(
    Surv(OS_months, OS_event) ~ gprc5a_expr + subtype_moffitt,
    data = analysis_df
  ),
  "Fully_Adjusted" = coxph(
    Surv(OS_months, OS_event) ~ gprc5a_expr + received_gemcitabine +
      subtype_moffitt + age_years + stage_binary,
    data = analysis_df
  ),
  "Interaction_GemxGPRC5A" = coxph(
    Surv(OS_months, OS_event) ~ gprc5a_expr * received_gemcitabine,
    data = analysis_df
  )
)

# Extract and compile results
cox_results <- dplyr::bind_rows(
  lapply(names(cox_models), function(nm) {
    res <- broom::tidy(cox_models[[nm]], exponentiate = TRUE, conf.int = TRUE)
    data.frame(
      Model    = nm,
      term     = res$term,
      HR       = round(res$estimate,  3),
      CI_lower = round(res$conf.low,  3),
      CI_upper = round(res$conf.high, 3),
      p_value  = format.pval(res$p.value, digits = 3),
      stringsAsFactors = FALSE
    )
  })
)

write.csv(cox_results,
          file.path(TABLES_DIR, "aim2_multivariable_cox.csv"),
          row.names = FALSE)

# ── 8. Forest Plot: Multivariable Cox Results ─────────────────────────────────
message("Generating forest plot...")

# Extract GPRC5A coefficient across all models
gprc5a_forest <- cox_results %>%
  filter(term == "gprc5a_expr") %>%
  mutate(
    model_label = recode(Model,
      "Univariate_GPRC5A"     = "Unadjusted",
      "Adjusted_Treatment"    = "Adjusted: Gemcitabine",
      "Adjusted_Subtype"      = "Adjusted: Subtype",
      "Fully_Adjusted"        = "Fully adjusted\n(Gem + Subtype + Age + Stage)",
      "Interaction_GemxGPRC5A"= "Interaction model\n(GPRC5A × Gemcitabine)"
    ),
    sig = as.numeric(p_value) < 0.05
  )

p_forest <- ggplot(gprc5a_forest,
                   aes(x = HR, y = fct_rev(fct_inorder(model_label)),
                       xmin = CI_lower, xmax = CI_upper,
                       color = sig, shape = sig)) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "grey50", linewidth = 0.9) +
  geom_errorbarh(height = 0.25, linewidth = 1.2) +
  geom_point(size = 4.5) +
  scale_color_manual(values = c("TRUE" = "#C0392B", "FALSE" = "#7F8C8D"),
                     labels = c("p ≥ 0.05", "p < 0.05")) +
  scale_shape_manual(values = c("TRUE" = 18, "FALSE" = 16),
                     labels = c("p ≥ 0.05", "p < 0.05")) +
  scale_x_log10(breaks = c(0.5, 0.75, 1.0, 1.25, 1.5, 2.0)) +
  # Annotate HR and CI
  geom_text(aes(label = sprintf("HR=%.2f [%.2f–%.2f]", HR, CI_lower, CI_upper)),
            hjust = -0.1, size = 3.2, color = "grey30") +
  labs(
    title    = "GPRC5A Hazard Ratio Before and After Treatment Adjustment",
    subtitle = "Key finding: does HR direction change after adjusting for gemcitabine?",
    x        = "Hazard Ratio for GPRC5A (log scale)\n← Protective | Harmful →",
    y        = NULL,
    color    = "Significance",
    shape    = "Significance",
    caption  = "HR > 1: higher GPRC5A = worse outcome (oncogenic)\nHR < 1: higher GPRC5A = better outcome (paradoxical)"
  ) +
  expand_limits(x = 2.5) +
  theme_bw(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey40", size = 10),
    legend.position = "bottom",
    panel.grid.major.y = element_blank()
  )

ggsave(file.path(FIGURES_DIR, "aim2_forest_multivariable.pdf"),
       plot = p_forest, width = 11, height = 7, device = "pdf")

# ── 9. Mediation Analysis: Does gemcitabine mediate GPRC5A → survival? ─────────
# Tests: gemcitabine treatment ~ mediator between GPRC5A expression and outcome
# Model: GPRC5A → Gemcitabine treatment → Survival
# If significant mediation: treatment is a true confounder on the causal path

message("Running mediation analysis...")

# Step 1: Does GPRC5A predict gemcitabine treatment receipt?
med_model <- glm(received_gemcitabine ~ gprc5a_expr + age_years + stage_binary,
                 data = analysis_df, family = binomial())

# Step 2: Outcome model with mediator
out_model <- glm(OS_event ~ gprc5a_expr + received_gemcitabine +
                   age_years + stage_binary,
                 data = analysis_df, family = binomial())

mediation_result <- tryCatch({
  med_out <- mediate(med_model, out_model,
                     treat = "gprc5a_expr", mediator = "received_gemcitabine",
                     treat.value = quantile(analysis_df$gprc5a_expr, 0.75),
                     control.value = quantile(analysis_df$gprc5a_expr, 0.25),
                     sims = 500, boot = TRUE)
  summary(med_out)
}, error = function(e) {
  message("  Mediation analysis skipped (insufficient data): ", e$message)
  NULL
})

# ── 10. Treatment-Stratified Cox (primary deconfounding result) ────────────────
message("Running treatment-stratified Cox models...")

strat_cox_results <- dplyr::bind_rows(lapply(names(strata), function(nm) {
  sub <- strata[[nm]]
  if (nrow(sub) < 20) return(NULL)
  fit <- coxph(Surv(OS_months, OS_event) ~ gprc5a_expr, data = sub)
  res <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  data.frame(
    Stratum  = nm,
    n        = nrow(sub),
    HR       = round(res$estimate,  3),
    CI_lower = round(res$conf.low,  3),
    CI_upper = round(res$conf.high, 3),
    p_value  = format.pval(res$p.value, digits = 3),
    stringsAsFactors = FALSE
  )
}))

write.csv(strat_cox_results,
          file.path(TABLES_DIR, "aim2_treatment_stratified_cox.csv"),
          row.names = FALSE)

# ── 11. Figure: Treatment-Stratified Forest ────────────────────────────────────
p_strat_forest <- ggplot(strat_cox_results,
                          aes(x = HR, y = fct_rev(fct_inorder(Stratum)),
                              xmin = CI_lower, xmax = CI_upper)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(height = 0.2, linewidth = 1.2, color = "#2E75B6") +
  geom_point(size = 4.5, color = "#2E75B6") +
  geom_text(aes(label = sprintf("n=%d | HR=%.2f [%.2f–%.2f] p=%s",
                                 n, HR, CI_lower, CI_upper, p_value)),
            hjust = -0.05, size = 3.3, color = "grey30") +
  scale_x_log10() +
  expand_limits(x = 3) +
  labs(
    title    = "GPRC5A Hazard Ratio Stratified by Treatment Status",
    subtitle = "Resolution of paradox if treatment-naive HR flips direction",
    x        = "Hazard Ratio (log scale)",
    y        = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.major.y = element_blank())

ggsave(file.path(FIGURES_DIR, "aim2_treatment_stratified_forest.pdf"),
       plot = p_strat_forest, width = 11, height = 5, device = "pdf")

# ── 12. Summary ───────────────────────────────────────────────────────────────
message("\n=== Aim 2 Complete ===")
message("Key result — GPRC5A HR by treatment stratum:")
print(strat_cox_results[, c("Stratum", "n", "HR", "p_value")])

message("\nInterpretation guide:")
message("  - If treatment-naive HR > 1 (p < 0.05): oncogenic signal confirmed")
message("  - If treatment-naive HR < 1 (p < 0.05): true suppressive signal")
message("  - If treatment-naive HR crosses 1 or p > 0.05: paradox persists")
message("  - Compare to fully-adjusted Cox: aim2_multivariable_cox.csv")

message("\nOutputs saved:")
message(paste0("  ", FIGURES_DIR, "/aim2_*.pdf"))
message(paste0("  ", TABLES_DIR,  "/aim2_*.csv"))
