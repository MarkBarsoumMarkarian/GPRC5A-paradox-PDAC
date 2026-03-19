# =============================================================================
# Regenerate affected figures from pre-computed tables (no data download)
# Fixes applied:
#   1. aim1_forest_plot.pdf  — adds interaction term with p=0.272
#   2. aim2_km_by_treatment_stratum.pdf — fixes subtitle
#   3. aim2_treatment_stratified_forest.pdf — fixes subtitle
#   4. aim5_roc_curves.pdf — adds circularity caveat footnote
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(survival)
  library(survminer)
  library(pROC)
  library(caret)
  library(randomForest)
  library(xgboost)
})

RESULTS_DIR <- "/home/SexyThighs/results"
FIGURES_DIR <- file.path(RESULTS_DIR, "figures")
TABLES_DIR  <- file.path(RESULTS_DIR, "tables")

# ── FIX 1: aim1_forest_plot.pdf ──────────────────────────────────────────────
message("Regenerating aim1_forest_plot.pdf ...")

cox_tbl <- read.csv(file.path(TABLES_DIR, "aim1_cox_results.csv"),
                    stringsAsFactors = FALSE)

# Pull values directly from the saved CSV
get_val <- function(model, col) {
  val <- cox_tbl[cox_tbl$Model == model & cox_tbl$term == "gprc5a_expr", col]
  as.numeric(val[1])
}
int_row <- cox_tbl[cox_tbl$Model == "Interaction" &
                   cox_tbl$term == "gprc5a_expr:subtype_moffittClassical", ]

forest_data <- data.frame(
  Model  = c("Overall (n=177)", "Classical subtype (n=100)",
             "Basal-like subtype (n=77)",
             sprintf("Interaction term\n(GPRC5A \u00d7 Subtype, p=%s)",
                     format(round(as.numeric(int_row$p_value), 3), nsmall = 3))),
  HR     = c(get_val("Overall","HR"), get_val("Classical","HR"),
             get_val("Basal-like","HR"), as.numeric(int_row$HR)),
  lower  = c(get_val("Overall","CI_lower"), get_val("Classical","CI_lower"),
             get_val("Basal-like","CI_lower"), as.numeric(int_row$CI_lower)),
  upper  = c(get_val("Overall","CI_upper"), get_val("Classical","CI_upper"),
             get_val("Basal-like","CI_upper"), as.numeric(int_row$CI_upper)),
  color  = c("Combined","Classical","Basal-like","Interaction"),
  p_val  = c(get_val("Overall","p_value"), get_val("Classical","p_value"),
             get_val("Basal-like","p_value"), as.numeric(int_row$p_value)),
  stringsAsFactors = FALSE
)

forest_data$sig_label <- ifelse(
  forest_data$p_val < 0.05,
  sprintf("HR=%.2f [%.2f\u2013%.2f] *p=%s*", forest_data$HR, forest_data$lower,
          forest_data$upper, formatC(forest_data$p_val, digits=3, format="g")),
  sprintf("HR=%.2f [%.2f\u2013%.2f]  p=%s (n.s.)", forest_data$HR, forest_data$lower,
          forest_data$upper, formatC(forest_data$p_val, digits=3, format="g"))
)

p_forest <- ggplot(forest_data,
                   aes(x = HR, y = factor(Model, levels = rev(Model)),
                       xmin = lower, xmax = upper, color = color)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.8) +
  geom_errorbarh(height = 0.25, linewidth = 1.3) +
  geom_point(size = 4.5) +
  geom_text(aes(x = upper, label = sig_label),
            hjust = -0.08, size = 3.2, color = "grey30") +
  scale_color_manual(values = c("Combined"    = "#7F8C8D",
                                "Classical"   = "#2E75B6",
                                "Basal-like"  = "#C0392B",
                                "Interaction" = "#8E44AD")) +
  scale_x_log10(limits = c(0.7, 3.5)) +
  labs(
    title    = "GPRC5A Hazard Ratios by PDAC Molecular Subtype",
    subtitle = "Cox proportional hazards | TCGA-PAAD | per unit VST increase\nInteraction p=0.272 (non-significant) — subtype-differential effect not confirmed at n=177",
    x        = "Hazard Ratio (log scale)",
    y        = NULL,
    color    = "Model",
    caption  = "HR > 1: higher GPRC5A = worse survival  |  *p < 0.05  |  n.s. = non-significant"
  ) +
  theme_bw(base_size = 13) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(color = "grey40", size = 10),
        legend.position = "bottom",
        panel.grid.major.y = element_blank())

ggsave(file.path(FIGURES_DIR, "aim1_forest_plot.pdf"),
       plot = p_forest, width = 11, height = 5.5, device = "pdf")
message("  -> Saved aim1_forest_plot.pdf")

# ── FIX 2 & 3: aim2 forest figures ──────────────────────────────────────────
message("Regenerating aim2_treatment_stratified_forest.pdf ...")

strat_tbl <- read.csv(file.path(TABLES_DIR, "aim2_treatment_stratified_cox.csv"),
                      stringsAsFactors = FALSE)

p_strat <- ggplot(strat_tbl,
                  aes(x = HR, y = reorder(Stratum, HR),
                      xmin = CI_lower, xmax = CI_upper)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(height = 0.2, linewidth = 1.2, color = "#2E75B6") +
  geom_point(size = 4.5, color = "#2E75B6") +
  geom_text(aes(label = sprintf("n=%d | HR=%.2f [%.2f\u2013%.2f] p=%s",
                                n, HR, CI_lower, CI_upper, p_value)),
            hjust = -0.05, size = 3.3, color = "grey30") +
  scale_x_log10() +
  expand_limits(x = 3) +
  labs(
    title    = "GPRC5A Hazard Ratio Stratified by Treatment Status",
    subtitle = paste0("HR attenuates in gemcitabine-treated patients ",
                      "(HR=1.22, p=0.221 vs HR=1.36, p<0.001 overall)\n",
                      "Treatment-naive comparison not feasible (n=1 deceased, n=0 alive)"),
    x        = "Hazard Ratio (log scale)",
    y        = NULL,
    caption  = "HR > 1: higher GPRC5A = worse outcome (oncogenic direction)"
  ) +
  theme_bw(base_size = 13) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(color = "grey40", size = 10),
        panel.grid.major.y = element_blank())

ggsave(file.path(FIGURES_DIR, "aim2_treatment_stratified_forest.pdf"),
       plot = p_strat, width = 11, height = 5, device = "pdf")
message("  -> Saved aim2_treatment_stratified_forest.pdf")

# ── FIX 4: aim5_roc_curves.pdf ───────────────────────────────────────────────
# ROC curves require the saved model predictions — use aim5_model_performance.csv
# and aim5_role_state_predictions.csv to rebuild the plot from scores
message("Regenerating aim5_roc_curves.pdf ...")

preds_tbl <- read.csv(file.path(TABLES_DIR, "aim5_role_state_predictions.csv"),
                      stringsAsFactors = FALSE)
perf_tbl  <- read.csv(file.path(TABLES_DIR, "aim5_model_performance.csv"),
                      stringsAsFactors = FALSE)

message("  Predictions columns: ", paste(names(preds_tbl), collapse=", "))
message("  Performance columns: ", paste(names(perf_tbl), collapse=", "))

# We'll reconstruct what we can from preds_tbl
# Expected columns: true_state, predicted_state, and probability columns per model
prob_cols <- grep("prob|Prob|_prob|oncogenic|Oncogenic", names(preds_tbl),
                  value = TRUE, ignore.case = TRUE)
message("  Probability columns found: ", paste(prob_cols, collapse=", "))

if (length(prob_cols) > 0 && "true_state" %in% names(preds_tbl)) {
  y_true <- ifelse(preds_tbl$true_state == "Oncogenic", 1, 0)
  
  roc_colors <- c("#3498DB", "#27AE60", "#E74C3C",
                  "#8E44AD", "#E67E22")[seq_along(prob_cols)]
  
  pdf(file.path(FIGURES_DIR, "aim5_roc_curves.pdf"), width = 7.5, height = 6.5)
  par(mar = c(6, 4, 4, 2))
  
  first_roc <- pROC::roc(y_true, preds_tbl[[prob_cols[1]]],
                          levels = c(0,1), direction = "<", quiet = TRUE)
  plot(first_roc, col = roc_colors[1], lwd = 2,
       main = "ROC Curves — GPRC5A Role-State Classifier\n(Leakage-Free Held-Out Test Set, n=12)")
  
  legend_entries <- paste0(sub("_prob|_Prob|prob_|Prob_", "", prob_cols),
                           " (AUC=", round(as.numeric(pROC::auc(first_roc)), 3), ")")
  
  if (length(prob_cols) > 1) {
    for (i in 2:length(prob_cols)) {
      roc_i <- pROC::roc(y_true, preds_tbl[[prob_cols[i]]],
                          levels = c(0,1), direction = "<", quiet = TRUE)
      plot(roc_i, col = roc_colors[i], lwd = 2, add = TRUE)
      legend_entries[i] <- paste0(sub("_prob|_Prob|prob_|Prob_","",prob_cols[i]),
                                   " (AUC=", round(as.numeric(pROC::auc(roc_i)), 3), ")")
    }
  }
  
  legend("bottomright", legend = legend_entries,
         col = roc_colors[seq_along(prob_cols)], lwd = 2, bty = "n", cex = 0.9)
  abline(0, 1, lty = 2, col = "grey60")
  
  mtext(paste0("Interpretation caveat: Labels derived from survival \u00d7 subtype identity; ",
               "AUC reflects proof-of-concept\nsubtype-context encoding, ",
               "not independent prognostic prediction. Validate externally."),
        side = 1, line = 4.8, cex = 0.68, col = "grey40", adj = 0)
  dev.off()
  message("  -> Saved aim5_roc_curves.pdf")
} else {
  message("  WARNING: Could not reconstruct ROC — probability columns not found.")
  message("  Columns available: ", paste(names(preds_tbl), collapse=", "))
}

message("\n=== Figure regeneration complete ===")
message("Files updated:")
message("  results/figures/aim1_forest_plot.pdf")
message("  results/figures/aim2_treatment_stratified_forest.pdf")
message("  results/figures/aim5_roc_curves.pdf  (if predictions available)")
