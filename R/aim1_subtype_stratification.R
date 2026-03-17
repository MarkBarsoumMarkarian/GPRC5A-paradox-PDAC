# =============================================================================
# GPRC5A Paradox in PDAC — Aim 1: Subtype-Stratified Analysis
# Project: Decoding the GPRC5A Paradox in Pancreatic Ductal Adenocarcinoma
# Author:  Mark Barsoum Markarian
# Date:    2026
#
# Description:
#   Re-examines GPRC5A expression in TCGA-PAAD stratified by PDAC molecular
#   subtype (classical vs. basal-like) using the Moffitt 2015 gene signature.
#   Tests the hypothesis that the paradoxical inverse expression-survival
#   association from the previous study reflects subtype mixing.
#
# Outputs (saved to results/):
#   - figures/aim1_gprc5a_subtype_boxplot.pdf
#   - figures/aim1_kaplan_meier_by_subtype.pdf
#   - figures/aim1_expression_heatmap.pdf
#   - tables/aim1_cox_results.csv
#   - tables/aim1_subtype_assignments.csv
# =============================================================================

# ── 0. Libraries ──────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(DESeq2)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(RColorBrewer)
  library(viridis)
  library(sva)
})

# ── 1. Configuration ──────────────────────────────────────────────────────────
set.seed(42)
RESULTS_DIR   <- here::here("results")
FIGURES_DIR   <- file.path(RESULTS_DIR, "figures")
TABLES_DIR    <- file.path(RESULTS_DIR, "tables")
DATA_DIR      <- here::here("data")
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR,    showWarnings = FALSE)

message("=== Aim 1: GPRC5A Subtype Stratification ===")

# ── 2. Moffitt 2015 Subtype Gene Signatures ───────────────────────────────────
# Moffitt et al. Nat Genet 2015 — 25-gene tumor subtype classifier
# Classical (epithelial) vs Basal-like (mesenchymal/aggressive)

moffitt_classical <- c(
  "BTNL8", "FAM83A", "CEACAM5", "CEACAM6", "LYZ", "TFF1", "TFF2", "TFF3",
  "LGALS4", "CYP3A7", "MYO1A", "CLRN3", "SLC40A1", "ANXA10", "CTSE",
  "AGR2", "ST6GALNAC1", "LOC400573", "VSIG2", "REG4", "LRRC26",
  "MSLN", "KRT20", "LTBP2", "CLDN18"
)

moffitt_basal <- c(
  "VGLL1", "UCA1", "S100A2", "LY6D", "LEMD1", "KRT15", "KRT17", "KRT19",
  "AREG", "SPRR1B", "SPRR2B", "CDKN2A", "ANLN", "FOXC2", "SCEL",
  "CSTA", "DSP", "EMP1", "LAMC2", "SERPINB13", "TNS4", "DYRK3",
  "FAT2", "TP63", "ST14"
)

# Bailey 2016 (Nature) subtypes — supplementary marker genes for cross-validation
bailey_markers <- list(
  ADEX        = c("NR5A2", "RBPJL", "PTF1A", "NR5A2", "MLXIPL"),
  Immunogenic = c("CD3E", "CD4", "CD8A", "TIGIT", "LAG3", "PDCD1"),
  Classical   = c("TFF1", "TFF2", "TFF3", "AGR2", "CEACAM5", "CEACAM6"),
  Squamous    = c("KRT5", "KRT14", "TP63", "SOX2", "EGFR", "S100A2")
)

# ── 3. Download & Cache TCGA-PAAD Data ────────────────────────────────────────
tcga_cache <- file.path(DATA_DIR, "tcga_paad_vst.rds")
clinical_cache <- file.path(DATA_DIR, "tcga_paad_clinical.rds")

if (!file.exists(tcga_cache)) {
  message("Downloading TCGA-PAAD RNA-seq data...")
  
  query <- GDCquery(
    project           = "TCGA-PAAD",
    data.category     = "Transcriptome Profiling",
    data.type         = "Gene Expression Quantification",
    workflow.type     = "STAR - Counts"
  )
  GDCdownload(query, method = "api", files.per.chunk = 10)
  se <- GDCprepare(query)
  
  # Filter to primary tumors only
  se <- se[, se$sample_type == "Primary Tumor"]
  
  # VST normalization
  dds <- DESeqDataSet(se, design = ~1)
  dds <- dds[rowSums(counts(dds)) > 10, ]
  vst_data <- vst(dds, blind = TRUE)
  vst_mat  <- assay(vst_data)
  
  # Use HGNC symbols
  rownames(vst_mat) <- rowData(vst_data)$gene_name
  vst_mat <- vst_mat[!duplicated(rownames(vst_mat)), ]
  
  saveRDS(vst_mat, tcga_cache)
  saveRDS(colData(se), clinical_cache)
  message("TCGA-PAAD data downloaded and cached.")
  
} else {
  message("Loading cached TCGA-PAAD data...")
  vst_mat  <- readRDS(tcga_cache)
  clinical <- readRDS(clinical_cache)
}

clinical <- readRDS(clinical_cache)

# ── 4. Clinical Data Preparation ──────────────────────────────────────────────
# Uses shared utility that handles TCGAbiolinks column name differences
source(here::here("R", "utils_clinical.R"))

message("Building clinical dataframe...")
clin_df <- build_clinical_df(clinical)

# ── 5. Molecular Subtype Classification (Moffitt SSP) ─────────────────────────
# Single-Sample Predictor: compute mean expression of classical vs. basal
# signature genes, assign subtype by higher z-score centroid distance

classify_moffitt_ssp <- function(expr_mat, classical_genes, basal_genes) {
  # Subset to available signature genes
  classical_avail <- intersect(classical_genes, rownames(expr_mat))
  basal_avail     <- intersect(basal_genes, rownames(expr_mat))
  
  message(sprintf("  Classical signature: %d/%d genes available",
                  length(classical_avail), length(classical_genes)))
  message(sprintf("  Basal signature: %d/%d genes available",
                  length(basal_avail), length(basal_genes)))
  
  # Z-score each gene across samples
  expr_z <- t(scale(t(expr_mat[c(classical_avail, basal_avail), ])))
  
  # Mean signature score per sample
  classical_score <- colMeans(expr_z[classical_avail, ], na.rm = TRUE)
  basal_score     <- colMeans(expr_z[basal_avail, ],     na.rm = TRUE)
  
  subtype <- ifelse(classical_score > basal_score, "Classical", "Basal-like")
  
  data.frame(
    barcode          = names(subtype),
    subtype_moffitt  = subtype,
    classical_score  = classical_score,
    basal_score      = basal_score,
    subtype_score    = classical_score - basal_score,
    stringsAsFactors = FALSE
  )
}

message("Classifying PDAC molecular subtypes...")
subtype_df <- classify_moffitt_ssp(vst_mat, moffitt_classical, moffitt_basal)

# Merge with clinical
analysis_df <- clin_df %>%
  inner_join(subtype_df, by = "barcode") %>%
  mutate(
    gprc5a_expr = as.numeric(vst_mat["GPRC5A", barcode]),
    gprc5a_high = ifelse(gprc5a_expr >= median(gprc5a_expr, na.rm = TRUE),
                         "High", "Low")
  ) %>%
  filter(!is.na(gprc5a_expr))

message(sprintf("Analysis cohort: n=%d (Classical=%d, Basal-like=%d)",
                nrow(analysis_df),
                sum(analysis_df$subtype_moffitt == "Classical"),
                sum(analysis_df$subtype_moffitt == "Basal-like")))

# Save subtype assignments
write.csv(analysis_df %>%
            select(barcode, subtype_moffitt, classical_score, basal_score,
                   subtype_score, vital_status, OS_months, gprc5a_expr),
          file.path(TABLES_DIR, "aim1_subtype_assignments.csv"),
          row.names = FALSE)

# ── 6. Figure 1: GPRC5A Expression by Subtype × Survival ─────────────────────
message("Generating Figure 1: GPRC5A expression by subtype x survival...")

subtype_palette <- c("Classical" = "#2E75B6", "Basal-like" = "#C0392B")
survival_palette <- c("Alive" = "#27AE60", "Deceased" = "#E74C3C")

p1 <- ggplot(analysis_df, aes(x = survival_group, y = gprc5a_expr,
                               fill = survival_group)) +
  geom_boxplot(alpha = 0.8, outlier.shape = 21, outlier.size = 1.5,
               outlier.fill = "grey60", width = 0.6) +
  geom_jitter(aes(color = survival_group), width = 0.15, alpha = 0.4, size = 1) +
  facet_wrap(~subtype_moffitt, ncol = 2) +
  scale_fill_manual(values  = survival_palette) +
  scale_color_manual(values = survival_palette) +
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     label.x = 1.3, size = 4) +
  labs(
    title    = "GPRC5A Expression by PDAC Molecular Subtype and Survival Status",
    subtitle = "Moffitt 2015 single-sample predictor classification | TCGA-PAAD",
    x        = "Survival Status",
    y        = "GPRC5A Expression (VST-normalized)",
    fill     = "Survival",
    color    = "Survival",
    caption  = "Wilcoxon rank-sum test p-values shown"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey40", size = 11),
    strip.background = element_rect(fill = "#1F3864"),
    strip.text       = element_text(color = "white", face = "bold", size = 12),
    legend.position  = "bottom"
  )

ggsave(file.path(FIGURES_DIR, "aim1_gprc5a_subtype_boxplot.pdf"),
       plot = p1, width = 10, height = 6, device = "pdf")

# ── 7. Figure 2: Kaplan-Meier by Subtype × GPRC5A High/Low ───────────────────
message("Generating Figure 2: Kaplan-Meier survival curves...")

km_plots <- list()

for (subtype in c("Classical", "Basal-like")) {
  sub_data <- analysis_df %>% filter(subtype_moffitt == subtype)
  
  fit <- survfit(Surv(OS_months, OS_event) ~ gprc5a_high, data = sub_data)
  
  km_plots[[subtype]] <- ggsurvplot(
    fit,
    data          = sub_data,
    pval          = TRUE,
    pval.method   = TRUE,
    conf.int      = TRUE,
    risk.table    = TRUE,
    risk.table.height = 0.28,
    palette       = c("High" = "#E74C3C", "Low" = "#2E75B6"),
    title         = paste0("GPRC5A Prognostic Association — ", subtype, " Subtype"),
    xlab          = "Time (months)",
    ylab          = "Overall Survival Probability",
    legend.labs   = c("GPRC5A High", "GPRC5A Low"),
    legend.title  = "GPRC5A\nExpression",
    ggtheme       = theme_bw(base_size = 12),
    font.title    = c(13, "bold"),
    surv.median.line = "hv"
  )
}

# Save KM plots
pdf(file.path(FIGURES_DIR, "aim1_kaplan_meier_by_subtype.pdf"),
    width = 12, height = 8)
  print(km_plots[["Classical"]])
  print(km_plots[["Basal-like"]])
dev.off()

# ── 8. Figure 3: Subtype Score vs GPRC5A Scatter ─────────────────────────────
message("Generating Figure 3: Subtype score vs GPRC5A correlation...")

p3 <- ggplot(analysis_df, aes(x = subtype_score, y = gprc5a_expr,
                               color = survival_group, shape = subtype_moffitt)) +
  geom_point(alpha = 0.7, size = 2.5) +
  geom_smooth(aes(group = subtype_moffitt, color = subtype_moffitt),
              method = "lm", se = TRUE, linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50") +
  scale_color_manual(values = c(survival_palette, subtype_palette)) +
  scale_shape_manual(values = c("Classical" = 16, "Basal-like" = 17)) +
  labs(
    title    = "GPRC5A Expression vs. Molecular Subtype Score",
    subtitle = "Positive score = Classical; Negative = Basal-like",
    x        = "Moffitt Subtype Score (Classical - Basal-like)",
    y        = "GPRC5A Expression (VST)",
    color    = "Group", shape = "Subtype"
  ) +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(FIGURES_DIR, "aim1_subtype_score_scatter.pdf"),
       plot = p3, width = 10, height = 7, device = "pdf")

# ── 9. Figure 4: Expression Heatmap of Signature + GPRC5A ────────────────────
message("Generating Figure 4: Signature gene heatmap...")

sig_genes  <- c("GPRC5A",
                intersect(moffitt_classical, rownames(vst_mat))[1:10],
                intersect(moffitt_basal,     rownames(vst_mat))[1:10])

heat_mat <- vst_mat[sig_genes, analysis_df$barcode]
heat_mat_z <- t(scale(t(heat_mat)))

annotation_col <- data.frame(
  Subtype   = analysis_df$subtype_moffitt,
  Survival  = analysis_df$survival_group,
  GPRC5A_Level = analysis_df$gprc5a_high,
  row.names = analysis_df$barcode
)

ann_colors <- list(
  Subtype      = subtype_palette,
  Survival     = survival_palette,
  GPRC5A_Level = c("High" = "#E74C3C", "Low" = "#2E75B6")
)

# Sort samples: Classical | Basal-like, then by survival
heat_order <- analysis_df %>%
  arrange(subtype_moffitt, desc(subtype_score)) %>%
  pull(barcode)

annotation_row <- data.frame(
  Gene_Group = c("Target", rep("Classical Sig", 10), rep("Basal Sig", 10)),
  row.names  = sig_genes
)
ann_row_colors <- list(Gene_Group = c("Target"       = "#F39C12",
                                       "Classical Sig" = "#2E75B6",
                                       "Basal Sig"     = "#C0392B"))

pdf(file.path(FIGURES_DIR, "aim1_expression_heatmap.pdf"), width = 14, height = 8)
pheatmap(
  heat_mat_z[, heat_order],
  annotation_col  = annotation_col[heat_order, ],
  annotation_row  = annotation_row,
  annotation_colors = c(ann_colors, ann_row_colors),
  cluster_cols    = FALSE,
  cluster_rows    = TRUE,
  show_colnames   = FALSE,
  color           = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  breaks          = seq(-3, 3, length.out = 101),
  fontsize_row    = 10,
  main            = "GPRC5A & Moffitt Signature Genes | TCGA-PAAD\nOrdered by Molecular Subtype",
  border_color    = NA
)
dev.off()

# ── 10. Cox Proportional Hazards Models ───────────────────────────────────────
message("Running Cox proportional hazards models...")

cox_results <- list()

# A. Overall (replicates original finding)
cox_overall <- coxph(Surv(OS_months, OS_event) ~ gprc5a_expr, data = analysis_df)
cox_results[["Overall"]] <- broom::tidy(cox_overall, exponentiate = TRUE, conf.int = TRUE)

# B. Within Classical subtype
cox_classical <- coxph(Surv(OS_months, OS_event) ~ gprc5a_expr,
                        data = filter(analysis_df, subtype_moffitt == "Classical"))
cox_results[["Classical"]] <- broom::tidy(cox_classical, exponentiate = TRUE, conf.int = TRUE)

# C. Within Basal-like subtype
cox_basal <- coxph(Surv(OS_months, OS_event) ~ gprc5a_expr,
                   data = filter(analysis_df, subtype_moffitt == "Basal-like"))
cox_results[["Basal-like"]] <- broom::tidy(cox_basal, exponentiate = TRUE, conf.int = TRUE)

# D. Interaction model: does subtype modify GPRC5A's effect?
cox_interaction <- coxph(Surv(OS_months, OS_event) ~ gprc5a_expr * subtype_moffitt,
                          data = analysis_df)
cox_results[["Interaction"]] <- broom::tidy(cox_interaction, exponentiate = TRUE, conf.int = TRUE)

# Compile results
cox_summary <- bind_rows(cox_results, .id = "Model") %>%
  select(Model, term, HR = estimate, CI_lower = conf.low, CI_upper = conf.high,
         p_value = p.value) %>%
  mutate(across(c(HR, CI_lower, CI_upper), ~round(.x, 3)),
         p_value   = format.pval(p_value, digits = 3),
         Significant = ifelse(as.numeric(p_value) < 0.05, "*", ""))

write.csv(cox_summary, file.path(TABLES_DIR, "aim1_cox_results.csv"), row.names = FALSE)

# ── 11. Figure 5: Forest Plot of Cox Results ─────────────────────────────────
message("Generating Figure 5: Forest plot...")

forest_data <- data.frame(
  Model   = c("Overall (all subtypes)", "Classical subtype", "Basal-like subtype"),
  HR      = c(coef(summary(cox_overall))[,"exp(coef)"],
              coef(summary(cox_classical))[,"exp(coef)"],
              coef(summary(cox_basal))[,"exp(coef)"]),
  lower   = c(exp(confint(cox_overall)),
              exp(confint(cox_classical)),
              exp(confint(cox_basal)))[c(1,3,5)],
  upper   = c(exp(confint(cox_overall)),
              exp(confint(cox_classical)),
              exp(confint(cox_basal)))[c(2,4,6)],
  color   = c("Combined", "Classical", "Basal-like")
)

p_forest <- ggplot(forest_data, aes(x = HR, y = reorder(Model, HR),
                                     xmin = lower, xmax = upper,
                                     color = color)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.8) +
  geom_errorbarh(height = 0.2, linewidth = 1.2) +
  geom_point(size = 4) +
  scale_color_manual(values = c("Combined"  = "#7F8C8D",
                                 "Classical"  = "#2E75B6",
                                 "Basal-like" = "#C0392B")) +
  scale_x_log10() +
  labs(
    title    = "GPRC5A Hazard Ratios by PDAC Molecular Subtype",
    subtitle = "Cox proportional hazards | TCGA-PAAD | GPRC5A per unit increase",
    x        = "Hazard Ratio (log scale)",
    y        = NULL,
    color    = "Subtype",
    caption  = "HR > 1: higher GPRC5A = worse prognosis (oncogenic direction)\nHR < 1: higher GPRC5A = better prognosis (paradoxical direction)"
  ) +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")

ggsave(file.path(FIGURES_DIR, "aim1_forest_plot.pdf"),
       plot = p_forest, width = 9, height = 5, device = "pdf")

# ── 12. Summary Statistics Table ─────────────────────────────────────────────
summary_stats <- analysis_df %>%
  group_by(subtype_moffitt, survival_group) %>%
  summarise(
    n             = n(),
    median_gprc5a = round(median(gprc5a_expr), 3),
    mean_gprc5a   = round(mean(gprc5a_expr),   3),
    sd_gprc5a     = round(sd(gprc5a_expr),     3),
    median_OS_mo  = round(median(OS_months),   1),
    .groups       = "drop"
  )
write.csv(summary_stats,
          file.path(TABLES_DIR, "aim1_gprc5a_summary_by_subtype.csv"),
          row.names = FALSE)

# ── 13. Session Summary ───────────────────────────────────────────────────────
message("\n=== Aim 1 Complete ===")
message(sprintf("Samples analyzed: n=%d", nrow(analysis_df)))
message(sprintf("Classical: n=%d | Basal-like: n=%d",
                sum(analysis_df$subtype_moffitt == "Classical"),
                sum(analysis_df$subtype_moffitt == "Basal-like")))
message(sprintf("Overall GPRC5A HR: %.3f (p=%s)",
                coef(summary(cox_overall))[,"exp(coef)"],
                format.pval(coef(summary(cox_overall))[,"Pr(>|z|)"], digits=3)))
message(sprintf("Classical GPRC5A HR: %.3f (p=%s)",
                coef(summary(cox_classical))[,"exp(coef)"],
                format.pval(coef(summary(cox_classical))[,"Pr(>|z|)"], digits=3)))
message(sprintf("Basal-like GPRC5A HR: %.3f (p=%s)",
                coef(summary(cox_basal))[,"exp(coef)"],
                format.pval(coef(summary(cox_basal))[,"Pr(>|z|)"], digits=3)))

message("\nOutputs saved:")
message(paste0("  ", FIGURES_DIR, "/aim1_*.pdf"))
message(paste0("  ", TABLES_DIR,  "/aim1_*.csv"))
