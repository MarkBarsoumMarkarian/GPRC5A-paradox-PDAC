# =============================================================================
# GPRC5A Paradox in PDAC — Aim 3: Protein-Level Validation via CPTAC-PAAD
# Project: Decoding the GPRC5A Paradox in Pancreatic Ductal Adenocarcinoma
# Author:  Mark Barsoum Markarian
# Date:    2026
#
# Description:
#   Uses the CPTAC-PAAD dataset (matched RNA + proteomics + phosphoproteomics)
#   to ask whether GPRC5A's RNA paradox is reflected at the protein level.
#
#   Three core questions:
#     Q1: Does GPRC5A protein level correlate with its RNA level?
#         A low correlation = post-transcriptional regulation is active,
#         meaning the RNA signal alone (from our original study) is incomplete.
#     Q2: Does GPRC5A protein level associate with survival differently
#         than RNA — i.e., does protein resolve the paradox?
#     Q3: Are there GPRC5A phosphorylation sites that differ between
#         patient groups, pointing to specific functional states?
#
#   Data source:
#     CPTAC-PAAD (Cao et al. Cancer Cell 2021) via the cptac R package.
#     Dataset: ~140 PDAC patients with matched tumor RNA-seq, protein LC-MS/MS,
#     and phosphoproteomics. Survival annotations available.
#
# Prerequisite: aim1_subtype_stratification.R must have been run first
#
# Outputs:
#   - figures/aim3_rna_protein_correlation.pdf
#   - figures/aim3_protein_survival_boxplot.pdf
#   - figures/aim3_km_protein_level.pdf
#   - figures/aim3_phosphosite_heatmap.pdf
#   - figures/aim3_rna_protein_scatter_by_subtype.pdf
#   - tables/aim3_rna_protein_correlation.csv
#   - tables/aim3_protein_cox_results.csv
#   - tables/aim3_phosphosite_differential.csv
# =============================================================================

# ── 0. Libraries ──────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(broom)
  library(pheatmap)
  library(RColorBrewer)
  library(ggpubr)
  library(patchwork)
  library(stringr)
  library(readr)
  library(limma)
})

# ── 1. Configuration ──────────────────────────────────────────────────────────
set.seed(42)
RESULTS_DIR <- here::here("results")
FIGURES_DIR <- file.path(RESULTS_DIR, "figures")
TABLES_DIR  <- file.path(RESULTS_DIR, "tables")
DATA_DIR    <- here::here("data")
CPTAC_DIR   <- file.path(DATA_DIR, "cptac")

dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES_DIR,  showWarnings = FALSE)
dir.create(CPTAC_DIR,   showWarnings = FALSE)

message("=== Aim 3: CPTAC Protein-Level Validation ===")

# ── 2. Load Aim 1 Outputs ─────────────────────────────────────────────────────
subtype_cache <- file.path(TABLES_DIR, "aim1_subtype_assignments.csv")
if (!file.exists(subtype_cache)) {
  stop("Run aim1_subtype_stratification.R first.")
}
subtype_df <- read.csv(subtype_cache, stringsAsFactors = FALSE)
message(sprintf("Loaded subtype assignments: n=%d", nrow(subtype_df)))

# ── 3. CPTAC-PDAC Data Acquisition ────────────────────────────────────────────
# Source: LinkedOmics CPTAC-PDAC (Cao et al. Cancer Cell 2021)
# URL: https://linkedomics.org/data_download/CPTAC-PDAC/
# All files are tab-separated .cct matrices: rows = genes, cols = samples
# n=140 tumors with matched RNA, proteomics, phosphoproteomics, and clinical data

cptac_rna_cache     <- file.path(CPTAC_DIR, "cptac_rna.rds")
cptac_prot_cache    <- file.path(CPTAC_DIR, "cptac_protein.rds")
cptac_phospho_cache <- file.path(CPTAC_DIR, "cptac_phospho.rds")
cptac_clin_cache    <- file.path(CPTAC_DIR, "cptac_clinical.rds")

# Correct LinkedOmics CPTAC-PDAC file URLs (verified February 2026)
LINKEDOMICS_BASE <- "https://linkedomics.org/data_download/CPTAC-PDAC"
CPTAC_FILES <- list(
  rna     = "mRNA_RSEM_UQ_log2_Tumor.cct",
  prot    = "proteomics_gene_level_MD_abundance_tumor.cct",
  phospho = "phosphoproteomics_site_level_MD_abundance_tumor.cct",
  clin    = "clinical_table_140.tsv"
)

# Helper: download one file with retries
download_cptac_file <- function(filename, dest, base = LINKEDOMICS_BASE) {
  url <- paste0(base, "/", filename)
  message(sprintf("  Downloading: %s", filename))
  result <- tryCatch({
    download.file(url, dest, mode = "wb", quiet = TRUE, method = "libcurl")
    TRUE
  }, warning = function(w) {
    message(sprintf("    Warning: %s", w$message))
    file.exists(dest) && file.size(dest) > 1000
  }, error = function(e) {
    message(sprintf("    Failed: %s", e$message))
    FALSE
  })
  result
}

# Helper: parse a .cct file into a genes × samples matrix
parse_cct <- function(path) {
  if (!file.exists(path) || file.size(path) < 100) return(NULL)
  tryCatch({
    df  <- read.table(path, sep = "\t", header = TRUE,
                      row.names = 1, check.names = FALSE,
                      comment.char = "", quote = "")
    mat <- as.matrix(df)
    storage.mode(mat) <- "numeric"
    message(sprintf("    Parsed: %d features × %d samples", nrow(mat), ncol(mat)))
    mat
  }, error = function(e) {
    message(sprintf("    Parse failed: %s", e$message))
    NULL
  })
}

if (!all(file.exists(c(cptac_rna_cache, cptac_prot_cache, cptac_clin_cache)))) {
  message("Downloading CPTAC-PDAC data from LinkedOmics...")
  
  # Download all files
  raw_paths <- lapply(names(CPTAC_FILES), function(nm) {
    ext  <- tools::file_ext(CPTAC_FILES[[nm]])
    dest <- file.path(CPTAC_DIR, paste0("cptac_raw_", nm, ".", ext))
    ok   <- download_cptac_file(CPTAC_FILES[[nm]], dest)
    if (ok) dest else NULL
  })
  names(raw_paths) <- names(CPTAC_FILES)
  
  # Parse and cache
  if (!is.null(raw_paths$rna)) {
    rna_mat <- parse_cct(raw_paths$rna)
    if (!is.null(rna_mat)) saveRDS(rna_mat, cptac_rna_cache)
  }
  if (!is.null(raw_paths$prot)) {
    prot_mat <- parse_cct(raw_paths$prot)
    if (!is.null(prot_mat)) saveRDS(prot_mat, cptac_prot_cache)
  }
  if (!is.null(raw_paths$phospho)) {
    phospho_mat <- parse_cct(raw_paths$phospho)
    if (!is.null(phospho_mat)) saveRDS(phospho_mat, cptac_phospho_cache)
  }
  if (!is.null(raw_paths$clin)) {
    clin_df_raw <- tryCatch(
      read.table(raw_paths$clin, sep = "\t", header = TRUE,
                 check.names = FALSE, quote = ""),
      error = function(e) NULL
    )
    if (!is.null(clin_df_raw)) saveRDS(clin_df_raw, cptac_clin_cache)
  }
  
  message("Downloads complete.")
  
} else {
  message("Loading cached CPTAC-PDAC data...")
}

# ── 4. Load and Validate CPTAC Data ───────────────────────────────────────────
load_cptac_matrix <- function(cache_path, name) {
  if (!file.exists(cache_path)) {
    message(sprintf("  %s not available — skipping", name))
    return(NULL)
  }
  mat <- readRDS(cache_path)
  
  # Standardize: ensure samples are columns, genes are rows
  if (is.data.frame(mat)) mat <- as.matrix(mat)
  
  # If samples appear to be rows (more rows than cols and rownames look like IDs)
  if (nrow(mat) > ncol(mat) && !grepl("^[A-Z]", rownames(mat)[1])) {
    mat <- t(mat)
  }
  
  message(sprintf("  %s: %d features × %d samples", name, nrow(mat), ncol(mat)))
  mat
}

cptac_rna    <- load_cptac_matrix(cptac_rna_cache,    "RNA")
cptac_prot   <- load_cptac_matrix(cptac_prot_cache,   "Protein")
cptac_phospho<- load_cptac_matrix(cptac_phospho_cache,"Phospho")
cptac_clin   <- if (file.exists(cptac_clin_cache)) readRDS(cptac_clin_cache) else NULL

if (is.null(cptac_prot)) {
  stop(paste(
    "CPTAC proteomics data could not be loaded.",
    "Please download manually from:",
    "https://cptac-data-portal.georgetown.edu/study-summary/S044",
    "or https://linkedomics.org/data_download/CPTAC-PAAD/",
    "Save protein abundance matrix to:", cptac_prot_cache,
    sep = "\n"
  ))
}

# ── 5. Align Samples Across Omics Layers ──────────────────────────────────────
message("Aligning samples across omics layers...")

# Standardize sample IDs — CPTAC uses C3L/C3N prefixes
clean_sample_ids <- function(ids) {
  # Remove trailing suffixes like _T (tumor), _N (normal), _NAT
  gsub("_T$|_N$|_NAT$|\\.T$|\\.N$", "", ids)
}

if (!is.null(cptac_rna))    colnames(cptac_rna)    <- clean_sample_ids(colnames(cptac_rna))
if (!is.null(cptac_prot))   colnames(cptac_prot)   <- clean_sample_ids(colnames(cptac_prot))
if (!is.null(cptac_phospho))colnames(cptac_phospho)<- clean_sample_ids(colnames(cptac_phospho))

# Find samples with both RNA and protein data
if (!is.null(cptac_rna) && !is.null(cptac_prot)) {
  shared_samples <- intersect(colnames(cptac_rna), colnames(cptac_prot))
  message(sprintf("  Samples with matched RNA + protein: n=%d", length(shared_samples)))
  
  cptac_rna_aligned  <- cptac_rna[,  shared_samples]
  cptac_prot_aligned <- cptac_prot[, shared_samples]
} else {
  shared_samples     <- colnames(cptac_prot)
  cptac_rna_aligned  <- NULL
  cptac_prot_aligned <- cptac_prot
}

# Add phospho if available
if (!is.null(cptac_phospho)) {
  phospho_samples <- intersect(shared_samples, colnames(cptac_phospho))
  cptac_phospho_aligned <- cptac_phospho[, phospho_samples]
  message(sprintf("  Samples with phosphoproteomics: n=%d", length(phospho_samples)))
}

# ── 6. Clinical Survival Data for CPTAC ───────────────────────────────────────
build_cptac_clinical <- function(clin_obj, sample_ids) {
  if (is.null(clin_obj)) {
    message("  No CPTAC clinical data — generating survival-unknown placeholder")
    return(data.frame(
      sample_id = sample_ids, OS_months = NA_real_,
      OS_event = NA_integer_, survival_group = NA_character_,
      stringsAsFactors = FALSE
    ))
  }
  
  df   <- as.data.frame(clin_obj)
  cols <- colnames(df)
  message(sprintf("  Clinical columns available: %s",
                  paste(head(cols, 15), collapse = ", ")))
  
  # CPTAC-PDAC clinical_table_140.tsv known column names
  # (from Cao et al. 2021 supplementary)
  id_candidates      <- c("case_id", "CASE_ID", "sample_id", "Patient_ID",
                           "Proteomics_Participant_ID", "id")
  vital_candidates   <- c("Vital.Status", "vital_status", "Overall_Survival_Status",
                           "OS_Status", "status")
  os_day_candidates  <- c("OS", "Overall_Survival", "Days.to.death",
                           "days_to_death", "OS_days", "Survival_Days")
  lastfu_candidates  <- c("Days.to.Last.Follow.up", "days_to_last_follow_up",
                           "Last_Contact_Days")
  
  find_col_clin <- function(candidates, cols) {
    hit <- candidates[candidates %in% cols]
    if (length(hit) > 0) return(hit[1])
    # Case-insensitive fallback
    hit2 <- cols[tolower(cols) %in% tolower(candidates)]
    if (length(hit2) > 0) return(hit2[1])
    NA_character_
  }
  
  col_id     <- find_col_clin(id_candidates,    cols)
  col_vital  <- find_col_clin(vital_candidates,  cols)
  col_os     <- find_col_clin(os_day_candidates, cols)
  col_lastfu <- find_col_clin(lastfu_candidates, cols)
  
  message(sprintf("  ID col: %s | Vital: %s | OS days: %s | Last FU: %s",
                  col_id, col_vital, col_os, col_lastfu))
  
  # Build sample IDs — use rownames if no ID column found
  ids <- if (!is.na(col_id)) df[[col_id]] else rownames(df)
  
  # Build survival variables
  days_os    <- if (!is.na(col_os))     suppressWarnings(as.numeric(df[[col_os]]))    else NA_real_
  days_lastfu<- if (!is.na(col_lastfu)) suppressWarnings(as.numeric(df[[col_lastfu]])) else NA_real_
  vital_raw  <- if (!is.na(col_vital))  df[[col_vital]] else rep(NA_character_, nrow(df))
  
  result <- data.frame(
    sample_id      = ids,
    OS_time        = ifelse(!is.na(days_os) & days_os > 0, days_os, days_lastfu),
    vital_raw      = as.character(vital_raw),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      OS_event       = ifelse(grepl("dead|deceased|yes|1", tolower(vital_raw)), 1L, 0L),
      OS_months      = OS_time / 30.44,
      survival_group = ifelse(OS_event == 1, "Deceased", "Alive")
    ) %>%
    dplyr::select(sample_id, OS_months, OS_event, survival_group)
  
  # Filter to samples in our analysis set
  result <- result[result$sample_id %in% sample_ids, ]
  message(sprintf("  Clinical matched: n=%d | Alive=%d | Deceased=%d | OS unknown=%d",
                  nrow(result),
                  sum(result$survival_group == "Alive",   na.rm = TRUE),
                  sum(result$survival_group == "Deceased", na.rm = TRUE),
                  sum(is.na(result$OS_months))))
  result
}

cptac_clin_df <- build_cptac_clinical(cptac_clin, shared_samples)

# ── 7. Extract GPRC5A Values Across Omics Layers ─────────────────────────────
message("Extracting GPRC5A across omics layers...")

get_gene <- function(mat, gene) {
  if (is.null(mat)) return(NULL)
  # Try exact match first
  if (gene %in% rownames(mat)) return(as.numeric(mat[gene, ]))
  # Try case-insensitive
  hit <- rownames(mat)[tolower(rownames(mat)) == tolower(gene)]
  if (length(hit) > 0) return(as.numeric(mat[hit[1], ]))
  message(sprintf("  %s not found in matrix", gene))
  NULL
}

gprc5a_rna_cptac   <- get_gene(cptac_rna_aligned,  "GPRC5A")
gprc5a_prot_cptac  <- get_gene(cptac_prot_aligned, "GPRC5A")

# Phospho: collect all GPRC5A phosphosites (rownames contain gene name)
gprc5a_phospho_df <- NULL
if (!is.null(cptac_phospho)) {
  phospho_rows <- rownames(cptac_phospho_aligned)[
    grepl("^GPRC5A", rownames(cptac_phospho_aligned), ignore.case = TRUE)
  ]
  message(sprintf("  GPRC5A phosphosites found: %d", length(phospho_rows)))
  if (length(phospho_rows) > 0) {
    gprc5a_phospho_df <- as.data.frame(
      t(cptac_phospho_aligned[phospho_rows, , drop = FALSE])
    )
    gprc5a_phospho_df$sample_id <- rownames(gprc5a_phospho_df)
  }
}

# ── 8. Build CPTAC Analysis Dataframe ─────────────────────────────────────────
cptac_analysis <- data.frame(
  sample_id       = shared_samples,
  gprc5a_protein  = gprc5a_prot_cptac,
  stringsAsFactors = FALSE
)

if (!is.null(gprc5a_rna_cptac)) {
  cptac_analysis$gprc5a_rna <- gprc5a_rna_cptac
}

# Add survival
cptac_analysis <- dplyr::left_join(
  cptac_analysis, cptac_clin_df, by = "sample_id"
) %>%
  dplyr::mutate(
    protein_high = ifelse(gprc5a_protein >= median(gprc5a_protein, na.rm = TRUE),
                          "Protein High", "Protein Low"),
    rna_high     = if ("gprc5a_rna" %in% names(.))
                     ifelse(gprc5a_rna >= median(gprc5a_rna, na.rm = TRUE),
                            "RNA High", "RNA Low")
                   else NA_character_
  )

message(sprintf("CPTAC analysis cohort: n=%d", nrow(cptac_analysis)))

# ── 9. Q1: RNA-Protein Correlation ────────────────────────────────────────────
message("Q1: Computing GPRC5A RNA-protein correlation...")

if (!is.null(gprc5a_rna_cptac)) {
  cor_test <- cor.test(cptac_analysis$gprc5a_rna,
                       cptac_analysis$gprc5a_protein,
                       method = "spearman", use = "complete.obs")
  
  r_val <- round(cor_test$estimate, 3)
  p_val <- format.pval(cor_test$p.value, digits = 3)
  
  message(sprintf("  GPRC5A RNA-Protein Spearman r = %s (p = %s)", r_val, p_val))
  
  # Interpretation
  if (abs(r_val) < 0.3) {
    message("  => LOW correlation: strong post-transcriptional regulation likely")
  } else if (abs(r_val) < 0.6) {
    message("  => MODERATE correlation: partial post-transcriptional regulation")
  } else {
    message("  => HIGH correlation: RNA is a good proxy for protein level")
  }
  
  # Save
  write.csv(
    data.frame(Gene = "GPRC5A", Spearman_r = r_val, p_value = p_val,
               n = sum(!is.na(cptac_analysis$gprc5a_rna) &
                         !is.na(cptac_analysis$gprc5a_protein)),
               Interpretation = ifelse(abs(r_val) < 0.3, "Low — post-transcriptional regulation likely",
                                ifelse(abs(r_val) < 0.6, "Moderate", "High — RNA proxies protein well"))),
    file.path(TABLES_DIR, "aim3_rna_protein_correlation.csv"), row.names = FALSE
  )
  
  # Figure 1: RNA vs Protein scatter
  p_corr <- ggplot(cptac_analysis, aes(x = gprc5a_rna, y = gprc5a_protein)) +
    geom_point(aes(color = survival_group), alpha = 0.7, size = 2.5) +
    geom_smooth(method = "lm", se = TRUE, color = "#2C3E50", linewidth = 1) +
    scale_color_manual(values = c("Alive" = "#27AE60", "Deceased" = "#E74C3C"),
                       na.value = "grey60") +
    annotate("text", x = -Inf, y = Inf,
             label = sprintf("Spearman r = %s\np = %s", r_val, p_val),
             hjust = -0.1, vjust = 1.3, size = 4.5, color = "#2C3E50",
             fontface = "bold") +
    labs(
      title    = "GPRC5A RNA vs. Protein Abundance — CPTAC-PAAD",
      subtitle = sprintf("Spearman r = %s | n = %d matched samples",
                         r_val, nrow(cptac_analysis)),
      x        = "GPRC5A mRNA Expression (CPTAC RNA-seq)",
      y        = "GPRC5A Protein Abundance (LC-MS/MS)",
      color    = "Survival"
    ) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))
  
  # Figure 1b: Scatter faceted by subtype (if subtype info available)
  # Map TCGA subtypes to CPTAC if CPTAC has its own subtype calls
  p_corr_save <- p_corr
  
  ggsave(file.path(FIGURES_DIR, "aim3_rna_protein_correlation.pdf"),
         plot = p_corr_save, width = 8, height = 6, device = "pdf")
  
} else {
  message("  RNA data not available — skipping RNA-protein correlation")
}

# ── 10. Genome-wide RNA-Protein Correlation (context) ─────────────────────────
# How does GPRC5A's RNA-protein correlation compare to all other genes?
# Places GPRC5A in the global landscape of post-transcriptional regulation

if (!is.null(cptac_rna_aligned) && !is.null(cptac_prot_aligned)) {
  message("Computing genome-wide RNA-protein correlations for context...")
  
  common_genes <- intersect(rownames(cptac_rna_aligned),
                             rownames(cptac_prot_aligned))
  message(sprintf("  Genes with matched RNA + protein: %d", length(common_genes)))
  
  # Compute Spearman r for each gene
  all_cors <- sapply(common_genes, function(g) {
    r <- tryCatch(
      cor(as.numeric(cptac_rna_aligned[g, shared_samples]),
          as.numeric(cptac_prot_aligned[g, shared_samples]),
          method = "spearman", use = "complete.obs"),
      error = function(e) NA_real_
    )
    r
  })
  
  cor_df <- data.frame(
    gene      = names(all_cors),
    r_rna_prot = all_cors,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(!is.na(r_rna_prot)) %>%
    dplyr::arrange(desc(r_rna_prot))
  
  gprc5a_r    <- cor_df$r_rna_prot[cor_df$gene == "GPRC5A"]
  gprc5a_pct  <- if (length(gprc5a_r) > 0)
    round(100 * mean(cor_df$r_rna_prot <= gprc5a_r, na.rm = TRUE), 1)
  else NA
  
  message(sprintf("  GPRC5A percentile in genome-wide distribution: %.1f%%", gprc5a_pct))
  
  # Figure 2: Genome-wide distribution with GPRC5A highlighted
  p_genome <- ggplot(cor_df, aes(x = r_rna_prot)) +
    geom_histogram(bins = 80, fill = "steelblue", alpha = 0.7, color = NA) +
    geom_vline(xintercept = ifelse(length(gprc5a_r) > 0, gprc5a_r, NA),
               color = "#E74C3C", linewidth = 1.2, linetype = "dashed") +
    annotate("text",
             x = ifelse(length(gprc5a_r) > 0, gprc5a_r + 0.02, 0),
             y = Inf,
             label = sprintf("GPRC5A\nr=%.2f\n(%s%%ile)",
                             ifelse(length(gprc5a_r) > 0, gprc5a_r, NA),
                             gprc5a_pct),
             hjust = 0, vjust = 1.5, color = "#E74C3C",
             size = 3.8, fontface = "bold") +
    labs(
      title    = "Genome-wide RNA-Protein Correlation — CPTAC-PAAD",
      subtitle = "GPRC5A highlighted: where does it fall in the global landscape?",
      x        = "Spearman r (RNA vs. Protein)",
      y        = "Number of genes",
      caption  = sprintf("n = %d genes with matched RNA + protein data", nrow(cor_df))
    ) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(file.path(FIGURES_DIR, "aim3_genome_wide_rna_protein_cor.pdf"),
         plot = p_genome, width = 9, height = 5, device = "pdf")
}

# ── 11. Q2: GPRC5A Protein Level and Survival ────────────────────────────────
message("Q2: GPRC5A protein level vs. survival...")

if (any(!is.na(cptac_analysis$OS_months))) {
  
  cptac_surv <- cptac_analysis %>%
    dplyr::filter(!is.na(OS_months), !is.na(OS_event), !is.na(gprc5a_protein))
  
  # Boxplot: protein by survival group
  p_box <- ggplot(cptac_surv, aes(x = survival_group, y = gprc5a_protein,
                                   fill = survival_group)) +
    geom_boxplot(alpha = 0.8, outlier.shape = 21, width = 0.5) +
    geom_jitter(aes(color = survival_group), width = 0.12, alpha = 0.5, size = 2) +
    scale_fill_manual(values  = c("Alive" = "#27AE60", "Deceased" = "#E74C3C")) +
    scale_color_manual(values = c("Alive" = "#27AE60", "Deceased" = "#E74C3C")) +
    stat_compare_means(method = "wilcox.test", label = "p.format",
                       label.x = 1.3, size = 4) +
    labs(
      title    = "GPRC5A Protein Abundance by Survival Status — CPTAC-PAAD",
      subtitle = "Key test: does protein-level direction match or resolve the RNA paradox?",
      x        = "Survival Status", y = "GPRC5A Protein (log2 abundance)",
      fill = "Survival", color = "Survival"
    ) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold"), legend.position = "none")
  
  ggsave(file.path(FIGURES_DIR, "aim3_protein_survival_boxplot.pdf"),
         plot = p_box, width = 7, height = 5, device = "pdf")
  
  # KM curve by protein level
  fit_prot <- survfit(Surv(OS_months, OS_event) ~ protein_high, data = cptac_surv)
  
  km_prot <- tryCatch(
    ggsurvplot(
      fit_prot,
      data              = cptac_surv,
      pval              = TRUE,
      pval.method       = TRUE,
      conf.int          = TRUE,
      risk.table        = TRUE,
      risk.table.height = 0.28,
      palette           = c("#E74C3C", "#2E75B6"),
      title             = "Survival by GPRC5A Protein Level — CPTAC-PAAD",
      subtitle          = "Compare direction to TCGA RNA-based KM (Aim 1)",
      xlab              = "Time (months)",
      ylab              = "Overall Survival",
      legend.labs       = c("Protein High", "Protein Low"),
      legend.title      = "GPRC5A Protein",
      ggtheme           = theme_bw(base_size = 12),
      surv.median.line  = "hv"
    ),
    error = function(e) {
      message("  KM plot failed: ", e$message); NULL
    }
  )
  
  if (!is.null(km_prot)) {
    pdf(file.path(FIGURES_DIR, "aim3_km_protein_level.pdf"), width = 10, height = 8)
      print(km_prot)
    dev.off()
  }
  
  # Cox model: protein level
  cox_prot <- coxph(Surv(OS_months, OS_event) ~ gprc5a_protein, data = cptac_surv)
  cox_prot_sum <- broom::tidy(cox_prot, exponentiate = TRUE, conf.int = TRUE)
  
  write.csv(
    data.frame(
      Model    = "GPRC5A_Protein_Univariate",
      term     = cox_prot_sum$term,
      HR       = round(cox_prot_sum$estimate,  3),
      CI_lower = round(cox_prot_sum$conf.low,  3),
      CI_upper = round(cox_prot_sum$conf.high, 3),
      p_value  = format.pval(cox_prot_sum$p.value, digits = 3),
      n        = nrow(cptac_surv),
      stringsAsFactors = FALSE
    ),
    file.path(TABLES_DIR, "aim3_protein_cox_results.csv"), row.names = FALSE
  )
  
  message(sprintf("  Protein Cox HR = %.3f (p = %s)",
                  cox_prot_sum$estimate,
                  format.pval(cox_prot_sum$p.value, digits = 3)))
  message("  Interpretation:")
  message("    RNA HR < 1 (paradoxical) + Protein HR > 1 => RNA paradox, protein tells truth")
  message("    Both HR < 1                               => True suppressive signal")
  message("    Both HR > 1                               => Consistent oncogenic signal")
  
} else {
  message("  No CPTAC survival data available — skipping survival analysis")
}

# ── 12. Q3: Phosphosite Differential Analysis ────────────────────────────────
if (!is.null(gprc5a_phospho_df) && any(!is.na(cptac_analysis$survival_group))) {
  message("Q3: GPRC5A phosphosite differential analysis...")
  
  phospho_surv <- gprc5a_phospho_df %>%
    dplyr::inner_join(
      cptac_analysis[, c("sample_id", "survival_group", "OS_event")],
      by = "sample_id"
    ) %>%
    dplyr::filter(!is.na(survival_group))
  
  phospho_sites <- setdiff(colnames(phospho_surv),
                            c("sample_id", "survival_group", "OS_event"))
  
  if (length(phospho_sites) > 0) {
    message(sprintf("  Testing %d GPRC5A phosphosites", length(phospho_sites)))
    
    # Wilcoxon test per phosphosite: Alive vs Deceased
    phospho_results <- lapply(phospho_sites, function(site) {
      vals_alive <- phospho_surv[[site]][phospho_surv$survival_group == "Alive"]
      vals_dead  <- phospho_surv[[site]][phospho_surv$survival_group == "Deceased"]
      
      vals_alive <- vals_alive[!is.na(vals_alive)]
      vals_dead  <- vals_dead[!is.na(vals_dead)]
      
      if (length(vals_alive) < 3 || length(vals_dead) < 3) return(NULL)
      
      wt <- wilcox.test(vals_alive, vals_dead, exact = FALSE)
      
      data.frame(
        phosphosite   = site,
        mean_alive    = round(mean(vals_alive), 3),
        mean_deceased = round(mean(vals_dead),  3),
        log2FC        = round(mean(vals_dead) - mean(vals_alive), 3),
        p_value       = wt$p.value,
        stringsAsFactors = FALSE
      )
    })
    
    phospho_df <- dplyr::bind_rows(phospho_results) %>%
      dplyr::mutate(
        p_adj     = p.adjust(p_value, method = "BH"),
        sig       = p_adj < 0.05,
        direction = ifelse(log2FC > 0, "Higher in Deceased", "Higher in Alive")
      ) %>%
      dplyr::arrange(p_value)
    
    write.csv(phospho_df,
              file.path(TABLES_DIR, "aim3_phosphosite_differential.csv"),
              row.names = FALSE)
    
    # Heatmap of phosphosites
    if (nrow(phospho_df) >= 2) {
      heat_phospho <- t(as.matrix(
        phospho_surv[, phospho_df$phosphosite, drop = FALSE]
      ))
      heat_phospho_z <- t(scale(t(heat_phospho)))
      
      ann_col <- data.frame(
        Survival = phospho_surv$survival_group,
        row.names = phospho_surv$sample_id
      )
      
      pdf(file.path(FIGURES_DIR, "aim3_phosphosite_heatmap.pdf"),
          width = 12, height = max(4, nrow(heat_phospho_z) * 0.5 + 3))
      pheatmap(
        heat_phospho_z,
        annotation_col   = ann_col,
        annotation_colors= list(Survival = c("Alive" = "#27AE60", "Deceased" = "#E74C3C")),
        cluster_cols     = TRUE,
        show_colnames    = FALSE,
        color            = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
        breaks           = seq(-3, 3, length.out = 101),
        main             = "GPRC5A Phosphosite Abundances | CPTAC-PAAD",
        border_color     = NA,
        fontsize_row     = 9
      )
      dev.off()
    }
    
    message(sprintf("  Significant phosphosites (FDR<0.05): %d",
                    sum(phospho_df$sig, na.rm = TRUE)))
  }
} else {
  message("  Phosphoproteomics not available — skipping phosphosite analysis")
}

# ── 13. Summary comparison: TCGA RNA vs. CPTAC Protein ──────────────────────
message("\n=== Aim 3 Summary ===")
message("Cross-study comparison:")
message("  TCGA RNA signal: GPRC5A LOWER in deceased (paradoxical)")
message("  CPTAC Protein  : check aim3_protein_cox_results.csv for direction")

if (file.exists(file.path(TABLES_DIR, "aim3_protein_cox_results.csv"))) {
  prot_cox <- read.csv(file.path(TABLES_DIR, "aim3_protein_cox_results.csv"))
  hr <- prot_cox$HR[1]
  pv <- prot_cox$p_value[1]
  message(sprintf("  Protein HR = %.3f (p = %s)", hr, pv))
  if (hr > 1) {
    message("  => Protein confirms oncogenic direction — RNA was misleading")
  } else {
    message("  => Protein also shows protective signal — paradox is real")
  }
}

message("\nOutputs saved:")
message(paste0("  ", FIGURES_DIR, "/aim3_*.pdf"))
message(paste0("  ", TABLES_DIR,  "/aim3_*.csv"))
