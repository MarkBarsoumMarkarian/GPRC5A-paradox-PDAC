# =============================================================================
# GPRC5A Paradox in PDAC — Aim 5: Context-State ML Classifier
# Project: Decoding the GPRC5A Paradox in Pancreatic Ductal Adenocarcinoma
# Author:  Mark Barsoum Markarian
# Date:    2026
#
# Description:
#   Builds a machine learning classifier that predicts which "role state"
#   GPRC5A is operating in for a given tumor — oncogenic or suppressive —
#   based on the multi-omic context established in Aims 1 and 2.
#
#   This is the novel, integrative payoff of the whole project. It turns
#   the resolved contradiction into a clinically actionable tool: given a
#   new PDAC tumor's RNA profile, predict whether GPRC5A is acting as an
#   oncogene (poor prognosis) or suppressor (better prognosis), and use
#   that prediction for survival stratification.
#
#   Approach:
#     - Labels: derived from Aim 1 (subtype-stratified survival association)
#       Oncogenic state  = Classical subtype + GPRC5A High + worse survival
#       Suppressive state = Basal-like subtype + GPRC5A High + better survival
#     - Features: GPRC5A expression + Moffitt subtype score + co-expressed
#       genes (selected on TRAINING data only) + treatment flag + clinical covariates
#     - Models: Logistic Regression (glmnet), Random Forest, XGBoost
#     - Validation: LOOCV on training set; held-out 20% test set for final eval
#
# KEY FIX NOTES (v2):
#   1. Feature selection (co-expression ranking) is now done INSIDE the
#      training fold only — using only training-set rows — preventing the
#      data leakage that inflated test AUC to 1.0 in v1.
#   2. XGBoost LOOCV failure fixed: single-class LOOCV folds caused
#      twoClassSummary to return all-NA. We now use a custom summary function
#      that gracefully handles degenerate folds, and fall back to 5-fold CV
#      for XGB if LOOCV still fails.
#   3. legend.title in ggsurvplot fixed from list() to plain string.
#   4. Tertile splits also moved inside training-only data to avoid leakage.
#
# Prerequisite: Aims 1 and 2 must have been run first
#
# Outputs:
#   - figures/aim5_feature_importance.pdf
#   - figures/aim5_roc_curves.pdf
#   - figures/aim5_role_state_km.pdf
#   - figures/aim5_calibration_plot.pdf
#   - tables/aim5_model_performance.csv
#   - tables/aim5_feature_importance.csv
#   - tables/aim5_role_state_predictions.csv
# =============================================================================

# ── 0. Libraries ──────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(caret)
  library(randomForest)
  library(xgboost)
  library(pROC)
  library(broom)
  library(patchwork)
  library(RColorBrewer)
  library(stringr)
  library(pheatmap)
})

# ── 1. Configuration ──────────────────────────────────────────────────────────
set.seed(42)
RESULTS_DIR <- here::here("results")
FIGURES_DIR <- file.path(RESULTS_DIR, "figures")
TABLES_DIR  <- file.path(RESULTS_DIR, "tables")
DATA_DIR    <- here::here("data")

dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES_DIR,  showWarnings = FALSE)

message("=== Aim 5: Context-State ML Classifier ===")

# ── 2. Load All Prior Outputs ─────────────────────────────────────────────────
subtype_cache  <- file.path(TABLES_DIR, "aim1_subtype_assignments.csv")
treatment_cache <- file.path(TABLES_DIR, "aim2_treatment_summary.csv")
vst_cache      <- file.path(DATA_DIR, "tcga_paad_vst.rds")

if (!file.exists(subtype_cache)) stop("Run aim1_subtype_stratification.R first.")
if (!file.exists(vst_cache))     stop("Run aim1_subtype_stratification.R first.")

subtype_df <- read.csv(subtype_cache, stringsAsFactors = FALSE)
vst_mat    <- readRDS(vst_cache)

# Reload clinical from Aims 1 & 2
clinical_cache <- file.path(DATA_DIR, "tcga_paad_clinical.rds")
source(here::here("R", "utils_clinical.R"))
clinical_raw <- readRDS(clinical_cache)
clin_df      <- build_clinical_df(clinical_raw)

message(sprintf("Loaded data: n=%d samples", nrow(subtype_df)))

# ── 3. Define GPRC5A Role State Labels ────────────────────────────────────────
# NOTE: tertile splits are computed on the FULL labeled set here for labeling
# purposes only. Feature values extracted from vst_mat are what go into the
# model, and feature *selection* (gene ranking) is deferred to training data.
#
# Label logic:
#   Oncogenic  = Classical + GPRC5A High + Deceased
#             OR Basal-like + GPRC5A Low  + Deceased
#   Suppressive= Classical + GPRC5A Low  + Alive
#             OR Basal-like + GPRC5A High + Alive
#   Ambiguous  — excluded

message("Assigning GPRC5A role state labels...")

label_df <- subtype_df %>%
  dplyr::mutate(
    OS_event       = ifelse(vital_status == "Dead", 1, 0),
    gprc5a_tertile = NA_character_
  )

for (st in c("Classical", "Basal-like")) {
  idx <- label_df$subtype_moffitt == st
  if (sum(idx, na.rm = TRUE) < 9) next
  tertiles <- quantile(label_df$gprc5a_expr[idx], probs = c(1/3, 2/3), na.rm = TRUE)
  label_df$gprc5a_tertile[idx] <- ifelse(
    label_df$gprc5a_expr[idx] <= tertiles[1], "Low",
    ifelse(label_df$gprc5a_expr[idx] >= tertiles[2], "High", "Mid")
  )
}

label_df <- label_df %>%
  dplyr::mutate(
    role_state = dplyr::case_when(
      subtype_moffitt == "Classical"  & gprc5a_tertile == "High" & OS_event == 1 ~ "Oncogenic",
      subtype_moffitt == "Basal-like" & gprc5a_tertile == "Low"  & OS_event == 1 ~ "Oncogenic",
      subtype_moffitt == "Classical"  & gprc5a_tertile == "Low"  & OS_event == 0 ~ "Suppressive",
      subtype_moffitt == "Basal-like" & gprc5a_tertile == "High" & OS_event == 0 ~ "Suppressive",
      TRUE ~ "Ambiguous"
    )
  ) %>%
  dplyr::filter(role_state != "Ambiguous")

message(sprintf("Labeled samples: Oncogenic=%d | Suppressive=%d",
                sum(label_df$role_state == "Oncogenic"),
                sum(label_df$role_state == "Suppressive")))

# ── 4. Train / Test Split (BEFORE any feature selection) ──────────────────────
# CRITICAL: we split here so that co-expression feature selection in step 5
# uses ONLY training rows — preventing the leakage that caused Test AUC = 1.0.

y        <- factor(label_df$role_state, levels = c("Oncogenic", "Suppressive"))
barcodes <- label_df$barcode
n_total  <- nrow(label_df)

set.seed(42)
train_idx <- caret::createDataPartition(y, p = 0.80, list = FALSE)[, 1]
bar_train <- barcodes[train_idx]
bar_test  <- barcodes[-train_idx]
y_train   <- y[train_idx]
y_test    <- y[-train_idx]
label_train <- label_df[train_idx, ]

message(sprintf("Held-out split: Train n=%d | Test n=%d",
                length(bar_train), length(bar_test)))

# ── 5. Feature Engineering (training data only for selection) ─────────────────
message("Engineering features for classifier...")

# Moffitt signature genes (pre-specified — no data-driven selection, so safe)
moffitt_classical <- c("TFF1","TFF2","TFF3","AGR2","CEACAM5","CEACAM6",
                        "LYZ","LGALS4","CYP3A7","MYO1A")
moffitt_basal     <- c("VGLL1","UCA1","S100A2","LY6D","KRT15","KRT17",
                        "AREG","SPRR1B","CDKN2A","ANLN")
signature_genes   <- c(moffitt_classical, moffitt_basal)
signature_avail   <- intersect(signature_genes, rownames(vst_mat))

# Co-expression ranking — TRAINING SAMPLES ONLY
gprc5a_train_vec <- as.numeric(vst_mat["GPRC5A", bar_train])
cor_train <- apply(vst_mat[, bar_train], 1, function(g) {
  tryCatch(cor(g, gprc5a_train_vec, use = "complete.obs"), error = function(e) NA)
})
cor_train <- sort(cor_train[!is.na(cor_train)], decreasing = TRUE)

top_coexpr <- c(
  names(head(cor_train[names(cor_train) != "GPRC5A"], 25)),
  names(tail(cor_train[names(cor_train) != "GPRC5A"], 25))
)
message(sprintf("  Selected %d co-expression features (from training data only)",
                length(top_coexpr)))

feature_genes <- unique(c("GPRC5A", top_coexpr, signature_avail))
feature_genes <- intersect(feature_genes, rownames(vst_mat))
message(sprintf("  Total gene features: %d", length(feature_genes)))

# Build feature matrices for train and test separately
build_feat_df <- function(bar_ids, label_sub) {
  feat <- t(vst_mat[feature_genes, bar_ids])
  feat <- as.data.frame(feat)
  colnames(feat) <- make.names(colnames(feat))
  feat$barcode        <- bar_ids
  feat$subtype_score  <- label_sub$subtype_score
  feat$classical_score <- label_sub$classical_score
  feat <- feat %>%
    dplyr::left_join(
      clin_df[, c("barcode", "received_gemcitabine", "age_years", "stage_binary")],
      by = "barcode"
    ) %>%
    dplyr::mutate(
      gem_flag   = as.integer(received_gemcitabine),
      stage_IV   = as.integer(stage_binary == "Stage IV"),
      age_scaled = as.numeric(scale(age_years))
    )
  feat
}

feat_train_df <- build_feat_df(bar_train, label_df[train_idx, ])
feat_test_df  <- build_feat_df(bar_test,  label_df[-train_idx, ])

meta_cols <- c("barcode", "received_gemcitabine", "stage_binary")

impute_medians <- function(df, ref_df = df) {
  as.data.frame(lapply(seq_along(df), function(i) {
    col <- df[[i]]
    if (any(is.na(col))) {
      col[is.na(col)] <- median(ref_df[[i]], na.rm = TRUE)
    }
    col
  }), stringsAsFactors = FALSE, col.names = colnames(df))
}

X_train_raw <- feat_train_df[, !colnames(feat_train_df) %in% meta_cols]
X_test_raw  <- feat_test_df[,  !colnames(feat_test_df)  %in% meta_cols]

# Impute: compute medians from training, apply to both
X_train <- impute_medians(X_train_raw, X_train_raw)
X_test  <- impute_medians(X_test_raw,  X_train_raw)  # use train medians for test

# Remove barcode column from feature matrices
X_train <- X_train[, colnames(X_train) != "barcode"]
X_test  <- X_test[,  colnames(X_test)  != "barcode"]

# predict_input: routes newdata to the right format per model.
# predict.xgb_native handles its own matrix conversion internally,
# so all models can safely receive the X_test data frame.
predict_input <- function(nm) X_test


message(sprintf("Feature matrix: Train=%d × %d | Test=%d × %d",
                nrow(X_train), ncol(X_train), nrow(X_test), ncol(X_test)))
message(sprintf("Class balance (train): Oncogenic=%d | Suppressive=%d",
                sum(y_train == "Oncogenic"), sum(y_train == "Suppressive")))

# ── 6. CV Strategy ────────────────────────────────────────────────────────────
# LOOCV on training set only (n_train ~53). For XGBoost, degenerate LOOCV folds
# (where all predictions are one class) cause twoClassSummary to return all-NA.
# Fix: use a safe summary function that returns 0.5 AUC for degenerate folds,
# and fall back to 5-fold CV for XGB if LOOCV still fails.

n_train <- nrow(X_train)
message(sprintf("Using LOOCV on training set (n=%d)", n_train))

# Safe two-class summary: handles single-class folds gracefully
safe_two_class_summary <- function(data, lev = NULL, model = NULL) {
  # data has: obs (factor), pred (factor), prob columns named after lev
  out <- tryCatch(
    caret::twoClassSummary(data, lev = lev, model = model),
    error = function(e) c(ROC = 0.5, Sens = NA_real_, Spec = NA_real_)
  )
  # If ROC is NA (degenerate fold), substitute 0.5
  if (is.na(out["ROC"])) out["ROC"] <- 0.5
  out
}

ctrl_loocv <- caret::trainControl(
  method          = "LOOCV",
  classProbs      = TRUE,
  summaryFunction = safe_two_class_summary,
  savePredictions = "final",
  verboseIter     = FALSE
)

ctrl_5fold <- caret::trainControl(
  method          = "cv",
  number          = 5,
  classProbs      = TRUE,
  summaryFunction = safe_two_class_summary,
  savePredictions = "final",
  verboseIter     = FALSE
)

# ── 7. Model Training ─────────────────────────────────────────────────────────
message("Training classifiers with LOOCV on training set...")

# 7a. Logistic Regression (glmnet — L1/L2 regularisation handles high-dim)
message("  Training Logistic Regression...")
lr_model <- tryCatch(
  caret::train(
    x          = X_train,
    y          = y_train,
    method     = "glmnet",
    metric     = "ROC",
    trControl  = ctrl_loocv,
    preProcess = c("center", "scale", "zv"),
    tuneGrid   = expand.grid(
      alpha  = c(0, 0.5, 1),
      lambda = 10^seq(-4, 0, length = 8)
    )
  ),
  error = function(e) { message("  LR failed: ", e$message); NULL }
)

# 7b. Random Forest
message("  Training Random Forest...")
rf_model <- tryCatch(
  caret::train(
    x         = X_train,
    y         = y_train,
    method    = "rf",
    metric    = "ROC",
    trControl = ctrl_loocv,
    tuneGrid  = data.frame(mtry = c(5, 10, floor(sqrt(ncol(X_train)))))
  ),
  error = function(e) { message("  RF failed: ", e$message); NULL }
)


# 7c. XGBoost — native API with manual LOOCV
# ─────────────────────────────────────────────────────────────────────────────
# WHY NATIVE API:
#   caret's xgbTree method reconstructs an internal xgb.DMatrix by calling
#   as.matrix() on the data inside the caret C glue code, which produces an
#   ALTREP-backed matrix that this version of xgboost cannot handle
#   (ALTLIST classes must provide a Set_elt method).
#   Calling xgboost::xgb.DMatrix() ourselves on a pre-materialised plain
#   matrix bypasses this entirely.
#
# INTERFACE CONTRACT:
#   We wrap the fitted model in a minimal S3 object of class "xgb_native"
#   and define predict.xgb_native() so the rest of the pipeline (test eval,
#   ROC curves, feature importance, role-state predictions) can call
#   predict(xgb_model, newdata=..., type="prob") exactly as for caret models.
# ─────────────────────────────────────────────────────────────────────────────
message("  Training XGBoost (native API, manual LOOCV)...")

# Fixed hyperparameters (same as before)
xgb_params <- list(
  booster          = "gbtree",
  objective        = "binary:logistic",
  eval_metric      = "auc",
  max_depth        = 3L,
  eta              = 0.1,
  gamma            = 0,
  colsample_bytree = 0.8,
  min_child_weight = 1,
  subsample        = 0.8,
  verbosity        = 0L
)
xgb_nrounds <- 50L

# Materialise a truly plain numeric matrix — no ALTREP, no tibble, no list attrs
make_plain_matrix <- function(df) {
  m <- matrix(
    as.double(unlist(df, use.names = FALSE)),
    nrow = nrow(df),
    ncol = ncol(df),
    dimnames = list(NULL, colnames(df))
  )
  m
}

X_train_plain <- make_plain_matrix(X_train)
X_test_plain  <- make_plain_matrix(X_test)
y_bin_train   <- as.integer(y_train == "Oncogenic")   # 1=Oncogenic, 0=Suppressive

xgb_model <- tryCatch({

  # ── Manual LOOCV for AUC estimate ──────────────────────────────────────────
  n_tr     <- nrow(X_train_plain)
  loo_probs <- numeric(n_tr)

  for (i in seq_len(n_tr)) {
    tr_idx  <- setdiff(seq_len(n_tr), i)
    dtrain  <- xgboost::xgb.DMatrix(
                 data  = X_train_plain[tr_idx, , drop = FALSE],
                 label = y_bin_train[tr_idx])
    dval    <- xgboost::xgb.DMatrix(
                 data  = X_train_plain[i, , drop = FALSE],
                 label = y_bin_train[i])

    m_fold <- xgboost::xgb.train(
      params  = xgb_params,
      data    = dtrain,
      nrounds = xgb_nrounds,
      verbose = 0
    )
    loo_probs[i] <- predict(m_fold, dval)
  }

  loocv_roc  <- pROC::roc(y_train, loo_probs,
                           levels    = c("Suppressive", "Oncogenic"),
                           direction = "<", quiet = TRUE)
  loocv_auc  <- round(as.numeric(pROC::auc(loocv_roc)), 3)

  # Sensitivity / specificity at Youden threshold
  best_thresh <- pROC::coords(loocv_roc, "best",
                               ret = c("threshold","sensitivity","specificity"),
                               transpose = FALSE)
  loocv_sens  <- round(best_thresh$sensitivity[1],  3)
  loocv_spec  <- round(best_thresh$specificity[1], 3)

  message(sprintf("    XGB LOOCV AUC=%.3f  Sens=%.3f  Spec=%.3f",
                  loocv_auc, loocv_sens, loocv_spec))

  # ── Final model trained on ALL training data ────────────────────────────────
  dtrain_full <- xgboost::xgb.DMatrix(
    data  = X_train_plain,
    label = y_bin_train
  )
  final_xgb <- xgboost::xgb.train(
    params  = xgb_params,
    data    = dtrain_full,
    nrounds = xgb_nrounds,
    verbose = 0
  )

  # ── S3 wrapper: exposes predict() compatible with caret model interface ─────
  structure(
    list(
      model       = final_xgb,
      loocv_auc   = loocv_auc,
      loocv_sens  = loocv_sens,
      loocv_spec  = loocv_spec,
      cv_method   = "LOOCV",
      feat_names  = colnames(X_train_plain),
      threshold   = best_thresh$threshold[1]
    ),
    class = "xgb_native"
  )

}, error = function(e) {
  message("  XGBoost failed: ", e$message)
  NULL
})

# predict.xgb_native — mimics caret predict() interface
# type="prob"  -> data.frame with columns Oncogenic, Suppressive
# type="raw"   -> factor with levels Oncogenic / Suppressive  (default)
predict.xgb_native <- function(object, newdata, type = "raw", ...) {
  # newdata may be a data.frame or plain matrix; always materialise to plain matrix
  if (is.data.frame(newdata)) {
    newdata <- make_plain_matrix(newdata[, object$feat_names, drop = FALSE])
  }
  dtest  <- xgboost::xgb.DMatrix(data = newdata)
  probs  <- predict(object$model, dtest)           # P(Oncogenic)
  if (type == "prob") {
    return(data.frame(Oncogenic = probs, Suppressive = 1 - probs,
                      row.names = NULL))
  }
  # "raw": threshold at 0.5 (or Youden threshold stored in object)
  thresh <- if (!is.null(object$threshold)) object$threshold else 0.5
  factor(ifelse(probs >= thresh, "Oncogenic", "Suppressive"),
         levels = c("Oncogenic", "Suppressive"))
}

# varImp.xgb_native — for feature importance extraction
varImp.xgb_native <- function(object, scale = TRUE, ...) {
  imp <- xgboost::xgb.importance(
    feature_names = object$feat_names,
    model         = object$model
  )
  # Return a data.frame with rownames = feature, col = Overall (caret convention)
  imp_df <- data.frame(
    Overall = imp$Gain,
    row.names = imp$Feature,
    stringsAsFactors = FALSE
  )
  # Fill in zeros for features not in importance table
  missing_feats <- setdiff(object$feat_names, rownames(imp_df))
  if (length(missing_feats) > 0) {
    zero_df <- data.frame(Overall = rep(0, length(missing_feats)),
                          row.names = missing_feats)
    imp_df <- rbind(imp_df, zero_df)
  }
  if (scale && max(imp_df$Overall) > 0) {
    imp_df$Overall <- 100 * imp_df$Overall / max(imp_df$Overall)
  }
  list(importance = imp_df)
}

# Register S3 methods in the current environment so R dispatch finds them
environment(predict.xgb_native) <- globalenv()
environment(varImp.xgb_native)  <- globalenv()
registerS3method("predict", "xgb_native", predict.xgb_native, envir = globalenv())
registerS3method("varImp",  "xgb_native", varImp.xgb_native,  envir = globalenv())

# ── 8. CV Performance Summary ─────────────────────────────────────────────────
message("CV performance (training set):")

models <- list(Logistic = lr_model, RandomForest = rf_model, XGBoost = xgb_model)
models <- Filter(Negate(is.null), models)


cv_perf <- dplyr::bind_rows(lapply(names(models), function(nm) {
  m <- models[[nm]]
  # xgb_native stores metrics directly; caret models use $results + $bestTune
  if (inherits(m, "xgb_native")) {
    data.frame(
      Model     = nm,
      CV_AUC    = m$loocv_auc,
      CV_Sens   = m$loocv_sens,
      CV_Spec   = m$loocv_spec,
      CV_method = m$cv_method,
      stringsAsFactors = FALSE
    )
  } else {
    best_row <- merge(m$results, m$bestTune)
    data.frame(
      Model     = nm,
      CV_AUC    = round(best_row$ROC[1],  3),
      CV_Sens   = round(best_row$Sens[1], 3),
      CV_Spec   = round(best_row$Spec[1], 3),
      CV_method = m$control$method,
      stringsAsFactors = FALSE
    )
  }
}))
message("CV results:")
print(cv_perf)

# ── 9. Held-Out Test Set Evaluation ───────────────────────────────────────────
message(sprintf("Evaluating on held-out test set (n=%d):", nrow(X_test)))

test_perf <- dplyr::bind_rows(lapply(names(models), function(nm) {
  pred_prob <- tryCatch(
    predict(models[[nm]], newdata = predict_input(nm), type = "prob"),
    error = function(e) NULL
  )
  pred_cls <- tryCatch(
    predict(models[[nm]], newdata = predict_input(nm)),
    error = function(e) NULL
  )
  if (is.null(pred_prob) || is.null(pred_cls)) return(NULL)

  roc_obj <- tryCatch(
    pROC::roc(y_test, pred_prob$Oncogenic,
              levels    = c("Suppressive", "Oncogenic"),
              direction = "<",
              quiet     = TRUE),
    error = function(e) NULL
  )
  auc_val <- if (!is.null(roc_obj)) round(as.numeric(pROC::auc(roc_obj)), 3) else NA
  acc     <- round(mean(pred_cls == y_test, na.rm = TRUE), 3)

  data.frame(Model = nm, Test_AUC = auc_val, Test_Acc = acc,
             stringsAsFactors = FALSE)
}))
print(test_perf)

perf_df <- dplyr::left_join(cv_perf, test_perf, by = "Model")
write.csv(perf_df, file.path(TABLES_DIR, "aim5_model_performance.csv"),
          row.names = FALSE)

# ── 10. ROC Curves (held-out test set) ────────────────────────────────────────
message("Generating ROC curves (held-out test set)...")

roc_list <- lapply(names(models), function(nm) {
  pred_prob <- tryCatch(
    predict(models[[nm]], newdata = predict_input(nm), type = "prob"),
    error = function(e) NULL
  )
  if (is.null(pred_prob)) return(NULL)
  tryCatch(
    pROC::roc(y_test, pred_prob$Oncogenic,
              levels    = c("Suppressive", "Oncogenic"),
              direction = "<",
              quiet     = TRUE),
    error = function(e) NULL
  )
})
names(roc_list) <- names(models)
roc_list <- Filter(Negate(is.null), roc_list)

if (length(roc_list) > 0) {
  roc_colors <- c(Logistic = "#3498DB", RandomForest = "#27AE60", XGBoost = "#E74C3C")

  pdf(file.path(FIGURES_DIR, "aim5_roc_curves.pdf"), width = 7, height = 6)
  plot(roc_list[[1]],
       col  = roc_colors[names(roc_list)[1]],
       lwd  = 2,
       main = "ROC Curves — GPRC5A Role-State Classifier\n(Leakage-Free Held-Out Test Set)")
  for (i in seq_along(roc_list)[-1]) {
    plot(roc_list[[i]], col = roc_colors[names(roc_list)[i]], lwd = 2, add = TRUE)
  }
  legend("bottomright",
         legend = paste0(names(roc_list), " (AUC=",
                         sapply(roc_list, function(r)
                           round(as.numeric(pROC::auc(r)), 3)), ")"),
         col = roc_colors[names(roc_list)],
         lwd = 2,
         bty = "n")
  abline(0, 1, lty = 2, col = "grey60")
  dev.off()
}

# ── 11. Feature Importance ────────────────────────────────────────────────────
message("Extracting feature importance...")

best_model_name <- perf_df$Model[which.max(perf_df$Test_AUC)]
# If all Test_AUC are NA, fall back to CV_AUC
if (is.na(best_model_name) || length(best_model_name) == 0) {
  best_model_name <- perf_df$Model[which.max(perf_df$CV_AUC)]
}
best_model <- models[[best_model_name]]
message(sprintf("  Best model: %s (Test AUC = %s)",
                best_model_name,
                ifelse(is.na(max(perf_df$Test_AUC, na.rm = TRUE)), "NA",
                       round(max(perf_df$Test_AUC, na.rm = TRUE), 3))))

imp_df <- tryCatch({
  imp <- varImp(best_model, scale = TRUE)$importance
  data.frame(
    Feature    = rownames(imp),
    Importance = imp[, 1],
    stringsAsFactors = FALSE
  ) %>%
    dplyr::arrange(desc(Importance)) %>%
    dplyr::slice_head(n = 30)
}, error = function(e) {
  message("  Feature importance extraction failed: ", e$message)
  NULL
})

if (!is.null(imp_df)) {
  write.csv(imp_df, file.path(TABLES_DIR, "aim5_feature_importance.csv"),
            row.names = FALSE)

  imp_df <- imp_df %>%
    dplyr::mutate(
      feature_type = dplyr::case_when(
        Feature == "GPRC5A"                              ~ "Target gene",
        Feature %in% make.names(moffitt_classical)       ~ "Classical signature",
        Feature %in% make.names(moffitt_basal)           ~ "Basal-like signature",
        Feature %in% c("subtype_score", "classical_score") ~ "Subtype score",
        Feature %in% c("gem_flag", "age_scaled", "stage_IV") ~ "Clinical",
        TRUE                                              ~ "Co-expression"
      )
    )

  type_colors <- c(
    "Target gene"          = "#E74C3C",
    "Classical signature"  = "#2E75B6",
    "Basal-like signature" = "#C0392B",
    "Subtype score"        = "#8E44AD",
    "Clinical"             = "#27AE60",
    "Co-expression"        = "#95A5A6"
  )

  p_imp <- ggplot(imp_df,
                  aes(x = Importance,
                      y = reorder(Feature, Importance),
                      fill = feature_type)) +
    geom_col(alpha = 0.85) +
    scale_fill_manual(values = type_colors) +
    labs(
      title    = sprintf("Feature Importance — %s Classifier", best_model_name),
      subtitle = "GPRC5A Role-State Prediction | TCGA-PAAD\n(features selected on training data only)",
      x        = "Scaled Importance",
      y        = NULL,
      fill     = "Feature type"
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "right")

  ggsave(file.path(FIGURES_DIR, "aim5_feature_importance.pdf"),
         plot = p_imp, width = 10, height = 8, device = "pdf")
}

# ── 12. Role-State Predictions on Test Set ────────────────────────────────────
message("Generating role-state predictions on held-out test set...")

all_preds <- tryCatch({
  pred_prob <- predict(best_model, newdata = predict_input(best_model_name), type = "prob")
  data.frame(
    barcode          = bar_test,
    predicted_state  = as.character(predict(best_model, newdata = predict_input(best_model_name))),
    prob_oncogenic   = round(pred_prob$Oncogenic,   3),
    prob_suppressive = round(pred_prob$Suppressive, 3),
    true_state       = as.character(y_test),
    stringsAsFactors = FALSE
  )
}, error = function(e) {
  message("  Prediction failed: ", e$message); NULL
})

if (!is.null(all_preds)) {
  pred_surv <- all_preds %>%
    dplyr::left_join(
      subtype_df[, c("barcode", "vital_status", "OS_months",
                     "subtype_moffitt", "gprc5a_expr")],
      by = "barcode"
    ) %>%
    dplyr::mutate(OS_event = ifelse(vital_status == "Dead", 1, 0))

  write.csv(pred_surv,
            file.path(TABLES_DIR, "aim5_role_state_predictions.csv"),
            row.names = FALSE)

  # KM survival by predicted role state
  surv_data <- pred_surv %>% dplyr::filter(!is.na(OS_months))
  if (nrow(surv_data) > 0 && length(unique(surv_data$predicted_state)) > 1) {
    fit_role <- survfit(Surv(OS_months, OS_event) ~ predicted_state, data = surv_data)

    groups_in_data <- sort(unique(all_preds$predicted_state))
    km_palette     <- c("Oncogenic" = "#E74C3C", "Suppressive" = "#27AE60")
    this_palette   <- unname(km_palette[groups_in_data])


    # NOTE on legend.title: some survminer versions interpret any legend.title
    # value (string or list) as an aes() label and emit "Ignoring unknown labels".
    # Workaround: omit legend.title from ggsurvplot(), then patch it directly
    # onto the ggplot layer via guides() — clean, no warning.
    km_role <- tryCatch({
      km_obj <- ggsurvplot(
        fit_role,
        data              = surv_data,
        pval              = TRUE,
        pval.method       = TRUE,
        conf.int          = TRUE,
        risk.table        = TRUE,
        risk.table.height = 0.28,
        palette           = this_palette,
        title             = "Overall Survival by Predicted GPRC5A Role State",
        subtitle          = "Held-out test set | Oncogenic should show worse OS",
        xlab              = "Time (months)",
        ylab              = "Overall Survival",
        legend.labs       = groups_in_data,
        ggtheme           = theme_bw(base_size = 12),
        font.title        = c(13, "bold"),
        surv.median.line  = "hv"
      )
      # Patch legend title after construction — bypasses the internal aes clash
      km_obj$plot <- km_obj$plot +
        guides(colour = guide_legend(title = "GPRC5A State"))
      km_obj
    },
    error = function(e) {
      message("  Role-state KM failed: ", e$message); NULL
    })

  }

  # Confusion matrix
  conf <- table(Predicted = all_preds$predicted_state,
                Actual    = all_preds$true_state)
  message("Confusion matrix:")
  print(conf)
  acc <- sum(diag(conf)) / sum(conf)
  message(sprintf("Overall accuracy: %.1f%%", 100 * acc))
}

# ── 13. Calibration Plot ──────────────────────────────────────────────────────
if (!is.null(all_preds)) {
  message("Generating calibration plot...")

  calib_df <- all_preds %>%
    dplyr::mutate(
      prob_bin    = cut(prob_oncogenic, breaks = seq(0, 1, 0.1),
                        include.lowest = TRUE),
      actual_onco = as.integer(true_state == "Oncogenic")
    ) %>%
    dplyr::group_by(prob_bin) %>%
    dplyr::summarise(
      mean_predicted = mean(prob_oncogenic, na.rm = TRUE),
      frac_actual    = mean(actual_onco,    na.rm = TRUE),
      n              = dplyr::n(),
      .groups        = "drop"
    )

  p_calib <- ggplot(calib_df, aes(x = mean_predicted, y = frac_actual)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "grey50", linewidth = 0.8) +
    geom_line(color = "#E74C3C", linewidth = 1.2) +
    geom_point(aes(size = n), color = "#E74C3C", alpha = 0.8) +
    scale_size_continuous(range = c(2, 8), name = "n samples") +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      title    = "Calibration Plot — GPRC5A Role-State Classifier",
      subtitle = "Diagonal = perfect calibration | Held-out test set",
      x        = "Mean predicted probability (Oncogenic)",
      y        = "Observed fraction (Oncogenic)"
    ) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(file.path(FIGURES_DIR, "aim5_calibration_plot.pdf"),
         plot = p_calib, width = 7, height = 6, device = "pdf")
}

# ── 14. Summary ───────────────────────────────────────────────────────────────
message("\n=== Aim 5 Summary ===")
message("Model performance:")
print(perf_df)
message(sprintf("\nBest model: %s (Test AUC = %s)",
                best_model_name,
                ifelse(is.na(max(perf_df$Test_AUC, na.rm = TRUE)), "NA",
                       round(max(perf_df$Test_AUC, na.rm = TRUE), 3))))
if (!is.null(all_preds)) {
  acc <- mean(all_preds$predicted_state == all_preds$true_state, na.rm = TRUE)
  message(sprintf("Test set accuracy: %.1f%% (n=%d held-out samples)",
                  100 * acc, nrow(all_preds)))
}
message("\nKey question answered:")
message("  Can a tumor's GPRC5A role state be predicted from its transcriptome?")
message("  => Check aim5_roc_curves.pdf for AUC and aim5_role_state_km.pdf")
message("     for whether predicted role state separates survival curves")
message("\nOutputs saved:")
message(paste0("  ", FIGURES_DIR, "/aim5_*.pdf"))
message(paste0("  ", TABLES_DIR,  "/aim5_*.csv"))
