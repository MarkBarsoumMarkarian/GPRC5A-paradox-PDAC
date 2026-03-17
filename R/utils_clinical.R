# =============================================================================
# utils_clinical.R — Shared clinical data utilities
# Project: Decoding the GPRC5A Paradox in PDAC
#
# Provides robust, version-agnostic TCGA clinical column detection.
# Source this file at the top of any aim script.
# =============================================================================

#' Find the best matching column from a list of candidates
#' Handles the fact that TCGAbiolinks changes clinical column names across versions
find_col <- function(candidates, cols, required = TRUE, label = NULL) {
  # Exact match first
  hit <- candidates[candidates %in% cols]
  if (length(hit) > 0) return(hit[1])
  
  # Partial / pattern match as fallback
  pattern <- paste(candidates, collapse = "|")
  partial <- cols[grep(pattern, cols, ignore.case = TRUE)]
  if (length(partial) > 0) {
    msg <- sprintf("  [find_col] Exact match not found for {%s}, using: %s",
                   paste(candidates, collapse="|"), partial[1])
    message(msg)
    return(partial[1])
  }
  
  if (required) {
    stop(sprintf(
      "Required clinical column not found.\nTried: %s\nAvailable cols (first 40): %s",
      paste(candidates, collapse = ", "),
      paste(head(cols, 40), collapse = ", ")
    ))
  }
  return(NA_character_)
}

#' Robustly build a survival-ready clinical dataframe from TCGA colData
#' Works with TCGAbiolinks 2.x and 3.x column naming conventions
build_clinical_df <- function(clinical_raw) {
  df    <- as.data.frame(clinical_raw)
  cols  <- colnames(df)
  
  message("  Detecting clinical columns...")
  
  # Core survival columns
  col_vital  <- find_col(c("vital_status"), cols)
  col_death  <- find_col(c("days_to_death"), cols)
  col_lastfu <- find_col(c("days_to_last_follow_up",
                            "days_to_last_followup",
                            "last_contact_days_to"), cols)
  
  # Optional columns — gracefully handled if absent
  col_stage  <- find_col(c("ajcc_pathologic_stage", "tumor_stage",
                            "pathologic_stage", "clinical_stage",
                            "paper_Tumor.grade"), cols, required = FALSE)
  col_age    <- find_col(c("age_at_diagnosis", "age_at_index",
                            "days_to_birth"), cols, required = FALSE)
  col_gender <- find_col(c("gender", "sex"), cols, required = FALSE)
  
  # Treatment columns — TCGA encodes these very inconsistently
  tx_cols <- grep("treatment|therapy|drug|pharmaceutical|regimen|chemotherapy",
                  cols, ignore.case = TRUE, value = TRUE)
  message(sprintf("  Found %d treatment-related columns: %s",
                  length(tx_cols),
                  paste(head(tx_cols, 5), collapse=", ")))
  
  result <- df %>%
    mutate(
      # Survival
      vital_status   = .data[[col_vital]],
      days_to_death  = suppressWarnings(as.numeric(.data[[col_death]])),
      days_to_lastfu = suppressWarnings(as.numeric(.data[[col_lastfu]])),
      
      # Stage
      tumor_stage = if (!is.na(col_stage)) .data[[col_stage]] else NA_character_,
      
      # Age — stored variously as years, days, or negative days
      age_raw = if (!is.na(col_age)) suppressWarnings(as.numeric(.data[[col_age]])) else NA_real_,
      age_years = case_when(
        !is.na(age_raw) & age_raw < 0     ~ abs(age_raw) / 365.25,  # negative days
        !is.na(age_raw) & age_raw > 200   ~ age_raw / 365.25,        # positive days
        !is.na(age_raw)                   ~ age_raw,                  # already years
        TRUE                              ~ NA_real_
      ),
      
      # Gender
      gender = if (!is.na(col_gender)) .data[[col_gender]] else NA_character_,
      
      # Derived survival variables
      OS_time   = ifelse(!is.na(days_to_death) & days_to_death > 0,
                         days_to_death, days_to_lastfu),
      OS_event  = ifelse(vital_status == "Dead", 1, 0),
      OS_months = OS_time / 30.44,
      survival_group = ifelse(vital_status == "Dead", "Deceased", "Alive"),
      
      stage_binary = case_when(
        grepl("IV",        tumor_stage, ignore.case = TRUE) ~ "Stage IV",
        grepl("I|II|III",  tumor_stage, ignore.case = TRUE) ~ "Stage I-III",
        TRUE ~ "Unknown"
      ),
      
      # Treatment: concatenate all tx fields into one searchable string
      tx_string = if (length(tx_cols) > 0) {
        apply(df[, tx_cols, drop = FALSE], 1,
              function(x) paste(tolower(as.character(x[!is.na(x)])), collapse = " "))
      } else {
        rep("", nrow(df))
      }
    ) %>%
    filter(!is.na(OS_time), OS_time > 0)
  
  # Treatment flags derived from tx_string
  result <- result %>%
    mutate(
      received_gemcitabine  = stringr::str_detect(tx_string,
        "gemcitabine|gemzar|\\bgem\\b"),
      received_fluorouracil = stringr::str_detect(tx_string,
        "fluorouracil|5-fu|folfirinox|capecitabine"),
      received_any_chemo    = received_gemcitabine | received_fluorouracil |
        stringr::str_detect(tx_string, "chemo|oxaliplatin|irinotecan|abraxane|paclitaxel"),
      received_radiation    = stringr::str_detect(tx_string,
        "radiation|radiotherapy|xrt|sbrt"),
      treatment_naive       = !received_any_chemo & !received_radiation,
      
      treatment_group = dplyr::case_when(
        received_gemcitabine  ~ "Gemcitabine",
        received_fluorouracil ~ "Fluorouracil/FOLFIRINOX",
        received_any_chemo    ~ "Other Chemotherapy",
        received_radiation    ~ "Radiation only",
        TRUE                  ~ "Treatment-naive / Unknown"
      )
    )
  
  message(sprintf("  Clinical df built: n=%d samples", nrow(result)))
  message(sprintf("  Alive=%d | Deceased=%d",
                  sum(result$survival_group == "Alive"),
                  sum(result$survival_group == "Deceased")))
  
  return(result)
}
