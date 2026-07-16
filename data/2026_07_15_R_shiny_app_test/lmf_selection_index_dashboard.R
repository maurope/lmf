# ==============================================================================
# LMF Selection Index Dashboard
# ==============================================================================
#
# A single-file R Shiny application that implements a Livestock and Methane
# Forage (LMF) Selection Index Dashboard. This tool enables researchers to:
#
#   1. Upload raw multi-trial accession data (CSV, TSV, Excel)
#   2. Map dataset columns to semantic roles (accession ID, trial ID, traits)
#   3. Aggregate replicated observations into accession-level best estimates
#   4. Configure classification thresholds and scoring weights
#   5. Calculate selection indices across four predefined selection approaches
#   6. Display ranked results in sortable, filterable tables
#   7. Visualize accession distribution via interactive scatter plots
#   8. Diagnose data quality issues
#   9. Export complete results as a multi-sheet Excel workbook
#
# The application replaces an Excel-based workflow, providing an interactive
# interface for the complete LMF selection index analysis pipeline.
#
# ==============================================================================

# --- Install Missing Packages -------------------------------------------------

required_packages <- c(
  "shiny", "bslib", "shinyjs", "dplyr", "tidyr", "readr", "readxl", "openxlsx2", 
  "DT", "ggplot2", "ggExtra", "patchwork", "plotly", "lme4", "emmeans"
)

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

# --- Package Loading ----------------------------------------------------------

library(shiny)
library(shinyjs)
library(bslib)
library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(openxlsx2)
library(DT)
library(ggplot2)
library(ggExtra)
library(patchwork)
library(plotly)
library(lme4)
library(emmeans)

# --- Global Options -----------------------------------------------------------

options(shiny.maxRequestSize = 50 * 1024^2)  # 50 MB upload limit

# --- Standalone Functions -----------------------------------------------------
# (Defined at the top level so they can be sourced and tested independently)

#' Parse an uploaded file and return a data.frame
#'
#' Reads a user-uploaded file (CSV, TSV, or Excel) and returns its contents
#' as a data.frame. The file format is detected from the extension of the
#' original file name. The first row is treated as column headers.
#'
#' @param file_path Character. Temporary file path from Shiny file input.
#' @param file_name Character. Original file name (used to detect format).
#' @return A data.frame with all columns as-is from the file, OR a list with
#'   element \code{error} containing an error message string.
#' @details Detects format from extension: .csv/.tsv → readr; .xlsx/.xls → readxl.
#'          First row is treated as headers. Returns error message on failure.
#'          Error conditions handled:
#'          - Unsupported file extension
#'          - File with zero data rows (headers only)
#'          - Parse failures (encoding errors, file corruption)
parse_upload <- function(file_path, file_name) {
  # Extract file extension (lowercase, without dot)
  ext <- tolower(tools::file_ext(file_name))

  # Validate supported formats
  supported_formats <- c("csv", "tsv", "xlsx", "xls")
  if (!(ext %in% supported_formats)) {
    return(list(error = "Unsupported file format. Please upload a CSV, TSV, .xlsx, or .xls file."))
  }

  # Attempt to read the file based on detected format
  df <- tryCatch(
    {
      if (ext == "csv") {
        readr::read_csv(file_path, col_names = TRUE, show_col_types = FALSE)
      } else if (ext == "tsv") {
        readr::read_tsv(file_path, col_names = TRUE, show_col_types = FALSE)
      } else if (ext %in% c("xlsx", "xls")) {
        readxl::read_excel(file_path, sheet = 1, col_names = TRUE)
      }
    },
    error = function(e) {
      return(NULL)
    }
  )

  # Handle parse failure

  if (is.null(df)) {
    return(list(error = "Unable to read the file. Please verify the file integrity and encoding."))
  }

  # Convert tibble to plain data.frame
  df <- as.data.frame(df)

  # Check for empty file (0 data rows)
  if (nrow(df) == 0) {
    return(list(error = "The uploaded file contains no data rows (only headers detected)."))
  }

  return(df)
}

#' Validate whether a column contains a sufficient proportion of numeric values
#'
#' Checks if at least 50% of non-missing values in a column can be coerced to
#' numeric. Missing values (NA, blank "", or whitespace-only strings) are excluded
#' from both the numerator and denominator of the proportion calculation.
#'
#' @param x A vector (typically character) representing a column from a data.frame.
#' @return Logical. TRUE if the proportion of successfully numeric-coercible
#'   non-missing values is >= 0.5, or if all values are missing. FALSE otherwise.
#' @details
#'   - NA, blank (""), and whitespace-only strings are treated as missing and
#'     excluded from the check (not counted as non-numeric failures).
#'   - Remaining non-missing values are tested via suppressWarnings(as.numeric(x)).
#'   - If all values are missing (count_non_missing == 0), returns TRUE.
#' @examples
#' validate_numeric_column(c("1", "2", "abc", NA, ""))
#' # 3 non-missing: "1", "2", "abc"; 2 numeric out of 3 => 0.67 >= 0.5 => TRUE
#'
#' validate_numeric_column(c("abc", "def", "ghi"))
#' # 3 non-missing, 0 numeric => 0/3 = 0 < 0.5 => FALSE
validate_numeric_column <- function(x) {
  # Convert to character to handle mixed-type inputs uniformly
  x <- as.character(x)

  # Identify missing values: NA, blank, or whitespace-only
  is_missing <- is.na(x) | trimws(x) == ""

  # Extract non-missing values
  non_missing <- x[!is_missing]

  # Edge case: if all values are missing, return TRUE (no non-numeric failures)
  count_non_missing <- length(non_missing)
  if (count_non_missing == 0) {
    return(TRUE)
  }

  # Attempt numeric coercion on non-missing values
  numeric_values <- suppressWarnings(as.numeric(non_missing))

  # Count how many successfully converted (non-NA after coercion)
  count_numeric <- sum(!is.na(numeric_values))

  # Compute proportion and check threshold
  proportion <- count_numeric / count_non_missing
  return(proportion >= 0.5)
}

#' Aggregate Accession-Level Estimates from Raw Multi-Trial Data
#'
#' Groups raw observations by accession and computes summary statistics for each
#' numerical trait and categorical grouping variable. Returns a best-estimate
#' dataset (one row per accession) and a quality-control dataset with variability
#' and replication metrics.
#'
#' @param raw_data data.frame. The uploaded raw data containing per-accession,
#'   per-trial observations.
#' @param accession_col Character. Name of the column containing the accession
#'   (genotype) identifier.
#' @param trial_col Character. Name of the column containing the trial/run
#'   identifier.
#' @param trait_cols Character vector. Names of numerical trait columns to
#'   aggregate (e.g., CH4_Intensity, TDDM).
#' @param group_cols Character vector. Names of categorical/grouping columns
#'   (e.g., Functional_Group, Taxonomic_Name). May be empty (character(0)).
#' @return A named list with two elements:
#'   \describe{
#'     \item{estimates}{data.frame with one row per unique accession. Contains
#'       the accession ID column, the arithmetic mean of each trait (NA when all
#'       observations are NA), and the modal value for each grouping variable.}
#'     \item{qc}{data.frame with one row per unique accession. Contains the
#'       accession ID column, sample standard deviation (sd_{trait}, NA when
#'       n=1), and count of non-missing observations (n_obs_{trait}) for each
#'       trait.}
#'   }
aggregate_accessions <- function(raw_data, accession_col, trial_col, trait_cols, group_cols) {
  # Get unique accessions preserving original order of appearance
  accessions <- unique(raw_data[[accession_col]])

  # Initialize output data.frames
  estimates <- data.frame(accession = accessions, stringsAsFactors = FALSE)
  names(estimates)[1] <- accession_col

  qc <- data.frame(accession = accessions, stringsAsFactors = FALSE)
  names(qc)[1] <- accession_col

  # Split raw data by accession for efficient grouped computation
  split_data <- split(raw_data, raw_data[[accession_col]])

  # --- Aggregate numerical trait columns ---
  for (trait in trait_cols) {
    trait_means <- vapply(accessions, function(acc) {
      values <- split_data[[as.character(acc)]][[trait]]
      non_na <- values[!is.na(values)]
      if (length(non_na) == 0L) {
        return(NA_real_)
      }
      mean(non_na)
    }, numeric(1))

    trait_sds <- vapply(accessions, function(acc) {
      values <- split_data[[as.character(acc)]][[trait]]
      non_na <- values[!is.na(values)]
      if (length(non_na) <= 1L) {
        return(NA_real_)
      }
      sd(non_na)
    }, numeric(1))

    trait_nobs <- vapply(accessions, function(acc) {
      values <- split_data[[as.character(acc)]][[trait]]
      sum(!is.na(values))
    }, integer(1))

    estimates[[trait]] <- trait_means
    qc[[paste0("sd_", trait)]] <- trait_sds
    qc[[paste0("n_obs_", trait)]] <- trait_nobs
  }

  # --- Aggregate categorical/grouping columns ---
  for (gcol in group_cols) {
    group_modes <- vapply(accessions, function(acc) {
      values <- split_data[[as.character(acc)]][[gcol]]
      # Remove NA values before computing mode
      values <- values[!is.na(values)]
      if (length(values) == 0L) {
        return(NA_character_)
      }
      # table() preserves order of first occurrence; sort decreasing picks
      # the most frequent value with first-encountered tie-break
      freq_table <- sort(table(values), decreasing = TRUE)
      names(freq_table)[1]
    }, character(1))

    estimates[[gcol]] <- group_modes
  }

  list(estimates = estimates, qc = qc)
}

#' Calculate Selection Index Across Four Selection Approaches
#'
#' Computes reference baseline, CH4 decrease percentages, classifies accessions
#' into CH4 and TDDM categories, calculates Z-scores and composite scores within
#' each selection approach subset, assigns competition rankings, and computes
#' additional reduction metrics (CH4_8h, CH4_24h, Benchmark vs Control).
#'
#' @param estimates data.frame. Accession-level estimates with columns for
#'        CH4_Intensity and TDDM at minimum. Must contain a column named
#'        "ch4_intensity" and "tddm" (or the mapped column names).
#' @param config list. Configuration parameters:
#'   - ch4_col: character. Name of CH4 Intensity column in estimates
#'   - tddm_col: character. Name of TDDM column in estimates
#'   - accession_col: character. Name of accession ID column
#'   - ch4_low_threshold: numeric (percentage, e.g. 35)
#'   - ch4_medium_threshold: numeric (percentage, e.g. 20)
#'   - tddm_high_cutoff: numeric (mg/g, e.g. 60)
#'   - reference_baseline_pct: numeric (percentage, e.g. 10)
#'   - reference_control_id: character (e.g. "ABC-18447-1")
#'   - weights: named list of 4 approaches, each with w1 and w2
#'     - cat_1a = list(w1 = 2, w2 = 1)
#'     - cat_1b = list(w1 = 2, w2 = 1)
#'     - cat_2  = list(w1 = 1, w2 = 0)
#'     - cat_3  = list(w1 = 1, w2 = 2)
#'   - ch4_8h_col: character or NULL (column name for CH4 8h data)
#'   - ch4_24h_col: character or NULL (column name for CH4 24h data)
#' @return A named list with:
#'   - all_estimates: data.frame with CH4_Decrease_Pct, categories, Z-scores,
#'                    composite scores for each approach
#'   - approach_results: named list of 4 data.frames (cat_1a, cat_1b, cat_2, cat_3),
#'                       each filtered and ranked by Overall_Rank
#'   - reference_baseline_avg: numeric
calculate_selection_index <- function(estimates, config) {
  # Extract configuration parameters
  ch4_col <- config$ch4_col
  tddm_col <- config$tddm_col
  accession_col <- config$accession_col
  ch4_low_threshold <- config$ch4_low_threshold
  ch4_medium_threshold <- config$ch4_medium_threshold
  tddm_high_cutoff <- config$tddm_high_cutoff
  reference_baseline_pct <- config$reference_baseline_pct
  reference_control_id <- config$reference_control_id
  weights <- config$weights
  ch4_8h_col <- config$ch4_8h_col
  ch4_24h_col <- config$ch4_24h_col

  # Work on a copy of estimates
  df <- estimates

  # --- Step 1: Compute Reference Baseline Average ---
  ch4_values <- df[[ch4_col]]
  n_total <- length(ch4_values)
  n_top <- ceiling(n_total * reference_baseline_pct / 100)
  # Sort descending and take top n_top values

  sorted_ch4 <- sort(ch4_values, decreasing = TRUE, na.last = TRUE)
  top_values <- sorted_ch4[seq_len(n_top)]
  reference_baseline_avg <- mean(top_values, na.rm = TRUE)

  # --- Step 2: Compute CH4 Decrease Percentage ---
  df$CH4_Decrease_Pct <- ((reference_baseline_avg - df[[ch4_col]]) / reference_baseline_avg) * 100

  # --- Step 3: Classify CH4 Category ---
  df$CH4_Category <- ifelse(
    df$CH4_Decrease_Pct >= ch4_low_threshold, "Low",
    ifelse(df$CH4_Decrease_Pct >= ch4_medium_threshold, "Medium", "High")
  )

  # --- Step 4: Classify TDDM Category ---
  df$TDDM_Category <- ifelse(df[[tddm_col]] >= tddm_high_cutoff, "High", "Not High")

  # --- Step 5: CH4_8h Reduction Percentage ---
  if (!is.null(ch4_8h_col) && ch4_8h_col %in% names(df)) {
    mean_ch4_8h <- mean(df[[ch4_8h_col]], na.rm = TRUE)
    df$CH4_8h_Reduction_Pct <- ((mean_ch4_8h - df[[ch4_8h_col]]) / mean_ch4_8h) * 100
  } else {
    df$CH4_8h_Reduction_Pct <- NA_real_
  }

  # --- Step 6: CH4_24h Reduction Percentage ---
  if (!is.null(ch4_24h_col) && ch4_24h_col %in% names(df)) {
    mean_ch4_24h <- mean(df[[ch4_24h_col]], na.rm = TRUE)
    df$CH4_24h_Reduction_Pct <- ((mean_ch4_24h - df[[ch4_24h_col]]) / mean_ch4_24h) * 100
  } else {
    df$CH4_24h_Reduction_Pct <- NA_real_
  }

  # --- Step 7: Benchmark vs Control ---
  control_rows <- df[[accession_col]] == reference_control_id
  if (any(control_rows, na.rm = TRUE)) {
    control_ch4 <- df[[ch4_col]][which(control_rows)[1]]
    df$Benchmark_vs_Control <- ((control_ch4 - df[[ch4_col]]) / control_ch4) * 100
  } else {
    df$Benchmark_vs_Control <- NA_real_
  }

  # --- Step 8: Define approach filtering criteria ---
  approach_filters <- list(
    cat_1a = function(d) d$CH4_Category == "Low" & d$TDDM_Category == "High",
    cat_1b = function(d) d$CH4_Category == "Medium" & d$TDDM_Category == "High",
    cat_2  = function(d) (d$CH4_Category == "Low" | d$CH4_Category == "Medium") & d$TDDM_Category == "Not High",
    cat_3  = function(d) d$CH4_Category == "High" & d$TDDM_Category == "High"
  )

  # --- Step 9: For each approach, filter, compute Z-scores, composite, rank ---
  approach_names <- c("cat_1a", "cat_1b", "cat_2", "cat_3")

  # Compute Z-scores using the WHOLE dataset (matching Excel formula)
  # Z_Score_TDDM = (TDDM - mean(TDDM)) / sd(TDDM)
  # Z_Score_CH4 = -1 * (CH4 - mean(CH4)) / sd(CH4)  (negated: lower CH4 is better)
  all_ch4 <- df[[ch4_col]]
  all_tddm <- df[[tddm_col]]

  mean_ch4_all <- mean(all_ch4, na.rm = TRUE)
  sd_ch4_all <- sd(all_ch4, na.rm = TRUE)
  mean_tddm_all <- mean(all_tddm, na.rm = TRUE)
  sd_tddm_all <- sd(all_tddm, na.rm = TRUE)

  # Compute global Z-scores for all accessions
  if (is.na(sd_ch4_all) || sd_ch4_all == 0) {
    df$Z_Score_CH4 <- rep(0, nrow(df))
  } else {
    df$Z_Score_CH4 <- -1 * (all_ch4 - mean_ch4_all) / sd_ch4_all
  }

  if (is.na(sd_tddm_all) || sd_tddm_all == 0) {
    df$Z_Score_TDDM <- rep(0, nrow(df))
  } else {
    df$Z_Score_TDDM <- (all_tddm - mean_tddm_all) / sd_tddm_all
  }

  # Initialize composite and rank columns for all_estimates
  for (app_name in approach_names) {
    df[[paste0("Composite_", app_name)]] <- NA_real_
    df[[paste0("Rank_", app_name)]] <- NA_integer_
  }

  approach_results <- list()

  for (app_name in approach_names) {
    # Get filter mask
    mask <- approach_filters[[app_name]](df)
    # Handle NA in mask (from NA categories)
    mask[is.na(mask)] <- FALSE
    subset_indices <- which(mask)
    subset_df <- df[subset_indices, , drop = FALSE]

    w1 <- weights[[app_name]]$w1
    w2 <- weights[[app_name]]$w2

    n_subset <- nrow(subset_df)

    if (n_subset < 2) {
      # Edge case: fewer than 2 accessions - composite = NA
      composite <- rep(NA_real_, n_subset)
      ranks <- rep(NA_integer_, n_subset)
    } else {
      # Use global Z-scores (already computed above) for composite
      z_ch4 <- subset_df$Z_Score_CH4
      z_tddm <- subset_df$Z_Score_TDDM

      # Composite score
      composite <- w1 * z_ch4 + w2 * z_tddm

      # Competition ranking (descending, highest = rank 1)
      ranks <- rank(-composite, ties.method = "min")
    }

    # Store composite and rank back into main df for all_estimates
    if (n_subset > 0) {
      df[[paste0("Composite_", app_name)]][subset_indices] <- composite
      df[[paste0("Rank_", app_name)]][subset_indices] <- ranks
    }

    # Build approach result data.frame
    result_df <- subset_df
    result_df$Z_Score_CH4 <- subset_df$Z_Score_CH4
    result_df$Z_Score_TDDM <- subset_df$Z_Score_TDDM
    result_df$Composite_Score <- composite
    result_df$Overall_Rank <- ranks

    # Sort by Overall_Rank ascending (best first)
    if (n_subset > 0 && !all(is.na(ranks))) {
      result_df <- result_df[order(result_df$Overall_Rank), , drop = FALSE]
    }

    approach_results[[app_name]] <- result_df
  }

  # --- Return results ---
  list(
    all_estimates = df,
    approach_results = approach_results,
    reference_baseline_avg = reference_baseline_avg
  )
}

#' Compute Data Quality Diagnostics
#'
#' Generates a comprehensive set of data quality metrics and reference control
#' verification results from the raw uploaded data and aggregated accession
#' estimates. Used to populate the diagnostic panel and the diagnostics sheet
#' in the exported workbook.
#'
#' @param raw_data data.frame. Original uploaded data containing per-accession,
#'   per-trial observations.
#' @param estimates data.frame. Aggregated accession-level estimates (one row
#'   per unique accession).
#' @param mapping list. Column mapping configuration with elements:
#'   \describe{
#'     \item{accession_col}{character. Name of the accession ID column.}
#'     \item{trial_col}{character. Name of the trial/run ID column.}
#'     \item{trait_cols}{character vector. Names of mapped numerical trait columns.}
#'   }
#' @param config list. Configuration parameters with element:
#'   \describe{
#'     \item{reference_control_id}{character. The accession identifier designated
#'       as the reference control (e.g., "ABC-18447-1").}
#'   }
#' @return A named list with the following elements:
#'   \describe{
#'     \item{n_records}{integer. Total number of rows in the raw data.}
#'     \item{n_accessions}{integer. Number of unique accessions.}
#'     \item{n_trials}{integer. Number of unique trials/runs.}
#'     \item{obs_per_accession}{list with elements min (integer), median (numeric),
#'       max (integer) summarising the distribution of observation counts per accession.}
#'     \item{missing_per_trait}{named integer vector. Count of missing (NA or blank)
#'       values per mapped trait in the raw data.}
#'     \item{sd_per_trait}{named numeric vector. Standard deviation of each trait
#'       across all accession estimates.}
#'     \item{control_in_trials}{named logical vector. For each trial, TRUE if the
#'       reference control is present, FALSE otherwise. NULL if control not found at all.}
#'     \item{control_missing_trials}{character vector. Trial names where the reference
#'       control accession is absent. NULL if control not found at all.}
#'     \item{low_replication_count}{integer. Number of accessions with fewer than
#'       2 total observations in the raw data.}
#'     \item{control_error}{character or NULL. Error message if the reference control
#'       accession is not found in any record of the raw data; NULL otherwise.}
#'   }
compute_diagnostics <- function(raw_data, estimates, mapping, config) {
  accession_col <- mapping$accession_col
  trial_col <- mapping$trial_col
  trait_cols <- mapping$trait_cols
  reference_control_id <- config$reference_control_id


  # --- Basic counts ---
  n_records <- nrow(raw_data)
  n_accessions <- length(unique(raw_data[[accession_col]]))
  n_trials <- length(unique(raw_data[[trial_col]]))

  # --- Observations per accession summary ---
  obs_counts <- table(raw_data[[accession_col]])
  obs_per_accession <- list(
    min = as.integer(min(obs_counts)),
    median = median(as.numeric(obs_counts)),
    max = as.integer(max(obs_counts))
  )

  # --- Missing values per trait (NA or blank/empty string) ---
  missing_per_trait <- vapply(trait_cols, function(trait) {
    values <- raw_data[[trait]]
    # Count NA values
    na_count <- sum(is.na(values))
    # Count blank/empty string values (only for character columns)
    if (is.character(values)) {
      blank_count <- sum(trimws(values) == "", na.rm = TRUE)
    } else {
      blank_count <- 0L
    }
    as.integer(na_count + blank_count)
  }, integer(1))
  names(missing_per_trait) <- trait_cols

  # --- SD per trait across accession estimates ---
  sd_per_trait <- vapply(trait_cols, function(trait) {
    if (trait %in% names(estimates)) {
      sd(estimates[[trait]], na.rm = TRUE)
    } else {
      NA_real_
    }
  }, numeric(1))
  names(sd_per_trait) <- trait_cols

  # --- Reference control verification ---
  control_error <- NULL
  control_in_trials <- NULL
  control_missing_trials <- NULL

  if (!is.null(reference_control_id) && nzchar(reference_control_id)) {
    # Check if control exists anywhere in the data
    control_present <- reference_control_id %in% raw_data[[accession_col]]

    if (!control_present) {
      control_error <- paste0(
        "Reference control '", reference_control_id,
        "' was not found in the dataset."
      )
    } else {
      # For each trial, check if the control accession appears
      trials <- unique(raw_data[[trial_col]])
      control_in_trials <- vapply(trials, function(trial) {
        trial_data <- raw_data[raw_data[[trial_col]] == trial, , drop = FALSE]
        reference_control_id %in% trial_data[[accession_col]]
      }, logical(1))
      names(control_in_trials) <- trials

      # Identify trials where control is missing
      control_missing_trials <- names(control_in_trials[!control_in_trials])
    }
  }

  # --- Low replication count (accessions with < 2 observations) ---
  low_replication_count <- as.integer(sum(obs_counts < 2))

  # --- Return diagnostics named list ---
  list(
    n_records = n_records,
    n_accessions = n_accessions,
    n_trials = n_trials,
    obs_per_accession = obs_per_accession,
    missing_per_trait = missing_per_trait,
    sd_per_trait = sd_per_trait,
    control_in_trials = control_in_trials,
    control_missing_trials = control_missing_trials,
    low_replication_count = low_replication_count,
    control_error = control_error
  )
}


#' Calculate Technical-Replicate CV% for Each Trait
#'
#' Computes the coefficient of variation (CV%) within each accession x trial
#' combination for each numerical trait. CV% = 100 * SD / abs(mean).
#'
#' @param raw_data data.frame. The uploaded raw data.
#' @param trait_cols Character vector. Names of numerical trait columns.
#' @param accession_col Character. Name of the accession ID column.
#' @param trial_col Character. Name of the trial/run ID column.
#' @return A data.frame with one row per trait and summary statistics:
#'   trait, groups, groups_with_cv, groups_lt_3_obs, mean_cv_percent,
#'   median_cv_percent, sd_cv_percent, q90_cv_percent.
calculate_trait_cv <- function(raw_data, trait_cols, accession_col, trial_col) {
  cv_summary <- lapply(trait_cols, function(trait) {
    cv_by_group <- raw_data %>%
      dplyr::group_by(
        id_lab = .data[[accession_col]],
        batch = .data[[trial_col]]
      ) %>%
      dplyr::summarise(
        n_samples = sum(!is.na(.data[[trait]])),
        mean_value = mean(.data[[trait]], na.rm = TRUE),
        sd_value = sd(.data[[trait]], na.rm = TRUE),
        cv_percent = ifelse(
          n_samples >= 2 & mean_value != 0,
          100 * sd_value / abs(mean_value),
          NA_real_
        ),
        .groups = "drop"
      )

    data.frame(
      Trait = trait,
      Groups = nrow(cv_by_group),
      Groups_with_CV = sum(!is.na(cv_by_group$cv_percent)),
      Groups_lt_3_obs = sum(cv_by_group$n_samples < 3),
      Mean_CV_Pct = round(mean(cv_by_group$cv_percent, na.rm = TRUE), 3),
      Median_CV_Pct = round(median(cv_by_group$cv_percent, na.rm = TRUE), 3),
      SD_CV_Pct = round(sd(cv_by_group$cv_percent, na.rm = TRUE), 3),
      Q90_CV_Pct = round(as.numeric(quantile(cv_by_group$cv_percent, 0.90, na.rm = TRUE)), 3),
      stringsAsFactors = FALSE
    )
  })

  dplyr::bind_rows(cv_summary)
}


#' Compute BLUEs, Variance Components, and Average LSD Using Mixed Models
#'
#' For each trait, fits a mixed model with accession as fixed effect and trial
#' as random effect. Extracts BLUEs (best linear unbiased estimates) via emmeans,
#' computes variance components from a random-effects-only model, and calculates
#' the average LSD for pairwise comparisons.
#'
#' @param raw_data data.frame. The uploaded raw data with multiple observations
#'   per accession across trials.
#' @param trait_cols Character vector. Names of numerical trait columns.
#' @param accession_col Character. Name of the accession ID column.
#' @param trial_col Character. Name of the trial/run ID column.
#' @return A named list with:
#'   - blues: data.frame with one row per accession, columns for each trait BLUE
#'   - variance_components: data.frame with variance sources per trait
#'   - avg_lsd: named numeric vector of average LSD per trait
#'   - errors: character vector of traits that failed model fitting
compute_blues <- function(raw_data, trait_cols, accession_col, trial_col) {
  # Ensure factors
  raw_data[[accession_col]] <- factor(raw_data[[accession_col]])
  raw_data[[trial_col]] <- factor(raw_data[[trial_col]])

  blues_list <- list()
  variance_list <- list()
  lsd_values <- c()
  errors <- character(0)

  for (trait in trait_cols) {
    tryCatch({
      # Skip if trait has too few non-NA values
      non_na_count <- sum(!is.na(raw_data[[trait]]))
      if (non_na_count < 5) {
        errors <- c(errors, paste0(trait, " (insufficient data)"))
        next
      }

      # Model 1: random-effects only for variance components
      formula_vc <- as.formula(paste0("`", trait, "` ~ (1 | `", accession_col, "`) + (1 | `", trial_col, "`)"))
      model_vc <- lme4::lmer(formula_vc, data = raw_data, REML = TRUE)

      vc <- as.data.frame(lme4::VarCorr(model_vc))
      total_var <- sum(vc$vcov)
      vc_df <- data.frame(
        Trait = trait,
        Source = vc$grp,
        Variance = round(vc$vcov, 4),
        Pct_of_Total = round(100 * vc$vcov / total_var, 2),
        stringsAsFactors = FALSE
      )
      variance_list[[trait]] <- vc_df

      # Model 2: accession fixed, trial random for BLUEs
      formula_blue <- as.formula(paste0("`", trait, "` ~ `", accession_col, "` + (1 | `", trial_col, "`)"))
      model_blue <- lme4::lmer(formula_blue, data = raw_data, REML = TRUE)

      # Extract BLUEs via emmeans
      blues_emm <- emmeans::emmeans(model_blue, as.formula(paste0("~ `", accession_col, "`")), lmer.df = "asymptotic")
      blues_df <- as.data.frame(summary(blues_emm))

      # Compute average LSD
      vcov_mat <- vcov(blues_emm)
      pairwise_se <- sqrt(outer(diag(vcov_mat), diag(vcov_mat), "+") - 2 * as.matrix(vcov_mat))
      pairwise_se_vals <- pairwise_se[upper.tri(pairwise_se)]
      df_val <- summary(blues_emm)$df[1]
      avg_lsd <- signif(mean(qt(0.975, df = df_val) * pairwise_se_vals, na.rm = TRUE), 4)
      lsd_values[trait] <- avg_lsd

      # Store BLUEs
      blue_result <- data.frame(
        accession = blues_df[[accession_col]],
        value = blues_df$emmean,
        se = blues_df$SE,
        stringsAsFactors = FALSE
      )
      names(blue_result) <- c(accession_col, trait, paste0(trait, "_SE"))
      blues_list[[trait]] <- blue_result

    }, error = function(e) {
      errors <<- c(errors, paste0(trait, " (", e$message, ")"))
    })
  }

  # Merge all BLUE estimates into one data.frame
  if (length(blues_list) > 0) {
    blues_merged <- Reduce(function(x, y) merge(x, y, by = accession_col, all = TRUE), blues_list)
  } else {
    blues_merged <- data.frame()
  }

  # Combine variance components
  variance_df <- if (length(variance_list) > 0) dplyr::bind_rows(variance_list) else data.frame()

  list(
    blues = blues_merged,
    variance_components = variance_df,
    avg_lsd = lsd_values,
    errors = errors
  )
}


#' Export Results as a Multi-Sheet Excel Workbook
#'
#' Generates a multi-sheet Excel (.xlsx) workbook containing selection approach
#' results, all accession estimates, configuration settings, and diagnostic
#' metrics. Each approach sheet is sorted by Overall_Rank ascending (rank 1 first).
#'
#' @param approach_results Named list of 4 data.frames (cat_1a, cat_1b, cat_2, cat_3).
#'   Each data.frame contains the filtered and ranked accessions for that approach.
#' @param all_estimates data.frame. Complete accession estimates with all computed fields.
#' @param config list. Configuration settings for the config sheet. Contains:
#'   - ch4_low_threshold, ch4_medium_threshold, tddm_high_cutoff,
#'   - reference_baseline_pct, reference_control_id,
#'   - weights (list of 4 approaches each with w1, w2)
#'   - reference_baseline_avg (computed value)
#' @param diagnostics list. Diagnostic summary for the diagnostics sheet. Contains:
#'   - n_records, n_accessions, n_trials,
#'   - obs_per_accession (list with min, median, max),
#'   - missing_per_trait, sd_per_trait,
#'   - control_in_trials, control_missing_trials, low_replication_count
#' @param mapping list. Column mapping with accession_col, group_cols, trait_cols.
#' @return Character. Path to the temporary .xlsx file.
export_workbook <- function(approach_results, all_estimates, config, diagnostics, mapping = NULL) {
  # Create a new workbook

  wb <- openxlsx2::wb_workbook()

  # Helper to reorder columns to match Excel dashboard layout
  reorder_cols <- function(app_df) {
    if (is.null(mapping)) return(app_df)
    desired_order <- c(
      mapping$accession_col,
      mapping$group_cols,
      mapping$trait_cols,
      "CH4_Decrease_Pct",
      "Benchmark_vs_Control",
      "CH4_Category",
      "TDDM_Category",
      "Z_Score_CH4",
      "Z_Score_TDDM",
      "Composite_Score",
      "Overall_Rank",
      "CH4_8h_Reduction_Pct",
      "CH4_24h_Reduction_Pct"
    )
    ordered_cols <- desired_order[desired_order %in% names(app_df)]
    remaining <- setdiff(names(app_df), ordered_cols)
    app_df[, c(ordered_cols, remaining), drop = FALSE]
  }

  # --- Approach result sheets ---
  # Map internal names to display sheet names
  sheet_names <- c(
    cat_1a = "Category_1a",
    cat_1b = "Category_1b",
    cat_2  = "Category_2",
    cat_3  = "Category_3"
  )

  for (app_name in names(sheet_names)) {
    sheet_label <- sheet_names[[app_name]]
    app_df <- approach_results[[app_name]]

    # Sort by Overall_Rank ascending (best rank first)
    if (!is.null(app_df) && nrow(app_df) > 0 && "Overall_Rank" %in% names(app_df)) {
      app_df <- app_df[order(app_df$Overall_Rank), , drop = FALSE]
    }

    # Remove intermediate per-approach columns
    drop_pattern <- "^(Composite_cat_|Rank_cat_|Z_CH4_cat_|Z_TDDM_cat_)"
    drop_cols <- grep(drop_pattern, names(app_df), value = TRUE)
    if (length(drop_cols) > 0) {
      app_df <- app_df[, !(names(app_df) %in% drop_cols), drop = FALSE]
    }

    # Reorder columns to match Excel dashboard layout
    app_df <- reorder_cols(app_df)

    wb$add_worksheet(sheet = sheet_label)
    wb$add_data(sheet = sheet_label, x = app_df)
  }

  # --- Configuration sheet ---
  # Build a 2-column data.frame summarising configuration parameters
  config_params <- data.frame(
    Parameter = character(0),
    Value = character(0),
    stringsAsFactors = FALSE
  )

  config_params <- rbind(config_params, data.frame(
    Parameter = "CH4 Low Threshold (%)",
    Value = as.character(config$ch4_low_threshold),
    stringsAsFactors = FALSE
  ))
  config_params <- rbind(config_params, data.frame(
    Parameter = "CH4 Medium Threshold (%)",
    Value = as.character(config$ch4_medium_threshold),
    stringsAsFactors = FALSE
  ))
  config_params <- rbind(config_params, data.frame(
    Parameter = "TDDM High Cutoff (mg/g)",
    Value = as.character(config$tddm_high_cutoff),
    stringsAsFactors = FALSE
  ))
  config_params <- rbind(config_params, data.frame(
    Parameter = "Reference Baseline (%)",
    Value = as.character(config$reference_baseline_pct),
    stringsAsFactors = FALSE
  ))
  config_params <- rbind(config_params, data.frame(
    Parameter = "Reference Control ID",
    Value = as.character(config$reference_control_id),
    stringsAsFactors = FALSE
  ))

  # Reference baseline average (computed value)
  if (!is.null(config$reference_baseline_avg)) {
    config_params <- rbind(config_params, data.frame(
      Parameter = "Reference Baseline Average (computed)",
      Value = as.character(round(config$reference_baseline_avg, 4)),
      stringsAsFactors = FALSE
    ))
  }

  # Weights per approach
  weight_names <- c("cat_1a", "cat_1b", "cat_2", "cat_3")
  for (wn in weight_names) {
    if (!is.null(config$weights[[wn]])) {
      config_params <- rbind(config_params, data.frame(
        Parameter = paste0("Weight ", wn, " W1 (CH4)"),
        Value = as.character(config$weights[[wn]]$w1),
        stringsAsFactors = FALSE
      ))
      config_params <- rbind(config_params, data.frame(
        Parameter = paste0("Weight ", wn, " W2 (TDDM)"),
        Value = as.character(config$weights[[wn]]$w2),
        stringsAsFactors = FALSE
      ))
    }
  }

  wb$add_worksheet(sheet = "Configuration")
  wb$add_data(sheet = "Configuration", x = config_params)

  # --- Diagnostics sheet ---
  # Build a 2-column data.frame summarising diagnostic metrics
  diag_params <- data.frame(
    Metric = character(0),
    Value = character(0),
    stringsAsFactors = FALSE
  )

  diag_params <- rbind(diag_params, data.frame(
    Metric = "Total Records (Raw Data)",
    Value = as.character(diagnostics$n_records),
    stringsAsFactors = FALSE
  ))
  diag_params <- rbind(diag_params, data.frame(
    Metric = "Unique Accessions",
    Value = as.character(diagnostics$n_accessions),
    stringsAsFactors = FALSE
  ))
  diag_params <- rbind(diag_params, data.frame(
    Metric = "Unique Trials/Runs",
    Value = as.character(diagnostics$n_trials),
    stringsAsFactors = FALSE
  ))

  # Observations per accession summary
  if (!is.null(diagnostics$obs_per_accession)) {
    diag_params <- rbind(diag_params, data.frame(
      Metric = "Obs per Accession (min)",
      Value = as.character(diagnostics$obs_per_accession$min),
      stringsAsFactors = FALSE
    ))
    diag_params <- rbind(diag_params, data.frame(
      Metric = "Obs per Accession (median)",
      Value = as.character(diagnostics$obs_per_accession$median),
      stringsAsFactors = FALSE
    ))
    diag_params <- rbind(diag_params, data.frame(
      Metric = "Obs per Accession (max)",
      Value = as.character(diagnostics$obs_per_accession$max),
      stringsAsFactors = FALSE
    ))
  }

  # Missing values per trait
  if (!is.null(diagnostics$missing_per_trait) && length(diagnostics$missing_per_trait) > 0) {
    for (trait_name in names(diagnostics$missing_per_trait)) {
      diag_params <- rbind(diag_params, data.frame(
        Metric = paste0("Missing Values: ", trait_name),
        Value = as.character(diagnostics$missing_per_trait[[trait_name]]),
        stringsAsFactors = FALSE
      ))
    }
  }

  # Standard deviation per trait
  if (!is.null(diagnostics$sd_per_trait) && length(diagnostics$sd_per_trait) > 0) {
    for (trait_name in names(diagnostics$sd_per_trait)) {
      diag_params <- rbind(diag_params, data.frame(
        Metric = paste0("SD (Estimates): ", trait_name),
        Value = as.character(round(diagnostics$sd_per_trait[[trait_name]], 4)),
        stringsAsFactors = FALSE
      ))
    }
  }

  # Control verification
  if (!is.null(diagnostics$control_missing_trials) && length(diagnostics$control_missing_trials) > 0) {
    diag_params <- rbind(diag_params, data.frame(
      Metric = "Control Missing from Trials",
      Value = paste(diagnostics$control_missing_trials, collapse = ", "),
      stringsAsFactors = FALSE
    ))
  } else {
    diag_params <- rbind(diag_params, data.frame(
      Metric = "Control Missing from Trials",
      Value = "None (present in all trials)",
      stringsAsFactors = FALSE
    ))
  }

  # Low replication count
  diag_params <- rbind(diag_params, data.frame(
    Metric = "Accessions with < 2 Observations",
    Value = as.character(diagnostics$low_replication_count),
    stringsAsFactors = FALSE
  ))

  wb$add_worksheet(sheet = "Diagnostics")
  wb$add_data(sheet = "Diagnostics", x = diag_params)

  # --- Save workbook to temporary file ---
  tmp_path <- tempfile(fileext = ".xlsx")
  wb$save(file = tmp_path)

  return(tmp_path)
}

#' Create Interactive Scatter Plot with Threshold Lines and Sector Annotations
#'
#' Generates a scatter plot of TDDM (x-axis) versus CH4_Intensity (y-axis) with
#' configurable threshold lines that divide the plot into six sectors. Computes
#' accession counts per sector and annotates them on the plot. Returns both an
#' interactive plotly version with tooltips and a static ggExtra version with
#' marginal histograms.
#'
#' @param estimates data.frame. Accession-level estimates containing at minimum the
#'   CH4_Intensity and TDDM columns.
#' @param config list. Configuration parameters:
#'   - ch4_col: character. Name of CH4 Intensity column
#'   - tddm_col: character. Name of TDDM column
#'   - accession_col: character. Name of accession ID column
#'   - tddm_high_cutoff: numeric (mg/g) - vertical threshold line
#'   - reference_baseline_avg: numeric - the computed reference baseline average
#'   - ch4_low_threshold: numeric (percentage) - for computing absolute CH4 cutoff
#'   - ch4_medium_threshold: numeric (percentage) - for computing absolute CH4 cutoff
#' @param colour_var Character or NULL. Column name for point colour mapping.
#' @param size_var Character or NULL. Column name for point size mapping.
#' @return A named list with:
#'   - plotly_scatter: a plotly object (interactive scatter with tooltip)
#'   - static_marginal: a ggplot/ggExtra object (static scatter with marginal histograms)
plot_scatter_with_thresholds <- function(estimates, config, colour_var = NULL, size_var = NULL) {
  # Extract configuration parameters
  ch4_col <- config$ch4_col
  tddm_col <- config$tddm_col
  accession_col <- config$accession_col
  tddm_high_cutoff <- config$tddm_high_cutoff
  reference_baseline_avg <- config$reference_baseline_avg
  ch4_low_threshold <- config$ch4_low_threshold
  ch4_medium_threshold <- config$ch4_medium_threshold

  # --- Step 1: Compute absolute CH4 cutoffs ---
  ch4_low_cutoff <- reference_baseline_avg * (1 - ch4_low_threshold / 100)
  ch4_medium_cutoff <- reference_baseline_avg * (1 - ch4_medium_threshold / 100)

  # --- Step 2: Prepare plot data ---
  plot_df <- data.frame(
    accession_id = estimates[[accession_col]],
    tddm = estimates[[tddm_col]],
    ch4_intensity = estimates[[ch4_col]],
    stringsAsFactors = FALSE
  )

  # Add colour variable if specified and present in data
  has_colour <- !is.null(colour_var) && colour_var %in% names(estimates)
  if (has_colour) {
    plot_df$colour_var <- estimates[[colour_var]]
  }

  # Add size variable if specified and present in data
  has_size <- !is.null(size_var) && size_var %in% names(estimates)
  if (has_size) {
    plot_df$size_var <- estimates[[size_var]]
  }

  # Build tooltip text
  tooltip_text <- paste0(
    "Accession: ", plot_df$accession_id, "\n",
    "TDDM: ", round(plot_df$tddm, 2), "\n",
    "CH4 Intensity: ", round(plot_df$ch4_intensity, 2)
  )
  if (has_colour) {
    tooltip_text <- paste0(tooltip_text, "\n", colour_var, ": ", plot_df$colour_var)
  }

  plot_df$tooltip_text <- tooltip_text

  # --- Step 3: Build base ggplot with appropriate aesthetics ---
  if (has_colour && has_size) {
    p <- ggplot(plot_df, aes(x = tddm, y = ch4_intensity, colour = colour_var,
                             size = size_var, text = tooltip_text))
  } else if (has_colour) {
    p <- ggplot(plot_df, aes(x = tddm, y = ch4_intensity, colour = colour_var,
                             text = tooltip_text))
  } else if (has_size) {
    p <- ggplot(plot_df, aes(x = tddm, y = ch4_intensity, size = size_var,
                             text = tooltip_text))
  } else {
    p <- ggplot(plot_df, aes(x = tddm, y = ch4_intensity, text = tooltip_text))
  }

  p <- p + geom_point(alpha = 0.7)

  # --- Step 4: Add threshold lines ---
  p <- p +
    geom_vline(xintercept = tddm_high_cutoff, linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = ch4_low_cutoff, linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = ch4_medium_cutoff, linetype = "dashed", colour = "grey40")

  # --- Step 5: Compute sector counts ---
  # 6 sectors created by 1 vertical + 2 horizontal lines:
  # Columns: TDDM < cutoff (left) vs TDDM >= cutoff (right)
  # Rows: CH4 >= medium_cutoff (top/high intensity),
  #        low_cutoff <= CH4 < medium_cutoff (middle),
  #        CH4 < low_cutoff (bottom/low intensity = high decrease = good)
  tddm_vals <- plot_df$tddm
  ch4_vals <- plot_df$ch4_intensity

  count_top_left <- sum(tddm_vals < tddm_high_cutoff & ch4_vals >= ch4_medium_cutoff, na.rm = TRUE)
  count_top_right <- sum(tddm_vals >= tddm_high_cutoff & ch4_vals >= ch4_medium_cutoff, na.rm = TRUE)
  count_mid_left <- sum(tddm_vals < tddm_high_cutoff & ch4_vals >= ch4_low_cutoff & ch4_vals < ch4_medium_cutoff, na.rm = TRUE)
  count_mid_right <- sum(tddm_vals >= tddm_high_cutoff & ch4_vals >= ch4_low_cutoff & ch4_vals < ch4_medium_cutoff, na.rm = TRUE)
  count_bot_left <- sum(tddm_vals < tddm_high_cutoff & ch4_vals < ch4_low_cutoff, na.rm = TRUE)
  count_bot_right <- sum(tddm_vals >= tddm_high_cutoff & ch4_vals < ch4_low_cutoff, na.rm = TRUE)

  # Compute annotation positions (center of each sector region)
  tddm_range <- range(tddm_vals, na.rm = TRUE)
  ch4_range <- range(ch4_vals, na.rm = TRUE)

  # X positions: center of left sector, center of right sector
  x_left <- mean(c(tddm_range[1], tddm_high_cutoff))
  x_right <- mean(c(tddm_high_cutoff, tddm_range[2]))

  # Y positions: center of top, middle, bottom sectors
  y_top <- mean(c(ch4_medium_cutoff, ch4_range[2]))
  y_mid <- mean(c(ch4_low_cutoff, ch4_medium_cutoff))
  y_bot <- mean(c(ch4_range[1], ch4_low_cutoff))

  # --- Step 6: Add sector count annotations ---
  p <- p +
    annotate("text", x = x_left, y = y_top, label = paste0("n=", count_top_left),
             size = 3.5, fontface = "bold", colour = "grey30") +
    annotate("text", x = x_right, y = y_top, label = paste0("n=", count_top_right),
             size = 3.5, fontface = "bold", colour = "grey30") +
    annotate("text", x = x_left, y = y_mid, label = paste0("n=", count_mid_left),
             size = 3.5, fontface = "bold", colour = "grey30") +
    annotate("text", x = x_right, y = y_mid, label = paste0("n=", count_mid_right),
             size = 3.5, fontface = "bold", colour = "grey30") +
    annotate("text", x = x_left, y = y_bot, label = paste0("n=", count_bot_left),
             size = 3.5, fontface = "bold", colour = "grey30") +
    annotate("text", x = x_right, y = y_bot, label = paste0("n=", count_bot_right),
             size = 3.5, fontface = "bold", colour = "grey30")

  # --- Step 7: Add labels and theme ---
  colour_label <- if (has_colour) colour_var else NULL
  size_label <- if (has_size) size_var else NULL

  p <- p +
    labs(
      x = "TDDM (mg/g)",
      y = "CH4 Intensity (ml CH4/g TDDM)",
      colour = colour_label,
      size = size_label
    ) +
    theme_minimal()

  # --- Step 8: Convert to plotly for interactive version ---
  plotly_scatter <- plotly::ggplotly(p, tooltip = "text")

  # --- Step 9: Create static version with marginal density plots ---
  if (has_colour && has_size) {
    p_main <- ggplot(plot_df, aes(x = tddm, y = ch4_intensity,
                                  colour = colour_var, size = size_var))
  } else if (has_colour) {
    p_main <- ggplot(plot_df, aes(x = tddm, y = ch4_intensity,
                                  colour = colour_var))
  } else if (has_size) {
    p_main <- ggplot(plot_df, aes(x = tddm, y = ch4_intensity,
                                  size = size_var))
  } else {
    p_main <- ggplot(plot_df, aes(x = tddm, y = ch4_intensity))
  }

  p_main <- p_main +
    geom_point(alpha = 0.7) +
    geom_vline(xintercept = tddm_high_cutoff, linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = ch4_low_cutoff, linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = ch4_medium_cutoff, linetype = "dashed", colour = "grey40") +
    annotate("text", x = x_left, y = y_top, label = paste0("n=", count_top_left),
             size = 3.5, fontface = "bold", colour = "grey30") +
    annotate("text", x = x_right, y = y_top, label = paste0("n=", count_top_right),
             size = 3.5, fontface = "bold", colour = "grey30") +
    annotate("text", x = x_left, y = y_mid, label = paste0("n=", count_mid_left),
             size = 3.5, fontface = "bold", colour = "grey30") +
    annotate("text", x = x_right, y = y_mid, label = paste0("n=", count_mid_right),
             size = 3.5, fontface = "bold", colour = "grey30") +
    annotate("text", x = x_left, y = y_bot, label = paste0("n=", count_bot_left),
             size = 3.5, fontface = "bold", colour = "grey30") +
    annotate("text", x = x_right, y = y_bot, label = paste0("n=", count_bot_right),
             size = 3.5, fontface = "bold", colour = "grey30") +
    labs(
      x = "TDDM (mg/g)",
      y = "CH4 Intensity (ml CH4/g TDDM)",
      colour = colour_label,
      size = size_label
    ) +
    theme_minimal()

  # Add marginal density plots
  static_plot <- ggExtra::ggMarginal(p_main, type = "density",
                                     fill = "steelblue", alpha = 0.4)

  # --- Return plot objects ---
  list(
    plotly_scatter = plotly_scatter,
    static_plot = static_plot
  )
}


# --- UI -----------------------------------------------------------------------

ui <- page_sidebar(
  title = "LMF Selection Index Dashboard - ver 0.1.0 (Beta)",
  theme = bs_theme(version = 5, bootswatch = "sandstone"),
  shinyjs::useShinyjs(),

  sidebar = sidebar(
    width = 350,

    # Tighten spacing in sidebar
    tags$style(HTML("
      .sidebar .form-group { margin-bottom: 8px; }
      .sidebar .row { margin-bottom: 0; }
      .sidebar hr { margin-top: 8px; margin-bottom: 8px; }
    ")),

    # === File Upload Section ===
    tags$h5("Data Upload"),
    fileInput("file_upload", "Upload Data File",
              accept = c(".csv", ".tsv", ".xlsx", ".xls")),
    tags$hr(),

    # === Column Mapping Section (shown after file upload, collapsible) ===
    conditionalPanel(
      condition = "output.file_uploaded",
      tags$div(
        id = "mapping_section",
        tags$h5(
          tags$a(
            tags$span(id = "mapping_icon", class = "me-1", "\u25BC"),
            "Column Mapping",
            href = "#",
            onclick = "Shiny.setInputValue('toggle_mapping', Math.random())",
            style = "text-decoration: none; color: inherit; cursor: pointer;"
          )
        ),
        tags$div(
          id = "mapping_panel",
          selectInput("accession_col", "Accession ID Column", choices = NULL),
          selectInput("trial_col", "Trial/Run ID Column", choices = NULL),
          selectInput("ch4_col_select", "CH4 Intensity Column", choices = NULL),
          selectInput("tddm_col_select", "TDDM Column", choices = NULL),
          selectInput("ch4_8h_col_select", "CH4 8h Column (optional)",
                      choices = c("(none)" = "")),
          selectInput("ch4_24h_col_select", "CH4 24h Column (optional)",
                      choices = c("(none)" = "")),
          selectizeInput("group_cols", "Grouping Columns (categorical)",
                         choices = NULL, multiple = TRUE),
          selectizeInput("aux_cols", "Auxiliary Columns (numeric, for viz only)",
                         choices = NULL, multiple = TRUE),
          actionButton("confirm_mapping", "Confirm Mapping",
                       class = "btn-info", width = "100%")
        )
      ),
      tags$hr()
    ),

    # === Action Buttons (right below column mapping) ===
    actionButton("run_analysis", "Run Analysis",
                 class = "btn-primary", width = "100%"),
    tags$div(style = "margin-top: 6px;",
      shinyjs::disabled(
        downloadButton("download_results", "Download Excel",
                       class = "btn-success", style = "width: 100%;")
      )
    ),

    tags$hr(),

    # === Configuration Section (collapsible) ===
    tags$h5(
      tags$a(
        tags$span(id = "config_icon", class = "me-1", "\u25B6"),
        "Configuration",
        href = "#",
        onclick = "Shiny.setInputValue('toggle_config', Math.random())",
        style = "text-decoration: none; color: inherit; cursor: pointer;"
      )
    ),
    tags$div(
      id = "config_panel",
      style = "display: none;",

      # Aggregation Method
      tags$strong("Aggregation Method"),
      radioButtons("agg_method", NULL,
                   choices = c("Arithmetic Mean" = "mean", "BLUEs (mixed model)" = "blues"),
                   selected = "blues", inline = TRUE),
      tags$hr(),

      # CH4 Thresholds
      tags$strong("CH4 Intensity Thresholds"),
      numericInput("ch4_low_threshold", "CH4 Low Threshold (%)",
                   value = 35, min = 0.1, max = 99.9, step = 0.1),
      numericInput("ch4_medium_threshold", "CH4 Medium Threshold (%)",
                   value = 20, min = 0.1, max = 99.9, step = 0.1),

      # TDDM Cutoff
      tags$strong("TDDM Classification"),
      numericInput("tddm_high_cutoff", "TDDM High Cutoff (mg/g)",
                   value = 60, min = 1, max = 999, step = 1),

      # Reference Baseline
      tags$strong("Reference Settings"),
      numericInput("reference_baseline_pct", "Reference Baseline (%)",
                   value = 10, min = 1, max = 100, step = 1),
      selectizeInput("reference_control_id", "Reference Control ID",
                choices = NULL, multiple = FALSE,
                options = list(placeholder = "Select after mapping columns")),

      # Computed values display
      tags$div(
        style = "background-color: #f8f9fa; padding: 8px; border-radius: 4px; margin-bottom: 10px;",
        tags$small(tags$strong("Computed Values:")),
        textOutput("computed_baseline_avg"),
        textOutput("computed_ch4_cutoffs")
      ),

      tags$hr(),

      # Weights Section
      tags$strong("Score Weights (per approach)"),

      # Category 1a weights
      tags$div(style = "margin-top: 6px; margin-bottom: 2px;",
        tags$em(tags$small("Category 1a (Low CH4 & High TDDM)"))
      ),
      fluidRow(
        column(6, numericInput("w1_cat1a", "W1 (CH4)", value = 2, min = 0, max = 10, step = 1)),
        column(6, numericInput("w2_cat1a", "W2 (TDDM)", value = 1, min = 0, max = 10, step = 1))
      ),

      # Category 1b weights
      tags$div(style = "margin-bottom: 2px;",
        tags$em(tags$small("Category 1b (Medium CH4 & High TDDM)"))
      ),
      fluidRow(
        column(6, numericInput("w1_cat1b", "W1 (CH4)", value = 2, min = 0, max = 10, step = 1)),
        column(6, numericInput("w2_cat1b", "W2 (TDDM)", value = 1, min = 0, max = 10, step = 1))
      ),

      # Category 2 weights
      tags$div(style = "margin-bottom: 2px;",
        tags$em(tags$small("Category 2 (Not High CH4 & Not High TDDM)"))
      ),
      fluidRow(
        column(6, numericInput("w1_cat2", "W1 (CH4)", value = 1, min = 0, max = 10, step = 1)),
        column(6, numericInput("w2_cat2", "W2 (TDDM)", value = 0, min = 0, max = 10, step = 1))
      ),

      # Category 3 weights
      tags$div(style = "margin-bottom: 2px;",
        tags$em(tags$small("Category 3 (High CH4 & High TDDM)"))
      ),
      fluidRow(
        column(6, numericInput("w1_cat3", "W1 (CH4)", value = 1, min = 0, max = 10, step = 1)),
        column(6, numericInput("w2_cat3", "W2 (TDDM)", value = 2, min = 0, max = 10, step = 1))
      )
    )
  ),

  # === Main Panel Content (tabbed layout) ===
  navset_card_pill(
    id = "main_tabs",

    # --- Data Preview Tab ---
    nav_panel(
      "Data Preview",
      card_body(
        textOutput("upload_status"),
        DT::dataTableOutput("data_preview")
      )
    ),

    # --- Selection Results Tab ---
    nav_panel(
      "Selection Results",
      navset_card_tab(
        nav_panel(
          "Category 1a",
          textOutput("count_cat1a"),
          DT::dataTableOutput("table_cat1a")
        ),
        nav_panel(
          "Category 1b",
          textOutput("count_cat1b"),
          DT::dataTableOutput("table_cat1b")
        ),
        nav_panel(
          "Category 2",
          textOutput("count_cat2"),
          DT::dataTableOutput("table_cat2")
        ),
        nav_panel(
          "Category 3",
          textOutput("count_cat3"),
          DT::dataTableOutput("table_cat3")
        )
      )
    ),

    # --- Scatter Plot Tab ---
    nav_panel(
      "Scatter Plot",
      card_body(
        fluidRow(
          column(4, selectInput("colour_var", "Colour Variable",
                                choices = c("(none)" = ""))),
          column(4, selectInput("size_var", "Size Variable",
                                choices = c("(none)" = ""))),
          column(4, tags$div(style = "margin-top: 25px;",
            downloadButton("download_plot", "Download Plot",
                           class = "btn-outline-primary", style = "width: 100%;")
          ))
        ),
        plotOutput("scatter_static", height = "600px")
      )
    ),

    # --- Diagnostics Tab ---
    nav_panel(
      "Diagnostics",
      card_body(
        uiOutput("diagnostics_panel")
      )
    )
  )
)

# --- Server -------------------------------------------------------------------

server <- function(input, output, session) {

  # Flag to indicate whether a file has been uploaded (used by conditionalPanel)
  output$file_uploaded <- reactive({
    !is.null(input$file_upload)
  })
  outputOptions(output, "file_uploaded", suspendWhenHidden = FALSE)

  # --- File Upload Reactive and Preview (Task 12.1) ---------------------------


  # ReactiveVal to store parsed data.frame or NULL
  raw_data <- reactiveVal(NULL)
  upload_message <- reactiveVal("")

  # Observe file input change, parse the uploaded file

  observeEvent(input$file_upload, {
    req(input$file_upload)
    file_info <- input$file_upload
    result <- parse_upload(file_info$datapath, file_info$name)

    if (is.list(result) && !is.null(result$error)) {
      raw_data(NULL)
      upload_message(result$error)
    } else {
      raw_data(result)
      # Build success message
      ext <- tolower(tools::file_ext(file_info$name))
      msg <- paste0("Successfully loaded: ", file_info$name,
                    " (", nrow(result), " rows, ", ncol(result), " columns)")
      if (ext %in% c("xlsx", "xls")) {
        msg <- paste0(msg, " \u2014 Sheet: ", readxl::excel_sheets(file_info$datapath)[1])
      }
      upload_message(msg)
    }
  })

  # Render upload status message
  output$upload_status <- renderText({
    upload_message()
  })

  # Render data preview table (first 10 rows)
  output$data_preview <- DT::renderDataTable({
    req(raw_data())
    DT::datatable(raw_data(),
                  options = list(scrollX = TRUE, pageLength = 10),
                  rownames = FALSE)
  })

  # Update column mapping dropdowns when raw_data changes
  observeEvent(raw_data(), {
    req(raw_data())
    cols <- names(raw_data())
    # All dropdowns start blank — user must explicitly select
    updateSelectInput(session, "accession_col",
                      choices = c("-- Select --" = "", cols), selected = "")
    updateSelectInput(session, "trial_col",
                      choices = c("-- Select --" = "", cols), selected = "")
    updateSelectInput(session, "ch4_col_select",
                      choices = c("-- Select --" = "", cols),
                      selected = grep("ch4.*int|ch4_int", cols, ignore.case = TRUE, value = TRUE)[1])
    updateSelectInput(session, "tddm_col_select",
                      choices = c("-- Select --" = "", cols),
                      selected = grep("tddm|tdm", cols, ignore.case = TRUE, value = TRUE)[1])
    updateSelectInput(session, "ch4_8h_col_select",
                      choices = c("-- Select --" = "", cols),
                      selected = grep("ch4.*8h|8h", cols, ignore.case = TRUE, value = TRUE)[1])
    updateSelectInput(session, "ch4_24h_col_select",
                      choices = c("-- Select --" = "", cols),
                      selected = grep("ch4.*24h|24h", cols, ignore.case = TRUE, value = TRUE)[1])
    updateSelectizeInput(session, "group_cols", choices = cols)
    updateSelectizeInput(session, "aux_cols", choices = cols)
  })

  # --- Column Mapping Reactive with Validation (Task 12.2) --------------------

  # ReactiveVal to store confirmed column mapping or NULL
  column_mapping <- reactiveVal(NULL)
  mapping_message <- reactiveVal("")

  # Observe "Confirm Mapping" button: validate and store the mapping
  observeEvent(input$confirm_mapping, {
    req(raw_data())

    # Get selected columns
    acc_col <- input$accession_col
    trial_col <- input$trial_col
    ch4_col <- input$ch4_col_select
    tddm_col <- input$tddm_col_select
    ch4_8h_col <- input$ch4_8h_col_select
    ch4_24h_col <- input$ch4_24h_col_select
    group_cols <- input$group_cols
    aux_cols <- input$aux_cols

    # Build trait_cols from the explicit selectors
    trait_cols <- c(ch4_col, tddm_col)
    if (!is.null(ch4_8h_col) && ch4_8h_col != "") trait_cols <- c(trait_cols, ch4_8h_col)
    if (!is.null(ch4_24h_col) && ch4_24h_col != "") trait_cols <- c(trait_cols, ch4_24h_col)

    # Validation: accession ID required
    if (is.null(acc_col) || acc_col == "") {
      mapping_message("Please map exactly one column as the Accession ID.")
      column_mapping(NULL)
      return()
    }

    # Validation: trial ID required
    if (is.null(trial_col) || trial_col == "") {
      mapping_message("Please map exactly one column as the Trial/Run ID.")
      column_mapping(NULL)
      return()
    }

    # Validation: CH4 and TDDM required
    if (is.null(ch4_col) || ch4_col == "") {
      mapping_message("Please select a CH4 Intensity column.")
      column_mapping(NULL)
      return()
    }
    if (is.null(tddm_col) || tddm_col == "") {
      mapping_message("Please select a TDDM column.")
      column_mapping(NULL)
      return()
    }

    # Validation: no duplicate assignments
    all_mapped <- c(acc_col, trial_col, trait_cols, group_cols, aux_cols)
    all_mapped <- all_mapped[all_mapped != ""]
    if (any(duplicated(all_mapped))) {
      dup_col <- all_mapped[duplicated(all_mapped)][1]
      mapping_message(paste0("Column '", dup_col, "' is already assigned to another role."))
      column_mapping(NULL)
      return()
    }

    # Validation: numeric proportion >= 50% for trait columns
    rejected_traits <- character(0)
    valid_traits <- character(0)
    for (trait in trait_cols) {
      if (validate_numeric_column(raw_data()[[trait]])) {
        valid_traits <- c(valid_traits, trait)
      } else {
        rejected_traits <- c(rejected_traits, trait)
      }
    }

    # Also validate auxiliary columns for numeric proportion
    valid_aux <- character(0)
    rejected_aux <- character(0)
    if (!is.null(aux_cols) && length(aux_cols) > 0) {
      for (acol in aux_cols) {
        if (validate_numeric_column(raw_data()[[acol]])) {
          valid_aux <- c(valid_aux, acol)
        } else {
          rejected_aux <- c(rejected_aux, acol)
        }
      }
    }

    # Build warning messages for rejected columns
    warnings <- character(0)
    if (length(rejected_traits) > 0) {
      for (rt in rejected_traits) {
        warnings <- c(warnings, paste0("Column '", rt, "' contains more than 50% non-numeric values and cannot be used as a trait."))
      }
    }
    if (length(rejected_aux) > 0) {
      for (ra in rejected_aux) {
        warnings <- c(warnings, paste0("Column '", ra, "' contains more than 50% non-numeric values and cannot be used as auxiliary."))
      }
    }

    # Check that at least CH4 and TDDM are valid
    if (!(ch4_col %in% valid_traits) || !(tddm_col %in% valid_traits)) {
      mapping_message(paste(c("CH4 Intensity and TDDM columns must contain valid numeric data.", warnings), collapse = "\n"))
      column_mapping(NULL)
      return()
    }

    # Store confirmed mapping
    column_mapping(list(
      accession_col = acc_col,
      trial_col = trial_col,
      trait_cols = valid_traits,
      group_cols = if (!is.null(group_cols) && length(group_cols) > 0) group_cols else character(0),
      aux_cols = valid_aux
    ))

    # Populate Reference Control ID dropdown with unique accession IDs
    accession_ids <- sort(unique(raw_data()[[acc_col]]))
    default_ctrl <- if ("ABC-18447-1" %in% accession_ids) "ABC-18447-1" else accession_ids[1]
    updateSelectizeInput(session, "reference_control_id",
                         choices = accession_ids,
                         selected = default_ctrl)

    # Success message (with warnings if any)
    msg <- paste0("Column mapping confirmed. ", length(valid_traits), " trait(s), ",
                  length(group_cols), " grouping, ", length(valid_aux), " auxiliary columns mapped.")
    if (length(warnings) > 0) {
      msg <- paste(c(msg, "Warnings:", warnings), collapse = "\n")
    }
    mapping_message(msg)

    # Collapse the mapping panel after successful confirmation
    shinyjs::hide("mapping_panel")
    shinyjs::runjs("document.getElementById('mapping_icon').textContent = '\u25B6';")
  })

  # Toggle mapping panel visibility on header click
  observeEvent(input$toggle_mapping, {
    shinyjs::toggle("mapping_panel")
    shinyjs::runjs("
      var icon = document.getElementById('mapping_icon');
      icon.textContent = icon.textContent === '\u25B6' ? '\u25BC' : '\u25B6';
    ")
  })

  # Toggle configuration panel visibility on header click
  observeEvent(input$toggle_config, {
    shinyjs::toggle("config_panel")
    shinyjs::runjs("
      var icon = document.getElementById('config_icon');
      icon.textContent = icon.textContent === '\u25B6' ? '\u25BC' : '\u25B6';
    ")
  })

  # Display mapping validation results as notification
  observe({
    req(mapping_message())
    msg <- mapping_message()
    if (nzchar(msg)) {
      showNotification(msg, type = if (is.null(column_mapping())) "error" else "message", duration = 8)
    }
  })

  # --- Analysis Trigger and Calculation Pipeline (Tasks 12.3, 12.4) -----------

  # Reactive values to store analysis results
  accession_estimates <- reactiveVal(NULL)

  selection_results <- reactiveVal(NULL)
  diagnostics_data <- reactiveVal(NULL)
  blues_diagnostics <- reactiveVal(NULL)
  cv_diagnostics <- reactiveVal(NULL)
  analysis_done <- reactiveVal(FALSE)

  # Helper: build config list from current inputs
  build_config <- reactive({
    mapping <- column_mapping()
    req(mapping)

    list(
      ch4_col = input$ch4_col_select,
      tddm_col = input$tddm_col_select,
      accession_col = mapping$accession_col,
      ch4_low_threshold = input$ch4_low_threshold,
      ch4_medium_threshold = input$ch4_medium_threshold,
      tddm_high_cutoff = input$tddm_high_cutoff,
      reference_baseline_pct = input$reference_baseline_pct,
      reference_control_id = input$reference_control_id,
      weights = list(
        cat_1a = list(w1 = input$w1_cat1a, w2 = input$w2_cat1a),
        cat_1b = list(w1 = input$w1_cat1b, w2 = input$w2_cat1b),
        cat_2  = list(w1 = input$w1_cat2, w2 = input$w2_cat2),
        cat_3  = list(w1 = input$w1_cat3, w2 = input$w2_cat3)
      ),
      ch4_8h_col = if (!is.null(input$ch4_8h_col_select) && input$ch4_8h_col_select != "") input$ch4_8h_col_select else NULL,
      ch4_24h_col = if (!is.null(input$ch4_24h_col_select) && input$ch4_24h_col_select != "") input$ch4_24h_col_select else NULL
    )
  })

  # Run analysis pipeline
  run_analysis <- function() {
    mapping <- column_mapping()
    data <- raw_data()
    req(mapping, data)

    # Build rename mapping: user column name -> standard name
    config <- build_config()
    rename_map <- c()
    rename_map[config$ch4_col] <- "CH4_Intensity"
    rename_map[config$tddm_col] <- "TDDM"
    if (!is.null(config$ch4_8h_col)) rename_map[config$ch4_8h_col] <- "CH4_8h"
    if (!is.null(config$ch4_24h_col)) rename_map[config$ch4_24h_col] <- "CH4_24h"
    rename_map[mapping$accession_col] <- "Accession_ID"
    rename_map[mapping$trial_col] <- "Trial_ID"

    # Rename columns in raw data before aggregation
    data_renamed <- data
    for (old_name in names(rename_map)) {
      if (old_name %in% names(data_renamed)) {
        names(data_renamed)[names(data_renamed) == old_name] <- rename_map[[old_name]]
      }
    }

    # Update mapping and config to use standard names
    std_mapping <- mapping
    std_mapping$accession_col <- "Accession_ID"
    std_mapping$trial_col <- "Trial_ID"
    std_mapping$trait_cols <- vapply(mapping$trait_cols, function(col) {
      if (col %in% names(rename_map)) rename_map[[col]] else col
    }, character(1), USE.NAMES = FALSE)

    std_config <- config
    std_config$ch4_col <- "CH4_Intensity"
    std_config$tddm_col <- "TDDM"
    std_config$accession_col <- "Accession_ID"
    if (!is.null(config$ch4_8h_col)) std_config$ch4_8h_col <- "CH4_8h"
    if (!is.null(config$ch4_24h_col)) std_config$ch4_24h_col <- "CH4_24h"

    # Step 1: Aggregate accessions (method depends on user selection)
    agg_method <- input$agg_method

    if (agg_method == "blues") {
      # Use BLUEs from mixed model
      showNotification("Computing BLUEs via mixed models... this may take a moment.", type = "message", duration = 3)
      blues_result <- compute_blues(
        raw_data = data_renamed,
        trait_cols = std_mapping$trait_cols,
        accession_col = std_mapping$accession_col,
        trial_col = std_mapping$trial_col
      )

      if (nrow(blues_result$blues) == 0) {
        showNotification("BLUEs computation failed for all traits. Falling back to arithmetic mean.", type = "warning")
        agg_method <- "mean"
      } else {
        # Build estimates from BLUEs (keep only trait columns without _SE)
        estimates <- blues_result$blues
        se_cols <- grep("_SE$", names(estimates), value = TRUE)
        estimates_clean <- estimates[, !(names(estimates) %in% se_cols), drop = FALSE]

        # Add grouping columns from arithmetic mean aggregation
        agg_result <- aggregate_accessions(
          raw_data = data_renamed,
          accession_col = std_mapping$accession_col,
          trial_col = std_mapping$trial_col,
          trait_cols = std_mapping$trait_cols,
          group_cols = std_mapping$group_cols
        )
        # Merge grouping columns into BLUE estimates
        if (length(std_mapping$group_cols) > 0) {
          group_df <- agg_result$estimates[, c(std_mapping$accession_col, std_mapping$group_cols), drop = FALSE]
          estimates_clean <- merge(estimates_clean, group_df, by = std_mapping$accession_col, all.x = TRUE)
        }

        accession_estimates(estimates_clean)

        # Store BLUE-specific diagnostics
        blues_diagnostics(blues_result)

        if (length(blues_result$errors) > 0) {
          showNotification(
            paste0("BLUEs warnings: ", paste(blues_result$errors, collapse = "; ")),
            type = "warning", duration = 8
          )
        }
      }
    }

    if (agg_method == "mean") {
      # Use arithmetic mean
      agg_result <- aggregate_accessions(
        raw_data = data_renamed,
        accession_col = std_mapping$accession_col,
        trial_col = std_mapping$trial_col,
        trait_cols = std_mapping$trait_cols,
        group_cols = std_mapping$group_cols
      )
      accession_estimates(agg_result$estimates)
      blues_diagnostics(NULL)
    }

    # Step 1b: Compute CV% (always, for diagnostics)
    n_trials <- length(unique(data_renamed[[std_mapping$trial_col]]))
    if (n_trials > 1) {
      cv_result <- calculate_trait_cv(
        raw_data = data_renamed,
        trait_cols = std_mapping$trait_cols,
        accession_col = std_mapping$accession_col,
        trial_col = std_mapping$trial_col
      )
      cv_diagnostics(cv_result)
    } else {
      cv_diagnostics(NULL)
    }

    # Step 2: Calculate selection index
    if (is.null(std_config$tddm_col)) {
      showNotification("At least 2 trait columns needed (CH4 Intensity and TDDM).", type = "error")
      return()
    }

    sel_result <- calculate_selection_index(accession_estimates(), std_config)
    selection_results(sel_result)

    # Step 3: Compute diagnostics
    diag_result <- compute_diagnostics(
      raw_data = data_renamed,
      estimates = accession_estimates(),
      mapping = std_mapping,
      config = std_config
    )
    diagnostics_data(diag_result)

    # Mark analysis as done and enable download
    analysis_done(TRUE)
    shinyjs::enable("download_results")

    showNotification("Analysis complete!", type = "message", duration = 5)
  }

  # Trigger analysis on button click
  observeEvent(input$run_analysis, {
    req(column_mapping(), raw_data())
    run_analysis()
  })

  # Recalculate on configuration changes (after initial analysis)
  observeEvent(list(
    input$ch4_col_select,
    input$tddm_col_select,
    input$ch4_8h_col_select,
    input$ch4_24h_col_select,
    input$ch4_low_threshold,
    input$ch4_medium_threshold,
    input$tddm_high_cutoff,
    input$reference_baseline_pct,
    input$reference_control_id,
    input$w1_cat1a, input$w2_cat1a,
    input$w1_cat1b, input$w2_cat1b,
    input$w1_cat2, input$w2_cat2,
    input$w1_cat3, input$w2_cat3
  ), {
    # Only recalculate if analysis has been done before
    if (isTRUE(analysis_done())) {
      run_analysis()
    }
  }, ignoreInit = TRUE)

  # Display computed reference baseline average
  output$computed_baseline_avg <- renderText({
    res <- selection_results()
    if (!is.null(res)) {
      paste0("Reference Baseline Avg: ", round(res$reference_baseline_avg, 2))
    } else {
      "Run analysis to compute"
    }
  })

  # Display computed CH4 absolute cutoffs
  output$computed_ch4_cutoffs <- renderText({
    res <- selection_results()
    if (!is.null(res)) {
      low_cutoff <- res$reference_baseline_avg * (1 - input$ch4_low_threshold / 100)
      med_cutoff <- res$reference_baseline_avg * (1 - input$ch4_medium_threshold / 100)
      paste0("Low cutoff: ", round(low_cutoff, 2), " | Medium cutoff: ", round(med_cutoff, 2))
    } else {
      ""
    }
  })

  # --- Results Tables (Task 12.5) ---------------------------------------------

  # Helper: reorder columns to match Excel dashboard layout
  reorder_result_columns <- function(app_df, mapping) {
    # Since columns are now standardized, use fixed standard names
    # Grouping columns keep their original names (not renamed)
    group_cols <- mapping$group_cols
    
    # Desired order with standard names:
    desired_order <- c(
      "Accession_ID",
      group_cols,
      "CH4_Intensity",
      "TDDM",
      "CH4_8h",
      "CH4_24h",
      "CH4_Decrease_Pct",
      "Benchmark_vs_Control",
      "CH4_Category",
      "TDDM_Category",
      "Z_Score_CH4",
      "Z_Score_TDDM",
      "Composite_Score",
      "Overall_Rank",
      "CH4_8h_Reduction_Pct",
      "CH4_24h_Reduction_Pct"
    )
    
    # Keep only columns that exist in app_df, in the desired order
    ordered_cols <- desired_order[desired_order %in% names(app_df)]
    # Append any remaining columns not in the desired list
    remaining <- setdiff(names(app_df), ordered_cols)
    final_order <- c(ordered_cols, remaining)
    
    app_df[, final_order, drop = FALSE]
  }

  # Helper to render an approach table
  render_approach_table <- function(app_name) {
    DT::renderDataTable({
      res <- selection_results()
      mapping <- column_mapping()
      req(res, mapping)
      app_df <- res$approach_results[[app_name]]
      
      if (is.null(app_df) || nrow(app_df) == 0) {
        return(DT::datatable(
          data.frame(Message = "No accessions matched the selection criteria for this approach."),
          options = list(dom = "t"), rownames = FALSE
        ))
      }
      
      # Remove intermediate per-approach columns (Composite_cat_*, Rank_cat_*)
      drop_pattern <- "^(Composite_cat_|Rank_cat_|Z_CH4_cat_|Z_TDDM_cat_)"
      drop_cols <- grep(drop_pattern, names(app_df), value = TRUE)
      if (length(drop_cols) > 0) {
        app_df <- app_df[, !(names(app_df) %in% drop_cols), drop = FALSE]
      }
      
      # Reorder columns to match Excel dashboard layout
      app_df <- reorder_result_columns(app_df, mapping)
      
      # Round numeric columns to 2 decimal places
      numeric_cols <- sapply(app_df, is.numeric)
      app_df[numeric_cols] <- lapply(app_df[numeric_cols], round, digits = 2)
      
      DT::datatable(app_df,
        options = list(
          order = list(list(which(names(app_df) == "Overall_Rank") - 1, "asc")),
          pageLength = 10,
          scrollX = TRUE
        ),
        rownames = FALSE
      )
    })
  }

  output$table_cat1a <- render_approach_table("cat_1a")
  output$table_cat1b <- render_approach_table("cat_1b")
  output$table_cat2 <- render_approach_table("cat_2")
  output$table_cat3 <- render_approach_table("cat_3")

  # Render accession counts per approach
  output$count_cat1a <- renderText({
    res <- selection_results()
    if (!is.null(res)) paste0("Selected accessions: ", nrow(res$approach_results$cat_1a)) else ""
  })
  output$count_cat1b <- renderText({
    res <- selection_results()
    if (!is.null(res)) paste0("Selected accessions: ", nrow(res$approach_results$cat_1b)) else ""
  })
  output$count_cat2 <- renderText({
    res <- selection_results()
    if (!is.null(res)) paste0("Selected accessions: ", nrow(res$approach_results$cat_2)) else ""
  })
  output$count_cat3 <- renderText({
    res <- selection_results()
    if (!is.null(res)) paste0("Selected accessions: ", nrow(res$approach_results$cat_3)) else ""
  })

  # --- Scatter Plot (Task 12.6) -----------------------------------------------

  # Update colour/size variable selectors after analysis
  observe({
    mapping <- column_mapping()
    estimates <- accession_estimates()
    req(mapping, estimates)
    
    # Available variables for colour: all categorical + all numeric traits + aux
    all_vars <- c(mapping$group_cols, mapping$trait_cols, mapping$aux_cols)
    
    # Default colour to none (so user sees the option clearly)
    updateSelectInput(session, "colour_var", 
                      choices = c("None" = "", all_vars),
                      selected = "")
    updateSelectInput(session, "size_var",
                      choices = c("None" = "", mapping$trait_cols, mapping$aux_cols),
                      selected = "")
  })

  # Render scatter plot with marginal density curves
  output$scatter_static <- renderPlot({
    res <- selection_results()
    estimates <- accession_estimates()
    req(res, estimates)
    
    config <- build_config()
    config$reference_baseline_avg <- res$reference_baseline_avg
    
    colour <- if (!is.null(input$colour_var) && input$colour_var != "") input$colour_var else NULL
    size <- if (!is.null(input$size_var) && input$size_var != "") input$size_var else NULL
    
    plot_result <- plot_scatter_with_thresholds(estimates, config, colour, size)
    plot_result$static_plot
  })

  # Download scatter plot as PNG
  output$download_plot <- downloadHandler(
    filename = function() {
      paste0("LMF_Scatter_Plot_", format(Sys.Date(), "%Y%m%d"), ".png")
    },
    content = function(file) {
      res <- selection_results()
      estimates <- accession_estimates()
      req(res, estimates)
      
      config <- build_config()
      config$reference_baseline_avg <- res$reference_baseline_avg
      
      colour <- if (!is.null(input$colour_var) && input$colour_var != "") input$colour_var else NULL
      size <- if (!is.null(input$size_var) && input$size_var != "") input$size_var else NULL
      
      plot_result <- plot_scatter_with_thresholds(estimates, config, colour, size)
      
      png(file, width = 10, height = 8, units = "in", res = 300)
      print(plot_result$static_plot)
      dev.off()
    }
  )

  # --- Diagnostics Panel (Task 12.7) ------------------------------------------

  output$diagnostics_panel <- renderUI({
    diag <- diagnostics_data()
    req(diag)
    
    cv_data <- cv_diagnostics()
    blues_data <- blues_diagnostics()
    
    # Build diagnostics display
    tagList(
      tags$h6("Dataset Summary"),
      tags$ul(
        tags$li(paste0("Total records: ", diag$n_records)),
        tags$li(paste0("Unique accessions: ", diag$n_accessions)),
        tags$li(paste0("Unique trials/runs: ", diag$n_trials)),
        tags$li(paste0("Observations per accession \u2014 Min: ", diag$obs_per_accession$min,
                       ", Median: ", diag$obs_per_accession$median,
                       ", Max: ", diag$obs_per_accession$max))
      ),
      
      tags$h6("Missing Values per Trait"),
      tags$ul(
        lapply(names(diag$missing_per_trait), function(trait) {
          tags$li(paste0(trait, ": ", diag$missing_per_trait[[trait]]))
        })
      ),
      
      tags$h6("Standard Deviation per Trait (across estimates)"),
      tags$ul(
        lapply(names(diag$sd_per_trait), function(trait) {
          tags$li(paste0(trait, ": ", round(diag$sd_per_trait[[trait]], 4)))
        })
      ),
      
      # CV% section
      if (!is.null(cv_data) && nrow(cv_data) > 0) {
        tagList(
          tags$hr(),
          tags$h6("Technical-Replicate CV% by Trait"),
          DT::renderDataTable(
            DT::datatable(cv_data, options = list(dom = "t", pageLength = 20, scrollX = TRUE), rownames = FALSE)
          )
        )
      },
      
      # Variance components and LSD (from BLUEs)
      if (!is.null(blues_data)) {
        tagList(
          tags$hr(),
          tags$h6("Variance Components (from mixed model)"),
          if (nrow(blues_data$variance_components) > 0) {
            # Simple text summary: Trait (Accession: X%, Trial: Y%, Residual: Z%)
            vc_df <- blues_data$variance_components
            traits_in_vc <- unique(vc_df$Trait)
            tagList(
              tags$ul(
                lapply(traits_in_vc, function(tr) {
                  tr_data <- vc_df[vc_df$Trait == tr, ]
                  parts <- paste0(tr_data$Source, ": ", tr_data$Pct_of_Total, "%")
                  tags$li(paste0(tr, " (", paste(parts, collapse = ", "), ")"))
                })
              ),
              # 100% stacked column chart (vertical)
              renderPlot({
                vc_plot_df <- vc_df
                ggplot(vc_plot_df, aes(x = Trait, y = Pct_of_Total, fill = Source)) +
                  geom_col(position = "stack", width = 0.6) +
                  scale_fill_brewer(palette = "Set2") +
                  labs(x = NULL, y = "% of Total Variance", fill = "Source") +
                  theme_minimal() +
                  theme(legend.position = "bottom",
                        axis.text.x = element_text(angle = 45, hjust = 1))
              }, height = 350)
            )
          },
          tags$h6("Average LSD per Trait"),
          tags$ul(
            lapply(names(blues_data$avg_lsd), function(trait) {
              tags$li(paste0(trait, ": ", blues_data$avg_lsd[[trait]]))
            })
          ),
          if (length(blues_data$errors) > 0) {
            tags$div(class = "alert alert-warning",
              tags$strong("Model fitting warnings:"),
              tags$ul(lapply(blues_data$errors, tags$li))
            )
          }
        )
      },
      
      tags$hr(),
      tags$h6("Reference Control Verification"),
      if (!is.null(diag$control_error)) {
        tags$div(class = "alert alert-danger", diag$control_error)
      } else if (!is.null(diag$control_missing_trials) && length(diag$control_missing_trials) > 0) {
        tags$div(class = "alert alert-warning",
          paste0("Reference control missing from trials: ",
                 paste(diag$control_missing_trials, collapse = ", "))
        )
      } else {
        tags$div(class = "alert alert-success", "Reference control present in all trials.")
      },
      
      tags$h6("Replication"),
      tags$p(paste0("Accessions with < 2 observations: ", diag$low_replication_count))
    )
  })

  # --- Excel Download Handler (Task 12.8) -------------------------------------

  output$download_results <- downloadHandler(
    filename = function() {
      paste0("LMF_Selection_Results_", format(Sys.Date(), "%Y%m%d"), ".xlsx")
    },
    content = function(file) {
      res <- selection_results()
      diag <- diagnostics_data()
      
      if (is.null(res) || is.null(diag)) {
        showNotification("Please run the analysis before exporting results.", type = "error")
        return(NULL)
      }
      
      config <- build_config()
      config$reference_baseline_avg <- res$reference_baseline_avg
      
      # Generate workbook
      wb_path <- export_workbook(
        approach_results = res$approach_results,
        all_estimates = res$all_estimates,
        config = config,
        diagnostics = diag,
        mapping = column_mapping()
      )
      
      # Copy temp file to the download location
      file.copy(wb_path, file)
      unlink(wb_path)
    }
  )

}

# --- Launch Application -------------------------------------------------------

shinyApp(ui, server, options = list(launch.browser = TRUE))
