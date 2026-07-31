#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

find_repo_root <- function() {
  candidates <- unique(normalizePath(c(getwd(), file.path(getwd(), "..")), mustWork = FALSE))
  for (candidate in candidates) {
    if (file.exists(file.path(candidate, "lai_study1_replication_vh", "estimators.R"))) {
      return(candidate)
    }
  }
  stop("Could not locate the repository root.")
}

repo_root <- find_repo_root()
source(file.path(repo_root, "lai_study1_replication_vh", "estimators.R"), local = TRUE)

rows <- tibble::tibble(
  method = c(
    "oracle_dual", "naive_dual_blup", "closed_form_dual", "lai_2spa",
    "lai_2spa", "msem", "msem", "fuller_closed_form", "oracle_dual"
  ),
  estimate = c(0.4, 0.4, NA_real_, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4),
  se = c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, NA_real_, 0.1),
  ci_low = rep(0.2, 9L),
  ci_high = rep(0.6, 9L),
  status_code = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 10L),
  analysis_eligible = c(TRUE, FALSE, NA, NA, NA, NA, NA, NA, TRUE),
  analysis_exclusion_reason = c(NA, "stage2_near_collinear", NA, NA, NA, NA, NA, NA, NA),
  mx_issue_class = c(NA, NA, NA, "not_converged", "ok", NA, NA, NA, NA),
  mx_info_definite = c(NA, NA, NA, TRUE, FALSE, NA, NA, NA, NA),
  mx_condition_number = c(NA, NA, NA, NA, NA, NA, NA, NA, NA),
  mplus_critical_warning = c(NA, NA, NA, NA, NA, TRUE, FALSE, NA, NA),
  mplus_target_parameter_count = c(NA, NA, NA, NA, NA, 1L, 1L, NA, NA)
)

classified <- add_lai_study1_vh_analysis_eligibility(rows)

stopifnot(
  identical(
    classified$vh_analysis_eligible,
    c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE)
  ),
  # Lai's historical rule only uses status code, so rows rejected by VH-only
  # diagnostics remain historically retained whenever code is zero.
  identical(
    classified$lai_original_eligible,
    c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)
  ),
  identical(classified$analysis_eligible, classified$vh_analysis_eligible),
  identical(
    classified$eligibility_comparison,
    c("both", "lai_original_only", "lai_original_only", "lai_original_only",
      "lai_original_only", "lai_original_only", "both", "lai_original_only",
      "neither")
  ),
  identical(
    classified$vh_analysis_exclusion_reason,
    c(
      NA_character_, "stage2_near_collinear", "estimation_unavailable",
      "openmx_issue", "openmx_information_not_definite",
      "mplus_critical_warning", NA_character_, "invalid_standard_error",
      "estimation_status_nonzero"
    )
  ),
  identical(
    classified$lai_original_exclusion_reason,
    c(rep(NA_character_, 8L), "status_code_nonzero")
  )
)

cat("Lai Study 1 VH eligibility tests ok\n")
