#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(dplyr))

source(file.path("R", "core_utils.R"), local = TRUE)
source(file.path("R", "sim_diagnostics.R"), local = TRUE)
source(file.path("vig_hallquist_2026", "vh_study3.R"), local = TRUE)

rows <- tibble::tibble(
  method = c(
    "oracle_dual", "naive_blup_on_blup", "point_only_diagnostic",
    "fuller_alpha_stepdown_closed_form", "lai_2spa", "sem"
  ),
  estimate = rep(0.4, 6L),
  se = c(0.1, 0.1, NA_real_, 0.1, 0.1, 0.1),
  ci_low = c(0.2, 0.2, NA_real_, 0.2, 0.2, 0.2),
  ci_high = c(0.6, 0.6, NA_real_, 0.6, 0.6, 0.6),
  status_code = rep(0L, 6L),
  analysis_eligible = c(TRUE, FALSE, NA, NA, NA, NA),
  analysis_exclusion_reason = c(NA, "stage2_near_collinear", NA, NA, NA, NA),
  fuller_auto_guard_pass = c(NA, NA, NA, FALSE, NA, NA),
  fuller_auto_guard_reason = c(NA, NA, NA, "condition_cap", NA, NA),
  mx_issue_class = c(NA, NA, NA, NA, "ok", NA),
  mx_info_definite = c(NA, NA, NA, NA, FALSE, NA),
  mx_condition_number = rep(NA_real_, 6L),
  mplus_critical_warning = c(NA, NA, NA, NA, NA, TRUE),
  mplus_target_parameter_count = c(NA_integer_, NA_integer_, NA_integer_, NA_integer_, NA_integer_, 1L)
)

classified <- add_study3_analysis_eligibility(rows)

stopifnot(
  identical(classified$point_eligible, c(TRUE, FALSE, TRUE, TRUE, FALSE, FALSE)),
  identical(classified$analysis_eligible, classified$point_eligible),
  identical(classified$interval_eligible, c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE)),
  identical(
    classified$point_exclusion_reason,
    c(
      NA_character_, "stage2_near_collinear", NA_character_,
      NA_character_, "openmx_information_not_definite",
      "mplus_critical_warning"
    )
  ),
  identical(
    classified$interval_exclusion_reason,
    c(
      NA_character_, "stage2_near_collinear", "invalid_standard_error",
      NA_character_, "openmx_information_not_definite",
      "mplus_critical_warning"
    )
  ),
  identical(
    add_study3_analysis_eligibility(classified)$interval_eligible,
    classified$interval_eligible
  )
)

cat("Study 3 analysis-eligibility tests ok\n")
