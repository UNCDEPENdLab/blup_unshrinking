#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(dplyr))

source(file.path("R", "mplus_helpers.R"), local = TRUE)
source(file.path("vig_hallquist_2026", "vh_study2.R"), local = TRUE)

critical_warning <- mplus_warning_diagnostics(
  list("THE MODEL ESTIMATION DID NOT TERMINATE NORMALLY DUE TO A SINGULAR MATRIX.")
)
benign_warning <- mplus_warning_diagnostics(list("A non-critical reporting note."))

stopifnot(
  critical_warning$mplus_warning_count == 1L,
  critical_warning$mplus_critical_warning,
  !benign_warning$mplus_critical_warning
)

rows <- tibble::tibble(
  method = c(
    "oracle_dual", "naive_dual_blup", "fuller_alpha_stepdown_closed_form",
    "lai_2spa", "msem", "msem"
  ),
  estimate = rep(0.4, 6L),
  se = rep(0.1, 6L),
  status_code = rep(0L, 6L),
  analysis_eligible = c(TRUE, FALSE, NA, NA, NA, NA),
  analysis_exclusion_reason = c(NA, "stage2_near_collinear", NA, NA, NA, NA),
  fuller_auto_guard_pass = c(NA, NA, FALSE, NA, NA, NA),
  fuller_auto_guard_reason = c(NA, NA, "condition_cap", NA, NA, NA),
  mx_issue_class = c(NA, NA, NA, "ok", NA, NA),
  mx_info_definite = c(NA, NA, NA, FALSE, NA, NA),
  mx_condition_number = rep(NA_real_, 6L),
  mplus_critical_warning = c(NA, NA, NA, NA, TRUE, FALSE),
  mplus_target_parameter_count = c(NA_integer_, NA_integer_, NA_integer_, NA_integer_, 1L, 1L)
)

classified <- add_study2_analysis_eligibility(rows)

stopifnot(
  identical(classified$analysis_eligible, c(TRUE, FALSE, FALSE, FALSE, FALSE, TRUE)),
  identical(
    classified$analysis_exclusion_reason,
    c(
      NA_character_, "stage2_near_collinear", "fuller_guard_condition_cap",
      "openmx_information_not_definite", "mplus_critical_warning", NA_character_
    )
  )
)

cat("Study 2 analysis-eligibility tests ok\n")
