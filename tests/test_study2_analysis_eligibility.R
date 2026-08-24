#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(dplyr))

source(file.path("R", "core_utils.R"), local = TRUE)
source(file.path("R", "sim_diagnostics.R"), local = TRUE)
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

# Endogenous Mplus latent variables are reported under Residual.Variances,
# whereas exogenous latent variables use Variances. The audit must recognize
# both headers when reconstructing intercept/slope covariance blocks.
synthetic_parameters <- tibble::tibble(
  paramHeader = c(
    "Variances", "Variances", "I1.WITH",
    "Variances", "Residual.Variances", "I2.WITH"
  ),
  param = c("I1", "S1", "S1", "I2", "S2", "S2"),
  est = c(0.8, 0.6, 0.1, 0.7, 0.5, -0.08)
)
latent_diagnostics <- populate_mplus_latent_diagnostics(
  mplus_diagnostics_template(),
  pars = synthetic_parameters,
  xvar = "S1",
  latent_covariance_blocks = list(
    predictor = c("I1", "S1"),
    outcome = c("I2", "S2")
  )
)
stopifnot(
  latent_diagnostics$mplus_predictor_latent_slope_variance == 0.6,
  latent_diagnostics$mplus_outcome_latent_slope_variance == 0.5,
  latent_diagnostics$mplus_outcome_latent_intercept_slope_covariance == -0.08,
  is.finite(latent_diagnostics$mplus_outcome_latent_covariance_min_eigenvalue),
  !latent_diagnostics$mplus_outcome_latent_covariance_boundary
)

rows <- tibble::tibble(
  method = c(
    "oracle_dual", "naive_dual_blup", "fuller_alpha_stepdown_closed_form",
    "lai_2spa", "msem", "msem"
  ),
  estimate = rep(0.4, 6L),
  se = rep(0.1, 6L),
  ci_low = rep(0.2, 6L),
  ci_high = rep(0.6, 6L),
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
  identical(classified$point_eligible, c(TRUE, FALSE, TRUE, FALSE, FALSE, TRUE)),
  identical(classified$analysis_eligible, classified$point_eligible),
  identical(classified$interval_eligible, classified$point_eligible),
  identical(
    classified$point_exclusion_reason,
    c(
      NA_character_, "stage2_near_collinear", NA_character_,
      "openmx_information_not_definite", "mplus_critical_warning", NA_character_
    )
  ),
  identical(classified$analysis_exclusion_reason, classified$point_exclusion_reason)
)

cat("Study 2 analysis-eligibility tests ok\n")
