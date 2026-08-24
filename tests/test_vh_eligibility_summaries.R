#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(dplyr))

source(file.path("R", "core_utils.R"), local = TRUE)
source(file.path("R", "sim_diagnostics.R"), local = TRUE)
source(file.path("vig_hallquist_2026", "vh_runner.R"), local = TRUE)

null_rows <- tibble::tibble(
  condition_id = 1L,
  study = "study2v2",
  method = "fixture",
  estimate = c(0, 0.2, -0.3, 99, 0.1),
  truth = 0,
  se = c(0.1, 0.1, 0.1, 0.1, NA_real_),
  ci_low = c(-0.2, 0.01, -0.5, 98, NA_real_),
  ci_high = c(0.2, 0.4, -0.01, 100, NA_real_),
  status_code = 0L,
  analysis_eligible = c(TRUE, TRUE, TRUE, FALSE, TRUE),
  analysis_exclusion_reason = c(
    NA, NA, NA, "stage2_rank_deficient", NA
  ),
  interval_eligible = c(TRUE, TRUE, TRUE, FALSE, FALSE),
  interval_exclusion_reason = c(
    NA, NA, NA, "stage2_rank_deficient", "invalid_standard_error"
  )
)

nonnull_rows <- tibble::tibble(
  condition_id = 2L,
  study = "study2v2",
  method = "fixture",
  estimate = c(0.4, 0.6, NA_real_),
  truth = 0.4,
  se = c(0.1, 0.1, NA_real_),
  ci_low = c(-0.1, 0.41, NA_real_),
  ci_high = c(0.6, 0.79, NA_real_),
  status_code = c(0L, 0L, NA_integer_),
  analysis_eligible = NA,
  analysis_exclusion_reason = NA_character_,
  interval_eligible = NA,
  interval_exclusion_reason = NA_character_
)

summary_rows <- summarize_results_df(dplyr::bind_rows(null_rows, nonnull_rows))
null_summary <- dplyr::filter(summary_rows, condition_id == 1L)
nonnull_summary <- dplyr::filter(summary_rows, condition_id == 2L)

expected_point_sd <- stats::sd(c(0, 0.2, -0.3, 0.1))
expected_interval_sd <- stats::sd(c(0, 0.2, -0.3))
stopifnot(
  null_summary$n_rep == 5L,
  null_summary$n_point_eligible == 4L,
  null_summary$n_interval_eligible == 3L,
  null_summary$n_success == null_summary$n_point_eligible,
  null_summary$n_analysis_eligible == null_summary$n_point_eligible,
  isTRUE(all.equal(null_summary$point_eligibility, 4 / 5)),
  isTRUE(all.equal(null_summary$interval_eligibility, 3 / 5)),
  isTRUE(all.equal(null_summary$mean_estimate, 0)),
  isTRUE(all.equal(null_summary$bias, 0)),
  isTRUE(all.equal(null_summary$empirical_sd, expected_point_sd)),
  isTRUE(all.equal(
    null_summary$empirical_sd_interval_subset,
    expected_interval_sd
  )),
  isTRUE(all.equal(
    null_summary$mean_se_to_empirical_sd,
    0.1 / expected_interval_sd
  )),
  isTRUE(all.equal(null_summary$coverage, 1 / 3)),
  isTRUE(all.equal(null_summary$success_and_cover, 1 / 5)),
  isTRUE(all.equal(null_summary$type1_error, 2 / 3)),
  isTRUE(all.equal(null_summary$type1_error_positive, 1 / 3)),
  isTRUE(all.equal(null_summary$type1_error_negative, 1 / 3)),
  isTRUE(all.equal(null_summary$operational_type1_error, 2 / 5)),
  is.na(null_summary$power),
  grepl("stage2_rank_deficient=1", null_summary$point_exclusion_reasons),
  grepl("invalid_standard_error=1", null_summary$interval_exclusion_reasons),
  is.finite(null_summary$coverage_mc_low),
  is.finite(null_summary$coverage_mc_high),
  null_summary$bias_mc_low <= null_summary$bias,
  null_summary$bias_mc_high >= null_summary$bias
)

stopifnot(
  nonnull_summary$n_point_eligible == 2L,
  nonnull_summary$n_interval_eligible == 2L,
  isTRUE(all.equal(nonnull_summary$coverage, 1 / 2)),
  isTRUE(all.equal(nonnull_summary$power, 1 / 2)),
  isTRUE(all.equal(nonnull_summary$operational_power, 1 / 3)),
  is.na(nonnull_summary$type1_error)
)

stage1_summary <- null_rows %>%
  dplyr::mutate(
    stage1_singular_problem = FALSE,
    stage1_lmer_singular = FALSE,
    stage1_re_corr = 0.2,
    stage1_eb_corr = 0.3,
    stage1_design_kappa = 2,
    stage1_problem_detail = NA_character_
  ) %>%
  summarize_stage1_problem_df()
stopifnot(
  stage1_summary$n_point_eligible == 4L,
  stage1_summary$n_interval_eligible == 3L,
  isTRUE(all.equal(stage1_summary$coverage, 1 / 3)),
  isTRUE(all.equal(stage1_summary$success_and_cover, 1 / 5))
)

# Explicit estimator interval exclusions are monotone under reclassification.
reclassified <- add_vh_analysis_eligibility(null_rows)
reclassified_twice <- add_vh_analysis_eligibility(reclassified)
stopifnot(
  identical(reclassified$point_eligible, reclassified_twice$point_eligible),
  identical(reclassified$interval_eligible, reclassified_twice$interval_eligible),
  identical(
    reclassified$interval_exclusion_reason,
    reclassified_twice$interval_exclusion_reason
  )
)

cat("VH eligibility and inferential-summary tests ok\n")
