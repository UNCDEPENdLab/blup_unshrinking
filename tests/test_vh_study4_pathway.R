#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
  library(MASS)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(foreach)
  library(doParallel)
  library(OpenMx)
  library(glmnet)
  library(sandwich)
  library(geigen)
  library(MplusAutomation)
  library(glue)
})

source(file.path(
  "vig_hallquist_2026",
  "random_effects_structural_simulation.R"
))

design <- select_design("4")
stopifnot(
  nrow(design) == 111L,
  sum(design$is_falsification_control) == 3L,
  setequal(
    unique(design$information_profile),
    c(
      "homogeneous",
      "moderate",
      "severe",
      "severe_information_matched"
    )
  ),
  all(abs(design$achieved_reliability - design$target_reliability) < 1e-8),
  all(design$outcome_residual_variance > 0),
  all(design$mean_clus_size == 10),
  all(design$reliability_sd[design$information_profile == "homogeneous"] < 1e-12),
  all(design$reliability_sd[design$information_profile == "severe"] > 0),
  all(c(
    "reliability_iqr",
    "population_lambda22_mean",
    "population_lambda22_min",
    "population_lambda22_max",
    "population_lambda_matrix_frobenius_rms_dispersion",
    "population_theta22_mean",
    "population_theta_matrix_frobenius_rms_dispersion",
    "population_ols_var22_mean",
    "population_ols_var22_min",
    "population_ols_var22_max",
    "population_ols_cov_matrix_frobenius_rms_dispersion"
  ) %in% names(design)),
  all(design$population_lambda_matrix_frobenius_rms_dispersion[
    design$information_profile == "homogeneous"
  ] < 1e-12),
  all(design$population_lambda_matrix_frobenius_rms_dispersion[
    design$information_profile == "severe"
  ] > 0),
  all(design$population_ols_cov_matrix_frobenius_rms_dispersion[
    design$information_profile == "severe"
  ] > 0)
)

severe_condition <- design %>%
  dplyr::filter(
    num_clus == 50L,
    information_profile == "severe",
    target_reliability == 0.5,
    marginal_rho == 0,
    standardized_beta_target == 0.4
  ) %>%
  dplyr::slice(1L)

set.seed(17)
severe_sim <- simulate_study4(severe_condition)
size_table <- table(severe_sim$lv2_true$cluster_size)
stopifnot(
  nrow(severe_sim$lv1) == 500L,
  nrow(severe_sim$lv2_true) == 50L,
  identical(as.integer(names(size_table)), c(3L, 17L)),
  all(as.integer(size_table) == c(25L, 25L)),
  abs(mean(severe_sim$cluster_reliability) - 0.5) < 1e-10,
  abs(min(severe_sim$cluster_reliability) - severe_condition$reliability_min) < 1e-10,
  abs(max(severe_sim$cluster_reliability) - severe_condition$reliability_max) < 1e-10,
  identical(severe_sim$lv1$cid_chr, as.character(severe_sim$lv1$cid))
)

control_condition <- design %>%
  dplyr::filter(
    num_clus == 50L,
    information_profile == "severe_information_matched"
  ) %>%
  dplyr::slice(1L)
set.seed(18)
control_sim <- simulate_study4(control_condition)
cluster_slope_information <- control_sim$lv1 %>%
  dplyr::group_by(cid_chr) %>%
  dplyr::summarise(slope_information = sum(x^2), .groups = "drop")
stopifnot(
  all(abs(cluster_slope_information$slope_information - 9) < 1e-10),
  stats::sd(control_sim$cluster_reliability) < 1e-10,
  abs(
    control_condition$population_lambda22_max -
      control_condition$population_lambda22_min
  ) < 1e-12,
  abs(
    control_condition$population_ols_var22_max -
      control_condition$population_ols_var22_min
  ) < 1e-12
)

average_measurement_fixture <- tibble::tibble(
  u0_eb = c(0.4, -0.2, 0.7),
  u1_eb = c(-0.3, 0.8, 0.1),
  lambda11 = c(0.80, 0.90, 0.85),
  lambda12 = c(0.04, 0.02, 0.03),
  lambda21 = c(-0.01, 0.01, 0.00),
  lambda22 = c(0.55, 0.75, 0.65),
  theta11 = c(0.08, 0.04, 0.06),
  theta12 = c(0.01, 0.03, 0.02),
  theta22 = c(0.18, 0.10, 0.14)
)
average_measurement_prepared <- prepare_fuller_average_measurement(
  average_measurement_fixture
)
expected_lambda_bar <- matrix(
  c(0.85, 0.03, 0.00, 0.65),
  nrow = 2L,
  byrow = TRUE
)
expected_theta_bar <- matrix(
  c(0.06, 0.02, 0.02, 0.14),
  nrow = 2L,
  byrow = TRUE
)
expected_lambda_inv <- solve(expected_lambda_bar)
expected_average_scores <-
  as.matrix(average_measurement_fixture[, c("u0_eb", "u1_eb")]) %*%
  t(expected_lambda_inv)
expected_average_cov <-
  expected_lambda_inv %*% expected_theta_bar %*% t(expected_lambda_inv)
stopifnot(
  isTRUE(all.equal(
    unname(average_measurement_prepared$lambda_bar),
    expected_lambda_bar,
    tolerance = 1e-12
  )),
  isTRUE(all.equal(
    unname(average_measurement_prepared$theta_bar),
    expected_theta_bar,
    tolerance = 1e-12
  )),
  isTRUE(all.equal(
    unname(as.matrix(average_measurement_prepared$data[, c(
      "fuller_average_u0", "fuller_average_u1"
    )])),
    expected_average_scores,
    tolerance = 1e-12
  )),
  isTRUE(all.equal(
    unname(average_measurement_prepared$measurement_covariance_bar),
    expected_average_cov,
    tolerance = 1e-12
  )),
  length(unique(average_measurement_prepared$data$fuller_average_meas22)) == 1L
)
singular_average_fixture <- average_measurement_fixture[
  rep(seq_len(nrow(average_measurement_fixture)), length.out = 10L),
  ,
  drop = FALSE
]
singular_average_fixture$z <- seq_len(nrow(singular_average_fixture)) / 10
singular_average_fixture$lambda12 <- 0
singular_average_fixture$lambda21 <- 0
singular_average_fixture$lambda22 <- 0
singular_average_fit <- fit_fuller_average_measurement(
  singular_average_fixture,
  outcome = "z"
)
stopifnot(
  singular_average_fit$status_code == 3L,
  all(fuller_dual_result_columns() %in% names(singular_average_fit)),
  singular_average_fit$mx_issue_class ==
    "fuller_average_measurement_transform_failed",
  grepl("not invertible", singular_average_fit$mx_issue_detail, fixed = TRUE)
)

estimation_condition <- design %>%
  dplyr::filter(
    num_clus == 50L,
    information_profile == "moderate",
    target_reliability == 0.8,
    marginal_rho == 0,
    standardized_beta_target == 0.4
  ) %>%
  dplyr::slice(1L)
set.seed(1707)
results <- run_study_rep(estimation_condition)

locally_required_finite <- c(
  "oracle_dual",
  "naive_dual_blup",
  "closed_form_dual",
  "fuller_closed_form",
  "fuller_average_measurement",
  "fuller_alpha_stepdown_closed_form",
  "lai_2spa",
  "lai_2spaa"
)
stopifnot(
  nrow(results) == length(study4_methods()),
  setequal(results$method, study4_methods()),
  all(results$truth == estimation_condition$standardized_beta_target),
  all(!is.na(results$method_role)),
  all(is.finite(results$estimate[results$method %in% locally_required_finite])),
  results$realized_reliability_sd[[1]] > 0,
  results$realized_reliability_iqr[[1]] > 0,
  is.finite(results$lambda22_mean[[1]]),
  results$lambda22_sd[[1]] > 0,
  results$lambda22_min[[1]] < results$lambda22_max[[1]],
  results$lambda_matrix_frobenius_rms_dispersion[[1]] > 0,
  results$theta_matrix_frobenius_rms_dispersion[[1]] > 0,
  results$ols_var22_sd[[1]] > 0,
  results$ols_cov_matrix_frobenius_rms_dispersion[[1]] > 0,
  all(is.finite(results$blup_slope_rmse_small)),
  all(is.finite(results$blup_slope_rmse_large)),
  all(is.finite(results$corrected_slope_rmse_small)),
  all(is.finite(results$corrected_slope_rmse_large)),
  all(is.finite(results$average_measurement_slope_rmse_small)),
  all(is.finite(results$average_measurement_slope_rmse_large))
)

fuller_results <- results %>%
  dplyr::filter(method %in% c(
    "fuller_closed_form",
    "fuller_average_measurement",
    "fuller_alpha_stepdown_closed_form"
  ))
stopifnot(
  all(fuller_dual_result_columns() %in% names(results)),
  all(fuller_results$fuller_measurement_weight_used == 1),
  all(is.finite(fuller_results$fuller_alpha_step1_used)),
  all(is.finite(fuller_results$fuller_alpha_step3_used)),
  all(is.finite(fuller_results$fuller_correction1)),
  all(is.finite(fuller_results$fuller_sx1_star_min_eigen)),
  all(is.finite(fuller_results$fuller_sx_star_min_eigen))
)

summary <- summarize_results_df(
  dplyr::bind_cols(
    results,
    estimation_condition[rep(1L, nrow(results)), , drop = FALSE] %>%
      dplyr::select(-study),
    tibble::tibble(rep = 1L)
  )
)
stopifnot(
  nrow(summary) == length(study4_methods()),
  all(c(
    "information_profile",
    "reliability_sd",
    "reliability_min",
    "reliability_max",
    "median_absolute_error",
    "p95_absolute_error",
    "mean_realized_reliability_iqr",
    "mean_lambda_matrix_frobenius_rms_dispersion",
    "mean_theta_matrix_frobenius_rms_dispersion",
    "mean_ols_cov_matrix_frobenius_rms_dispersion",
    "mean_blup_slope_rmse_small",
    "mean_corrected_slope_rmse_large",
    "mean_average_measurement_slope_rmse_small",
    "mean_fuller_alpha_step1_used",
    "mean_fuller_sx_star_min_eigen"
  ) %in% names(summary)),
  all(c(
    "study4_profile_spec",
    "study4_time_design",
    "make_study4_cluster_sizes",
    "study4_weighted_quantile",
    "study4_matrix_rms_dispersion",
    "study4_measurement_matrix_summary",
    "simulate_data_study4",
    "study4_measurement_diagnostics",
    "prepare_fuller_average_measurement",
    "fit_fuller_average_measurement"
  ) %in% vig_hallquist_parallel_exports())
)

absolute_error_fixture <- tibble::tibble(
  condition_id = 9999L,
  study = "study4",
  method = "fixture",
  estimate = c(0, 1, 2, NA_real_),
  truth = 0,
  ci_low = c(-0.1, 0.9, 1.9, NA_real_),
  ci_high = c(0.1, 1.1, 2.1, NA_real_),
  status_code = c(0L, 0L, 0L, 3L)
)
absolute_error_summary <- summarize_results_df(absolute_error_fixture)
stopifnot(
  absolute_error_summary$median_absolute_error == 1,
  isTRUE(all.equal(
    absolute_error_summary$p95_absolute_error,
    unname(stats::quantile(c(0, 1, 2), 0.95)),
    tolerance = 1e-12
  ))
)

cat("VH Study 4 pathway test ok\n")
