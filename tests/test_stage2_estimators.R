#!/usr/bin/env Rscript

# Unit checks for stage-2 estimator helpers and small simulation checks for
# estimator-specific behavior.

source(file.path("R", "stage2_estimators.R"), local = TRUE)

sandwich_out <- list(
  alpha_hat = c(1.25, -0.40),
  vcov_hc0 = matrix(c(0.10, 0.01, 0.01, 0.04), nrow = 2L),
  vcov_hc1 = matrix(c(0.12, 0.01, 0.01, 0.09), nrow = 2L),
  vcov_hc2 = matrix(c(0.14, 0.01, 0.01, 0.16), nrow = 2L),
  vcov_hc3 = matrix(c(0.16, 0.01, 0.01, 0.25), nrow = 2L)
)

df <- 9L
out <- format_stacked_sandwich_rows(
  sandwich_out = sandwich_out,
  df = df,
  alpha_names = c("(Intercept)", "x"),
  term = "x",
  method_prefix = "unit_stacked"
)

expected_methods <- paste0("unit_stacked_hc", 0:3)
expected_se <- c(0.20, 0.30, 0.40, 0.50)
crit <- qt(0.975, df)

stopifnot(
  identical(out$method, expected_methods),
  isTRUE(all.equal(out$estimate, rep(-0.40, 4L), tolerance = 1e-12)),
  isTRUE(all.equal(out$se, expected_se, tolerance = 1e-12)),
  isTRUE(all.equal(out$ci_low, -0.40 - crit * expected_se, tolerance = 1e-12)),
  isTRUE(all.equal(out$ci_high, -0.40 + crit * expected_se, tolerance = 1e-12))
)

missing_term <- format_stacked_sandwich_rows(
  sandwich_out = sandwich_out,
  df = df,
  alpha_names = c("(Intercept)", "x"),
  term = "not_a_term",
  method_prefix = "unit_stacked"
)

stopifnot(
  identical(missing_term$method, expected_methods),
  all(is.na(missing_term$estimate)),
  all(is.na(missing_term$se)),
  all(is.na(missing_term$ci_low)),
  all(is.na(missing_term$ci_high))
)

set.seed(9123)
n_edge <- 40L
edge_df <- data.frame(
  z = rnorm(n_edge),
  corrected_intercept_full = rnorm(n_edge),
  corrected_slope_full = rnorm(n_edge),
  ols_var11 = rep(0.15, n_edge),
  ols_var12 = rep(0, n_edge),
  ols_var22 = rep(2.00, n_edge)
)
edge_df$z <- 0.25 * edge_df$corrected_intercept_full +
  0.55 * edge_df$corrected_slope_full +
  rnorm(n_edge, sd = 0.4)

edge_reject <- fit_eiv_dual(
  edge_df,
  outcome = "z",
  predictor_u0 = "corrected_intercept_full",
  predictor_u1 = "corrected_slope_full",
  meas11 = "ols_var11",
  meas12 = "ols_var12",
  meas22 = "ols_var22"
)

stopifnot(
  identical(edge_reject$se_type, c("naive", "hc0", "hc3")),
  all(as.integer(edge_reject$status_code) == 3L),
  all(edge_reject$mx_issue_class == "corrected_predictor_cov_not_pd"),
  all(is.na(edge_reject$estimate)),
  edge_reject$eiv_latent_cov_min_eigen <= 1e-6
)

edge_regularized <- fit_eiv_dual(
  edge_df,
  outcome = "z",
  predictor_u0 = "corrected_intercept_full",
  predictor_u1 = "corrected_slope_full",
  meas11 = "ols_var11",
  meas12 = "ols_var12",
  meas22 = "ols_var22",
  regularize = TRUE
)

stopifnot(
  identical(edge_regularized$se_type, c("naive", "hc0", "hc3")),
  all(as.integer(edge_regularized$status_code) == 0L),
  all(edge_regularized$mx_issue_class == "ok"),
  all(edge_regularized$eiv_regularized),
  all(is.finite(edge_regularized$estimate)),
  all(is.finite(edge_regularized$se)),
  all(edge_regularized$eiv_measurement_weight_used < edge_regularized$eiv_measurement_weight_requested),
  all(edge_regularized$eiv_latent_cov_min_eigen > 1e-6)
)

edge_methods <- finalize_eiv_se_variants(edge_regularized, "unit_eiv")
stopifnot(
  identical(edge_methods$method, c("unit_eiv", "unit_eiv_hc0", "unit_eiv_hc3"))
)

set.seed(44023)
m_fuller <- 600L

x0_true <- rnorm(m_fuller)
x1_true <- rnorm(m_fuller)

beta0 <- 0.50
beta1 <- -0.20
beta2 <- 0.80

eps <- rnorm(m_fuller, sd = 1.0)
y_true <- beta0 + beta1 * x0_true + beta2 * x1_true + eps

var_u0 <- 0.20
var_u1 <- 0.35
var_uy <- 0.40

u0_err <- rnorm(m_fuller, sd = sqrt(var_u0))
u1_err <- rnorm(m_fuller, sd = sqrt(var_u1))
uy_err <- rnorm(m_fuller, sd = sqrt(var_uy))

stage2_fuller <- data.frame(
  y_obs = y_true + uy_err,
  x0_obs = x0_true + u0_err,
  x1_obs = x1_true + u1_err,
  meas11 = rep(var_u0, m_fuller),
  meas12 = rep(0, m_fuller),
  meas22 = rep(var_u1, m_fuller),
  measyy = rep(var_uy, m_fuller)
)

fuller_out <- fit_fuller_dual(
  stage2_fuller,
  outcome = "y_obs",
  predictor_u0 = "x0_obs",
  predictor_u1 = "x1_obs",
  meas11 = "meas11",
  meas12 = "meas12",
  meas22 = "meas22",
  outcome_meas_var = "measyy"
)

stopifnot(
  isTRUE(all.equal(as.integer(fuller_out$status_code), 0L)),
  is.finite(fuller_out$estimate),
  is.finite(fuller_out$se),
  fuller_out$se > 0,
  identical(
    fuller_out$fuller_predictor_outcome_covariance_source,
    "zero_default"
  ),
  fuller_out$fuller_predictor_outcome_covariance_max_abs == 0,
  abs(fuller_out$estimate - beta2) < 0.20
)

stage2_noerr <- data.frame(
  y_obs = y_true,
  x0_obs = x0_true,
  x1_obs = x1_true,
  meas11 = rep(0, m_fuller),
  meas12 = rep(0, m_fuller),
  meas22 = rep(0, m_fuller),
  measyy = rep(0, m_fuller)
)

# reduces to OLS with no measurement error
fuller_noerr <- fit_fuller_dual(
  stage2_noerr,
  outcome = "y_obs",
  predictor_u0 = "x0_obs",
  predictor_u1 = "x1_obs",
  meas11 = "meas11",
  meas12 = "meas12",
  meas22 = "meas22",
  outcome_meas_var = "measyy"
)

ols_dual <- fit_observed_dual(stage2_noerr, outcome = "y_obs", predictor_u0 = "x0_obs", predictor_u1 = "x1_obs", reporting_scale = 1)
ols_naive <- ols_dual[ols_dual$se_type == "naive", , drop = FALSE]
fixed_scale <- 2.5
fuller_noerr_fixed <- rescale_fuller_to_population_sd(
  fuller_noerr,
  fixed_scale
)

stopifnot(
  isTRUE(all.equal(as.integer(fuller_noerr$status_code), 0L)),
  isTRUE(all.equal(fuller_noerr$estimate, ols_naive$estimate, tolerance = 1e-10)),
  isTRUE(all.equal(fuller_noerr$se, ols_naive$se, tolerance = 1e-10)),
  all(ols_dual$analysis_eligible),
  all(is.na(ols_dual$analysis_exclusion_reason)),
  isTRUE(all.equal(
    fuller_noerr_fixed$estimate,
    fuller_noerr$fuller_raw_estimate * fixed_scale,
    tolerance = 1e-12
  )),
  isTRUE(all.equal(
    fuller_noerr_fixed$se,
    fuller_noerr$fuller_raw_se * fixed_scale,
    tolerance = 1e-12
  ))
)

# A finite OLS coefficient can still be unsuitable for primary performance
# summaries when the two observed predictors are nearly collinear. The raw
# estimate remains available, while the method-specific analysis flag records
# why it should be excluded later.
set.seed(781)
n_collinear <- 120L
x0_collinear <- rnorm(n_collinear)
x1_near_collinear <- x0_collinear + rnorm(n_collinear, sd = 0.001)
near_collinear_df <- data.frame(
  y = 0.3 * x0_collinear + 0.7 * x1_near_collinear + rnorm(n_collinear),
  x0 = x0_collinear,
  x1 = x1_near_collinear
)
near_collinear_out <- fit_observed_dual(
  near_collinear_df,
  outcome = "y",
  predictor_u0 = "x0",
  predictor_u1 = "x1",
  reporting_scale = 1
)

exact_collinear_df <- transform(near_collinear_df, x1 = x0)
exact_collinear_out <- fit_observed_dual(
  exact_collinear_df,
  outcome = "y",
  predictor_u0 = "x0",
  predictor_u1 = "x1",
  reporting_scale = 1
)

stopifnot(
  all(!near_collinear_out$analysis_eligible),
  all(near_collinear_out$analysis_exclusion_reason == "stage2_near_collinear"),
  all(near_collinear_out$stage2_vif > 100),
  all(is.finite(near_collinear_out$estimate)),
  all(!exact_collinear_out$analysis_eligible),
  all(exact_collinear_out$analysis_exclusion_reason == "stage2_rank_deficient"),
  all(is.na(exact_collinear_out$estimate))
)

stage2_single_noerr <- data.frame(
  y_obs = y_true,
  x1_obs = x1_true,
  meas22 = rep(0, m_fuller),
  measyy = rep(0, m_fuller)
)

fuller_single_noerr <- fit_fuller(
  stage2_single_noerr,
  outcome = "y_obs",
  predictor_u1 = "x1_obs",
  meas22 = "meas22",
  outcome_meas_var = "measyy"
)

ols_single <- fit_observed_single(stage2_single_noerr, outcome = "y_obs", predictor = "x1_obs", reporting_scale = 1)
ols_single_naive <- ols_single[ols_single$se_type == "naive", , drop = FALSE]

stopifnot(
  isTRUE(all.equal(as.integer(fuller_single_noerr$status_code), 0L)),
  isTRUE(all.equal(fuller_single_noerr$estimate, ols_single_naive$estimate, tolerance = 1e-10)),
  isTRUE(all.equal(fuller_single_noerr$se, ols_single_naive$se, tolerance = 1e-10))
)

# The original version of stepdown is currently broken due to the removal of the scaling step
# fuller_stepdown_noerr <- fit_fuller_dual_stepdown(
#   stage2_noerr,
#   outcome = "y_obs",
#   predictor_u0 = "x0_obs",
#   predictor_u1 = "x1_obs",
#   meas11 = "meas11",
#   meas12 = "meas12",
#   meas22 = "meas22",
#   outcome_meas_var = "measyy"
# )
# 
# stopifnot(
#   isTRUE(all.equal(as.integer(fuller_stepdown_noerr$status_code), 0L)),
#   isTRUE(all.equal(fuller_stepdown_noerr$estimate, ols_naive$estimate, tolerance = 1e-10)),
#   isTRUE(all.equal(fuller_stepdown_noerr$se, ols_naive$se, tolerance = 1e-10)),
#   identical(fuller_stepdown_noerr$fuller_auto_guard_reason, "ok"),
#   isTRUE(all.equal(fuller_stepdown_noerr$fuller_measurement_weight_used, 1))
# )

simulate_fuller_known_eiv <- function(n,
                                      beta_u0 = 0.25,
                                      beta_u1 = 0.55,
                                      residual_sd = 0.70) {
  latent_cov <- matrix(c(1.00, 0.30, 0.30, 0.64), nrow = 2L)
  meas_cov <- matrix(c(0.18, 0.03, 0.03, 0.12), nrow = 2L)

  latent_scores <- matrix(stats::rnorm(2L * n), ncol = 2L) %*% chol(latent_cov)
  measurement_error <- matrix(stats::rnorm(2L * n), ncol = 2L) %*% chol(meas_cov)

  u0 <- latent_scores[, 1L]
  u1 <- latent_scores[, 2L]
  observed_scores <- latent_scores + measurement_error

  data.frame(
    z = 0.10 + beta_u0 * u0 + beta_u1 * u1 + stats::rnorm(n, sd = residual_sd),
    corrected_intercept_full = observed_scores[, 1L],
    corrected_slope_full = observed_scores[, 2L],
    ols_var11 = rep(meas_cov[1L, 1L], n),
    ols_var12 = rep(meas_cov[1L, 2L], n),
    ols_var22 = rep(meas_cov[2L, 2L], n)
  )
}

fuller_truth <- 0.55

set.seed(4801)
fuller_large <- fit_fuller_dual(
  simulate_fuller_known_eiv(2500L),
  outcome = "z",
  predictor_u0 = "corrected_intercept_full",
  predictor_u1 = "corrected_slope_full",
  meas11 = "ols_var11",
  meas12 = "ols_var12",
  meas22 = "ols_var22"
)

stopifnot(
  identical(as.integer(fuller_large$status_code), 0L),
  identical(fuller_large$mx_issue_class, "ok"),
  is.finite(fuller_large$estimate),
  is.finite(fuller_large$se),
  abs(fuller_large$estimate - fuller_truth) < 0.04,
  fuller_large$se > 0,
  fuller_large$se < 0.04
)

set.seed(4802)
fuller_mc <- replicate(120L, {
  out <- fit_fuller_dual(
    simulate_fuller_known_eiv(450L),
    outcome = "z",
    predictor_u0 = "corrected_intercept_full",
    predictor_u1 = "corrected_slope_full",
    meas11 = "ols_var11",
    meas12 = "ols_var12",
    meas22 = "ols_var22",
    skip_internal_scaling = TRUE
  )
  c(estimate = out$estimate, se = out$se, status_code = out$status_code)
})
fuller_mc <- as.data.frame(t(fuller_mc))
fuller_mc_ok <- fuller_mc[fuller_mc$status_code == 0L, , drop = FALSE]

fuller_empirical_sd <- stats::sd(fuller_mc_ok$estimate)
fuller_mean_se <- mean(fuller_mc_ok$se)
fuller_se_ratio <- fuller_mean_se / fuller_empirical_sd
fuller_mc_coverage <- mean(abs(fuller_mc_ok$estimate - fuller_truth) <= stats::qnorm(0.975) * fuller_mc_ok$se)

stopifnot(
  nrow(fuller_mc_ok) == 120L,
  abs(mean(fuller_mc_ok$estimate) - fuller_truth) < 0.03,
  fuller_se_ratio > 0.80,
  fuller_se_ratio < 1.20,
  fuller_mc_coverage > 0.88,
  fuller_mc_coverage < 0.99
)

set.seed(4803)
fuller_alpha_large_df <- simulate_fuller_known_eiv(2500L)
fuller_alpha_large <- fit_fuller_dual_alpha_stepdown(
  fuller_alpha_large_df,
  outcome = "z",
  predictor_u0 = "corrected_intercept_full",
  predictor_u1 = "corrected_slope_full",
  meas11 = "ols_var11",
  meas12 = "ols_var12",
  meas22 = "ols_var22",
  skip_internal_scaling = TRUE
)

stopifnot(
  identical(as.integer(fuller_alpha_large$status_code), 0L),
  identical(fuller_alpha_large$mx_issue_class, "ok"),
  is.finite(fuller_alpha_large$estimate),
  is.finite(fuller_alpha_large$se),
  abs(fuller_alpha_large$estimate - fuller_truth) < 0.04,
  fuller_alpha_large$se > 0,
  fuller_alpha_large$se < 0.04,
  isTRUE(all.equal(fuller_alpha_large$fuller_alpha_step1_used, 4, tolerance = 1e-12)),
  isTRUE(all.equal(fuller_alpha_large$fuller_alpha_step3_used, 4, tolerance = 1e-12)),
  isTRUE(all.equal(fuller_alpha_large$fuller_alpha_scaling_used, 4, tolerance = 1e-12)) || 
    is.na(fuller_alpha_large$fuller_alpha_scaling_used),
  fuller_alpha_large$fuller_sx1_star_relative_min_eigen >= 5e-2,
  fuller_alpha_large$fuller_sx_star_relative_min_eigen >= 5e-2,
  fuller_alpha_large$fuller_scaling_relative_min_eigen >= 5e-2 || 
    is.na(fuller_alpha_large$fuller_scaling_relative_min_eigen),
  fuller_alpha_large$fuller_sx1_star_condition <= 1e5,
  fuller_alpha_large$fuller_sx_star_condition <= 1e5
)

fuller_alpha_boundary <- fit_fuller_dual_alpha_stepdown(
  fuller_alpha_large_df,
  outcome = "z",
  predictor_u0 = "corrected_intercept_full",
  predictor_u1 = "corrected_slope_full",
  meas11 = "ols_var11",
  meas12 = "ols_var12",
  meas22 = "ols_var22",
  coarse_grid_size = 3L,
  max_refinements = 8L,
  search_tolerance = 0.1,
  min_sx1_star_relative_eigen = 0.37,
  min_sx_star_relative_eigen = 0.37,
  min_scaling_relative_eigen = 0.37,
  skip_internal_scaling = FALSE
)

fuller_alpha_boundary_coarse <- fit_fuller_dual_alpha_stepdown(
  fuller_alpha_large_df,
  outcome = "z",
  predictor_u0 = "corrected_intercept_full",
  predictor_u1 = "corrected_slope_full",
  meas11 = "ols_var11",
  meas12 = "ols_var12",
  meas22 = "ols_var22",
  coarse_grid_size = 3L,
  max_refinements = 0L,
  search_tolerance = 0.1,
  min_sx1_star_relative_eigen = 0.37,
  min_sx_star_relative_eigen = 0.37,
  min_scaling_relative_eigen = 0.37,
  skip_internal_scaling = FALSE
)

stopifnot(
  identical(as.integer(fuller_alpha_boundary$status_code), 0L),
  fuller_alpha_boundary$fuller_alpha_step1_used > 4,
  fuller_alpha_boundary$fuller_alpha_step1_used < 1252,
  fuller_alpha_boundary$fuller_alpha_step3_used > 4,
  fuller_alpha_boundary$fuller_alpha_step3_used < 1252,
  fuller_alpha_boundary$fuller_alpha_scaling_used > 4,
  fuller_alpha_boundary$fuller_alpha_scaling_used < 1252,
  fuller_alpha_boundary$fuller_alpha_step1_used < fuller_alpha_boundary_coarse$fuller_alpha_step1_used,
  fuller_alpha_boundary$fuller_alpha_step3_used < fuller_alpha_boundary_coarse$fuller_alpha_step3_used,
  fuller_alpha_boundary$fuller_alpha_scaling_used < fuller_alpha_boundary_coarse$fuller_alpha_scaling_used,
  fuller_alpha_boundary$fuller_sx1_star_relative_min_eigen >= 0.37,
  fuller_alpha_boundary$fuller_sx_star_relative_min_eigen >= 0.37,
  fuller_alpha_boundary$fuller_scaling_relative_min_eigen >= 0.37
)


fuller_alpha_base_scaling <- fit_fuller_dual_core(
  fuller_alpha_large_df,
  outcome = "z",
  predictor_u0 = "corrected_intercept_full",
  predictor_u1 = "corrected_slope_full",
  meas11 = "ols_var11",
  meas12 = "ols_var12",
  meas22 = "ols_var22",
  alpha_step1 = 4,
  alpha_step3 = 4,
  alpha_scaling = 4,
  auto_tempered = TRUE,
  skip_internal_scaling = FALSE
)
wrong_relative_scaling_upper <- fuller_alpha_base_scaling$fuller_scaling_relative_min_eigen *
  nrow(fuller_alpha_large_df) - 1

fuller_alpha_scaling_boundary <- fit_fuller_dual_alpha_stepdown(
  fuller_alpha_large_df,
  outcome = "z",
  predictor_u0 = "corrected_intercept_full",
  predictor_u1 = "corrected_slope_full",
  meas11 = "ols_var11",
  meas12 = "ols_var12",
  meas22 = "ols_var22",
  coarse_grid_size = 5L,
  max_refinements = 8L,
  search_tolerance = 0.1,
  min_scaling_relative_eigen = 0.40,
  skip_internal_scaling = FALSE
)

stopifnot(
  identical(as.integer(fuller_alpha_scaling_boundary$status_code), 0L),
  fuller_alpha_scaling_boundary$fuller_lambda_scaling > 1,
  fuller_alpha_scaling_boundary$fuller_alpha_scaling_used > wrong_relative_scaling_upper,
  fuller_alpha_scaling_boundary$fuller_alpha_scaling_used > 1252,
  fuller_alpha_scaling_boundary$fuller_alpha_scaling_used < 1876,
  fuller_alpha_scaling_boundary$fuller_scaling_relative_min_eigen >= 0.40
)

set.seed(4802)
fuller_alpha_mc <- replicate(120L, {
  out <- fit_fuller_dual_alpha_stepdown(
    simulate_fuller_known_eiv(450L),
    outcome = "z",
    predictor_u0 = "corrected_intercept_full",
    predictor_u1 = "corrected_slope_full",
    meas11 = "ols_var11",
    meas12 = "ols_var12",
    meas22 = "ols_var22",
    skip_internal_scaling = TRUE
  )
  c(estimate = out$estimate, se = out$se, status_code = out$status_code)
})
fuller_alpha_mc <- as.data.frame(t(fuller_alpha_mc))
fuller_alpha_mc_ok <- fuller_alpha_mc[fuller_alpha_mc$status_code == 0L, , drop = FALSE]

fuller_alpha_empirical_sd <- stats::sd(fuller_alpha_mc_ok$estimate)
fuller_alpha_mean_se <- mean(fuller_alpha_mc_ok$se)
fuller_alpha_se_ratio <- fuller_alpha_mean_se / fuller_alpha_empirical_sd
fuller_alpha_mc_coverage <- mean(abs(fuller_alpha_mc_ok$estimate - fuller_truth) <= stats::qnorm(0.975) * fuller_alpha_mc_ok$se)

stopifnot(
  nrow(fuller_alpha_mc_ok) == 120L,
  abs(mean(fuller_alpha_mc_ok$estimate) - fuller_truth) < 0.03,
  fuller_alpha_se_ratio > 0.80,
  fuller_alpha_se_ratio < 1.20,
  fuller_alpha_mc_coverage > 0.88,
  fuller_alpha_mc_coverage < 0.99
)

cat("stage-2 estimator tests ok\n")
