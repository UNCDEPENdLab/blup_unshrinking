#!/usr/bin/env Rscript

# Unit checks for formatting stacked-sandwich outputs into simulation estimator
# rows.

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

cat("stage-2 estimator formatter tests ok\n")
