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

target_scaled <- beta2 * stats::sd(x1_true)
stopifnot(
  isTRUE(all.equal(as.integer(fuller_out$status_code), 0L)),
  is.finite(fuller_out$estimate),
  is.finite(fuller_out$se),
  fuller_out$se > 0,
  abs(fuller_out$estimate - target_scaled) < 0.20 # not sure where 0.2 came from
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

ols_dual <- fit_observed_dual(stage2_noerr, outcome = "y_obs", predictor_u0 = "x0_obs", predictor_u1 = "x1_obs")
ols_naive <- ols_dual[ols_dual$se_type == "naive", , drop = FALSE]

stopifnot(
  isTRUE(all.equal(as.integer(fuller_noerr$status_code), 0L)),
  isTRUE(all.equal(fuller_noerr$estimate, ols_naive$estimate, tolerance = 1e-10)),
  isTRUE(all.equal(fuller_noerr$se, ols_naive$se, tolerance = 1e-10))
)

cat("stage-2 estimator formatter tests ok\n")
