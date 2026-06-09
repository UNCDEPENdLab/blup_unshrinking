#!/usr/bin/env Rscript

source(file.path("R", "reliability_calibration.R"), local = TRUE)

balanced <- calibrate_random_slope_condition(
  target_reliability = 0.25,
  structural_r_squared = 0.36,
  marginal_rho = 0.50,
  tau0 = 0.9,
  mean_n_trial = 25L,
  sigma = 1,
  balance_mode = "balanced"
)

stopifnot(
  abs(balanced$achieved_reliability - 0.25) < 1e-8,
  abs(balanced$standardized_beta - 0.60) < 1e-12,
  abs(
    balanced$gamma_x_on_slope^2 /
      balanced$slope_variance_marginal -
      0.36
  ) < 1e-10,
  balanced$rho_residual != balanced$marginal_rho,
  min(eigen(
    balanced$G_residual,
    symmetric = TRUE,
    only.values = TRUE
  )$values) > 0
)

unbalanced_a <- calibrate_random_slope_condition(
  target_reliability = 0.50,
  structural_r_squared = 0.16,
  marginal_rho = 0.30,
  tau0 = 0.9,
  mean_n_trial = 8L,
  sigma = 1.25,
  balance_mode = "unbalanced",
  min_n_trial = 2L,
  r_spec = list(structure = "ar1", rho = 0.30),
  n_reference = 1001L
)
unbalanced_b <- calibrate_random_slope_condition(
  target_reliability = 0.50,
  structural_r_squared = 0.16,
  marginal_rho = 0.30,
  tau0 = 0.9,
  mean_n_trial = 8L,
  sigma = 1.25,
  balance_mode = "unbalanced",
  min_n_trial = 2L,
  r_spec = list(structure = "ar1", rho = 0.30),
  n_reference = 1001L
)

stopifnot(
  identical(unbalanced_a$gamma_x_on_slope, unbalanced_b$gamma_x_on_slope),
  identical(unbalanced_a$tau1_residual, unbalanced_b$tau1_residual),
  identical(unbalanced_a$rho_residual, unbalanced_b$rho_residual),
  abs(unbalanced_a$achieved_reliability - 0.50) < 1e-8,
  unbalanced_a$reference_min_n_trial == 5L,
  unbalanced_a$reference_max_n_trial == 11L
)

cat("Reliability calibration tests ok\n")
