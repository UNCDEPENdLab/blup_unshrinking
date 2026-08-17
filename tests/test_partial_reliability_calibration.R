#!/usr/bin/env Rscript

source(file.path("R", "reliability_calibration.R"), local = TRUE)

G <- make_random_effect_covariance(
  intercept_variance = 0.81,
  slope_variance = 0.81,
  intercept_slope_correlation = 0.5
)
Z_list <- list(make_reliability_time_design(25L))
R_shape_list <- list(diag(25L))

marginal <- calibrate_residual_scale(
  target_reliability = 0.25,
  G = G,
  Z_list = Z_list,
  R_shape_list = R_shape_list,
  reliability_measure = "marginal_slope"
)
partial <- calibrate_residual_scale(
  target_reliability = 0.25,
  G = G,
  Z_list = Z_list,
  R_shape_list = R_shape_list,
  reliability_measure = "residualized_slope"
)
icc_anchor <- calibrate_shape_preserving_measurement(
  target_value = 0.10,
  Z_list = list(make_reliability_time_design(10L)),
  R_shape_list = list(diag(10L)),
  intercept_variance = 0.81,
  slope_intercept_variance_ratio = 1,
  intercept_slope_correlation = 0,
  calibration_target = "intercept_icc"
)

stopifnot(
  abs(marginal$achieved_marginal_slope_reliability - 0.25) < 1e-8,
  abs(partial$achieved_residualized_slope_reliability - 0.25) < 1e-8,
  marginal$achieved_residualized_slope_reliability < 0.25,
  partial$achieved_marginal_slope_reliability > 0.25,
  marginal$sigma > partial$sigma,
  abs(icc_anchor$intercept_icc - 0.10) < 1e-12,
  abs(icc_anchor$achieved_marginal_slope_reliability - 0.50) < 1e-10,
  abs(icc_anchor$slope_intercept_variance_ratio - 1) < 1e-12,
  abs(
    expected_contrast_reliability(
      G,
      Z_list,
      partial$residual_covariances,
      contrast = c(-G[1, 2] / G[1, 1], 1)
    ) - 0.25
  ) < 1e-8
)

cat("Partial reliability calibration test ok\n")
