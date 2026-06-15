#!/usr/bin/env Rscript

source(file.path("R", "reliability_calibration.R"), local = TRUE)

balanced <- calibrate_random_slope_condition(
  target_reliability = 0.25,
  standardized_beta = 0.60,
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
  standardized_beta = 0.40,
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
  standardized_beta = 0.40,
  marginal_rho = 0.30,
  tau0 = 0.9,
  mean_n_trial = 8L,
  sigma = 1.25,
  balance_mode = "unbalanced",
  min_n_trial = 2L,
  r_spec = list(structure = "ar1", rho = 0.30),
  n_reference = 1001L
)

predictor_conditions <- expand.grid(
  standardized_beta = c(0, 0.2, 0.4, 0.6),
  structural_target = c("slope_only", "intercept_slope"),
  marginal_rho = c(-0.5, 0, 0.5),
  stringsAsFactors = FALSE
)
predictor_calibrated <- lapply(seq_len(nrow(predictor_conditions)), function(i) {
  cell <- predictor_conditions[i, ]
  calibrate_blup_predictor_condition(
    target_reliability = 0.50,
    standardized_slope_beta = cell$standardized_beta,
    structural_target = cell$structural_target,
    marginal_rho = cell$marginal_rho,
    tau0 = 0.9,
    mean_n_trial = 5L,
    sigma = 1,
    nuisance_intercept_beta = 0.4,
    outcome_variance = 1
  )
})

G_y <- make_random_effect_covariance(
  intercept_variance = 0.9^2,
  slope_variance = 0.30,
  intercept_slope_correlation = 0.5
)
G_q <- make_random_effect_covariance(
  intercept_variance = 0.9^2,
  slope_variance = 0.20,
  intercept_slope_correlation = 0.5
)
dual_conditions <- expand.grid(
  standardized_beta = c(0, 0.2, 0.5),
  structural_target = c("slope_only", "intercept_slope"),
  stringsAsFactors = FALSE
)
dual_calibrated <- lapply(seq_len(nrow(dual_conditions)), function(i) {
  calibrate_dual_process_effect(
    G_predictor = G_y,
    G_outcome = G_q,
    standardized_slope_beta =
      dual_conditions$standardized_beta[[i]],
    structural_target = dual_conditions$structural_target[[i]],
    nuisance_intercept_standardized_beta = 0.4
  )
})

stopifnot(
  identical(unbalanced_a$gamma_x_on_slope, unbalanced_b$gamma_x_on_slope),
  identical(unbalanced_a$tau1_residual, unbalanced_b$tau1_residual),
  identical(unbalanced_a$rho_residual, unbalanced_b$rho_residual),
  abs(unbalanced_a$achieved_reliability - 0.50) < 1e-8,
  unbalanced_a$reference_min_n_trial == 5L,
  unbalanced_a$reference_max_n_trial == 11L,
  all(vapply(
    predictor_calibrated,
    function(x) abs(x$achieved_reliability - 0.50) < 1e-8,
    logical(1)
  )),
  all(vapply(
    predictor_calibrated,
    function(x) x$outcome_residual_variance > 0,
    logical(1)
  )),
  all(vapply(seq_along(predictor_calibrated), function(i) {
    abs(
      predictor_calibrated[[i]]$standardized_slope_beta -
        predictor_conditions$standardized_beta[[i]]
    ) < 1e-12
  }, logical(1))),
  all(vapply(seq_along(predictor_calibrated), function(i) {
    expected_beta1 <- if (predictor_conditions$structural_target[[i]] == "slope_only") 0 else 0.4
    abs(predictor_calibrated[[i]]$beta1_intercept - expected_beta1) < 1e-12
  }, logical(1))),
  all(vapply(which(predictor_conditions$standardized_beta == 0), function(i) {
    abs(predictor_calibrated[[i]]$beta2_slope) < 1e-12
  }, logical(1))),
  all(vapply(seq_along(dual_calibrated), function(i) {
    calibrated <- dual_calibrated[[i]]
    abs(
      calibrated$theta1_slope *
        sqrt(G_y[2, 2] / G_q[2, 2]) -
        dual_conditions$standardized_beta[[i]]
    ) < 1e-12
  }, logical(1))),
  all(vapply(seq_along(dual_calibrated), function(i) {
    calibrated <- dual_calibrated[[i]]
    max(abs(calibrated$G_joint_marginal[3:4, 3:4] - G_q)) < 1e-12
  }, logical(1))),
  all(vapply(dual_calibrated, function(x) {
    min(eigen(
      x$G_outcome_residual,
      symmetric = TRUE,
      only.values = TRUE
    )$values) > 0 &&
      min(eigen(
        x$G_joint_marginal,
        symmetric = TRUE,
        only.values = TRUE
      )$values) > 0
  }, logical(1))),
  all(vapply(which(dual_conditions$standardized_beta == 0), function(i) {
    abs(dual_calibrated[[i]]$theta1_slope) < 1e-12
  }, logical(1)))
)

cat("Reliability calibration tests ok\n")
