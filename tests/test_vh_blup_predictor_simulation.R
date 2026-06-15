#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

source(file.path("R", "source_helpers.R"), local = TRUE)
source_project_helpers(".")
source(file.path("vig_hallquist_2026", "vh_study_common.R"), local = TRUE)

fixed_params <- list(
  z_variance = 1,
  tau0 = 0.9,
  beta0z = 1.5,
  beta1z = 0.4,
  gamma0_predictor = 0,
  gamma1_predictor = 0.5
)

balance_mode_to_sim_arg <- function(balance_mode) {
  switch(
    as.character(balance_mode[[1]]),
    balanced = "balanced",
    unbalanced = FALSE,
    informative_unbalanced = "highly_unbalanced",
    stop("Unsupported balance mode: ", balance_mode)
  )
}

calibrated <- calibrate_blup_predictor_condition(
  target_reliability = 0.50,
  standardized_slope_beta = 0.20,
  structural_target = "intercept_slope",
  marginal_rho = 0.50,
  tau0 = fixed_params$tau0,
  mean_n_trial = 5L,
  sigma = 1,
  nuisance_intercept_beta = fixed_params$beta1z,
  outcome_variance = fixed_params$z_variance
)

condition <- tibble(
  num_clus = 10000L,
  mean_clus_size = 5L,
  marginal_rho = 0.50,
  tau1 = calibrated$tau1_marginal,
  sigma = 1,
  beta1z = calibrated$beta1_intercept,
  beta2z = calibrated$beta2_slope,
  outcome_residual_variance = calibrated$outcome_residual_variance,
  balance_mode = "balanced",
  min_clus_size = 2L,
  r_structure = "iid",
  r_rho = NA_real_
)

set.seed(20260615)
sim <- simulate_data_blup_as_predictor(condition)
oracle_fit <- stats::lm(z ~ true_u0 + true_u1, data = sim$lv2_true)
oracle_estimates <- unname(stats::coef(oracle_fit)[2:3])

stopifnot(
  nrow(sim$lv2_true) == condition$num_clus,
  nrow(sim$lv1) == condition$num_clus * condition$mean_clus_size,
  max(abs(
    oracle_estimates -
      c(condition$beta1z, condition$beta2z)
  )) < 0.04,
  abs(stats::var(sim$lv2_true$z) - fixed_params$z_variance) < 0.04
)

cat("VH BLUP-as-predictor simulation test ok\n")
