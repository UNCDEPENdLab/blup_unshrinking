#!/usr/bin/env Rscript

# Unit checks for Mplus fitting functions

suppressPackageStartupMessages({
  library(MplusAutomation)
  library(dplyr)
})

source(file.path(
  "vig_hallquist_2026",
  "random_effects_structural_simulation.R"
))

first_stage_grid <- tidyr::tibble(
    mean_clus_size = 50L,
    target_reliability = 0.8,
    marginal_rho = 0.5,
    balance_mode = "balanced",
    min_clus_size = 2L,
    highly_unbalanced_min_clus_size = 2L,
    highly_unbalanced_power = 3,
    r_structure = "iid",
    r_rho = NA_real_,
    sigma = 1.0
  )

calibrated_first_stage <- lapply(seq_len(nrow(first_stage_grid)), function(i) {
  condition <- first_stage_grid[i, , drop = FALSE]
  reference <- make_reliability_reference_design(
    mean_n_trial = condition$mean_clus_size[[1]],
    sigma = condition$sigma[[1]],
    balance_mode = condition$balance_mode[[1]],
    min_n_trial = condition$min_clus_size[[1]],
    highly_unbalanced_min_n_trial = condition$highly_unbalanced_min_clus_size[[1]],
    highly_unbalanced_power = condition$highly_unbalanced_power[[1]],
    r_spec = condition_to_r_spec(condition),
    n_reference = 1001L
  )
  reliability <- calibrate_slope_variance(
    target_reliability = condition$target_reliability[[1]],
    Z_list = reference$Z_list,
    R_list = reference$R_list,
    weights = reference$count_weights,
    intercept_variance = fixed_params$tau0^2,
    intercept_slope_correlation = condition$marginal_rho[[1]]
  )
  dplyr::bind_cols(
    condition,
    tibble::tibble(
      calibration_tau0 = fixed_params$tau0,
      achieved_reliability = reliability$achieved_reliability,
      slope_variance_marginal = reliability$slope_variance_marginal,
      tau1 = sqrt(reliability$slope_variance_marginal),
      rho = condition$marginal_rho[[1]],
      G_marginal = list(reliability$G_marginal),
      reference_mean_clus_size = reference$mean_trial_count,
      reference_min_clus_size = reference$min_trial_count,
      reference_max_clus_size = reference$max_trial_count,
      calibration_reference_n = length(reference$trial_counts)
    )
  )
}) %>%
  dplyr::bind_rows()

design <- tidyr::tibble(
  calibrated_first_stage,
  study = "study2",
  num_clus = 2500L,
  standardized_beta_target = 0.4,
  structural_target = "intercept_slope", # only one structural target for now
  study_label = "BLUP as Predictor",
  study_structure = "z"
)

structural_rows <- lapply(seq_len(nrow(design)), function(i) {
  condition <- design[i, , drop = FALSE]
  structural <- calibrate_blup_predictor_effect(
    G_marginal = condition$G_marginal[[1]],
    standardized_slope_beta = condition$standardized_beta_target[[1]],
    structural_target = condition$structural_target[[1]],
    nuisance_intercept_beta = fixed_params$beta1z,
    outcome_variance = fixed_params$z_variance
  )
  tibble::tibble(
    standardized_beta = structural$standardized_slope_beta,
    beta1z = structural$beta1_intercept,
    beta2z = structural$beta2_slope,
    structural_r2 = structural$total_structural_r_squared,
    focal_unique_r2 = structural$focal_unique_r_squared,
    outcome_residual_variance = structural$outcome_residual_variance
  )
}) %>%
  dplyr::bind_rows()

condition <- dplyr::bind_cols(
  dplyr::select(design, -G_marginal),
  structural_rows
)

sim <- simulate_study2(condition)

fit <- fit_mplus_blup_predictor(
  level1_data = sim$lv1,
  level2_data = sim$lv2_true,
  outcome_variable = "y",
  within_component = "x",
  between_component = "z",
  cluster_id = "cid",
  reporting_scale = 1
) 

stopifnot(
  is.data.frame(fit),
  all(c("estimate", "se", "ci_low", "ci_high", "status_code", "mx_issue_class", "mx_issue_detail") %in% names(fit)),
  fit$status_code == 0L,
  fit$mx_issue_class == "ok",
  fit$mx_issue_detail == "ok",
  abs(fit$estimate - condition$beta2z) < 0.1
)


first_stage_grid <- tidyr::crossing(
  mean_clus_size = 50L,
  target_reliability = 0.8,
  marginal_rho = 0.5,
  balance_mode = "balanced",
  min_clus_size = 2L,
  highly_unbalanced_min_clus_size = 2L,
  highly_unbalanced_power = 3,
  r_structure = "iid",
  r_rho = NA_real_,
  sigma = 1.0
)

calibrated_first_stage <- lapply(seq_len(nrow(first_stage_grid)), function(i) {
  condition <- first_stage_grid[i, , drop = FALSE]
  reference <- make_reliability_reference_design(
    mean_n_trial = condition$mean_clus_size[[1]],
    sigma = condition$sigma[[1]],
    balance_mode = condition$balance_mode[[1]],
    min_n_trial = condition$min_clus_size[[1]],
    highly_unbalanced_min_n_trial =
      condition$highly_unbalanced_min_clus_size[[1]],
    highly_unbalanced_power = condition$highly_unbalanced_power[[1]],
    r_spec = condition_to_r_spec(condition),
    n_reference = 1001L
  )
  reliability <- calibrate_slope_variance(
    target_reliability = condition$target_reliability[[1]],
    Z_list = reference$Z_list,
    R_list = reference$R_list,
    weights = reference$count_weights,
    intercept_variance = fixed_params$tau0^2,
    intercept_slope_correlation = condition$marginal_rho[[1]]
  )
  dplyr::bind_cols(
    condition,
    tibble::tibble(
      achieved_reliability = reliability$achieved_reliability,
      slope_variance_marginal = reliability$slope_variance_marginal,
      tau1 = sqrt(reliability$slope_variance_marginal),
      G_marginal = list(reliability$G_marginal),
      reference_mean_clus_size = reference$mean_trial_count,
      reference_min_clus_size = reference$min_trial_count,
      reference_max_clus_size = reference$max_trial_count,
      calibration_reference_n = length(reference$trial_counts)
    )
  )
}) %>%
  dplyr::bind_rows()

design <- tidyr::crossing(
  study = "study3",
  num_clus = 2500L,
  mean_clus_size_y = 50L,
  mean_clus_size_q = 50L,
  target_reliability_y = 0.8,
  target_reliability_q = 0.8,
  marginal_rho = 0.5,
  standardized_beta_target = 0.5,
  structural_target = "intercept_slope", # only one structural target for now
  study_label = "BLUP as Predictor and Outcome",
  study_structure = "dual_process"
)

y_calibration <- calibrated_first_stage %>%
  dplyr::transmute(
    mean_clus_size_y = mean_clus_size,
    target_reliability_y = target_reliability,
    marginal_rho,
    achieved_reliability_y = achieved_reliability,
    slope_variance_marginal_y = slope_variance_marginal,
    tau1_y = tau1,
    G_y_marginal = G_marginal,
    reference_mean_clus_size_y = reference_mean_clus_size,
    reference_min_clus_size_y = reference_min_clus_size,
    reference_max_clus_size_y = reference_max_clus_size
  )
q_calibration <- calibrated_first_stage %>%
  dplyr::transmute(
    mean_clus_size_q = mean_clus_size,
    target_reliability_q = target_reliability,
    marginal_rho,
    achieved_reliability_q = achieved_reliability,
    slope_variance_marginal_q = slope_variance_marginal,
    tau1_q = tau1,
    G_q_marginal = G_marginal,
    reference_mean_clus_size_q = reference_mean_clus_size,
    reference_min_clus_size_q = reference_min_clus_size,
    reference_max_clus_size_q = reference_max_clus_size
  )

design <- design %>%
  dplyr::left_join(
    y_calibration,
    by = c("mean_clus_size_y", "target_reliability_y", "marginal_rho")
  ) %>%
  dplyr::left_join(
    q_calibration,
    by = c("mean_clus_size_q", "target_reliability_q", "marginal_rho")
  )

structural_rows <- lapply(seq_len(nrow(design)), function(i) {
  condition <- design[i, , drop = FALSE]
  structural <- calibrate_dual_process_effect(
    G_predictor = condition$G_y_marginal[[1]],
    G_outcome = condition$G_q_marginal[[1]],
    standardized_slope_beta = condition$standardized_beta_target[[1]],
    structural_target = condition$structural_target[[1]],
    nuisance_intercept_standardized_beta =
      fixed_params$theta0_standardized
  )
  tibble::tibble(
    standardized_beta = structural$standardized_slope_beta,
    standardized_theta0 = structural$standardized_intercept_beta,
    theta0 = structural$theta0_intercept,
    theta1 = structural$theta1_slope,
    structural_r2 = structural$total_structural_r_squared,
    focal_unique_r2 = structural$focal_unique_r_squared,
    slope_variance_residual_q =
      structural$outcome_slope_residual_variance,
    tau1_residual_q =
      sqrt(structural$outcome_slope_residual_variance),
    rho_residual_q = structural$outcome_residual_correlation
  )
}) %>%
  dplyr::bind_rows()

condition <- dplyr::bind_cols(
  dplyr::select(design, -G_y_marginal, -G_q_marginal),
  structural_rows
) %>%
  dplyr::mutate(
    calibration_tau0 = fixed_params$tau0,
    sigma_y = 1.0,
    sigma_q = 1.0,
    balance_mode = "balanced",
    min_clus_size = 2L,
    highly_unbalanced_min_clus_size = 2L,
    highly_unbalanced_power = 3,
    r_structure = "iid",
    r_rho = NA_real_
  )

sim <- simulate_study3(condition)

fit <- fit_mplus_dual_process(
  proc1_data = sim$lv1_y,
  proc2_data = sim$lv1_q,
  outcome1_var = "y", 
  outcome2_var = "q",
  cluster_id = "cid",
  time_index_var = "trial_index",
  time_value_var = "x",
  reporting_scale = 1
)

stopifnot(
  is.data.frame(fit),
  all(c("estimate", "se", "ci_low", "ci_high", "status_code", "mx_issue_class", "mx_issue_detail") %in% names(fit)),
  fit$status_code == 0L,
  fit$mx_issue_class == "ok",
  fit$mx_issue_detail == "ok",
  abs(fit$estimate - condition$theta1) < 0.1
)

cat("Mplus fitting function test ok\n")
