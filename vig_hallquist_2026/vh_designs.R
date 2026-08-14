#' Study-specific condition grids for Vig-Hallquist (2026)

fixed_params <- list(
  x_mean = 0,
  x_variance = 1,
  z_variance = 1,
  tau0 = 0.9,
  beta0z = 1.5,
  beta1z = 0.4,
  theta0_standardized = 0.4,
  gamma0_outcome = 1.0,
  gamma1_outcome = 0.6,
  gamma0_predictor = 0.0,
  gamma1_predictor = 0.5,
  gamma0_process_y = 0.0,
  gamma1_process_y = 0.5,
  gamma0_process_q = 0.5,
  gamma1_process_q = 0.4
)

cluster_balance_modes <- function() {
  c("balanced", "unbalanced", "informative_unbalanced")
}

balance_mode_to_sim_arg <- function(balance_mode) {
  switch(
    as.character(balance_mode[[1]]),
    balanced = "balanced",
    unbalanced = FALSE,
    informative_unbalanced = "highly_unbalanced",
    stop("Unsupported balance mode: ", balance_mode)
  )
}

residual_structures <- function() {
  tibble::tribble(
    ~r_structure, ~r_rho,
    "iid", NA_real_,
    "ar1", 0.3,
    "ar1", 0.6
  )
}

#' Add fixed posterior-reliability calibration parameters to a condition grid.
#'
#' Calibration is performed once while the grid is constructed. The resulting
#' population parameters are stored in the manifest and reused unchanged across
#' Monte Carlo replications.
#' 
#' NOTE: BLUP as outcome study called cluster size `n_trial` so there is some
#' mismatch in this function in terms of naming to be aware of.
calibrate_reliability_design <- function(
    design,
    tau0 = 0.9,
    calibration_reference_n = 1001L) {
  required <- c(
    "target_reliability", "standardized_beta_target", "marginal_rho",
    "mean_clus_size", "sigma", "balance_mode", "min_clus_size",
    "highly_unbalanced_min_clus_size", "highly_unbalanced_power",
    "r_structure", "r_rho"
  )
  missing <- setdiff(required, names(design))
  if (length(missing) > 0L) {
    stop("Reliability design is missing columns: ", paste(missing, collapse = ", "))
  }
  if (!exists("calibrate_random_slope_condition", mode = "function")) {
    stop(
      "`calibrate_random_slope_condition()` is unavailable. ",
      "Source `R/reliability_calibration.R` before building this grid."
    )
  }
  
  calibrated_rows <- lapply(seq_len(nrow(design)), function(i) {
    condition <- design[i, , drop = FALSE]
    calibrated <- calibrate_random_slope_condition(
      target_reliability = condition$target_reliability[[1]],
      standardized_beta = condition$standardized_beta_target[[1]],
      marginal_rho = condition$marginal_rho[[1]],
      tau0 = tau0,
      mean_n_trial = condition$mean_clus_size[[1]],
      sigma = condition$sigma[[1]],
      balance_mode = condition$balance_mode[[1]],
      min_n_trial = condition$min_clus_size[[1]],
      highly_unbalanced_min_n_trial = condition$highly_unbalanced_min_clus_size[[1]],
      highly_unbalanced_power = condition$highly_unbalanced_power[[1]],
      r_spec = condition_to_r_spec(condition),
      n_reference = calibration_reference_n
    )
    
    data.frame(
      calibration_tau0 = tau0,
      achieved_reliability = calibrated$achieved_reliability,
      standardized_beta = calibrated$standardized_beta,
      structural_r2 = calibrated$structural_r_squared,
      gamma_x_on_slope = calibrated$gamma_x_on_slope,
      slope_variance_marginal = calibrated$slope_variance_marginal,
      slope_variance_residual = calibrated$slope_variance_residual,
      tau1_residual = calibrated$tau1_residual,
      rho_residual = calibrated$rho_residual,
      reference_mean_clus_size = calibrated$reference_mean_n_trial,
      reference_min_clus_size = calibrated$reference_min_n_trial,
      reference_max_clus_size = calibrated$reference_max_n_trial,
      calibration_reference_n = calibrated$calibration_reference_n
    )
  })
  
  calibrated <- dplyr::bind_rows(calibrated_rows)
  dplyr::bind_cols(design, calibrated) %>%
    # `rho` remains the planned marginal correlation for reporting. The runner
    # explicitly passes `rho_residual` to simulate_dataset().
    dplyr::mutate(
      rho = .data$marginal_rho,
      tau1 = .data$tau1_residual
    )
}

# this is a smoke version of study 1 for testing
make_study0_design <- function() {
  tidyr::tibble(
    study = "study0",
    num_clus = 30L,
    mean_clus_size = 5L,
    target_reliability = 0.5,
    marginal_rho = 0.0,
    standardized_beta_target = 0.6,
    balance_mode = "balanced",
    min_clus_size = 2L,
    highly_unbalanced_min_clus_size = 2L,
    highly_unbalanced_power = 3,
    r_structure = "iid",
    r_rho = NA_real_,
    sigma = 1.0,
    study_label = "BLUP as Outcome (Smoke)",
    study_structure = "w"
  ) %>%
    calibrate_reliability_design(
      tau0 = fixed_params$tau0,
      calibration_reference_n = 1001L
    ) %>%
    rename(
      beta1w = gamma_x_on_slope # level 2 structural slope is now "w"
    )
}

make_study1_design <- function() {
  tidyr::crossing(
    study = "study1",
    num_clus = c(30L, 50L, 100L, 150L, 300L),
    mean_clus_size = c(3L, 5L, 10L, 25L),
    target_reliability = c(0.25, 0.5, 0.8),
    marginal_rho = c(-0.50, 0, 0.50),
    standardized_beta_target = c(0, 0.2, 0.4, 0.6),
    balance_mode = "balanced",
    min_clus_size = 2L,
    highly_unbalanced_min_clus_size = 2L,
    highly_unbalanced_power = 3,
    r_structure = "iid",
    r_rho = NA_real_,
    sigma = 1.0,
    study_label = "BLUP as Outcome",
    study_structure = "w"
    ) %>%
    calibrate_reliability_design(
      tau0 = fixed_params$tau0,
      calibration_reference_n = 1001L
    ) %>%
    rename(
      beta1w = gamma_x_on_slope # level 2 structural slope is now "w"
    )
}

make_study2_design <- function() {
  first_stage_grid <- tidyr::crossing(
    mean_clus_size = c(3L, 5L, 10L, 25L),
    target_reliability = c(0.25, 0.5, 0.8),
    marginal_rho = c(-0.50, 0, 0.50),
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

  design <- tidyr::crossing(
    calibrated_first_stage,
    study = "study2",
    num_clus = c(30L, 50L, 100L, 150L, 300L),
    standardized_beta_target = c(0, 0.2, 0.4, 0.6),
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

  dplyr::bind_cols(
    dplyr::select(design, -G_marginal),
    structural_rows
  )
}

make_study3_design <- function() {
  first_stage_grid <- tidyr::crossing(
    mean_clus_size = c(3L, 5L, 10L),
    target_reliability = c(0.25, 0.8),
    marginal_rho = c(0, 0.5),
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
    num_clus = c(50L, 100L, 150L, 300L),
    mean_clus_size_y = c(3L, 5L, 10L),
    mean_clus_size_q = c(3L, 5L, 10L),
    target_reliability_y = c(0.25, 0.8),
    target_reliability_q = c(0.25, 0.8),
    marginal_rho = c(0, 0.5),
    standardized_beta_target = c(0, 0.2, 0.5),
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

  dplyr::bind_cols(
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
}

#' Return the cluster-information specification for one Study 4 profile.
#'
#' Primary profiles all have mean cluster size 10. The information-matched
#' profile is a focused falsification control: its cluster sizes vary, but its
#' slope columns are rescaled during calibration and generation so that
#' `sum(x^2) = 9` for every cluster.
study4_profile_spec <- function(information_profile) {
  information_profile <- as.character(information_profile[[1]])
  switch(
    information_profile,
    homogeneous = list(
      sizes = 10L,
      weights = 1,
      information_matched = FALSE
    ),
    moderate = list(
      sizes = c(5L, 15L),
      weights = c(0.5, 0.5),
      information_matched = FALSE
    ),
    severe = list(
      sizes = c(3L, 17L),
      weights = c(0.5, 0.5),
      information_matched = FALSE
    ),
    severe_information_matched = list(
      sizes = c(3L, 17L),
      weights = c(0.5, 0.5),
      information_matched = TRUE
    ),
    stop("Unsupported Study 4 information profile: ", information_profile)
  )
}

#' Construct the random-effect design for one Study 4 cluster type.
study4_time_design <- function(cluster_size, information_matched = FALSE) {
  Z <- make_reliability_time_design(cluster_size)
  if (isTRUE(information_matched)) {
    Z[, "slope"] <- Z[, "slope"] * sqrt(9 / (cluster_size - 1))
  }
  Z
}

#' Compute weighted quantiles for a discrete Study 4 profile.
study4_weighted_quantile <- function(x, weights, probs) {
  keep <- is.finite(x) & is.finite(weights) & weights > 0
  x <- as.numeric(x[keep])
  weights <- as.numeric(weights[keep])
  if (length(x) == 0L) {
    return(rep(NA_real_, length(probs)))
  }
  ord <- order(x)
  x <- x[ord]
  cumulative_weight <- cumsum(weights[ord] / sum(weights))
  vapply(
    probs,
    function(prob) {
      x[[which(cumulative_weight >= min(1, max(0, prob)))[[1L]]]]
    },
    numeric(1L)
  )
}

#' Weighted RMS Frobenius distance of matrices from their elementwise mean.
study4_matrix_rms_dispersion <- function(matrices, weights = NULL) {
  if (is.null(weights)) {
    weights <- rep(1, length(matrices))
  }
  valid <- vapply(
    matrices,
    function(mat) is.matrix(mat) && all(dim(mat) == c(2L, 2L)) && all(is.finite(mat)),
    logical(1L)
  ) & is.finite(weights) & weights > 0
  if (!any(valid)) {
    return(NA_real_)
  }
  matrices <- matrices[valid]
  weights <- weights[valid] / sum(weights[valid])
  mean_matrix <- Reduce(
    `+`,
    Map(function(mat, weight) mat * weight, matrices, weights)
  )
  sqrt(sum(vapply(
    seq_along(matrices),
    function(i) weights[[i]] * sum((matrices[[i]] - mean_matrix)^2),
    numeric(1L)
  )))
}

#' Summarize population or fitted Study 4 measurement matrices.
study4_measurement_matrix_summary <- function(
    lambda_matrices,
    theta_matrices,
    score_error_covariances,
    weights = NULL) {
  if (is.null(weights)) {
    weights <- rep(1, length(lambda_matrices))
  }
  weights <- weights / sum(weights)
  lambda22 <- vapply(lambda_matrices, function(mat) mat[2L, 2L], numeric(1L))
  theta22 <- vapply(theta_matrices, function(mat) mat[2L, 2L], numeric(1L))
  ols_var22 <- vapply(
    score_error_covariances,
    function(mat) mat[2L, 2L],
    numeric(1L)
  )

  list(
    lambda22_mean = stats::weighted.mean(lambda22, weights),
    lambda22_min = min(lambda22),
    lambda22_max = max(lambda22),
    lambda_matrix_frobenius_rms_dispersion =
      study4_matrix_rms_dispersion(lambda_matrices, weights),
    theta22_mean = stats::weighted.mean(theta22, weights),
    theta22_min = min(theta22),
    theta22_max = max(theta22),
    theta_matrix_frobenius_rms_dispersion =
      study4_matrix_rms_dispersion(theta_matrices, weights),
    ols_var22_mean = stats::weighted.mean(ols_var22, weights),
    ols_var22_min = min(ols_var22),
    ols_var22_max = max(ols_var22),
    ols_cov_matrix_frobenius_rms_dispersion =
      study4_matrix_rms_dispersion(score_error_covariances, weights)
  )
}

#' Calibrate a Study 4 first-stage information profile.
calibrate_study4_profile <- function(
    information_profile,
    target_reliability,
    marginal_rho,
    sigma = 1.0) {
  spec <- study4_profile_spec(information_profile)
  Z_list <- lapply(
    spec$sizes,
    study4_time_design,
    information_matched = spec$information_matched
  )
  R_list <- lapply(
    spec$sizes,
    make_reliability_residual_covariance,
    sigma = sigma,
    r_spec = list(structure = "iid")
  )
  calibration <- calibrate_slope_variance(
    target_reliability = target_reliability,
    Z_list = Z_list,
    R_list = R_list,
    weights = spec$weights,
    intercept_variance = fixed_params$tau0^2,
    intercept_slope_correlation = marginal_rho
  )
  type_reliability <- vapply(
    calibration$posterior_covariances,
    function(V_i) {
      1 - V_i[2L, 2L] / calibration$G_marginal[2L, 2L]
    },
    numeric(1L)
  )
  weighted_mean <- stats::weighted.mean(type_reliability, spec$weights)
  weighted_variance <- stats::weighted.mean(
    (type_reliability - weighted_mean)^2,
    spec$weights
  )
  reliability_quartiles <- study4_weighted_quantile(
    type_reliability,
    spec$weights,
    c(0.25, 0.75)
  )
  population_measurement <- Map(function(Z_i, R_i) {
    sigma_i <- Z_i %*% calibration$G_marginal %*% t(Z_i) + R_i
    A_i <- t(solve(sigma_i, Z_i %*% calibration$G_marginal))
    R_inv_Z <- solve(R_i, Z_i)
    list(
      lambda = A_i %*% Z_i,
      theta = A_i %*% R_i %*% t(A_i),
      score_error_covariance = solve(crossprod(Z_i, R_inv_Z))
    )
  }, Z_list, R_list)
  measurement_diagnostics <- study4_measurement_matrix_summary(
    lambda_matrices = lapply(population_measurement, `[[`, "lambda"),
    theta_matrices = lapply(population_measurement, `[[`, "theta"),
    score_error_covariances = lapply(
      population_measurement,
      `[[`,
      "score_error_covariance"
    ),
    weights = spec$weights
  )

  list(
    calibration = calibration,
    profile = spec,
    type_reliability = type_reliability,
    reliability_mean = weighted_mean,
    reliability_sd = sqrt(weighted_variance),
    reliability_iqr = diff(reliability_quartiles),
    reliability_min = min(type_reliability),
    reliability_max = max(type_reliability),
    measurement_diagnostics = measurement_diagnostics
  )
}

make_study4_design <- function() {
  primary <- tidyr::crossing(
    study = "study4",
    num_clus = c(50L, 150L, 300L),
    target_reliability = c(0.25, 0.5, 0.8),
    information_profile = c("homogeneous", "moderate", "severe"),
    marginal_rho = c(0, 0.5),
    standardized_beta_target = c(0, 0.4),
    structural_target = "intercept_slope",
    is_falsification_control = FALSE
  )
  information_matched_control <- tidyr::crossing(
    study = "study4",
    num_clus = c(50L, 150L, 300L),
    target_reliability = 0.5,
    information_profile = "severe_information_matched",
    marginal_rho = 0,
    standardized_beta_target = 0.4,
    structural_target = "intercept_slope",
    is_falsification_control = TRUE
  )
  design <- dplyr::bind_rows(primary, information_matched_control)

  calibration_rows <- lapply(seq_len(nrow(design)), function(i) {
    condition <- design[i, , drop = FALSE]
    calibrated <- calibrate_study4_profile(
      information_profile = condition$information_profile[[1]],
      target_reliability = condition$target_reliability[[1]],
      marginal_rho = condition$marginal_rho[[1]],
      sigma = 1.0
    )
    cal <- calibrated$calibration
    spec <- calibrated$profile
    reliability <- calibrated$type_reliability
    measurement <- calibrated$measurement_diagnostics
    structural <- calibrate_blup_predictor_effect(
      G_marginal = cal$G_marginal,
      standardized_slope_beta = condition$standardized_beta_target[[1]],
      structural_target = condition$structural_target[[1]],
      nuisance_intercept_beta = fixed_params$beta1z,
      outcome_variance = fixed_params$z_variance
    )

    tibble::tibble(
      mean_clus_size = stats::weighted.mean(spec$sizes, spec$weights),
      profile_min_clus_size = min(spec$sizes),
      profile_max_clus_size = max(spec$sizes),
      profile_small_clus_size = min(spec$sizes),
      profile_large_clus_size = max(spec$sizes),
      profile_small_weight = spec$weights[[1L]],
      profile_large_weight = if (length(spec$weights) > 1L) {
        spec$weights[[length(spec$weights)]]
      } else {
        0
      },
      information_matched = spec$information_matched,
      calibration_tau0 = fixed_params$tau0,
      achieved_reliability = cal$achieved_reliability,
      reliability_sd = calibrated$reliability_sd,
      reliability_iqr = calibrated$reliability_iqr,
      reliability_min = calibrated$reliability_min,
      reliability_max = calibrated$reliability_max,
      reliability_small = reliability[[1L]],
      reliability_large = reliability[[length(reliability)]],
      population_lambda22_mean = measurement$lambda22_mean,
      population_lambda22_min = measurement$lambda22_min,
      population_lambda22_max = measurement$lambda22_max,
      population_lambda_matrix_frobenius_rms_dispersion =
        measurement$lambda_matrix_frobenius_rms_dispersion,
      population_theta22_mean = measurement$theta22_mean,
      population_theta22_min = measurement$theta22_min,
      population_theta22_max = measurement$theta22_max,
      population_theta_matrix_frobenius_rms_dispersion =
        measurement$theta_matrix_frobenius_rms_dispersion,
      population_ols_var22_mean = measurement$ols_var22_mean,
      population_ols_var22_min = measurement$ols_var22_min,
      population_ols_var22_max = measurement$ols_var22_max,
      population_ols_cov_matrix_frobenius_rms_dispersion =
        measurement$ols_cov_matrix_frobenius_rms_dispersion,
      slope_variance_marginal = cal$slope_variance_marginal,
      tau1 = sqrt(cal$slope_variance_marginal),
      rho = condition$marginal_rho[[1]],
      standardized_beta = structural$standardized_slope_beta,
      beta1z = structural$beta1_intercept,
      beta2z = structural$beta2_slope,
      structural_r2 = structural$total_structural_r_squared,
      focal_unique_r2 = structural$focal_unique_r_squared,
      outcome_residual_variance = structural$outcome_residual_variance
    )
  }) %>%
    dplyr::bind_rows()

  dplyr::bind_cols(design, calibration_rows) %>%
    dplyr::mutate(
      balance_mode = "heterogeneous_profile",
      min_clus_size = .data$profile_min_clus_size,
      highly_unbalanced_min_clus_size = .data$profile_min_clus_size,
      highly_unbalanced_power = NA_real_,
      r_structure = "iid",
      r_rho = NA_real_,
      sigma = 1.0,
      study_label = "Heterogeneous Cluster Information",
      study_structure = "z"
    )
}

#' Build Study 5: a matched-reliability calibration bridge.
#'
#' @details
#' All three arms target reliability .25, but they reach it through different
#' one-dimensional calibrations:
#'
#' - `current_g22` reproduces Study 2 by fixing `G00 = .81` and `sigma = 1`
#'   and solving for `G22` using marginal slope reliability;
#' - `shape_preserving_marginal` fixes `G22/G00 = 1` and solves for `sigma`
#'   using the same marginal reliability;
#' - `shape_preserving_partial` fixes the same covariance shape and solves for
#'   `sigma` using reliability of the slope residualized on the intercept.
#'
#' Holding the estimator bundle, structural effect, sample size, and target
#' information constant makes the arm contrasts diagnostic of covariance
#' geometry versus the reliability definition.
make_study5_design <- function() {
  target <- 0.25
  shape_variance <- fixed_params$tau0^2
  first_stage_grid <- tidyr::crossing(
    calibration_arm = c(
      "current_g22",
      "shape_preserving_marginal",
      "shape_preserving_partial"
    ),
    mean_clus_size = c(10L, 25L),
    marginal_rho = c(0, 0.5),
    target_calibration_reliability = target,
    target_reliability = target,
    balance_mode = "balanced",
    min_clus_size = 2L,
    highly_unbalanced_min_clus_size = 2L,
    highly_unbalanced_power = 3,
    r_structure = "iid",
    r_rho = NA_real_
  )

  calibrated_first_stage <- lapply(seq_len(nrow(first_stage_grid)), function(i) {
    condition <- first_stage_grid[i, , drop = FALSE]
    reference <- make_reliability_reference_design(
      mean_n_trial = condition$mean_clus_size[[1]],
      sigma = 1,
      balance_mode = condition$balance_mode[[1]],
      min_n_trial = condition$min_clus_size[[1]],
      highly_unbalanced_min_n_trial =
        condition$highly_unbalanced_min_clus_size[[1]],
      highly_unbalanced_power = condition$highly_unbalanced_power[[1]],
      r_spec = condition_to_r_spec(condition),
      n_reference = 1001L
    )

    arm <- condition$calibration_arm[[1]]
    if (identical(arm, "current_g22")) {
      calibration_metric <- "marginal_slope"
      sigma <- 1
      marginal_calibration <- calibrate_slope_variance(
        target_reliability = target,
        Z_list = reference$Z_list,
        R_list = reference$R_list,
        weights = reference$count_weights,
        intercept_variance = shape_variance,
        intercept_slope_correlation = condition$marginal_rho[[1]]
      )
      G <- marginal_calibration$G_marginal
      achieved_marginal <- marginal_calibration$achieved_reliability
      achieved_partial <- expected_residualized_slope_reliability(
        G,
        reference$Z_list,
        reference$R_list,
        weights = reference$count_weights
      )
    } else {
      calibration_metric <- if (identical(
        arm, "shape_preserving_marginal"
      )) "marginal_slope" else "residualized_slope"
      G <- make_random_effect_covariance(
        intercept_variance = shape_variance,
        slope_variance = shape_variance,
        intercept_slope_correlation = condition$marginal_rho[[1]]
      )
      residual_calibration <- calibrate_residual_scale(
        target_reliability = target,
        G = G,
        Z_list = reference$Z_list,
        R_shape_list = reference$R_list,
        weights = reference$count_weights,
        reliability_measure = calibration_metric
      )
      sigma <- residual_calibration$sigma
      achieved_marginal <-
        residual_calibration$achieved_marginal_slope_reliability
      achieved_partial <-
        residual_calibration$achieved_residualized_slope_reliability
    }

    residualized_variance <-
      G[2L, 2L] - G[1L, 2L]^2 / G[1L, 1L]
    achieved_calibration <- if (identical(
      calibration_metric, "marginal_slope"
    )) achieved_marginal else achieved_partial

    dplyr::bind_cols(
      condition,
      tibble::tibble(
        calibration_metric = calibration_metric,
        calibration_tau0 = fixed_params$tau0,
        sigma = sigma,
        achieved_calibration_reliability = achieved_calibration,
        achieved_reliability = achieved_marginal,
        achieved_partial_reliability = achieved_partial,
        slope_variance_marginal = G[2L, 2L],
        residualized_slope_variance = residualized_variance,
        slope_intercept_variance_ratio = G[2L, 2L] / G[1L, 1L],
        G_condition_number = kappa(G, exact = TRUE),
        conditional_slope_icc =
          residualized_variance / (residualized_variance + sigma^2),
        tau1 = sqrt(G[2L, 2L]),
        rho = condition$marginal_rho[[1]],
        G_marginal = list(G),
        reference_mean_clus_size = reference$mean_trial_count,
        reference_min_clus_size = reference$min_trial_count,
        reference_max_clus_size = reference$max_trial_count,
        calibration_reference_n = length(reference$trial_counts)
      )
    )
  }) %>%
    dplyr::bind_rows()

  design <- tidyr::crossing(
    calibrated_first_stage,
    study = "study5",
    num_clus = 100L,
    standardized_beta_target = c(0, 0.2, 0.4),
    structural_target = "intercept_slope",
    study_label = "Matched Reliability Calibration Bridge",
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

  dplyr::bind_cols(
    dplyr::select(design, -G_marginal),
    structural_rows
  )
}

study_condition_counts <- function() {
  c(
    study1 = 5L * 4L * 3L * 3L * 4L,
    study2 = 4L * 3L * 3L * 5L * 4L,
    study3 = 4L * 3L * 3L * 2L * 2L * 2L * 3L,
    study4 = (3L * 3L * 3L * 2L * 2L) + 3L,
    study0 = 1L,
    study5 = 3L * 2L * 2L * 3L
  )
}

select_design <- function(study_arg = "all", max_conditions = NA_integer_) {
  study_arg <- tolower(study_arg)
  requested <- if (identical(study_arg, "all")) {
    1:5
  } else {
    as.integer(unlist(strsplit(gsub("study", "", study_arg, fixed = TRUE), ",")))
  }
  if (anyNA(requested) || any(!(requested %in% 0:5))) {
    stop("No study conditions selected. Use `all`, `0`, `1`, `2`, `3`, `4`, `5`, or a comma-separated combination like `1,2`.")
  }

  builders <- list(
    make_study0_design,
    make_study1_design,
    make_study2_design,
    make_study3_design,
    make_study4_design,
    make_study5_design
  )
  counts <- study_condition_counts()
  offsets <- c(0L, cumsum(counts)[-length(counts)])
  names(offsets) <- names(counts)
  built <- lapply(requested, function(i) {
    study_name <- paste0("study", i)
    study_design <- builders[[i + 1]]()
    if (is.null(study_design) || nrow(study_design) == 0L) {
      return(tibble::tibble())
    }
    if (nrow(study_design) != counts[[study_name]]) {
      stop(
        "Canonical condition count is out of sync for ", study_name,
        ": expected ", counts[[study_name]], ", got ", nrow(study_design), "."
      )
    }
    study_design %>%
      dplyr::mutate(
        condition_id = offsets[[study_name]] + seq_len(dplyr::n())
      ) %>%
      dplyr::relocate(condition_id)
  })
  design <- dplyr::bind_rows(built)

  if (!is.na(max_conditions) && max_conditions > 0L) {
    design <- design %>% dplyr::slice_head(n = max_conditions)
  }

  if (nrow(design) == 0L) {
    stop("No study conditions selected. Use `all`, `0`, `1`, `2`, `3`, `4`, `5`, or a comma-separated combination like `1,2`.")
  }

  design
}
