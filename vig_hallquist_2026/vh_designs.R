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
    structural_target = c("slope_only", "intercept_slope"),
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
    structural_target = c("slope_only", "intercept_slope"),
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

make_study4_design <- function() {
  # ...
}

study_condition_counts <- function() {
  c(
    study1 = 5L * 4L * 3L * 3L * 4L,
    study2 = 4L * 3L * 3L * 5L * 4L * 2L,
    study3 = 4L * 3L * 3L * 2L * 2L * 2L * 3L * 2L,
    study4 = 0L
  )
}

select_design <- function(study_arg = "all", max_conditions = NA_integer_) {
  study_arg <- tolower(study_arg)
  requested <- if (identical(study_arg, "all")) {
    1:4
  } else {
    as.integer(unlist(strsplit(gsub("study", "", study_arg, fixed = TRUE), ",")))
  }
  if (anyNA(requested) || any(!(requested %in% 1:4))) {
    stop("No study conditions selected. Use `all`, `1`, `2`, `3`, `4`, or a comma-separated combination like `1,2`.")
  }

  builders <- list(
    make_study1_design,
    make_study2_design,
    make_study3_design,
    make_study4_design
  )
  counts <- study_condition_counts()
  offsets <- c(0L, cumsum(counts)[-length(counts)])
  names(offsets) <- names(counts)
  built <- lapply(requested, function(i) {
    study_name <- paste0("study", i)
    study_design <- builders[[i]]()
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
    stop("No study conditions selected. Use `all`, `1`, `2`, `3`, `4`, or a comma-separated combination like `1,2`.")
  }

  design
}
