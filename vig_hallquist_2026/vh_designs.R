#' Study-specific condition grids for Vig-Hallquist (2026)

fixed_params <- list(
  x_mean = 0,
  x_variance = 1,
  z_variance = 1,
  tau0 = 0.9,
  beta0z = 1.5,
  beta1z = 0.4,
  gamma0_outcome = 1.0,
  gamma1_outcome = 0.6,
  gamma0_predictor = 0.0,
  gamma1_predictor = 0.5
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
    "target_reliability", "structural_r2", "marginal_rho",
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
      structural_r_squared = condition$structural_r2[[1]],
      marginal_rho = condition$marginal_rho[[1]],
      tau0 = tau0,
      mean_n_trial = condition$mean_clus_size[[1]],
      sigma = condition$sigma[[1]],
      balance_mode = condition$balance_mode[[1]],
      min_n_trial = condition$min_clus_size[[1]],
      highly_unbalanced_min_n_trial = condition$highly_unbalanced_min_clus_size[[1]],
      highly_unbalanced_power = condition$highly_unbalanced_power[[1]],
      r_spec = condition_to_r_spec(condition),
      n_reference = calibration_reference_n,
      study_structure = condition$study_structure[[1]]
    )
    
    data.frame(
      calibration_tau0 = tau0,
      achieved_reliability = calibrated$achieved_reliability,
      standardized_beta = calibrated$standardized_beta,
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
    structural_r2 = c(0, 0.04, 0.16, 0.36),
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
  tidyr::crossing(
    study = "study2",
    num_clus = c(30L, 50L, 100L, 150L, 300L),
    mean_clus_size = c(3L, 5L, 10L, 25L),
    target_reliability = c(0.25, 0.5, 0.8),
    marginal_rho = c(-0.50, 0, 0.50),
    structural_r2 = c(0, 0.04, 0.16, 0.36),
    structural_target = c("slope_only", "intercept_slope"),
    balance_mode = "balanced",
    min_clus_size = 2L,
    highly_unbalanced_min_clus_size = 2L,
    highly_unbalanced_power = 3,
    r_structure = "iid",
    r_rho = NA_real_,
    sigma = 1.0,
    study_label = "BLUP as Predictor",
    study_structure = "z"
  ) %>%
    calibrate_reliability_design(
      tau0 = fixed_params$tau0,
      calibration_reference_n = 1001L
    ) %>%
    rename(
      beta2z = gamma_x_on_slope # level 2 structural slope is now "w"
    )
}

make_study3_design <- function() {
  # ...
}

make_study4_design <- function() {
  # ...
}

all_designs <- dplyr::bind_rows(
  make_study1_design(),
  make_study2_design(),
  make_study3_design(),
  make_study4_design()
)

select_design <- function(study_arg = "all", max_conditions = NA_integer_) {
  study_arg <- tolower(study_arg)
  design <- if (identical(study_arg, "all")) {
    all_designs
  } else {
    all_designs %>%
      dplyr::filter(study %in% paste0("study", unlist(strsplit(study_arg, ","))))
  }

  if (!is.na(max_conditions) && max_conditions > 0L) {
    design <- design %>% dplyr::slice_head(n = max_conditions)
  }

  if (nrow(design) == 0L) {
    stop("No study conditions selected. Use `all`, `1`, `2`, `3`, `4`, or a comma-separated combination like `1,2`.")
  }

  design %>%
    dplyr::mutate(condition_id = seq_len(dplyr::n())) %>%
    dplyr::relocate(condition_id)
}
