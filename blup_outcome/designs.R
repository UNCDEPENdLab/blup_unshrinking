#' Condition grids for BLUP-as-stage-2-outcome simulations.

blup_outcome_balance_modes <- function() {
  c("balanced", "unbalanced", "informative_unbalanced")
}

blup_outcome_residual_structures <- function() {
  tibble::tribble(
    ~r_structure, ~r_rho,
    "iid", NA_real_,
    "ar1", 0.3,
    "ar1", 0.6
  )
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

condition_to_r_spec <- function(condition) {
  r_structure <- if ("r_structure" %in% names(condition)) {
    as.character(condition$r_structure[[1]])
  } else {
    "iid"
  }

  switch(
    r_structure,
    iid = list(structure = "iid"),
    ar1 = list(structure = "ar1", rho = as.numeric(condition$r_rho[[1]])),
    stop("Unsupported residual structure: ", r_structure)
  )
}

condition_to_nlme_correlation <- function(condition) {
  r_structure <- if ("r_structure" %in% names(condition)) {
    as.character(condition$r_structure[[1]])
  } else {
    "iid"
  }

  switch(
    r_structure,
    iid = NULL,
    ar1 = {
      if (!requireNamespace("nlme", quietly = TRUE)) {
        stop("The `nlme` package is required for AR(1) Stage-1 residual covariance fits.")
      }
      nlme::corAR1(form = ~trial_index | id)
    },
    stop("Unsupported residual structure for nlme Stage-1 fit: ", r_structure)
  )
}

condition_uses_non_iid_R <- function(condition) {
  r_structure <- if ("r_structure" %in% names(condition)) {
    as.character(condition$r_structure[[1]])
  } else {
    "iid"
  }
  !identical(r_structure, "iid")
}

#' Add fixed posterior-reliability calibration parameters to a condition grid.
#'
#' Calibration is performed once while the grid is constructed. The resulting
#' population parameters are stored in the manifest and reused unchanged across
#' Monte Carlo replications.
calibrate_blup_outcome_reliability_design <- function(
    design,
    tau0 = 0.9,
    calibration_reference_n = 1001L) {
  required <- c(
    "target_reliability", "structural_r2", "marginal_rho",
    "mean_n_trial", "sigma", "balance_mode", "min_n_trial",
    "highly_unbalanced_min_n_trial", "highly_unbalanced_power",
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
      mean_n_trial = condition$mean_n_trial[[1]],
      sigma = condition$sigma[[1]],
      balance_mode = condition$balance_mode[[1]],
      min_n_trial = condition$min_n_trial[[1]],
      highly_unbalanced_min_n_trial = condition$highly_unbalanced_min_n_trial[[1]],
      highly_unbalanced_power = condition$highly_unbalanced_power[[1]],
      r_spec = condition_to_r_spec(condition),
      n_reference = calibration_reference_n
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
      reference_mean_n_trial = calibrated$reference_mean_n_trial,
      reference_min_n_trial = calibrated$reference_min_n_trial,
      reference_max_n_trial = calibrated$reference_max_n_trial,
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

make_blup_outcome_design <- function(
    grid_mode = "base",
    max_conditions = NA_integer_,
    reliability_tau0 = 0.9,
    calibration_reference_n = 1001L) {
  grid_mode <- as.character(grid_mode[[1]])

  design <- switch(
    grid_mode,
    posterior_reliability_smoke = tidyr::crossing(
      n_id = 40L,
      mean_n_trial = 8L,
      target_reliability = 0.50,
      structural_r2 = 0.16,
      marginal_rho = 0.50,
      balance_mode = "balanced",
      min_n_trial = 2L,
      highly_unbalanced_min_n_trial = 2L,
      highly_unbalanced_power = 3,
      sigma = 1.0,
      r_structure = "iid",
      r_rho = NA_real_
    ),
    posterior_reliability = tidyr::crossing(
      n_id = c(30L, 50L, 100L, 150L, 300L),
      mean_n_trial = c(3L, 5L, 10L, 25L),
      target_reliability = c(0.25, 0.50, 0.80),
      structural_r2 = c(0, 0.04, 0.16, 0.36),
      marginal_rho = c(-0.50, 0, 0.50),
      balance_mode = "balanced",
      min_n_trial = 2L,
      highly_unbalanced_min_n_trial = 2L,
      highly_unbalanced_power = 3,
      sigma = 1.0,
      r_structure = "iid",
      r_rho = NA_real_
    ),
    smoke = tidyr::crossing(
      n_id = 40L,
      mean_n_trial = 8L,
      gamma_x_on_slope = 0.3,
      rho = 0.3,
      balance_mode = "unbalanced",
      min_n_trial = 4L,
      highly_unbalanced_min_n_trial = 3L,
      highly_unbalanced_power = 3,
      tau1 = 0.7,
      sigma = 1.0
    ),
    residual_ar1 = tidyr::crossing(
      n_id = c(50L, 200L),
      mean_n_trial = c(8L, 20L, 50L),
      gamma_x_on_slope = c(0.1, 0.5),
      rho = c(0.0, 0.5),
      balance_mode = c("balanced", "unbalanced"),
      min_n_trial = 2L,
      highly_unbalanced_min_n_trial = 2L,
      highly_unbalanced_power = 3,
      tau1 = 0.7,
      sigma = 1.0
    ) %>%
      tidyr::crossing(blup_outcome_residual_structures()),
    base = tidyr::crossing(
      n_id = c(50L, 200L),
      mean_n_trial = c(4L, 8L, 20L, 50L),
      gamma_x_on_slope = c(0.0, 0.1, 0.5),
      rho = c(-0.5, 0.0, 0.5),
      balance_mode = blup_outcome_balance_modes(),
      min_n_trial = 2L,
      highly_unbalanced_min_n_trial = 2L,
      highly_unbalanced_power = 3,
      tau1 = c(0.3, 1.2),
      sigma = c(0.5, 1.5)
    ),
    heteroscedastic = tidyr::crossing(
      n_id = c(50L, 100L, 200L),
      mean_n_trial = c(8L, 20L, 50L),
      gamma_x_on_slope = c(0.0, 0.5),
      rho = c(-0.5, 0.5),
      balance_mode = c("balanced", "informative_unbalanced"),
      min_n_trial = 2L,
      highly_unbalanced_min_n_trial = 2L,
      highly_unbalanced_power = c(2, 4),
      tau1 = c(0.7, 1.2),
      sigma = c(1.0, 1.5)
    ),
    heteroscedastic_sparse = tidyr::crossing(
      n_id = c(50L, 100L, 200L),
      mean_n_trial = c(3L, 4L, 6L, 8L, 12L, 20L),
      gamma_x_on_slope = c(0.0, 0.1, 0.5),
      rho = c(-0.5, 0.0, 0.5),
      balance_mode = c("balanced", "unbalanced", "informative_unbalanced"),
      min_n_trial = 2L,
      highly_unbalanced_min_n_trial = 2L,
      highly_unbalanced_power = 3,
      tau1 = c(0.3, 0.7, 1.2),
      sigma = c(0.5, 1.0, 1.5)
    ),
    expanded = tidyr::crossing(
      n_id = c(30L, 50L, 100L, 200L, 500L),
      mean_n_trial = c(3L, 4L, 6L, 8L, 12L, 20L, 50L, 100L),
      gamma_x_on_slope = c(0.0, 0.1, 0.3, 0.5),
      rho = c(-0.5, 0.0, 0.5, 0.8),
      balance_mode = blup_outcome_balance_modes(),
      min_n_trial = 2L,
      highly_unbalanced_min_n_trial = 2L,
      highly_unbalanced_power = c(2, 3, 4),
      tau1 = c(0.3, 0.7, 1.2),
      sigma = c(0.5, 1.0, 1.5)
    ),
    stop(
      "`grid_mode` must be one of: posterior_reliability_smoke, ",
      "posterior_reliability, smoke, residual_ar1, base, ",
      "heteroscedastic, heteroscedastic_sparse, expanded."
    )
  )

  design <- design %>%
    dplyr::mutate(
      r_structure = if ("r_structure" %in% names(.)) .data$r_structure else "iid",
      r_rho = if ("r_rho" %in% names(.)) .data$r_rho else NA_real_
    )

  if (grid_mode %in% c("posterior_reliability_smoke", "posterior_reliability")) {
    design <- calibrate_blup_outcome_reliability_design(
      design,
      tau0 = reliability_tau0,
      calibration_reference_n = calibration_reference_n
    )
  }

  design <- design %>%
    dplyr::mutate(
      condition_id = dplyr::row_number(),
      balanced = balance_mode,
      design_source = grid_mode
    ) %>%
    dplyr::relocate(condition_id)

  if (!is.na(max_conditions)) {
    design <- design %>% dplyr::slice(seq_len(min(max_conditions, nrow(design))))
  }

  design
}
