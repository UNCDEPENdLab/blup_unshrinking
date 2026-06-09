# Demonstrate calibration from posterior slope reliability and structural R-squared
# to the residual random-effects parameters expected by simulate_dataset().

repo_root <- dirname(rstudioapi::getActiveDocumentContext()$path)
source(file.path(repo_root, "R", "reliability_calibration.R"), local = TRUE)

calibration_row <- function(calibrated) {
  data.frame(
    achieved_reliability = calibrated$achieved_reliability,
    standardized_beta = calibrated$standardized_beta,
    gamma_x_on_slope = calibrated$gamma_x_on_slope,
    slope_variance_marginal = calibrated$slope_variance_marginal,
    slope_variance_residual = calibrated$slope_variance_residual,
    tau1_residual = calibrated$tau1_residual,
    rho_residual = calibrated$rho_residual
  )
}

run_examples <- function() {
  cat("\nExample 1: amended balanced simulation grid\n")
  grid <- expand.grid(
    cluster_size = c(3L, 5L, 10L, 25L),
    target_reliability = c(0.25, 0.50, 0.80),
    structural_r2 = c(0, 0.04, 0.16, 0.36),
    marginal_rho = c(-0.50, 0, 0.50)
  )

  calibrated_grid <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
    cell <- grid[i, , drop = FALSE]
    calibrated <- calibrate_random_slope_condition(
      target_reliability = cell$target_reliability,
      structural_r_squared = cell$structural_r2,
      marginal_rho = cell$marginal_rho,
      tau0 = 0.9,
      mean_n_trial = cell$cluster_size,
      sigma = 1,
      balance_mode = "balanced"
    )
    cbind(cell, calibration_row(calibrated))
  }))

  print(
    calibrated_grid[
      calibrated_grid$structural_r2 == 0.16 &
        calibrated_grid$marginal_rho == 0.50,
    ],
    row.names = FALSE,
    digits = 4
  )

  cat("\nExample 2: deterministic unbalanced AR(1) calibration\n")
  unbalanced <- calibrate_random_slope_condition(
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
  print(
    cbind(
      target_reliability = 0.50,
      structural_r2 = 0.16,
      marginal_rho = 0.30,
      calibration_row(unbalanced)
    ),
    row.names = FALSE,
    digits = 4
  )
  cat(
    "Reference count profile:",
    sprintf(
      "mean=%0.3f min=%d max=%d n=%d\n",
      unbalanced$reference_mean_n_trial,
      unbalanced$reference_min_n_trial,
      unbalanced$reference_max_n_trial,
      unbalanced$calibration_reference_n
    )
  )

  cat("\nValues passed to simulate_dataset() for this condition\n")
  print(list(
    gamma_x_on_slope = unbalanced$gamma_x_on_slope,
    tau1 = unbalanced$tau1_residual,
    rho = unbalanced$rho_residual,
    sigma = 1.25
  ))

  invisible(list(
    amended_grid = calibrated_grid,
    unbalanced = unbalanced
  ))
}

run_examples()

