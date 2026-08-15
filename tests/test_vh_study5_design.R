#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
})

source(file.path(
  "vig_hallquist_2026",
  "random_effects_structural_simulation.R"
))

design <- select_design("5")
first_stage <- design %>%
  distinct(
    calibration_arm,
    calibration_metric,
    mean_clus_size,
    marginal_rho,
    sigma,
    achieved_calibration_reliability,
    achieved_reliability,
    achieved_partial_reliability,
    slope_variance_marginal,
    slope_intercept_variance_ratio,
    G_condition_number
  )

stopifnot(
  nrow(design) == 36L,
  nrow(first_stage) == 12L,
  identical(sort(unique(design$condition_id)), 2417:2452),
  setequal(
    unique(design$calibration_arm),
    c(
      "current_g22",
      "shape_preserving_marginal",
      "shape_preserving_partial"
    )
  ),
  all(abs(design$achieved_calibration_reliability - 0.25) < 1e-8),
  all(design$outcome_residual_variance > 0),
  all(design$num_clus == 100L),
  setequal(unique(design$standardized_beta_target), c(0, 0.2, 0.4)),
  setequal(study5_methods(), study2_methods())
)

shape_rows <- first_stage %>%
  filter(grepl("^shape_preserving", calibration_arm))
current_correlated <- first_stage %>%
  filter(calibration_arm == "current_g22", marginal_rho == 0.5)
rho_zero <- first_stage %>%
  filter(marginal_rho == 0)

stopifnot(
  all(abs(shape_rows$slope_intercept_variance_ratio - 1) < 1e-12),
  all(shape_rows$G_condition_number < 4),
  all(current_correlated$achieved_partial_reliability < 0.05),
  all(current_correlated$G_condition_number > 200),
  all(abs(
    rho_zero$achieved_reliability - rho_zero$achieved_partial_reliability
  ) < 1e-10)
)

partial_condition <- design %>%
  filter(
    calibration_arm == "shape_preserving_partial",
    mean_clus_size == 10L,
    marginal_rho == 0.5,
    standardized_beta_target == 0.2
  ) %>%
  slice(1L)
set.seed(20260814)
sim <- simulate_study5(partial_condition)
stopifnot(
  nrow(sim$lv2_true) == 100L,
  nrow(sim$lv1) == 1000L,
  all(c("true_u0", "true_u1", "z") %in% names(sim$lv2_true))
)

cat("VH Study 5 design test ok\n")
