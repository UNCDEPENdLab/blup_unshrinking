#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
  library(tidyr)
})

source(file.path("R", "reliability_calibration.R"), local = TRUE)
source(file.path("R", "sim_helpers.R"), local = TRUE)
source(file.path("blup_outcome", "designs.R"), local = TRUE)
source(file.path("blup_outcome", "study_common.R"), local = TRUE)

smoke_design <- make_blup_outcome_design("smoke")
base_design <- make_blup_outcome_design("base")
ar1_design <- make_blup_outcome_design("residual_ar1")
reliability_smoke_design <- make_blup_outcome_design("posterior_reliability_smoke")
reliability_design <- make_blup_outcome_design("posterior_reliability")

stopifnot(
  nrow(smoke_design) == 1L,
  nrow(ar1_design) == 144L,
  all(c("balanced", "unbalanced", "informative_unbalanced") %in% unique(base_design$balance_mode)),
  all(c("iid", "ar1") %in% unique(ar1_design$r_structure)),
  all(c(0.3, 0.6) %in% unique(stats::na.omit(ar1_design$r_rho))),
  sum(ar1_design$r_structure == "iid") == 48L,
  sum(ar1_design$r_structure == "ar1" & ar1_design$r_rho == 0.3) == 48L,
  sum(ar1_design$r_structure == "ar1" & ar1_design$r_rho == 0.6) == 48L,
  identical(condition_to_r_spec(ar1_design[ar1_design$r_structure == "iid", ][1, ])$structure, "iid"),
  identical(condition_to_r_spec(ar1_design[ar1_design$r_rho == 0.3, ][1, ]), list(structure = "ar1", rho = 0.3)),
  0 %in% unique(base_design$gamma_x_on_slope),
  any(unique(base_design$rho) < 0),
  any(base_design$mean_n_trial <= 4L),
  identical(balance_mode_to_sim_arg("balanced"), "balanced"),
  identical(balance_mode_to_sim_arg("unbalanced"), FALSE),
  identical(balance_mode_to_sim_arg("informative_unbalanced"), "highly_unbalanced"),
  nrow(reliability_smoke_design) == 1L,
  nrow(reliability_design) == 720L,
  all(c(0.25, 0.50, 0.80) %in% unique(reliability_design$target_reliability)),
  all(c(0, 0.04, 0.16, 0.36) %in% unique(reliability_design$structural_r2)),
  max(abs(
    reliability_design$target_reliability -
      reliability_design$achieved_reliability
  )) < 1e-8,
  all(reliability_design$tau1 == reliability_design$tau1_residual),
  all(reliability_design$rho == reliability_design$marginal_rho),
  any(abs(reliability_design$rho_residual - reliability_design$marginal_rho) > 1e-8)
)

tau_sigma_pairs <- base_design %>%
  distinct(tau1, sigma)
expected_pairs <- tidyr::crossing(
  tau1 = sort(unique(base_design$tau1)),
  sigma = sort(unique(base_design$sigma))
)

stopifnot(
  nrow(tau_sigma_pairs) == nrow(expected_pairs),
  all(dplyr::anti_join(expected_pairs, tau_sigma_pairs, by = c("tau1", "sigma")) %>% nrow() == 0L)
)

stage2_schema <- ensure_blup_outcome_stage2_columns(data.frame(id = "a", x = 0))
stopifnot(
  all(c(
    "u0_eb", "u1_eb", "postvar11", "postvar12", "postvar22",
    "lambda11", "lambda12", "lambda21", "lambda22",
    "theta11", "theta12", "theta22",
    "corrected_z_var", "corrected_z_diag_var", "ols_var22"
  ) %in% names(stage2_schema)),
  all(is.na(stage2_schema[, c(
    "u0_eb", "u1_eb", "postvar11", "postvar12", "postvar22",
    "lambda11", "lambda12", "lambda21", "lambda22",
    "theta11", "theta12", "theta22",
    "corrected_z_var", "corrected_z_diag_var", "ols_var22"
  )]))
)

legacy_resolved <- resolve_blup_outcome_simulation_parameters(
  smoke_design,
  params = list(beta_0 = 1, beta_z = 0.6, tau0 = 0.9)
)
calibrated_resolved <- resolve_blup_outcome_simulation_parameters(
  reliability_smoke_design,
  params = list(beta_0 = 1, beta_z = 0.6, tau0 = 0.9)
)

stopifnot(
  !legacy_resolved$calibrated,
  identical(legacy_resolved$tau1, smoke_design$tau1[[1]]),
  identical(legacy_resolved$sim_params$rho, smoke_design$rho[[1]]),
  calibrated_resolved$calibrated,
  identical(
    calibrated_resolved$tau1,
    reliability_smoke_design$tau1_residual[[1]]
  ),
  identical(
    calibrated_resolved$sim_params$rho,
    reliability_smoke_design$rho_residual[[1]]
  )
)

set.seed(20260609)
calibrated_sim <- simulate_dataset(
  n_id = 1000L,
  mean_n_trial = reliability_smoke_design$mean_n_trial[[1]],
  params = calibrated_resolved$sim_params,
  tau1 = calibrated_resolved$tau1,
  sigma = reliability_smoke_design$sigma[[1]],
  has_random_slope = TRUE,
  balanced = balance_mode_to_sim_arg(reliability_smoke_design$balance_mode[[1]]),
  min_n_trial = reliability_smoke_design$min_n_trial[[1]],
  highly_unbalanced_min_n_trial =
    reliability_smoke_design$highly_unbalanced_min_n_trial[[1]],
  highly_unbalanced_power =
    reliability_smoke_design$highly_unbalanced_power[[1]],
  r_spec = condition_to_r_spec(reliability_smoke_design)
)
calibrated_diagnostics <- make_calibrated_dgp_diagnostics(
  condition = reliability_smoke_design,
  sim = calibrated_sim
)
legacy_diagnostics <- make_calibrated_dgp_diagnostics(
  condition = smoke_design,
  sim = calibrated_sim
)

stopifnot(
  all(is.finite(unlist(calibrated_diagnostics))),
  calibrated_diagnostics$empirical_g_error < 0.10,
  abs(
    calibrated_diagnostics$empirical_structural_r2 -
      reliability_smoke_design$structural_r2[[1]]
  ) < 0.05,
  abs(
    calibrated_diagnostics$empirical_reliability -
      reliability_smoke_design$target_reliability[[1]]
  ) < 1e-8,
  calibrated_diagnostics$residual_g_min_eigenvalue > 0,
  all(is.na(unlist(legacy_diagnostics)))
)

cat("BLUP-outcome design tests ok\n")
