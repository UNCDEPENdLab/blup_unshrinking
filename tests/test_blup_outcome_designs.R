#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

source(file.path("blup_outcome", "designs.R"), local = TRUE)

smoke_design <- make_blup_outcome_design("smoke")
base_design <- make_blup_outcome_design("base")
ar1_design <- make_blup_outcome_design("residual_ar1")

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
  identical(balance_mode_to_sim_arg("informative_unbalanced"), "highly_unbalanced")
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

cat("BLUP-outcome design tests ok\n")
