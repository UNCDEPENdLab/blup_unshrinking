#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

source(file.path("blup_outcome", "designs.R"), local = TRUE)

smoke_design <- make_blup_outcome_design("smoke")
base_design <- make_blup_outcome_design("base")

stopifnot(
  nrow(smoke_design) == 1L,
  all(c("balanced", "unbalanced", "informative_unbalanced") %in% unique(base_design$balance_mode)),
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
