#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
})

source(file.path("R", "sim_helpers.R"), local = TRUE)

params <- list(
  beta_0 = 0.2,
  beta_z = 0.5,
  gamma_x_on_slope = 0.1,
  tau0 = 0.6,
  tau1 = 0.4,
  tau = 0.6,
  rho = 0.2,
  sigma = 0.7,
  beta_1 = 0.3
)

set.seed(101)
iid_sim <- simulate_dataset(
  n_id = 3L,
  mean_n_trial = 4L,
  params = params,
  has_random_slope = TRUE,
  balanced = TRUE
)

stopifnot(
  length(iid_sim$R_list) == 3L,
  identical(names(iid_sim$R_list), as.character(seq_len(3L))),
  all(vapply(iid_sim$R_list, function(x) isTRUE(all.equal(x, params$sigma^2 * diag(4L))), logical(1))),
  "trial_index" %in% names(iid_sim$dat),
  all(iid_sim$dat %>% group_by(id) %>% summarise(ok = identical(trial_index, seq_len(n())), .groups = "drop") %>% pull(ok))
)

set.seed(102)
ar1_sim <- simulate_dataset(
  n_id = 2L,
  mean_n_trial = 5L,
  params = params,
  sigma = 1.2,
  has_random_slope = TRUE,
  balanced = TRUE,
  r_spec = list(structure = "ar1", rho = 0.5)
)
expected_ar1 <- 1.2^2 * outer(seq_len(5L), seq_len(5L), function(a, b) 0.5^abs(a - b))
stopifnot(
  identical(ar1_sim$r_spec$structure, "ar1"),
  isTRUE(all.equal(ar1_sim$R_list[[1]], expected_ar1))
)

set.seed(103)
toeplitz_sim <- simulate_dataset(
  n_id = 2L,
  mean_n_trial = 5L,
  params = params,
  sigma = 0.9,
  has_random_slope = TRUE,
  balanced = TRUE,
  r_spec = list(structure = "toeplitz", correlations = c(0.4, 0.1))
)
expected_toeplitz <- 0.9^2 * stats::toeplitz(c(1, 0.4, 0.1, 0, 0))
stopifnot(
  identical(toeplitz_sim$r_spec$structure, "toeplitz"),
  isTRUE(all.equal(toeplitz_sim$R_list[[1]], expected_toeplitz))
)

unsupported_alias_failed <- tryCatch({
  simulate_dataset(
    n_id = 2L,
    mean_n_trial = 4L,
    params = params,
    has_random_slope = TRUE,
    balanced = TRUE,
    r_spect = list(structure = "ar1", rho = 0.2)
  )
  FALSE
}, error = function(e) TRUE)
stopifnot(unsupported_alias_failed)

cat("Simulation residual covariance specification tests ok\n")
