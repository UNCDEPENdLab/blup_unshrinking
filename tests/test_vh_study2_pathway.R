#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
  library(MASS)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(foreach)
  library(doParallel)
  library(OpenMx)
  library(glmnet)
  library(sandwich)
  library(geigen)
})

source(file.path(
  "vig_hallquist_2026",
  "random_effects_structural_simulation.R"
))

design <- select_design("2")
stopifnot(
  nrow(design) == 1440L,
  all(design$outcome_residual_variance > 0),
  all(abs(design$achieved_reliability - design$target_reliability) < 1e-8),
  all(design$beta1z[design$structural_target == "slope_only"] == 0),
  all(
    design$beta1z[design$structural_target == "intercept_slope"] ==
      fixed_params$beta1z
  )
)

condition <- design %>%
  dplyr::filter(
    num_clus == 300L,
    mean_clus_size == 25L,
    target_reliability == 0.8,
    marginal_rho == 0.5,
    standardized_beta_target == 0.4,
    structural_target == "intercept_slope"
  ) %>%
  dplyr::slice(1L)

set.seed(20260615)
results <- run_study2_rep(condition)
primary_dual_methods <- c(
  "oracle_dual",
  "naive_dual_blup",
  "closed_form_dual",
  "fuller_closed_form",
  "fuller_alpha_stepdown_closed_form",
  "lai_2spa",
  "lai_2spaa"
)

stopifnot(
  setequal(results$method, study2_methods()),
  all(results$truth == condition$standardized_beta_target),
  all(primary_dual_methods %in% results$method),
  all(
    results$method_role[results$method %in% primary_dual_methods[-1]] ==
      "primary_dual"
  ),
  all(
    results$method_role[grepl("_slope", results$method)] ==
      "slope_only_diagnostic"
  ),
  abs(
    results$estimate[results$method == "oracle_dual"] -
      condition$standardized_beta_target
  ) < 0.08,
  results$estimate[results$method == "naive_slope_blup"] >
    results$estimate[results$method == "naive_dual_blup"]
)

failed_results <- make_failed_result(
  condition,
  study2_methods(),
  condition$standardized_beta_target[[1]]
) %>%
  add_study2_method_roles()
stopifnot(
  !anyNA(failed_results$method_role),
  setequal(failed_results$method, results$method)
)

summary <- summarize_results_df(
  dplyr::bind_cols(
    results,
    condition[rep(1L, nrow(results)), , drop = FALSE] %>%
      dplyr::select(-study),
    tibble::tibble(rep = 1L)
  )
)
stopifnot(
  nrow(summary) == length(study2_methods()),
  all(c("mean_clus_size", "standardized_beta_target", "structural_target") %in%
    names(summary))
)

cat("VH Study 2 pathway test ok\n")
