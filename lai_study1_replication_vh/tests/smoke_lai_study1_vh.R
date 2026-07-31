#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(MASS)
  library(lme4)
  library(OpenMx)
  library(sandwich)
  library(geigen)
  library(readr)
})

find_repo_root <- function() {
  candidates <- unique(normalizePath(c(getwd(), file.path(getwd(), "..")), mustWork = FALSE))
  for (candidate in candidates) {
    if (file.exists(file.path(candidate, "R", "source_helpers.R")) &&
        file.exists(file.path(candidate, "lai_study1_replication_vh", "runner.R"))) {
      return(candidate)
    }
  }
  stop("Could not locate the repository root.")
}

repo_root <- find_repo_root()
refresh_dir <- file.path(repo_root, "lai_study1_replication_vh")

source(file.path(repo_root, "R", "source_helpers.R"), local = TRUE)
source_project_helpers(repo_root)
source(file.path(refresh_dir, "designs.R"), local = TRUE)
source(file.path(refresh_dir, "study1.R"), local = TRUE)
source(file.path(refresh_dir, "estimators.R"), local = TRUE)
source(file.path(refresh_dir, "runner.R"), local = TRUE)

design <- make_lai_study1_vh_design()
stopifnot(
  nrow(design) == 486L,
  all(is.finite(design$dgm_posterior_slope_reliability)),
  all(design$dgm_posterior_slope_reliability > 0 &
    design$dgm_posterior_slope_reliability < 1)
)
condition <- design %>% filter(beta_zu1 == 0.4) %>% slice(1L)

# The DGM reliability must use the full covariance matrix.  With independent
# random effects, Lai's RMS x scaling gives the closed-form scalar result.
independent_condition <- design %>%
  filter(num_clus == 30L, clus_size == 10L, icc == 0.2,
         vr_u1_u0 == 1, cor_u0_u1 == 0, beta_zu1 == 0.4) %>%
  slice(1L)
expected_independent_reliability <- with(
  independent_condition,
  var_u1 * clus_size / (sigma2 + var_u1 * clus_size)
)
stopifnot(all.equal(
  independent_condition$dgm_posterior_slope_reliability,
  expected_independent_reliability,
  tolerance = 1e-12
))

# Verify the correlated-cell reliability against an independent direct matrix
# calculation, rather than a scalar reliability shortcut.
correlated_condition <- design %>%
  filter(num_clus == 30L, clus_size == 10L, icc == 0.2,
         vr_u1_u0 == 2, cor_u0_u1 == 0.5, beta_zu1 == 0.4) %>%
  slice(1L)
Z_correlated <- lai_study1_vh_within_design(correlated_condition$clus_size)
G_correlated <- lai_study1_vh_dgm_covariance(
  correlated_condition$icc,
  correlated_condition$var_u1,
  correlated_condition$cor_u0_u1
)
V_correlated <- solve(
  solve(G_correlated) +
    crossprod(Z_correlated) / correlated_condition$sigma2
)
expected_correlated_reliability <- 1 -
  V_correlated[2, 2] / G_correlated[2, 2]
stopifnot(all.equal(
  correlated_condition$dgm_posterior_slope_reliability,
  expected_correlated_reliability,
  tolerance = 1e-12
))

set.seed(20260730)
sim <- simulate_lai_study1_vh(condition)
stopifnot(
  nrow(sim$lv1) == condition$num_clus * condition$clus_size,
  nrow(sim$lv2_true) == condition$num_clus,
  abs(sum(sim$lv1$x[seq_len(condition$clus_size)]^2) - condition$clus_size) < 1e-12
)

results <- run_lai_study1_vh_condition(
  condition,
  n_sim = 1L,
  methods = "oracle_dual",
  n_cores = 1L
)
stopifnot(
  identical(sort(unique(results$reporting_scale)), c("latent_sd", "raw")),
  all(results$method == "oracle_dual"),
  all(results$truth[results$reporting_scale == "raw"] == condition$beta_zu1),
  all(results$truth[results$reporting_scale == "latent_sd"] ==
    condition$beta_zu1 * sqrt(condition$var_u1)),
  all(results$dgm_population_slope_sd == sqrt(condition$var_u1)),
  all(results$dgm_posterior_slope_reliability ==
    condition$dgm_posterior_slope_reliability)
)

# Array-task resume checks require an exact current-pipeline file: every
# requested method/replication needs one raw and one common-latent-SD row.
array_test_dir <- tempfile("lai_study1_vh_array_")
dir.create(file.path(array_test_dir, "conditions"), recursive = TRUE)
array_result_file <- file.path(
  array_test_dir, "conditions",
  sprintf("condition_%04d_replications.csv.gz", condition$condition_id)
)
readr::write_csv(results, array_result_file)
stopifnot(lai_study1_vh_condition_output_is_complete(
  array_result_file, condition, n_sim = 1L, methods = "oracle_dual"
))
stale_results <- dplyr::mutate(results, pipeline_version = "obsolete")
readr::write_csv(stale_results, array_result_file)
stopifnot(!lai_study1_vh_condition_output_is_complete(
  array_result_file, condition, n_sim = 1L, methods = "oracle_dual"
))
readr::write_csv(results, array_result_file)
aggregate_test <- rebuild_lai_study1_vh_summary(
  n_sim = 1L,
  out_dir = array_test_dir,
  methods = "oracle_dual",
  condition_ids = condition$condition_id
)
stopifnot(
  nrow(aggregate_test) == 2L,
  file.exists(file.path(array_test_dir, "lai_study1_vh_summary.csv")),
  identical(
    select_lai_study1_vh_conditions(design, c(1L, 3L))$condition_id,
    c(1L, 3L)
  )
)

raw <- results %>% filter(reporting_scale == "raw") %>% arrange(method, rep)
latent <- results %>% filter(reporting_scale == "latent_sd") %>% arrange(method, rep)
ok <- is.finite(raw$estimate) & is.finite(latent$estimate)
stopifnot(all.equal(
  latent$estimate[ok],
  raw$estimate[ok] * latent$reporting_multiplier[ok],
  tolerance = 1e-12
))

# Replication-level fitted reliability must be computed from the same posterior
# slope variances and fitted G used by the Stage-1 component extractor.
diag_raw <- raw %>% slice(1L)
if (is.finite(diag_raw$stage1_fitted_slope_sd) &&
    is.finite(diag_raw$stage1_mean_posterior_slope_variance)) {
  expected_fitted_reliability <- 1 -
    diag_raw$stage1_mean_posterior_slope_variance /
    diag_raw$stage1_fitted_slope_sd^2
  stopifnot(all.equal(
    diag_raw$stage1_fitted_posterior_slope_reliability,
    expected_fitted_reliability,
    tolerance = 1e-12
  ))
}
if (is.finite(diag_raw$stage1_eb_slope_sd)) {
  stopifnot(all.equal(
    diag_raw$stage1_eb_to_population_slope_sd,
    diag_raw$stage1_eb_slope_sd / diag_raw$dgm_population_slope_sd,
    tolerance = 1e-12
  ))
}

# Check the diagnostic helper itself against the fitted Stage-1 components on
# a well-powered correlated cell.  This independently verifies the fitted SD,
# residual SD, posterior reliability, and EB-SD ratio fields.
diagnostic_condition <- design %>%
  filter(num_clus == 100L, clus_size == 10L, icc == 0.2,
         vr_u1_u0 == 1, cor_u0_u1 == 0.5, beta_zu1 == 0.4) %>%
  slice(1L)
set.seed(517)
diagnostic_sim <- simulate_lai_study1_vh(diagnostic_condition)
diagnostic_fit <- lai_study1_vh_fit_stage1(diagnostic_sim$lv1)
stopifnot(!is.null(diagnostic_fit))
diagnostic_stage1 <- get_stage1_eb_components(
  diagnostic_fit, diagnostic_sim$lv1, "cid", "y", "x"
)
diagnostic_corrected <- get_closed_form_corrected_scores(
  diagnostic_fit, diagnostic_sim$lv1, "cid", "y", "x"
)
diagnostic_stage2 <- diagnostic_sim$lv2_true %>%
  left_join(diagnostic_stage1, by = "id") %>%
  left_join(
    dplyr::select(
      diagnostic_corrected, id, corrected_intercept_full, corrected_slope_full,
      ols_var11, ols_var12, ols_var22
    ),
    by = "id"
  )
diagnostic_values <- lai_study1_vh_slope_diagnostics(
  diagnostic_fit, diagnostic_sim$lv1, diagnostic_stage2, diagnostic_condition
)
diagnostic_components <- extract_stage1_components(
  diagnostic_fit, diagnostic_sim$lv1, "cid", "x"
)
manual_fitted_variance <- diagnostic_components$G_hat[2, 2]
manual_mean_postvar22 <- mean(diagnostic_stage2$postvar22)
manual_residual_sd <- sqrt(mean(unlist(lapply(diagnostic_components$R_list, diag))))
stopifnot(
  all.equal(
    diagnostic_values$stage1_realized_true_slope_sd,
    sd(diagnostic_stage2$true_u1), tolerance = 1e-12
  ),
  all.equal(
    diagnostic_values$stage1_eb_slope_sd,
    sd(diagnostic_stage2$u1_eb), tolerance = 1e-12
  ),
  all.equal(
    diagnostic_values$stage1_fitted_slope_sd,
    sqrt(manual_fitted_variance), tolerance = 1e-12
  ),
  all.equal(
    diagnostic_values$stage1_fitted_residual_sd,
    manual_residual_sd, tolerance = 1e-12
  ),
  all.equal(
    diagnostic_values$stage1_fitted_posterior_slope_reliability,
    1 - manual_mean_postvar22 / manual_fitted_variance, tolerance = 1e-12
  ),
  all.equal(
    diagnostic_values$stage1_eb_to_population_slope_sd,
    sd(diagnostic_stage2$u1_eb) / diagnostic_condition$dgm_population_slope_sd,
    tolerance = 1e-12
  )
)

diagnostic_summary <- summarize_lai_study1_vh_results(results)
stopifnot(
  all(c(
    "dgm_population_slope_sd", "dgm_posterior_slope_reliability",
    "mean_eb_slope_sd", "mean_fitted_posterior_slope_reliability",
    "mean_lai_original_naive_eb_slope_sd",
    "mean_lai_original_method_multiplier"
  ) %in% names(diagnostic_summary)),
  all(diagnostic_summary$dgm_population_slope_sd == sqrt(condition$var_u1)),
  all(diagnostic_summary$dgm_posterior_slope_reliability ==
    condition$dgm_posterior_slope_reliability)
)

# Lai's historical view has a method-specific multiplier and is therefore not
# simply another copy of the common latent-SD scale.  The naïve EB multiplier
# follows raw-data ML (denominator J), while 2S-PA and MSEM use fitted latent
# slope SDs.
scale_fixture <- tibble::tibble(
  method = c("naive_dual_blup", "lai_2spa", "msem", "oracle_dual"),
  estimate = c(1, 2, 3, 4),
  se = rep(0.1, 4L),
  ci_low = c(0.8, 1.8, 2.8, 3.8),
  ci_high = c(1.2, 2.2, 3.2, 4.2),
  lai_fitted_latent_slope_sd = c(NA_real_, 2, NA_real_, NA_real_),
  mplus_fitted_latent_slope_sd = c(NA_real_, NA_real_, 3, NA_real_)
)
scale_fixture <- add_lai_study1_vh_original_scale_diagnostics(
  scale_fixture,
  tibble::tibble(u1_eb = c(1, 2, 3))
)
stopifnot(
  all.equal(
    scale_fixture$lai_original_naive_eb_slope_sd[[1]],
    sqrt(2 / 3), tolerance = 1e-12
  ),
  all.equal(scale_fixture$lai_original_naive_eb_slope_sample_sd[[1]], 1, tolerance = 1e-12)
)
scale_views <- add_lai_study1_vh_reporting_scales(scale_fixture, condition)
historical_views <- scale_views %>% filter(reporting_scale == "lai_original_standardized")
stopifnot(
  nrow(scale_views) == 2L * nrow(scale_fixture) + 3L,
  setequal(historical_views$method, c("naive_dual_blup", "lai_2spa", "msem")),
  all.equal(
    historical_views$estimate,
    c(1, 2, 3) * c(sqrt(2 / 3), 2, 3), tolerance = 1e-12
  ),
  all(historical_views$truth == condition$beta_zu1 * sqrt(condition$var_u1)),
  !any(scale_views$method == "oracle_dual" &
    scale_views$reporting_scale == "lai_original_standardized")
)

# The full non-Mplus bundle retains raw and common-latent-SD rows for every
# method.  Its optional historical rows are limited to methods with a direct
# Lai counterpart and a finite historical multiplier.
non_mplus_methods <- setdiff(lai_study1_vh_methods(), "msem")
bundle_results <- run_lai_study1_vh_condition(
  condition,
  n_sim = 1L,
  methods = non_mplus_methods,
  n_cores = 1L
)
stopifnot(
  setequal(unique(bundle_results$method), non_mplus_methods),
  nrow(bundle_results[bundle_results$reporting_scale != "lai_original_standardized", ]) ==
    2L * length(non_mplus_methods),
  all(bundle_results$method[bundle_results$reporting_scale == "lai_original_standardized"] %in%
    c("naive_dual_blup", "lai_2spa"))
)

cat("lai_study1_replication_vh smoke test ok\n")
