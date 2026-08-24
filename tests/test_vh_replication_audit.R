#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

source(file.path(
  "vig_hallquist_2026",
  "random_effects_structural_simulation.R"
))

# Boundary/invalid Stage-1 covariance estimates must be recorded without the
# `cov2cor()` warning that previously flooded array logs.
valid_covariance <- assess_stage1_random_effect_covariance(
  matrix(c(2, 0.5, 0.5, 1), nrow = 2L, byrow = TRUE)
)
invalid_covariance <- assess_stage1_random_effect_covariance(
  matrix(c(NA_real_, 0, 0, -1), nrow = 2L, byrow = TRUE)
)
stopifnot(
  isTRUE(all.equal(valid_covariance$correlation, 0.5 / sqrt(2))),
  !valid_covariance$invalid,
  !valid_covariance$singular,
  is.na(invalid_covariance$correlation),
  invalid_covariance$invalid,
  invalid_covariance$singular
)

# A small hand-checkable Stage-1 example verifies that all nonrecoverable
# matrix entries and the two fitted reliability definitions are retained.
G_hat <- matrix(c(1, 0.2, 0.2, 0.5), nrow = 2L, byrow = TRUE)
scores <- tibble::tibble(
  true_u1 = c(-1, 0, 1),
  u1_eb = c(-0.8, 0.1, 0.9),
  corrected_slope_full = c(-1.2, 0.2, 1.1),
  lambda11 = c(0.7, 0.8, 0.9),
  lambda12 = c(0.01, 0.02, 0.03),
  lambda21 = c(0.04, 0.05, 0.06),
  lambda22 = c(0.6, 0.7, 0.8),
  theta11 = c(0.1, 0.2, 0.3),
  theta12 = c(-0.02, 0, 0.02),
  theta22 = c(0.2, 0.3, 0.4),
  postvar11 = c(0.2, 0.2, 0.2),
  postvar12 = c(0.02, 0.02, 0.02),
  postvar22 = c(0.1, 0.15, 0.2),
  ols_var11 = c(0.4, 0.5, 0.6),
  ols_var12 = c(-0.04, 0, 0.04),
  ols_var22 = c(0.7, 0.8, 0.9)
)
audit <- summarize_stage1_measurement_diagnostics(
  stage1_scores = scores,
  G_hat = G_hat,
  true_slope_col = "true_u1"
)

residualization <- G_hat[1L, 2L] / G_hat[1L, 1L]
residualized_variance <- G_hat[2L, 2L] -
  G_hat[1L, 2L]^2 / G_hat[1L, 1L]
residualized_postvar <- scores$postvar22 -
  2 * residualization * scores$postvar12 +
  residualization^2 * scores$postvar11
required_stage1_fields <- c(
  paste0("stage1_lambda", c("11", "12", "21", "22"), "_mean"),
  paste0("stage1_theta", c("11", "12", "22"), "_mean"),
  paste0("stage1_posterior_variance", c("11", "12", "22"), "_mean"),
  paste0(
    "stage1_corrected_score_error_covariance",
    c("11", "12", "22"),
    "_mean"
  ),
  "stage1_fitted_marginal_slope_posterior_reliability_mean",
  "stage1_fitted_residualized_slope_posterior_reliability_mean",
  "stage1_true_blup_slope_correlation",
  "stage1_true_blup_slope_r_squared",
  "stage1_blup_slope_bias",
  "stage1_blup_slope_rmse",
  "stage1_corrected_score_slope_bias",
  "stage1_corrected_score_slope_rmse"
)
stopifnot(
  all(required_stage1_fields %in% names(audit)),
  isTRUE(all.equal(
    audit$stage1_fitted_marginal_slope_posterior_reliability_mean,
    mean(1 - scores$postvar22 / G_hat[2L, 2L]),
    tolerance = 1e-12
  )),
  isTRUE(all.equal(
    audit$stage1_fitted_residualized_slope_posterior_reliability_mean,
    mean(1 - residualized_postvar / residualized_variance),
    tolerance = 1e-12
  )),
  audit$stage1_lambda12_mean == mean(scores$lambda12),
  audit$stage1_lambda21_mean == mean(scores$lambda21),
  audit$stage1_corrected_score_error_covariance12_mean ==
    mean(scores$ols_var12)
)

# Provenance and the deterministic crosswalk are separate, versioned artifacts.
bridge_chunk <- select_design("iccbridge", max_conditions = 2L) %>%
  slice_design_chunk(chunk_index = 1L, chunk_size = 1L)
provenance <- collect_vh_run_provenance(
  design = bridge_chunk,
  n_sim = 17L,
  study_arg = "iccbridge",
  chunk_meta = attr(bridge_chunk, "chunk_meta"),
  n_cores = 2L
)
required_provenance_fields <- c(
  "git_commit", "git_dirty", "pipeline_version", "requested_n_sim",
  "study_selector", "chunk_index", "chunk_size", "r_version",
  "package_openmx_version", "package_mplusautomation_version",
  "package_lme4_version", "package_nlme_version", "package_geigen_version",
  "mplus_version", "slurm_job_id", "slurm_array_job_id",
  "slurm_array_task_id", "run_started_at", "hostname"
)
stopifnot(
  nrow(provenance) == 1L,
  all(required_provenance_fields %in% names(provenance)),
  provenance$requested_n_sim == 17L,
  provenance$study_selector == "iccbridge"
)

audit_dir <- tempfile("vh-replication-audit-")
dir.create(audit_dir)
on.exit(unlink(audit_dir, recursive = TRUE, force = TRUE), add = TRUE)
crosswalk_path <- save_vh_icc_crosswalk(audit_dir)
crosswalk <- utils::read.csv(crosswalk_path)
stopifnot(
  basename(crosswalk_path) ==
    "icc_posterior_reliability_crosswalk_v1.csv",
  nrow(crosswalk) == 189L,
  all(c("crosswalk_version", "pipeline_version") %in% names(crosswalk)),
  all(crosswalk$crosswalk_version ==
    "icc_posterior_reliability_crosswalk_v1")
)

cat("VH replication diagnostics and provenance tests ok\n")
