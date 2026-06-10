#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
})

find_repo_root <- function() {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  script_path <- if (length(script_arg) > 0L) {
    normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = FALSE)
  } else {
    normalizePath(file.path("lai_replication", "tests", "smoke_lai_study4.R"), mustWork = FALSE)
  }

  candidates <- unique(normalizePath(c(
    getwd(),
    file.path(getwd(), ".."),
    file.path(dirname(script_path), "..", "..")
  ), mustWork = FALSE))

  for (candidate in candidates) {
    if (file.exists(file.path(candidate, "lai_replication", "mlm_random_slope_lai_apples_to_apples_sim.R"))) {
      return(candidate)
    }
  }

  getwd()
}

repo_root <- find_repo_root()
driver_path <- file.path(repo_root, "lai_replication", "mlm_random_slope_lai_apples_to_apples_sim.R")

if (!file.exists(driver_path)) {
  stop("Could not find Lai apples-to-apples driver from repo root: ", repo_root)
}

source(driver_path, local = TRUE)

study4_design <- select_design("4", max_conditions = 2L)
required_design_cols <- c("study", "r_structure", "r_rho", "condition_id")
missing_design_cols <- setdiff(required_design_cols, names(study4_design))
if (length(missing_design_cols) > 0L) {
  stop("Missing expected Study 4 design columns: ", paste(missing_design_cols, collapse = ", "))
}

if (!identical(unique(study4_design$study), "study4")) {
  stop("Study 4 design selection returned the wrong study label.")
}

if (!setequal(unique(study4_design$r_structure), "ar1")) {
  stop("Study 4 design should only use the ar1 residual structure.")
}

if (!setequal(sort(unique(study4_design$r_rho)), c(0.3, 0.6))) {
  stop("Study 4 design should expose both Study 4 residual autocorrelations.")
}

study4_condition <- study4_design[1, ]
sim <- simulate_study4(study4_condition)

r_spec <- lai_condition_to_r_spec(study4_condition)
expected_r_spec <- list(structure = "ar1", rho = as.numeric(study4_condition$r_rho))
if (!identical(r_spec, expected_r_spec)) {
  stop("Study 4 condition did not normalize to the expected AR(1) residual spec.")
}

expected_r_mat <- make_R_matrix(
  as.integer(study4_condition$clus_size),
  sigma = sqrt(study4_condition$sigma2),
  r_spec = r_spec
)
expected_cor_mat <- outer(
  seq_len(as.integer(study4_condition$clus_size)),
  seq_len(as.integer(study4_condition$clus_size)),
  function(a, b) r_spec$rho^abs(a - b)
)
if (!isTRUE(all.equal(expected_r_mat, study4_condition$sigma2 * expected_cor_mat))) {
  stop("Study 4 residual covariance matrix does not match the expected AR(1) form.")
}

resid_lookup <- match(as.character(sim$lv1$cid), sim$lv2_true$id)
resid_draws <- sim$lv1$y - fixed_params$gamma0 - fixed_params$gamma1 * sim$lv1$x -
  sim$lv2_true$true_u0[resid_lookup] - sim$lv2_true$true_u1[resid_lookup] * sim$lv1$x
resid_matrix <- matrix(resid_draws, nrow = as.integer(study4_condition$num_clus), byrow = TRUE)
emp_cov <- stats::cov(resid_matrix)
emp_cor <- stats::cor(resid_matrix)
expected_resid_cor <- outer(
  seq_len(as.integer(study4_condition$clus_size)),
  seq_len(as.integer(study4_condition$clus_size)),
  function(a, b) as.numeric(study4_condition$r_rho)^abs(a - b)
)

if (abs(mean(diag(emp_cov)) - study4_condition$sigma2) > 0.25 * study4_condition$sigma2) {
  stop("Empirical Study 4 residual variance is too far from the AR(1) target.")
}

emp_lag1 <- mean(c(emp_cor[1L, 2L], emp_cor[2L, 3L]))
emp_lag2 <- emp_cor[1L, 3L]

if (abs(emp_lag1 - as.numeric(study4_condition$r_rho)) > 0.30) {
  stop("Empirical Study 4 lag-1 correlation is too far from the AR(1) target.")
}

if (abs(emp_lag2 - as.numeric(study4_condition$r_rho)^2) > 0.30) {
  stop("Empirical Study 4 lag-2 correlation is too far from the AR(1) target.")
}

if (max(abs(emp_cor - expected_resid_cor)) > 0.30) {
  stop("Empirical Study 4 residual correlation matrix is too far from the AR(1) target.")
}

expected_lv1_n <- as.integer(study4_condition$num_clus) * as.integer(study4_condition$clus_size)
if (nrow(sim$lv1) != expected_lv1_n) {
  stop("Study 4 lv1 row count does not match the selected condition.")
}

if (nrow(sim$lv2_true) != as.integer(study4_condition$num_clus)) {
  stop("Study 4 lv2_true row count does not match the selected condition.")
}

required_lv1_cols <- c("cid", "cid_chr", "trial_index", "x", "y")
required_lv2_cols <- c("id", "z", "true_u0", "true_u1")
missing_lv1_cols <- setdiff(required_lv1_cols, names(sim$lv1))
missing_lv2_cols <- setdiff(required_lv2_cols, names(sim$lv2_true))
if (length(missing_lv1_cols) > 0L) {
  stop("Missing expected lv1 columns: ", paste(missing_lv1_cols, collapse = ", "))
}
if (length(missing_lv2_cols) > 0L) {
  stop("Missing expected lv2_true columns: ", paste(missing_lv2_cols, collapse = ", "))
}

if (!identical(sim$lv1$cid_chr, as.character(sim$lv1$cid))) {
  stop("Study 4 lv1 cluster labels are not internally consistent.")
}

rep_results <- run_study4_rep(study4_condition)

required_result_cols <- c("study", "method", "estimate", "se", "ci_low", "ci_high", "status_code", "truth")
missing_result_cols <- setdiff(required_result_cols, names(rep_results))
if (length(missing_result_cols) > 0L) {
  stop("Missing expected Study 4 result columns: ", paste(missing_result_cols, collapse = ", "))
}

expected_methods <- matched_study_methods(include_tempered_eiv = FALSE)
if (!setequal(unique(rep_results$method), expected_methods)) {
  stop("Study 4 replication methods do not match the expected matched-study set.")
}

if (nrow(rep_results) != length(expected_methods)) {
  stop("Study 4 replication results should contain one row per expected method.")
}

if (!identical(unique(rep_results$study), "study4")) {
  stop("Study 4 replication results should be labeled study4.")
}

if (!isTRUE(all.equal(unique(rep_results$truth), lai_truth(study4_condition)))) {
  stop("Study 4 replication truth does not match the design-implied target.")
}

cat("lai study 4 smoke test ok\n")