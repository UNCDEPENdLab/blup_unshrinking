library(data.table)
library(here)
library(tibble)
library(fst)

replication_select <- function() {
  c(
    "method",
    "estimate",
    "se",
    "ci_low",
    "ci_high",
    "status_code",
    "fuller_raw_estimate",
    "fuller_raw_se",
    "mx_issue_class",
    "mx_issue_detail",
    "fuller_lambda1",
    "fuller_lambda2",
    "fuller_sigma2",
    "fuller_weight_min",
    "fuller_weight_max",
    "fuller_correction_c",
    "fuller_lambda_scaling",
    "fuller_correction1",
    "fuller_correction_scaling",
    "fuller_measurement_weight_requested",
    "fuller_measurement_weight_used",
    "fuller_alpha_step1_requested",
    "fuller_alpha_step1_used",
    "fuller_alpha_step3_requested",
    "fuller_alpha_step3_used",
    "fuller_alpha_scaling_requested",
    "fuller_alpha_scaling_used",
    "fuller_auto_tempered",
    "fuller_sx1_star_condition",
    "fuller_sx1_star_min_eigen",
    "fuller_sx1_observed_max_eigen",
    "fuller_sx1_star_relative_min_eigen",
    "fuller_sx_star_condition",
    "fuller_sx_star_min_eigen",
    "fuller_sx_observed_max_eigen",
    "fuller_sx_star_relative_min_eigen",
    "fuller_scaling_condition",
    "fuller_scaling_min_eigen",
    "fuller_scaling_observed_max_eigen",
    "fuller_scaling_relative_min_eigen",
    "fuller_reference_se",
    "fuller_se_ratio",
    "fuller_auto_guard_pass",
    "fuller_auto_guard_reason",
    "fuller_auto_guard_score",
    "fuller_auto_full_weight_guard_pass",
    "fuller_auto_full_weight_guard_reason",
    "fuller_auto_full_weight_se_ratio",
    "fuller_auto_search_evaluations",
    "fuller_auto_search_nonmonotone",
    "mx_status_msg",
    "mx_info_definite",
    "mx_condition_number",
    "stage1_singular_problem",
    "stage1_problem_detail",
    "stage1_lmer_singular",
    "stage1_re_corr",
    "stage1_eb_corr",
    "stage1_design_kappa",
    "study",
    "truth",
    "condition_id",
    "num_clus",
    "mean_clus_size",
    "target_reliability",
    "marginal_rho",
    "standardized_beta_target",
    "balance_mode",
    "min_clus_size",
    "highly_unbalanced_min_clus_size",
    "highly_unbalanced_power",
    "r_structure",
    "r_rho",
    "sigma",
    "study_label",
    "study_structure",
    "calibration_tau0",
    "achieved_reliability",
    "standardized_beta",
    "structural_r2",
    "beta1w",
    "slope_variance_marginal",
    "slope_variance_residual",
    "tau1_residual",
    "rho_residual",
    "reference_mean_clus_size",
    "reference_min_clus_size",
    "reference_max_clus_size",
    "calibration_reference_n",
    "rho",
    "tau1",
    "rep",
    "pipeline_version"
  )
}

condition_summary_select <- function() {
  c(
    "condition_id",
    "study",
    "method",
    "num_clus",
    "mean_clus_size",
    "target_reliability",
    "achieved_reliability",
    "marginal_rho",
    "standardized_beta_target",
    "structural_r2",
    "r_rho",
    "sigma",
    "truth",
    "n_rep",
    "convergence",
    "mean_estimate",
    "mc_se_mean",
    "bias",
    "coverage",
    "rmse",
    "n_success",
    "n_status10_fail",
    "prop_status10_fail"
  )
}

manifest_select <- function () {
  c(
  "condition_id",
  "study",
  "num_clus",
  "mean_clus_size",
  "target_reliability",
  "marginal_rho",
  "standardized_beta_target",
  "balance_mode",
  "min_clus_size",
  "highly_unbalanced_min_clus_size",
  "highly_unbalanced_power",
  "r_structure",
  "r_rho",
  "sigma",
  "study_label",
  "study_structure",
  "calibration_tau0",
  "achieved_reliability",
  "standardized_beta",
  "structural_r2",
  "beta1w",
  "slope_variance_marginal",
  "slope_variance_residual",
  "tau1_residual",
  "rho_residual",
  "reference_mean_clus_size",
  "reference_min_clus_size",
  "reference_max_clus_size",
  "calibration_reference_n",
  "rho",
  "tau1",
  "structural_target",
  "beta1z",
  "beta2z",
  "focal_unique_r2",
  "outcome_residual_variance",
  "mean_clus_size_y",
  "mean_clus_size_q",
  "target_reliability_y",
  "target_reliability_q",
  "achieved_reliability_y",
  "slope_variance_marginal_y",
  "tau1_y",
  "reference_mean_clus_size_y",
  "reference_min_clus_size_y",
  "reference_max_clus_size_y",
  "achieved_reliability_q",
  "slope_variance_marginal_q",
  "tau1_q",
  "reference_mean_clus_size_q",
  "reference_min_clus_size_q",
  "reference_max_clus_size_q",
  "standardized_theta0",
  "theta0",
  "theta1",
  "slope_variance_residual_q",
  "tau1_residual_q",
  "rho_residual_q",
  "sigma_y",
  "sigma_q"
)
}

read_selected <- function(path, select, filter = NULL) {
  available <- names(fread(path, nrows = 0L, showProgress = FALSE))
  fread(path, select = intersect(select, available), showProgress = FALSE) %>%
    as_tibble() %>%
    {if (!is.null(filter)) filter(., !!filter) else .}
}

read_many_selected <- function(paths, select, missing_ok = FALSE, filter = NULL) {
  paths <- paths[file.exists(paths)]
  if (length(paths) == 0L) {
    if (missing_ok) {
      return(tibble())
    }
    stop("No input files found.", call. = FALSE)
  }
  rbindlist(
    lapply(paths, read_selected, select = select, filter = filter),
    fill = TRUE,
    use.names = TRUE
  )
}

project_dir <- function() {
  here()
}

output_dir <- function(this_run) {
  file.path(project_dir(), "outputs", this_run)
}

analysis_dir <- function() {
  file.path(project_dir(), "vig_hallquist_2026", "analysis")
}

conditions_dir <- function(this_run) {
  file.path(output_dir(this_run), "conditions")
}

condition_files <- function(this_run) {
  list.files(conditions_dir(this_run), pattern = "condition_[0-9]{4}_replication_results\\.csv\\.gz$", full.names = TRUE)
}

get_replication_results <- function(this_run, filter = NULL) {
  read_many_selected(condition_files(this_run), replication_select(), filter = filter) %>% as_tibble()
}

condition_summary_files <- function(this_run) {
  list.files(conditions_dir(this_run), pattern = "condition_[0-9]{4}_summary\\.csv$", full.names = TRUE)
}

get_condition_summary <- function(this_run, filter = NULL) {
  read_many_selected(condition_summary_files(this_run), condition_summary_select(), filter = filter) %>% as_tibble()
}

get_manifest <- function(this_run, study_arg = "all", filter = NULL) {
  manifest_files <- list.files(output_dir(this_run), pattern = sprintf("%s_chunk_[0-9]{3}_conditions_[0-9]{4}_[0-9]{4}_manifest\\.csv$", study_arg), full.names = TRUE)
  read_many_selected(manifest_files, manifest_select(), filter = filter) %>% as_tibble()
}

safe_mean <- function(x) {
  if (length(x) == 0L || all(is.na(x))) {
    return(NA_real_)
  }
  mean(x, na.rm = TRUE)
}

partition_cache_dir <- function(this_run = "vig_hallquist_slurm") {
  file.path(output_dir(this_run), "study_cache")
}

partition_cache_file <- function(study, this_run = "vig_hallquist_slurm") {
  file.path(partition_cache_dir(this_run), paste0(study, ".fst"))
}

read_study_cache <- function(study, this_run = "vig_hallquist_slurm") {
  cache_file <- partition_cache_file(study, this_run = this_run)
  if (!file.exists(cache_file)) {
    stop("Missing study cache: ", cache_file, call. = FALSE)
  }
  as_tibble(read_fst(cache_file))
}
