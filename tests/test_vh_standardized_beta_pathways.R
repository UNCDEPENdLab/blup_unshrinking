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

study1_design <- select_design("1")
study2_design <- select_design("2")
study3_design <- select_design("3")
all_design <- select_design("all")

stopifnot(
  nrow(study1_design) == 720L,
  all(c(0, 0.2, 0.4, 0.6) %in%
    unique(study1_design$standardized_beta_target)),
  max(abs(
    study1_design$structural_r2 -
      study1_design$standardized_beta_target^2
  )) < 1e-12,
  max(abs(
    study1_design$beta1w /
      sqrt(study1_design$slope_variance_marginal) -
      study1_design$standardized_beta_target
  )) < 1e-12,
  nrow(study3_design) == 1728L,
  all(abs(
    study3_design$achieved_reliability_y -
      study3_design$target_reliability_y
  ) < 1e-8),
  all(abs(
    study3_design$achieved_reliability_q -
      study3_design$target_reliability_q
  ) < 1e-8),
  max(abs(
    study3_design$theta1 *
      study3_design$tau1_y / study3_design$tau1_q -
      study3_design$standardized_beta_target
  )) < 1e-12,
  all(study3_design$slope_variance_residual_q > 0),
  all(study3_design$tau1_residual_q > 0),
  all(study3_design$theta1[
    study3_design$standardized_beta_target == 0
  ] == 0),
  all(study3_design$standardized_theta0[
    study3_design$structural_target == "slope_only"
  ] == 0),
  all(study3_design$standardized_theta0[
    study3_design$structural_target == "intercept_slope"
  ] == fixed_params$theta0_standardized),
  identical(
    study1_design$condition_id,
    all_design$condition_id[all_design$study == "study1"]
  ),
  identical(
    study2_design$condition_id,
    all_design$condition_id[all_design$study == "study2"]
  ),
  identical(
    study3_design$condition_id,
    all_design$condition_id[all_design$study == "study3"]
  ),
  length(intersect(study1_design$condition_id, study2_design$condition_id)) == 0L,
  length(intersect(study2_design$condition_id, study3_design$condition_id)) == 0L
)

resume_path <- tempfile(fileext = ".csv.gz")
resume_methods <- c("method_a", "method_b")
resume_rows <- tidyr::crossing(
  method = resume_methods,
  rep = 1:3
) %>%
  dplyr::mutate(pipeline_version = vh_pipeline_version())
data.table::fwrite(resume_rows, resume_path)
stopifnot(
  condition_output_is_complete(resume_path, resume_methods, 3L),
  !condition_output_is_complete(resume_path, resume_methods, 4L),
  !condition_output_is_complete(
    resume_path,
    resume_methods,
    3L,
    expected_pipeline_version = "older_pipeline"
  )
)
unlink(resume_path)

study3_condition <- study3_design %>%
  dplyr::filter(
    mean_clus_size_y == 10L,
    mean_clus_size_q == 10L,
    target_reliability_y == 0.8,
    target_reliability_q == 0.8,
    marginal_rho == 0.5,
    standardized_beta_target == 0.5,
    structural_target == "intercept_slope"
  ) %>%
  dplyr::slice(1L) %>%
  dplyr::mutate(num_clus = 10000L)

set.seed(20260615)
study3_sim <- simulate_study3(study3_condition)
oracle_fit <- stats::lm(
  true_q1 ~ true_y0 + true_y1,
  data = study3_sim$lv2_true
)
empirical_beta <- unname(stats::coef(oracle_fit)[["true_y1"]]) *
  study3_condition$tau1_y[[1]] / study3_condition$tau1_q[[1]]
empirical_q_cov <- stats::cov(
  study3_sim$lv2_true[, c("true_q0", "true_q1")]
)
target_q_cov <- make_random_effect_covariance(
  intercept_variance = fixed_params$tau0^2,
  slope_variance = study3_condition$slope_variance_marginal_q[[1]],
  intercept_slope_correlation = study3_condition$marginal_rho[[1]]
)

stopifnot(
  abs(empirical_beta - study3_condition$standardized_beta_target[[1]]) < 0.03,
  max(abs(empirical_q_cov - target_q_cov)) < 0.03
)

study1_condition <- study1_design %>%
  dplyr::filter(
    num_clus == 300L,
    mean_clus_size == 25L,
    target_reliability == 0.8,
    marginal_rho == 0.5,
    standardized_beta_target == 0.4
  ) %>%
  dplyr::slice(1L)
set.seed(20260616)
study1_results <- run_study1_rep(study1_condition)

study3_runner_condition <- study3_design %>%
  dplyr::filter(
    num_clus == 300L,
    mean_clus_size_y == 10L,
    mean_clus_size_q == 10L,
    target_reliability_y == 0.8,
    target_reliability_q == 0.8,
    marginal_rho == 0.5,
    standardized_beta_target == 0.5,
    structural_target == "intercept_slope"
  ) %>%
  dplyr::slice(1L)
set.seed(20260617)
study3_results <- run_study3_rep(study3_runner_condition)

stopifnot(
  setequal(study1_results$method, study1_methods()),
  all(study1_results$truth == study1_condition$standardized_beta_target),
  abs(
    study1_results$estimate[study1_results$method == "oracle"] -
      study1_condition$standardized_beta_target
  ) < 0.10,
  setequal(study3_results$method, study3_methods()),
  all(study3_results$truth ==
    study3_runner_condition$standardized_beta_target),
  all(c(
    "stage1_y_singular_problem",
    "stage1_q_singular_problem",
    "stage1_y_re_corr",
    "stage1_q_re_corr"
  ) %in% names(study3_results)),
  abs(
    study3_results$estimate[study3_results$method == "oracle_dual"] -
      study3_runner_condition$standardized_beta_target
  ) < 0.10
)

cat("VH standardized-beta pathway tests ok\n")
