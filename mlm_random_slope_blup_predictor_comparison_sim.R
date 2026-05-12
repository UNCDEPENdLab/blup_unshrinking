#!/usr/bin/env Rscript

#' ---
#' title: "BLUP Correction: Predicting Secondary Outcomes"
#' description: "Simulates a scenario where the true latent intercept and slope
#'               predict a set of 'downstream' or secondary outcomes. It demonstrates
#'               that correcting the BLUPs (using the full 2x2 unweighting) eliminates
#'               attenuation bias when the BLUPs are used as independent variables in a
#'               secondary regression."
#' ---

suppressPackageStartupMessages({
  library(lme4)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(foreach)
  library(doParallel)
  library(OpenMx)
})

# Load helpers
# Load helpers relative to script directory
script_dir <- if (basename(getwd()) == "mlm_blups") "." else "mlm_blups"
source(file.path(script_dir, "R", "sim_helpers.R"), local = TRUE)
source(file.path(script_dir, "R", "stats_helpers.R"), local = TRUE)
source(file.path(script_dir, "R", "blup_helpers.R"), local = TRUE)

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
OpenMx::mxOption(NULL, "Number of Threads", 1L)

# Argument handling
args <- commandArgs(trailingOnly = TRUE)
n_sim <- if (length(args) >= 1) as.integer(args[[1]]) else 100L
out_dir <- if (length(args) >= 2) args[[2]] else file.path(script_dir, "outputs", "blup_predictor_comparison")
n_cores <- if (length(args) >= 3) as.integer(args[[3]]) else 1L

if (is.na(n_sim) || n_sim < 1L) {
  stop("`n_sim` must be a positive integer.")
}

if (is.na(n_cores) || n_cores < 1L) {
  stop("`n_cores` must be a positive integer.")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Grid setup
sim_grid <- tidyr::crossing(
  n_id = c(50L, 100L, 200L),
  mean_n_trial = c(8L, 20L, 50L, 100L),
  tibble::tibble(
    tau1 = c(0.7, 0.3, 1.2),
    sigma = c(1.0, 1.5, 0.5)
  )
)

params <- list(
  beta_0 = 1.0,
  beta_z = 0.6,
  gamma_x_on_slope = 0.5,
  tau0 = 0.9,
  rho = 0.5,
  alpha_single = 0.4,
  beta_single_u1 = 0.5,
  alpha_dual = 0.4,
  beta_dual_u0 = 0.3,
  beta_dual_u1 = 0.5
)

#' Execute a single replication of the BLUP predictor comparison simulation.
#' In this simulation, the BLUPs are used as predictors in a secondary downstream regression
#' for both a single-predictor model (slope only) and a dual-predictor model (intercept + slope).
#' @param n_id Number of subjects.
#' @param mean_n_trial Target mean trials per subject.
#' @param tau1 Standard deviation of the random slope.
#' @param sigma Residual standard deviation.
#' @param params List of fixed base parameters.
#' @return A tibble with extraction results for Oracle, Naive BLUP, and Corrected BLUP in
#'         both single and dual predictor settings.
run_one_rep <- function(n_id, mean_n_trial, tau1, sigma, params) {
  sim <- simulate_dataset(
    n_id = n_id,
    mean_n_trial = mean_n_trial,
    params = params,
    tau1 = tau1,
    sigma = sigma,
    has_random_slope = TRUE,
    balanced = FALSE
  )

  # Data generation for this specific simulation includes secondary outcomes
  var_u0 <- params$tau0^2
  var_u1_total <- params$gamma_x_on_slope^2 + tau1^2
  cov_u0_u1 <- params$rho * params$tau0 * tau1

  resid_var_single <- max(0.1, 1 - params$beta_single_u1^2 * var_u1_total)
  resid_var_dual <- max(0.1, 1 - (params$beta_dual_u0^2 * var_u0 + params$beta_dual_u1^2 * var_u1_total + 2 * params$beta_dual_u0 * params$beta_dual_u1 * cov_u0_u1))

  id_df <- sim$id_df %>%
    mutate(
      outcome_single = params$alpha_single + params$beta_single_u1 * true_slope_dev + stats::rnorm(n(), mean = 0, sd = sqrt(resid_var_single)),
      outcome_dual = params$alpha_dual + params$beta_dual_u0 * true_intercept_dev + params$beta_dual_u1 * true_slope_dev + stats::rnorm(n(), mean = 0, sd = sqrt(resid_var_dual))
    )

  fit_null <- suppressWarnings(suppressMessages(
    lmer(y ~ 1 + z + (1 + z | id), data = sim$dat, REML = FALSE)
  ))

  score_df <- id_df %>%
    left_join(get_corrected_scores(fit_null), by = "id")

  # Regression 1: Single predictor (slope)
  res_single_oracle <- lm(outcome_single ~ true_slope_dev, data = score_df)
  res_single_blup <- lm(outcome_single ~ blup_z, data = score_df)
  res_single_corr <- lm(outcome_single ~ corrected_z, data = score_df)

  # Regression 2: Dual predictors (intercept and slope)
  res_dual_oracle <- lm(outcome_dual ~ true_intercept_dev + true_slope_dev, data = score_df)
  res_dual_blup <- lm(outcome_dual ~ blup_intercept + blup_z, data = score_df)
  res_dual_corr <- lm(outcome_dual ~ corrected_intercept + corrected_z, data = score_df)

  bind_rows(
    single_oracle = extract_lm_stats(res_single_oracle, "true_slope_dev"),
    single_blup = extract_lm_stats(res_single_blup, "blup_z"),
    single_corr = extract_lm_stats(res_single_corr, "corrected_z"),
    dual_oracle_u0 = extract_lm_stats(res_dual_oracle, "true_intercept_dev"),
    dual_oracle_u1 = extract_lm_stats(res_dual_oracle, "true_slope_dev"),
    dual_blup_u0 = extract_lm_stats(res_dual_blup, "blup_intercept"),
    dual_blup_u1 = extract_lm_stats(res_dual_blup, "blup_z"),
    dual_corr_u0 = extract_lm_stats(res_dual_corr, "corrected_intercept"),
    dual_corr_u1 = extract_lm_stats(res_dual_corr, "corrected_z"),
    .id = "method"
  )
}

# Execution
if (n_cores > 1L) {
  registerDoParallel(cores = n_cores)
}

set.seed(20260411)
results <- foreach(i = seq_len(nrow(sim_grid)), .combine = bind_rows, .packages = c("lme4", "dplyr", "tidyr", "purrr", "tibble")) %dopar% {
  n_id_i <- sim_grid$n_id[[i]]
  mean_n_trial_i <- sim_grid$mean_n_trial[[i]]
  tau1_i <- sim_grid$tau1[[i]]
  sigma_i <- sim_grid$sigma[[i]]

  map_dfr(seq_len(n_sim), function(rep_id) {
    out <- run_one_rep(n_id_i, mean_n_trial_i, tau1_i, sigma_i, params)
    mutate(out, rep = rep_id, n_id = n_id_i, mean_n_trial = mean_n_trial_i, tau1 = tau1_i, sigma = sigma_i)
  })
}

if (n_cores > 1L) {
  stopImplicitCluster()
}

# Analysis
results <- results %>%
  mutate(
    truth = case_when(
      method %in% c("single_oracle", "single_blup", "single_corr") ~ params$beta_single_u1,
      grepl("u0$", method) ~ params$beta_dual_u0,
      grepl("u1$", method) ~ params$beta_dual_u1,
      TRUE ~ NA_real_
    ),
    bias = estimate - truth,
    covered = ci_low <= truth & ci_high >= truth,
    sq_error = (estimate - truth)^2
  )

summary_df <- results %>%
  group_by(method, n_id, mean_n_trial, tau1, sigma) %>%
  summarise(
    mean_estimate = mean(estimate, na.rm = TRUE),
    mc_se_mean = sd(estimate, na.rm = TRUE) / sqrt(sum(!is.na(estimate))),
    bias = mean(bias, na.rm = TRUE),
    coverage = mean(covered, na.rm = TRUE),
    rmse = sqrt(mean(sq_error, na.rm = TRUE)),
    .groups = "drop"
  )

# Plotting helpers
summary_df <- summary_df %>%
  mutate(
    method_type = case_when(
      grepl("^single", method) ~ "Single Predictor",
      grepl("_u0$", method) ~ "Dual: Intercept",
      grepl("_u1$", method) ~ "Dual: Slope"
    ),
    estimator = case_when(
      grepl("oracle", method) ~ "Oracle",
      grepl("blup", method) ~ "Naive BLUP",
      grepl("corr", method) ~ "Corrected"
    ),
    estimator = factor(estimator, levels = c("Oracle", "Naive BLUP", "Corrected"))
  )

p_bias <- ggplot(summary_df, aes(x = mean_n_trial, y = bias, color = estimator)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line() +
  geom_point() +
  facet_grid(method_type ~ n_id, scales = "free_y") +
  theme_bw() +
  labs(title = "Bias in secondary regression", x = "Mean trials per subject", y = "Bias")

p_cov <- ggplot(summary_df, aes(x = mean_n_trial, y = coverage, color = estimator)) +
  geom_hline(yintercept = 0.95, linetype = 2) +
  geom_line() +
  geom_point() +
  facet_grid(method_type ~ n_id) +
  theme_bw() +
  labs(title = "Coverage of 95% CI", x = "Mean trials per subject", y = "Coverage")

p_se <- ggplot(summary_df, aes(x = mean_n_trial, y = mean_estimate, color = estimator)) +
  geom_line() +
  geom_point() +
  facet_grid(method_type ~ n_id, scales = "free_y") +
  theme_bw() +
  labs(title = "Mean Estimate", x = "Mean trials per subject", y = "Estimate")

write.csv(results, file = file.path(out_dir, "mlm_random_slope_blup_predictor_replication_results.csv"), row.names = FALSE)
write.csv(summary_df, file = file.path(out_dir, "mlm_random_slope_blup_predictor_summary.csv"), row.names = FALSE)
ggsave(file.path(out_dir, "mlm_random_slope_blup_predictor_bias.png"), p_bias, width = 10, height = 6.5, units = "in", dpi = 300)
ggsave(file.path(out_dir, "mlm_random_slope_blup_predictor_coverage.png"), p_cov, width = 10, height = 6.5, units = "in", dpi = 300)
ggsave(file.path(out_dir, "mlm_random_slope_blup_predictor_mean_se.png"), p_se, width = 10, height = 6.5, units = "in", dpi = 300)

message("Saved outputs to: ", normalizePath(out_dir))
