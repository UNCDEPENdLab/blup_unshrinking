#!/usr/bin/env Rscript

#' ---
#' title: "Random Slope BLUP Correction Simulation"
#' description: "Archived bivariate demonstration of why random-slope BLUP correction
#'               must use the full intercept-slope posterior covariance matrix."
#' ---
#'
#' # Purpose
#'
#' This script is a compact conceptual simulation for the two-dimensional BLUP
#' correction problem. It generates a random-intercept/random-slope multilevel
#' model in which a level-2 predictor `x` shifts the true subject-specific
#' random slope for the level-1 predictor `z`. It then compares five estimators
#' for the population `x -> random slope` association:
#'
#' 1. an oracle regression using the true simulated random slopes;
#' 2. a naive regression using ordinary random-slope BLUPs;
#' 3. a diagonal-only corrected-score regression that unweights each random
#'    effect separately;
#' 4. a full-matrix corrected-score regression that uses the joint
#'    intercept-slope posterior covariance; and
#' 5. a direct mixed model estimating the same `x:z` interaction in one model.
#'
#' The intended lesson is that the scalar correction does not generalize by
#' correcting each random-effect component independently. When random intercepts
#' and random slopes are correlated, the posterior covariance couples the two
#' empirical-Bayes scores. Ignoring that covariance creates a diagnostic
#' diagonal-only correction that can remain biased even when the full matrix
#' prior-unweighting recovers the target association.
#'
#' # Relationship to the current pipeline
#'
#' This file predates the more robust `blup_outcome/` pipeline and is now kept
#' only as an archived conceptual reference. The active pipeline includes this
#' script's core comparison, but with a broader design grid and additional
#' estimators: diagonal versus full-matrix corrected outcomes, closed-form
#' within-cluster scores, single-subject OLS slopes, Lai 2S-PA/2S-PAA,
#' stacked-sandwich variants, heteroscedastic/informative imbalance conditions,
#' chunked execution, resumable condition files, and first-stage diagnostics.
#' For substantive simulation work, use
#' `blup_outcome/mlm_random_slope_blup_outcome_sim.R`.
#'
#' The remaining virtue of this script is clarity. It isolates the specific
#' bivariate algebraic point that random-slope score correction is a matrix
#' operation, not a set of unrelated scalar corrections.

suppressPackageStartupMessages({
  library(lme4)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
})

# Load helpers relative to the repository root. This script is archived, so it
# cannot assume it is run from the historical `mlm_blups` directory.
locate_repo_root <- function() {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  script_path <- if (length(script_arg) > 0L) {
    normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = FALSE)
  } else {
    NA_character_
  }

  candidates <- unique(normalizePath(c(
    getwd(),
    if (!is.na(script_path)) dirname(script_path) else character(),
    if (!is.na(script_path)) file.path(dirname(script_path), "..") else character()
  ), mustWork = FALSE))

  roots <- candidates[file.exists(file.path(candidates, "R", "sim_helpers.R"))]
  if (length(roots) == 0L) {
    stop("Could not locate repository root containing R/sim_helpers.R.")
  }

  roots[[1]]
}

repo_root <- locate_repo_root()
source(file.path(repo_root, "R", "sim_helpers.R"), local = TRUE)
source(file.path(repo_root, "R", "stats_helpers.R"), local = TRUE)
source(file.path(repo_root, "R", "blup_helpers.R"), local = TRUE)

# Argument handling
args <- commandArgs(trailingOnly = TRUE)
n_sim <- if (length(args) >= 1) as.integer(args[[1]]) else 250L
out_dir <- if (length(args) >= 2) args[[2]] else file.path(repo_root, "outputs", "blup_random_slope")

if (is.na(n_sim) || n_sim < 1L) {
  stop("`n_sim` must be a positive integer.")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Grid setup
sim_grid <- tidyr::crossing(
  n_id = c(50L, 200L),
  mean_n_trial = c(8L, 20L, 50L),
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
  rho = 0.5
)

#' Execute a single replication of the random slope BLUP correction simulation.
#' @param n_id Number of subjects.
#' @param mean_n_trial Target mean trials per subject.
#' @param tau1 Standard deviation of the random slope.
#' @param sigma Residual standard deviation.
#' @param params List of fixed base parameters.
#' @return A tibble with extraction results for Oracle, Naive BLUP, Diagonal-only, Full-matrix, and Direct models.
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

  fit_null <- suppressWarnings(suppressMessages(
    lmer(y ~ 1 + z + (1 + z | id), data = sim$dat, REML = FALSE)
  ))

  fit_direct <- suppressWarnings(suppressMessages(
    lmer(y ~ x + z + x:z + (1 + z | id), data = sim$dat, REML = FALSE)
  ))

  score_df <- sim$id_df %>%
    left_join(
      get_stage1_eb_components(
        fit_obj = fit_null,
        data = sim$dat,
        cluster_var = "id",
        outcome_var = "y",
        within_var = "z"
      ),
      by = "id"
    )

  oracle_fit <- lm(true_slope_dev ~ x, data = score_df)
  naive_fit <- lm(blup_z ~ x, data = score_df)
  diag_fit <- lm(corrected_z_diag ~ x, data = filter(score_df, !is.na(corrected_z_diag)))
  full_fit <- lm(corrected_z ~ x, data = filter(score_df, !is.na(corrected_z)))

  bind_rows(
    oracle = extract_lm_stats(oracle_fit),
    naive_blup = extract_lm_stats(naive_fit),
    diag_only_correction = extract_lm_stats(diag_fit),
    corrected_full_matrix = extract_lm_stats(full_fit),
    direct_mlm = extract_lmer_stats(fit_direct, term = "x:z"),
    .id = "method"
  ) %>%
    mutate(mean_realized_trials = sim$mean_realized_trials)
}

set.seed(20260410)

results <- purrr::map_dfr(seq_len(nrow(sim_grid)), function(i) {
  n_id_i <- sim_grid$n_id[[i]]
  mean_n_trial_i <- sim_grid$mean_n_trial[[i]]
  tau1_i <- sim_grid$tau1[[i]]
  sigma_i <- sim_grid$sigma[[i]]

  purrr::map_dfr(seq_len(n_sim), function(rep_id) {
    out <- run_one_rep(
      n_id = n_id_i,
      mean_n_trial = mean_n_trial_i,
      tau1 = tau1_i,
      sigma = sigma_i,
      params = params
    )
    mutate(out, rep = rep_id, n_id = n_id_i, mean_n_trial = mean_n_trial_i, tau1 = tau1_i, sigma = sigma_i)
  })
})

# Analysis & Plotting
results <- results %>%
  mutate(
    truth = params$gamma_x_on_slope,
    bias = estimate - truth,
    sq_error = (estimate - truth)^2,
    covered = ci_low <= truth & ci_high >= truth
  )

summary_df <- results %>%
  group_by(method, n_id, mean_n_trial, tau1, sigma) %>%
  summarise(
    mean_estimate = mean(estimate, na.rm = TRUE),
    mc_se_mean = sd(estimate, na.rm = TRUE) / sqrt(sum(!is.na(estimate))),
    bias = mean(bias, na.rm = TRUE),
    rmse = sqrt(mean(sq_error, na.rm = TRUE)),
    coverage = mean(covered, na.rm = TRUE),
    realized_trials = mean(mean_realized_trials, na.rm = TRUE),
    n_success = sum(!is.na(estimate)),
    .groups = "drop"
  ) %>%
  mutate(
    method = factor(
      method,
      levels = c("oracle", "naive_blup", "diag_only_correction", "corrected_full_matrix", "direct_mlm"),
      labels = c("Oracle latent slope", "Naive slope BLUP regression", "Diagonal-only correction", "Full matrix correction", "Direct mixed model")
    )
  )

method_palette <- c(
  "Oracle latent slope" = "#111111",
  "Naive slope BLUP regression" = "#b33c2e",
  "Diagonal-only correction" = "#a66a00",
  "Full matrix correction" = "#1f6f8b",
  "Direct mixed model" = "#2a9d55"
)

p_mean <- ggplot(summary_df, aes(x = mean_n_trial, y = mean_estimate, color = method)) +
  geom_hline(yintercept = params$gamma_x_on_slope, linetype = 2, color = "grey40") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_estimate - 1.96 * mc_se_mean, ymax = mean_estimate + 1.96 * mc_se_mean), width = 0.2, linewidth = 0.5) +
  facet_grid(tau1 + sigma ~ n_id, labeller = label_both) +
  scale_x_continuous(breaks = sort(unique(summary_df$mean_n_trial))) +
  scale_color_manual(values = method_palette) +
  labs(
    title = "Recovering the external association from random-slope BLUPs",
    subtitle = "Full matrix correction uses the joint intercept-slope posterior covariance",
    x = "Target mean trials per subject", y = "Estimated slope of x on the subject-specific random slope", color = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

p_bias <- ggplot(summary_df, aes(x = mean_n_trial, y = bias, color = method)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_grid(tau1 + sigma ~ n_id, labeller = label_both) +
  scale_x_continuous(breaks = sort(unique(summary_df$mean_n_trial))) +
  scale_color_manual(values = method_palette) +
  labs(title = "Bias in estimating the x to random-slope association", x = "Target mean trials per subject", y = "Bias", color = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

write.csv(results, file = file.path(out_dir, "mlm_random_slope_blup_replication_results.csv"), row.names = FALSE)
write.csv(summary_df, file = file.path(out_dir, "mlm_random_slope_blup_summary.csv"), row.names = FALSE)
ggsave(filename = file.path(out_dir, "mlm_random_slope_blup_mean_estimates.png"), plot = p_mean, width = 9, height = 5.5, units = "in", dpi = 300)
ggsave(filename = file.path(out_dir, "mlm_random_slope_blup_bias.png"), plot = p_bias, width = 9, height = 5.5, units = "in", dpi = 300)

message("Saved outputs to: ", normalizePath(out_dir))
print(summary_df)
