#!/usr/bin/env Rscript

#' ---
#' title: "Basic MLM BLUP Correction Simulation"
#' description: "Archived scalar demonstration of BLUP shrinkage and prior unweighting
#'               in a random-intercept-only multilevel model."
#' ---
#'
#' # Purpose
#'
#' This script is a compact, pedagogical simulation for the one-dimensional
#' version of the BLUP correction problem. It generates a random-intercept-only
#' multilevel model where a level-2 predictor `x` affects the true latent
#' subject intercept. It then compares four estimators for that population
#' association:
#'
#' 1. an oracle regression using the true simulated latent intercept;
#' 2. a naive regression using the empirical-Bayes/BLUP intercept;
#' 3. a corrected-score regression that removes the Gaussian prior contribution
#'    from the BLUP/posterior mean using the scalar version of the matrix
#'    unweighting algebra; and
#' 4. a direct mixed model with `x` included in the level-1 model.
#'
#' The intended lesson is narrow: ordinary BLUPs are shrunken posterior means,
#' so using them as stage-2 outcomes attenuates the association with `x`; the
#' likelihood-only corrected score recovers the target slope in this simple
#' Gaussian setup.
#'
#' # Relationship to the current pipeline
#'
#' This file predates the more robust `blup_outcome/` pipeline and is now kept
#' only as an archived conceptual reference. It does not implement the current
#' production simulation machinery: no random-slope design grid, no diagonal
#' versus full-matrix random-slope comparison, no closed-form score benchmark,
#' no Lai 2S-PA/2S-PAA comparison, no stacked-sandwich variants, no chunked
#' execution, and no first-stage diagnostics. For substantive simulation work,
#' use `blup_outcome/mlm_random_slope_blup_outcome_sim.R`.
#'
#' The remaining virtue of this script is that the scalar random-intercept case
#' is easier to inspect than the full bivariate random-intercept/random-slope
#' pipeline. It is useful when explaining why BLUP-as-outcome analyses need
#' de-shrinking before moving to the full random-slope machinery.

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
out_dir <- if (length(args) >= 2) args[[2]] else file.path(repo_root, "outputs", "blup_mlm")

if (is.na(n_sim) || n_sim < 1L) {
  stop("`n_sim` must be a positive integer.")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Grid setup
sim_grid <- tidyr::crossing(
  n_id = c(50L, 200L),
  n_trial = c(5L, 10L, 25L, 100L)
)

params <- list(
  beta_0 = 1,
  beta_1 = 0.5,
  tau = 1,
  sigma = 1
)

#' Execute a single replication of the basic random intercept simulation.
#' @param n_id Number of subjects.
#' @param n_trial Number of trials per subject.
#' @param params List of fixed base parameters.
#' @return A tibble with extraction results for Oracle, Naive BLUP, Corrected BLUP, and Direct models.
run_one_rep <- function(n_id, n_trial, params) {
  sim <- simulate_dataset(
    n_id = n_id,
    mean_n_trial = n_trial,
    params = params,
    has_random_slope = FALSE,
    balanced = TRUE
  )

  fit_null <- suppressWarnings(suppressMessages(
    lmer(y ~ 1 + (1 | id), data = sim$dat, REML = TRUE)
  ))

  fit_direct <- suppressWarnings(suppressMessages(
    lmer(y ~ x + (1 | id), data = sim$dat, REML = TRUE)
  ))

  score_df <- sim$id_df %>%
    left_join(
      get_stage1_eb_components(
        fit_obj = fit_null,
        data = sim$dat,
        cluster_var = "id",
        outcome_var = "y",
        within_var = NULL
      ),
      by = "id"
    )

  oracle_fit <- lm(eta ~ x, data = score_df)
  naive_fit <- lm(blup_intercept ~ x, data = score_df)
  corrected_fit <- lm(corrected_intercept ~ x, data = filter(score_df, !is.na(corrected_intercept)))

  bind_rows(
    oracle = extract_lm_stats(oracle_fit),
    naive_blup = extract_lm_stats(naive_fit),
    corrected_blup = extract_lm_stats(corrected_fit),
    direct_mlm = extract_lmer_stats(fit_direct, term = "x"),
    .id = "method"
  )
}

set.seed(20260410)

results <- purrr::map_dfr(seq_len(nrow(sim_grid)), function(i) {
  n_id_i <- sim_grid$n_id[[i]]
  n_trial_i <- sim_grid$n_trial[[i]]

  purrr::map_dfr(seq_len(n_sim), function(rep_id) {
    out <- run_one_rep(
      n_id = n_id_i,
      n_trial = n_trial_i,
      params = params
    )
    mutate(out, rep = rep_id, n_id = n_id_i, n_trial = n_trial_i)
  })
})

# Analysis & Plotting (mostly unchanged logic)
results <- results %>%
  mutate(
    truth = params$beta_1,
    bias = estimate - truth,
    sq_error = (estimate - truth)^2,
    covered = ci_low <= truth & ci_high >= truth
  )

summary_df <- results %>%
  group_by(method, n_id, n_trial) %>%
  summarise(
    mean_estimate = mean(estimate, na.rm = TRUE),
    mc_se_mean = sd(estimate, na.rm = TRUE) / sqrt(sum(!is.na(estimate))),
    bias = mean(bias, na.rm = TRUE),
    rmse = sqrt(mean(sq_error, na.rm = TRUE)),
    coverage = mean(covered, na.rm = TRUE),
    n_success = sum(!is.na(estimate)),
    .groups = "drop"
  ) %>%
  mutate(
    method = factor(
      method,
      levels = c("oracle", "naive_blup", "corrected_blup", "direct_mlm"),
      labels = c("Oracle latent score", "Naive BLUP regression", "Corrected BLUP regression", "Direct mixed model")
    )
  )

method_palette <- c(
  "Oracle latent score" = "#111111",
  "Naive BLUP regression" = "#b33c2e",
  "Corrected BLUP regression" = "#1f6f8b",
  "Direct mixed model" = "#2a9d55"
)

p_mean <- ggplot(summary_df, aes(x = n_trial, y = mean_estimate, color = method)) +
  geom_hline(yintercept = params$beta_1, linetype = 2, color = "grey40") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(ymin = mean_estimate - 1.96 * mc_se_mean, ymax = mean_estimate + 1.96 * mc_se_mean),
    width = 0.12,
    linewidth = 0.5
  ) +
  facet_wrap(~n_id, labeller = label_both) +
  scale_x_continuous(breaks = sort(unique(summary_df$n_trial))) +
  scale_color_manual(values = method_palette) +
  labs(
    title = "Recovering the trait effect from BLUPs in a random-intercept MLM",
    subtitle = "Dashed line is the true population slope; points are Monte Carlo means",
    x = "Trials per subject",
    y = "Estimated slope",
    color = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

p_bias <- ggplot(summary_df, aes(x = n_trial, y = bias, color = method)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~n_id, labeller = label_both) +
  scale_x_continuous(breaks = sort(unique(summary_df$n_trial))) +
  scale_color_manual(values = method_palette) +
  labs(
    title = "Bias relative to the true population slope",
    x = "Trials per subject",
    y = "Bias",
    color = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

write.csv(results, file = file.path(out_dir, "mlm_blup_correction_replication_results.csv"), row.names = FALSE)
write.csv(summary_df, file = file.path(out_dir, "mlm_blup_correction_summary.csv"), row.names = FALSE)
ggsave(filename = file.path(out_dir, "mlm_blup_correction_mean_estimates.png"), plot = p_mean, width = 9, height = 5.5, units = "in", dpi = 300)
ggsave(filename = file.path(out_dir, "mlm_blup_correction_bias.png"), plot = p_bias, width = 9, height = 5.5, units = "in", dpi = 300)

message("Saved outputs to: ", normalizePath(out_dir))
print(summary_df)
