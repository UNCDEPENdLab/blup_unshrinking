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

locate_repo_root <- function() {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  script_path <- if (length(script_arg) > 0L) {
    normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = FALSE)
  } else {
    NA_character_
  }

  candidates <- unique(normalizePath(c(
    getwd(),
    if (!is.na(script_path)) dirname(script_path) else character()
  ), mustWork = FALSE))

  roots <- candidates[file.exists(file.path(candidates, "R", "source_helpers.R"))]
  if (length(roots) == 0L) {
    stop("Could not locate repository root containing shared R helpers.")
  }
  roots[[1]]
}

repo_root <- locate_repo_root()
source(file.path(repo_root, "R", "source_helpers.R"), local = TRUE)
source_project_helpers(repo_root)

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
OpenMx::mxOption(NULL, "Number of Threads", 1L)

# Argument handling
args <- commandArgs(trailingOnly = TRUE)
n_sim <- if (length(args) >= 1) as.integer(args[[1]]) else 100L
out_dir <- if (length(args) >= 2) args[[2]] else file.path(repo_root, "outputs", "blup_predictor_comparison")
n_cores <- if (length(args) >= 3) as.integer(args[[3]]) else 1L
grid_mode <- if (length(args) >= 4) args[[4]] else "base"
max_conditions <- if (length(args) >= 5) as.integer(args[[5]]) else NA_integer_

if (is.na(n_sim) || n_sim < 1L) {
  stop("`n_sim` must be a positive integer.")
}

if (is.na(n_cores) || n_cores < 1L) {
  stop("`n_cores` must be a positive integer.")
}

if (!is.na(max_conditions) && max_conditions < 1L) {
  stop("`max_conditions` must be a positive integer when supplied.")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

predictor_residual_structures <- function() {
  tibble::tribble(
    ~r_structure, ~r_rho,
    "iid", NA_real_,
    "ar1", 0.3,
    "ar1", 0.6
  )
}

predictor_condition_to_r_spec <- function(condition) {
  r_structure <- if ("r_structure" %in% names(condition)) {
    as.character(condition$r_structure[[1]])
  } else {
    "iid"
  }

  switch(
    r_structure,
    iid = list(structure = "iid"),
    ar1 = list(structure = "ar1", rho = as.numeric(condition$r_rho[[1]])),
    stop("Unsupported residual structure: ", r_structure)
  )
}

predictor_condition_uses_non_iid_R <- function(condition) {
  !identical(predictor_condition_to_r_spec(condition)$structure, "iid")
}

predictor_condition_to_nlme_correlation <- function(condition) {
  r_spec <- predictor_condition_to_r_spec(condition)
  switch(
    r_spec$structure,
    iid = NULL,
    ar1 = {
      if (!requireNamespace("nlme", quietly = TRUE)) {
        stop("The `nlme` package is required for AR(1) Stage-1 residual covariance fits.")
      }
      nlme::corAR1(form = ~trial_index | id)
    },
    stop("Unsupported residual structure for nlme Stage-1 fit: ", r_spec$structure)
  )
}

make_predictor_sim_grid <- function(grid_mode = "base", max_conditions = NA_integer_) {
  grid_mode <- as.character(grid_mode[[1]])

  base_grid <- switch(
    grid_mode,
    smoke = tidyr::crossing(
      n_id = 40L,
      mean_n_trial = 8L,
      tibble::tibble(tau1 = 0.7, sigma = 1.0)
    ),
    residual_ar1 = tidyr::crossing(
      n_id = c(50L, 200L),
      mean_n_trial = c(8L, 20L, 50L),
      tibble::tibble(
        tau1 = c(0.7, 0.3, 1.2),
        sigma = c(1.0, 1.5, 0.5)
      )
    ) %>%
      tidyr::crossing(predictor_residual_structures()),
    base = tidyr::crossing(
      n_id = c(50L, 100L, 200L),
      mean_n_trial = c(8L, 20L, 50L, 100L),
      tibble::tibble(
        tau1 = c(0.7, 0.3, 1.2),
        sigma = c(1.0, 1.5, 0.5)
      )
    ),
    stop("`grid_mode` must be one of: smoke, residual_ar1, base.")
  )

  out <- base_grid %>%
    dplyr::mutate(
      r_structure = if ("r_structure" %in% names(.)) .data$r_structure else "iid",
      r_rho = if ("r_rho" %in% names(.)) .data$r_rho else NA_real_,
      condition_id = dplyr::row_number(),
      design_source = grid_mode
    ) %>%
    dplyr::relocate(condition_id)

  if (!is.na(max_conditions)) {
    out <- out %>% dplyr::slice_head(n = min(max_conditions, nrow(out)))
  }

  out
}

fit_predictor_stage1 <- function(condition, data) {
  if (predictor_condition_uses_non_iid_R(condition)) {
    if (!requireNamespace("nlme", quietly = TRUE)) {
      return(NULL)
    }
    return(safe_lme(
      fixed = y ~ z,
      random = ~1 + z | id,
      data = data,
      correlation = predictor_condition_to_nlme_correlation(condition),
      method = "ML",
      control = nlme::lmeControl(returnObject = TRUE, msMaxIter = 100L, opt = "optim")
    ))
  }

  safe_lmer(y ~ 1 + z + (1 + z | id), data = data, REML = FALSE)
}

empty_predictor_results <- function() {
  tibble::tibble(
    method = c(
      "single_oracle", "single_blup", "single_matrix_corr", "single_diag_corr", "single_closed_form",
      "dual_oracle_u0", "dual_oracle_u1",
      "dual_blup_u0", "dual_blup_u1",
      "dual_matrix_corr_u0", "dual_matrix_corr_u1",
      "dual_diag_corr_u0", "dual_diag_corr_u1",
      "dual_closed_form_u0", "dual_closed_form_u1"
    ),
    estimate = NA_real_,
    se = NA_real_,
    ci_low = NA_real_,
    ci_high = NA_real_
  )
}

fit_predictor_lm_stats <- function(formula, data, term) {
  fit <- tryCatch(stats::lm(formula, data = data), error = function(e) NULL)
  if (is.null(fit)) {
    return(tibble::tibble(estimate = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
  }
  extract_lm_stats(fit, term)
}

sim_grid <- make_predictor_sim_grid(grid_mode, max_conditions = max_conditions)
write.csv(sim_grid, file = file.path(out_dir, "mlm_random_slope_blup_predictor_manifest.csv"), row.names = FALSE)

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
#' @param condition One-row simulation condition.
#' @param params List of fixed base parameters.
#' @return A tibble with oracle, naive BLUP, matrix-corrected, diagonal-
#' corrected, and closed-form GLS predictor results.
run_one_rep <- function(condition, params) {
  r_spec <- predictor_condition_to_r_spec(condition)
  sim <- simulate_dataset(
    n_id = condition$n_id[[1]],
    mean_n_trial = condition$mean_n_trial[[1]],
    params = params,
    tau1 = condition$tau1[[1]],
    sigma = condition$sigma[[1]],
    has_random_slope = TRUE,
    balanced = FALSE,
    r_spec = r_spec
  )

  # Data generation for this specific simulation includes secondary outcomes
  var_u0 <- params$tau0^2
  var_u1_total <- params$gamma_x_on_slope^2 + condition$tau1[[1]]^2
  cov_u0_u1 <- params$rho * params$tau0 * condition$tau1[[1]]

  resid_var_single <- max(0.1, 1 - params$beta_single_u1^2 * var_u1_total)
  resid_var_dual <- max(0.1, 1 - (params$beta_dual_u0^2 * var_u0 + params$beta_dual_u1^2 * var_u1_total + 2 * params$beta_dual_u0 * params$beta_dual_u1 * cov_u0_u1))

  id_df <- sim$id_df %>%
    mutate(
      outcome_single = params$alpha_single + params$beta_single_u1 * true_slope_dev + stats::rnorm(n(), mean = 0, sd = sqrt(resid_var_single)),
      outcome_dual = params$alpha_dual + params$beta_dual_u0 * true_intercept_dev + params$beta_dual_u1 * true_slope_dev + stats::rnorm(n(), mean = 0, sd = sqrt(resid_var_dual))
    )

  fit_null <- fit_predictor_stage1(condition, sim$dat)
  if (is.null(fit_null)) {
    return(empty_predictor_results())
  }

  stage1_components <- tryCatch(
      get_stage1_eb_components(
        fit_obj = fit_null,
        data = sim$dat,
        cluster_var = "id",
        outcome_var = "y",
        within_var = "z"
      ),
    error = function(e) tibble::tibble(id = as.character(id_df$id))
  )
  closed_form_scores <- tryCatch(
    get_closed_form_corrected_scores(
      fit_obj = fit_null,
      data = sim$dat,
      cluster_var = "id",
      outcome_var = "y",
      within_var = "z"
    ),
    error = function(e) tibble::tibble(id = as.character(id_df$id))
  )
  for (col in c("corrected_intercept_full", "corrected_slope_full")) {
    if (!(col %in% names(closed_form_scores))) {
      closed_form_scores[[col]] <- NA_real_
    }
  }

  score_df <- id_df %>%
    left_join(stage1_components, by = "id") %>%
    left_join(
      closed_form_scores %>%
        dplyr::select(id, closed_form_intercept = corrected_intercept_full, closed_form_slope = corrected_slope_full),
      by = "id"
    )

  # Regression 1: Single predictor (slope)
  bind_rows(
    single_oracle = fit_predictor_lm_stats(outcome_single ~ true_slope_dev, score_df, "true_slope_dev"),
    single_blup = fit_predictor_lm_stats(outcome_single ~ blup_z, score_df, "blup_z"),
    single_matrix_corr = fit_predictor_lm_stats(outcome_single ~ corrected_z, score_df, "corrected_z"),
    single_diag_corr = fit_predictor_lm_stats(outcome_single ~ corrected_z_diag, score_df, "corrected_z_diag"),
    single_closed_form = fit_predictor_lm_stats(outcome_single ~ closed_form_slope, score_df, "closed_form_slope"),
    # Regression 2: Dual predictors (intercept and slope)
    dual_oracle_u0 = fit_predictor_lm_stats(outcome_dual ~ true_intercept_dev + true_slope_dev, score_df, "true_intercept_dev"),
    dual_oracle_u1 = fit_predictor_lm_stats(outcome_dual ~ true_intercept_dev + true_slope_dev, score_df, "true_slope_dev"),
    dual_blup_u0 = fit_predictor_lm_stats(outcome_dual ~ blup_intercept + blup_z, score_df, "blup_intercept"),
    dual_blup_u1 = fit_predictor_lm_stats(outcome_dual ~ blup_intercept + blup_z, score_df, "blup_z"),
    dual_matrix_corr_u0 = fit_predictor_lm_stats(outcome_dual ~ corrected_intercept + corrected_z, score_df, "corrected_intercept"),
    dual_matrix_corr_u1 = fit_predictor_lm_stats(outcome_dual ~ corrected_intercept + corrected_z, score_df, "corrected_z"),
    dual_diag_corr_u0 = fit_predictor_lm_stats(outcome_dual ~ corrected_intercept_diag + corrected_z_diag, score_df, "corrected_intercept_diag"),
    dual_diag_corr_u1 = fit_predictor_lm_stats(outcome_dual ~ corrected_intercept_diag + corrected_z_diag, score_df, "corrected_z_diag"),
    dual_closed_form_u0 = fit_predictor_lm_stats(outcome_dual ~ closed_form_intercept + closed_form_slope, score_df, "closed_form_intercept"),
    dual_closed_form_u1 = fit_predictor_lm_stats(outcome_dual ~ closed_form_intercept + closed_form_slope, score_df, "closed_form_slope"),
    .id = "method"
  )
}

# Execution
if (n_cores > 1L) {
  registerDoParallel(cores = n_cores)
}

set.seed(20260411)
results <- foreach(i = seq_len(nrow(sim_grid)), .combine = bind_rows, .packages = c("lme4", "dplyr", "tidyr", "purrr", "tibble")) %dopar% {
  condition_i <- sim_grid[i, , drop = FALSE]

  map_dfr(seq_len(n_sim), function(rep_id) {
    out <- run_one_rep(condition_i, params)
    dplyr::bind_cols(
      out,
      condition_i[rep(1L, nrow(out)), , drop = FALSE] %>% dplyr::select(-condition_id),
      tibble::tibble(rep = rep_id, condition_id = condition_i$condition_id[[1]])
    )
  })
}

if (n_cores > 1L) {
  stopImplicitCluster()
}

# Analysis
results <- results %>%
  mutate(
    truth = case_when(
      grepl("^single", method) ~ params$beta_single_u1,
      grepl("u0$", method) ~ params$beta_dual_u0,
      grepl("u1$", method) ~ params$beta_dual_u1,
      TRUE ~ NA_real_
    ),
    bias = estimate - truth,
    covered = ci_low <= truth & ci_high >= truth,
    sq_error = (estimate - truth)^2
  )

summary_df <- results %>%
  group_by(method, n_id, mean_n_trial, tau1, sigma, r_structure, r_rho, design_source) %>%
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
      grepl("matrix_corr", method) ~ "Matrix corrected",
      grepl("diag_corr", method) ~ "Diagonal corrected",
      grepl("closed_form", method) ~ "Closed-form GLS"
    ),
    estimator = factor(estimator, levels = c("Oracle", "Naive BLUP", "Matrix corrected", "Diagonal corrected", "Closed-form GLS"))
  )

p_bias <- ggplot(summary_df, aes(x = mean_n_trial, y = bias, color = estimator)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_line() +
  geom_point() +
  facet_grid(method_type + r_structure + r_rho ~ n_id, scales = "free_y") +
  theme_bw() +
  labs(title = "Bias in secondary regression", x = "Mean trials per subject", y = "Bias")

p_cov <- ggplot(summary_df, aes(x = mean_n_trial, y = coverage, color = estimator)) +
  geom_hline(yintercept = 0.95, linetype = 2) +
  geom_line() +
  geom_point() +
  facet_grid(method_type + r_structure + r_rho ~ n_id) +
  theme_bw() +
  labs(title = "Coverage of 95% CI", x = "Mean trials per subject", y = "Coverage")

p_se <- ggplot(summary_df, aes(x = mean_n_trial, y = mean_estimate, color = estimator)) +
  geom_line() +
  geom_point() +
  facet_grid(method_type + r_structure + r_rho ~ n_id, scales = "free_y") +
  theme_bw() +
  labs(title = "Mean Estimate", x = "Mean trials per subject", y = "Estimate")

write.csv(results, file = file.path(out_dir, "mlm_random_slope_blup_predictor_replication_results.csv"), row.names = FALSE)
write.csv(summary_df, file = file.path(out_dir, "mlm_random_slope_blup_predictor_summary.csv"), row.names = FALSE)
ggsave(file.path(out_dir, "mlm_random_slope_blup_predictor_bias.png"), p_bias, width = 10, height = 6.5, units = "in", dpi = 300)
ggsave(file.path(out_dir, "mlm_random_slope_blup_predictor_coverage.png"), p_cov, width = 10, height = 6.5, units = "in", dpi = 300)
ggsave(file.path(out_dir, "mlm_random_slope_blup_predictor_mean_se.png"), p_se, width = 10, height = 6.5, units = "in", dpi = 300)

message("Saved outputs to: ", normalizePath(out_dir))
