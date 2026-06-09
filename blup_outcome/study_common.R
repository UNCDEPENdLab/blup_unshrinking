#' Shared BLUP-as-outcome replication machinery.
#'
#' This module contains the per-replication logic for simulations where
#' recovered random-slope scores are used as second-stage outcomes. It bridges
#' the generic simulation helpers, BLUP/corrected-score extraction helpers,
#' Lai/OpenMx wrappers, and stacked-sandwich machinery into one standardized
#' result table per Monte Carlo replication.

#' Return the BLUP-outcome method codes for an analysis mode.
#'
#' @details
#' The BLUP-outcome runner has two estimator sets. `screen` mode returns the
#' core, relatively cheap estimators used for fast design checks. `full` mode
#' appends the computationally heavier Lai 2S-PA/PAA rows and stacked-sandwich
#' corrected-score rows. The returned method vector is also used for failure
#' templates, so every replication has a stable row shape even when Stage 1
#' fitting fails.
#'
#' @param analysis_mode Character scalar. `"full"` includes Lai and
#'   stacked-sandwich estimators; any other value returns the core method set.
#'
#' @return Character vector of method codes in the expected reporting order.
blup_outcome_methods <- function(analysis_mode = "full") {
  base_methods <- c(
    "oracle",
    "naive_blup", "naive_blup_hc3",
    "diag_corrected", "diag_corrected_hc3",
    "matrix_corrected", "matrix_corrected_hc3",
    "closed_form", "closed_form_hc3",
    "single_subject_ols", "single_subject_ols_hc3",
    "direct_mlm"
  )

  if (identical(analysis_mode, "full")) {
    c(
      base_methods,
      "lai_2spa", "lai_2spaa",
      "fuller_blup", "fuller_diag_corrected", "fuller_matrix_corrected", "fuller_closed_form",
      "fuller_alpha_blup", "fuller_alpha_diag_corrected", "fuller_alpha_matrix_corrected", "fuller_alpha_closed_form",
      paste0("closed_form_stacked_hc", 0:3)
    )
  } else {
    base_methods
  }
}

#' Construct an all-missing BLUP-outcome result table.
#'
#' @details
#' Simulation loops should record failed replications as rows rather than
#' throwing away the replication. This helper creates one row per requested
#' method with missing estimates and intervals while preserving the simulation
#' truth.
#'
#' @param methods Character vector of method codes to include.
#' @param truth Numeric scalar true second-stage coefficient.
#'
#' @return Tibble with standardized estimator columns.
empty_blup_outcome_result <- function(methods, truth) {
  tibble::tibble(
    method = methods,
    estimate = NA_real_,
    se = NA_real_,
    ci_low = NA_real_,
    ci_high = NA_real_,
    status_code = NA_integer_,
    truth = truth
  )
}

#' Construct an all-missing BLUP-outcome diagnostics row.
#'
#' @details
#' Diagnostics are bound onto every estimator row. This helper defines the full
#' diagnostic schema for cases where Stage 1 fitting fails or diagnostic values
#' cannot be computed, while reusing the shared Stage 1 diagnostic schema.
#'
#' @return One-row tibble with missing diagnostic values.
empty_blup_outcome_diagnostics <- function() {
  dplyr::bind_cols(
    empty_stage1_diagnostics(),
    tibble::tibble(
      mean_postvar_u1 = NA_real_,
      mean_theta_u1 = NA_real_,
      mean_lambda_u1 = NA_real_,
      blup_variance_ratio = NA_real_,
      diag_corrected_variance_ratio = NA_real_,
      matrix_corrected_variance_ratio = NA_real_,
      closed_form_variance_ratio = NA_real_,
      blup_true_cor = NA_real_,
      matrix_corrected_true_cor = NA_real_,
      closed_form_true_cor = NA_real_,
      diag_corrected_failure_rate = NA_real_,
      matrix_corrected_failure_rate = NA_real_,
      closed_form_failure_rate = NA_real_,
      cluster_size_x_cor = NA_real_,
      empirical_g_error = NA_real_,
      empirical_structural_r2 = NA_real_,
      empirical_reliability = NA_real_,
      residual_g_min_eigenvalue = NA_real_
    )
  )
}

#' Attach empty diagnostics and realized-trial context to failed rows.
#'
#' @details
#' When the null mixed model fails, no score diagnostics are available. If a
#' simulated dataset exists, its realized trial-count descriptors are still
#' useful for diagnosing failures in sparse or unbalanced designs, so this
#' helper preserves them while leaving score/model diagnostics missing.
#'
#' @param results Estimator result tibble, usually from
#'   `empty_blup_outcome_result()`.
#' @param sim Optional simulated dataset list returned by `simulate_dataset()`.
#' @param condition Optional one-row design condition used to retain calibrated
#'   DGP diagnostics even when Stage 1 fitting fails.
#'
#' @return `results` with diagnostic and realized-trial columns appended.
add_empty_blup_outcome_context <- function(results, sim = NULL, condition = NULL) {
  if (is.null(sim)) {
    return(results %>%
      dplyr::bind_cols(empty_blup_outcome_diagnostics()[rep(1L, nrow(results)), , drop = FALSE]) %>%
      dplyr::mutate(
        mean_realized_trials = NA_real_,
        min_realized_trials = NA_real_,
        prop_ids_leq_2_trials = NA_real_,
        prop_ids_leq_3_trials = NA_real_
      ))
  }

  diagnostics <- empty_blup_outcome_diagnostics()
  calibrated_dgp <- make_calibrated_dgp_diagnostics(condition = condition, sim = sim)
  diagnostics[names(calibrated_dgp)] <- calibrated_dgp

  results %>%
    dplyr::bind_cols(diagnostics[rep(1L, nrow(results)), , drop = FALSE]) %>%
    dplyr::mutate(
      mean_realized_trials = sim$mean_realized_trials,
      min_realized_trials = sim$min_realized_trials,
      prop_ids_leq_2_trials = sim$prop_ids_leq_2_trials,
      prop_ids_leq_3_trials = sim$prop_ids_leq_3_trials
    )
}

#' Ensure all downstream Stage 2 score columns are present.
#'
#' @details
#' Different score-recovery helpers can fail independently. Rather than branch
#' every downstream estimator on column existence, this function adds missing
#' score, EB-measurement, and OLS-slope columns as `NA_real_`. Estimators then
#' fail through their normal complete-case checks.
#'
#' @param stage2_df Cluster-level Stage 2 data frame assembled from true scores,
#'   EB/Lai inputs, corrected scores, and OLS slopes.
#'
#' @return `stage2_df` with the required columns present.
ensure_blup_outcome_stage2_columns <- function(stage2_df) {
  needed <- c(
    "x", "true_slope_dev", "u0_eb", "u1_eb",
    "postvar11", "postvar12", "postvar22",
    "lambda11", "lambda12", "lambda21", "lambda22",
    "theta11", "theta12", "theta22",
    "corrected_z", "corrected_z_var",
    "corrected_z_diag", "corrected_z_diag_var",
    "corrected_slope_full", "ols_var22", "ols_slope"
  )
  for (col in setdiff(needed, names(stage2_df))) {
    stage2_df[[col]] <- NA_real_
  }
  stage2_df
}

#' Standardize estimator rows to the common result schema.
#'
#' @details
#' Individual estimator helpers return slightly different subsets of columns.
#' This helper row-binds a zero-row template first so the output always has the
#' expected columns and fills missing `truth` values with the replication truth.
#'
#' @param x Estimator result tibble.
#' @param truth Numeric scalar true second-stage coefficient.
#'
#' @return Tibble with standardized estimator columns.
standardize_estimator_rows <- function(x, truth) {
  required <- tibble::tibble(
    method = character(),
    estimate = numeric(),
    se = numeric(),
    ci_low = numeric(),
    ci_high = numeric(),
    status_code = integer(),
    truth = numeric()
  )

  out <- dplyr::bind_rows(required, x)
  if (!("truth" %in% names(x))) {
    out$truth <- truth
  } else {
    out$truth[is.na(out$truth)] <- truth
  }
  out
}

#' Fit OLS regression of a recovered score outcome on the true level-2 predictor.
#'
#' @details
#' The BLUP-outcome simulation evaluates many recovered random-slope outcomes
#' with the same Stage 2 model: `outcome ~ x`. This wrapper delegates fitting to
#' the shared observed-score OLS helper, maps the usual and HC3 standard-error
#' variants to method names, and optionally drops the HC3 row for oracle rows
#' where only one inferential variant is reported.
#'
#' @param stage2_df Cluster-level Stage 2 data frame.
#' @param outcome Character scalar naming the recovered score outcome column.
#' @param base_method Character scalar base method code.
#' @param include_hc3 Logical; if `TRUE`, retain the HC3 robust row.
#'
#' @return Tibble with one or two standardized estimator rows.
fit_score_outcome_ols <- function(stage2_df, outcome, base_method, include_hc3 = TRUE) {
  fit_tbl <- fit_observed_single(stage2_df, outcome = outcome, predictor = "x")
  out <- finalize_ols_se_variants(fit_tbl, base_method)
  if (!isTRUE(include_hc3)) {
    out <- out %>% dplyr::filter(.data$method == base_method)
  }

  out %>%
    dplyr::mutate(
      status_code = dplyr::if_else(is.finite(.data$estimate) & is.finite(.data$se), 0L, NA_integer_)
    )
}

#' Construct missing stacked-sandwich rows.
#'
#' @details
#' Stacked-sandwich estimation can fail because derivative, likelihood, or
#' covariance calculations are numerically fragile. This helper returns the
#' expected HC0-HC3 rows with missing values so full-mode result shape remains
#' stable.
#'
#' @param method_prefix Character scalar prefix for stacked-sandwich method
#'   codes.
#'
#' @return Four-row tibble for HC0-HC3 variants.
empty_stacked_rows <- function(method_prefix = "closed_form_stacked") {
  purrr::map_dfr(paste0("hc", 0:3), function(variant) {
    tibble::tibble(
      method = paste0(method_prefix, "_", variant),
      estimate = NA_real_,
      se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      status_code = NA_integer_
    )
  })
}

#' Extract fixed-effect statistics from a supported direct mixed model.
#'
#' @details
#' `extract_lmer_stats()` handles `lme4` fits, but the non-iid residual
#' conditions use `nlme::lme()` so that the direct benchmark is fit under the
#' same R-side covariance structure as the score extractors. This helper keeps
#' the interaction-term matching identical across both fit classes.
#'
#' @param fit Fitted `lme4` or `nlme` direct model.
#' @param term Fixed-effect term to extract. Interaction terms are matched in
#' either order.
#' @param use_t Logical; use a t critical value if `TRUE`.
#' @param df Degrees of freedom for the confidence interval.
#'
#' @return One-row tibble with `estimate`, `se`, `ci_low`, and `ci_high`.
extract_direct_mlm_stats <- function(fit, term = "x:z", use_t = TRUE, df = NULL) {
  coef_tab <- tryCatch({
    if (inherits(fit, "lme")) {
      summary(fit)$tTable
    } else {
      coef(summary(fit))
    }
  }, error = function(e) NULL)

  if (is.null(coef_tab)) {
    return(tibble::tibble(estimate = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
  }

  term_name <- if (term %in% rownames(coef_tab)) {
    term
  } else {
    matched <- grep(term, rownames(coef_tab), value = TRUE)
    if (length(matched) == 0L && identical(term, "x:z")) {
      matched <- grep("z:x", rownames(coef_tab), value = TRUE)
    }
    if (length(matched) == 1L) matched else NA_character_
  }

  if (is.na(term_name)) {
    return(tibble::tibble(estimate = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
  }

  estimate_col <- if ("Estimate" %in% colnames(coef_tab)) "Estimate" else "Value"
  se_col <- if ("Std. Error" %in% colnames(coef_tab)) "Std. Error" else "Std.Error"
  est <- unname(coef_tab[term_name, estimate_col])
  se <- unname(coef_tab[term_name, se_col])
  crit <- if (isTRUE(use_t)) stats::qt(0.975, df = df) else stats::qnorm(0.975)

  tibble::tibble(
    estimate = est,
    se = se,
    ci_low = est - crit * se,
    ci_high = est + crit * se
  )
}

#' Extract the direct mixed-model benchmark row.
#'
#' @details
#' The direct MLM benchmark estimates the target interaction in one model,
#' `y ~ x + z + x:z + (1 + z | id)`. The reported row uses a t critical value
#' with approximate residual degrees of freedom based on the number of clusters
#' minus the number of fixed effects. In non-iid residual conditions, callers
#' pass an `nlme::lme` fit so this benchmark uses the same R-side covariance
#' structure as the R-aware score extraction path.
#'
#' @param fit_direct Fitted direct mixed model or `NULL` if fitting failed.
#' @param n_id Integer number of subjects/clusters.
#'
#' @return One-row tibble for the `direct_mlm` method.
fit_direct_mlm_row <- function(fit_direct, n_id) {
  if (is.null(fit_direct)) {
    return(tibble::tibble(
      method = "direct_mlm",
      estimate = NA_real_,
      se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      status_code = NA_integer_
    ))
  }

  extract_direct_mlm_stats(fit_direct, term = "x:z", use_t = TRUE, df = n_id - length(stage1_fixef(fit_direct))) %>%
    dplyr::mutate(
      method = "direct_mlm",
      status_code = dplyr::if_else(is.finite(.data$estimate) & is.finite(.data$se), 0L, NA_integer_)
    ) %>%
    dplyr::select("method", dplyr::everything())
}

#' Compute Stage 1 and recovered-score diagnostics for one replication.
#'
#' @details
#' Diagnostics serve two purposes: they explain Stage 1 instability and they
#' characterize how each score recovery method changes the empirical random
#' slope signal. The variance ratios compare recovered slope-score variance
#' with the true latent slope variance; correlations compare recovered scores
#' against the true latent slope; failure rates record non-finite corrected
#' scores; and `cluster_size_x_cor` flags informative imbalance in realized
#' trial counts.
#'
#' @param stage1_score_fit Fitted null mixed model used to extract
#'   BLUP/corrected scores. In non-iid residual conditions this should be the
#'   R-aware `nlme::lme` fit, not the iid `lme4` starting-value fit.
#' @param stage2_df Cluster-level Stage 2 data frame.
#' @param sim Simulated dataset list from `simulate_dataset()`.
#' @param condition Optional one-row design condition. Calibrated conditions
#'   supply the fixed marginal/residual covariance targets needed for empirical
#'   DGP diagnostics; legacy conditions return `NA` for those fields.
#'
#' @return One-row tibble containing shared Stage 1 diagnostics plus
#'   BLUP-outcome-specific score diagnostics.
make_blup_outcome_diagnostics <- function(stage1_score_fit, stage2_df, sim, condition = NULL) {
  stage1_diag <- get_stage1_diagnostics(stage1_score_fit, stage2_df, predictor_u0 = "u0_eb", predictor_u1 = "u1_eb")

  safe_cor <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3L || stats::sd(x[ok]) <= sqrt(.Machine$double.eps) ||
      stats::sd(y[ok]) <= sqrt(.Machine$double.eps)) {
      return(NA_real_)
    }
    stats::cor(x[ok], y[ok])
  }

  variance_ratio <- function(measured, truth) {
    measured <- measured[is.finite(measured)]
    truth <- truth[is.finite(truth)]
    if (length(measured) < 2L || length(truth) < 2L || stats::var(truth) <= sqrt(.Machine$double.eps)) {
      return(NA_real_)
    }
    stats::var(measured) / stats::var(truth)
  }

  failure_rate <- function(x) {
    if (length(x) == 0L) {
      return(NA_real_)
    }
    mean(!is.finite(x))
  }

  # Realized trial counts can vary even within a nominal condition when the
  # design is unbalanced or highly unbalanced; correlate them with x to detect
  # informative imbalance induced by the simulation design.
  trial_counts <- sim$dat %>%
    dplyr::count(.data$id, name = "n_trial") %>%
    dplyr::mutate(id = as.character(.data$id)) %>%
    dplyr::left_join(sim$id_df %>% dplyr::select("id", "x"), by = "id")

  cluster_size_x_cor <- if (nrow(trial_counts) >= 3L &&
    stats::sd(trial_counts$n_trial) > sqrt(.Machine$double.eps) &&
    stats::sd(trial_counts$x) > sqrt(.Machine$double.eps)) {
    stats::cor(trial_counts$n_trial, trial_counts$x)
  } else {
    NA_real_
  }

  calibrated_dgp <- make_calibrated_dgp_diagnostics(condition = condition, sim = sim)

  dplyr::bind_cols(
    stage1_diag,
    tibble::tibble(
      mean_postvar_u1 = safe_mean(stage2_df$postvar22),
      mean_theta_u1 = safe_mean(stage2_df$theta22),
      mean_lambda_u1 = safe_mean(stage2_df$lambda22),
      blup_variance_ratio = variance_ratio(stage2_df$u1_eb, stage2_df$true_slope_dev),
      diag_corrected_variance_ratio = variance_ratio(stage2_df$corrected_z_diag, stage2_df$true_slope_dev),
      matrix_corrected_variance_ratio = variance_ratio(stage2_df$corrected_z, stage2_df$true_slope_dev),
      closed_form_variance_ratio = variance_ratio(stage2_df$corrected_slope_full, stage2_df$true_slope_dev),
      blup_true_cor = safe_cor(stage2_df$u1_eb, stage2_df$true_slope_dev),
      matrix_corrected_true_cor = safe_cor(stage2_df$corrected_z, stage2_df$true_slope_dev),
      closed_form_true_cor = safe_cor(stage2_df$corrected_slope_full, stage2_df$true_slope_dev),
      diag_corrected_failure_rate = failure_rate(stage2_df$corrected_z_diag),
      matrix_corrected_failure_rate = failure_rate(stage2_df$corrected_z),
      closed_form_failure_rate = failure_rate(stage2_df$corrected_slope_full),
      cluster_size_x_cor = cluster_size_x_cor
    ),
    calibrated_dgp
  )
}

#' Diagnose whether a calibrated replication realizes its intended DGP.
#'
#' @details
#' These diagnostics are computed from the true simulated subject effects and
#' the fixed condition-level calibration parameters, not from fitted Stage-1
#' variance-component estimates:
#'
#' - `empirical_g_error` is the largest absolute elementwise difference
#'   between the empirical covariance of `(true_intercept_dev, true_slope_dev)`
#'   and the target marginal G matrix.
#' - `empirical_structural_r2` is the sample R-squared from regressing the true
#'   total slope on the simulated standardized predictor `x`.
#' - `empirical_reliability` recomputes posterior reliability using the
#'   replication's realized `Z_i` and generated `R_i`, while holding the target
#'   marginal G fixed.
#' - `residual_g_min_eigenvalue` is the minimum eigenvalue of the fixed residual
#'   G matrix used to draw `(u0, u1)`.
#'
#' Finite-sample empirical covariance, R-squared, and design reliability should
#' fluctuate around their condition targets. The residual-G eigenvalue is fixed
#' across replications because population parameters are calibrated once per
#' condition.
#'
#' @param condition Optional one-row calibrated design condition.
#' @param sim Simulated dataset list from `simulate_dataset()`.
#'
#' @return One-row tibble with four calibrated-DGP diagnostics. Returns all
#'   missing values for legacy conditions.
make_calibrated_dgp_diagnostics <- function(condition, sim) {
  empty <- tibble::tibble(
    empirical_g_error = NA_real_,
    empirical_structural_r2 = NA_real_,
    empirical_reliability = NA_real_,
    residual_g_min_eigenvalue = NA_real_
  )

  required <- c(
    "calibration_tau0", "slope_variance_marginal",
    "slope_variance_residual", "marginal_rho", "rho_residual"
  )
  if (is.null(condition) || !all(required %in% names(condition)) ||
      is.null(sim) || is.null(sim$id_df) || is.null(sim$dat) ||
      is.null(sim$R_list)) {
    return(empty)
  }

  tau0 <- as.numeric(condition$calibration_tau0[[1]])
  slope_variance_marginal <- as.numeric(condition$slope_variance_marginal[[1]])
  slope_variance_residual <- as.numeric(condition$slope_variance_residual[[1]])
  marginal_rho <- as.numeric(condition$marginal_rho[[1]])
  residual_rho <- as.numeric(condition$rho_residual[[1]])

  G_marginal <- make_random_effect_covariance(
    intercept_variance = tau0^2,
    slope_variance = slope_variance_marginal,
    intercept_slope_correlation = marginal_rho
  )
  G_residual <- make_random_effect_covariance(
    intercept_variance = tau0^2,
    slope_variance = slope_variance_residual,
    intercept_slope_correlation = residual_rho
  )

  true_effects <- sim$id_df[, c("true_intercept_dev", "true_slope_dev"), drop = FALSE]
  complete_effects <- stats::complete.cases(true_effects)
  empirical_g_error <- if (sum(complete_effects) >= 2L) {
    empirical_G <- stats::cov(true_effects[complete_effects, , drop = FALSE])
    max(abs(empirical_G - G_marginal))
  } else {
    NA_real_
  }

  structural_complete <- stats::complete.cases(
    sim$id_df[, c("x", "true_slope_dev"), drop = FALSE]
  )
  empirical_structural_r2 <- if (sum(structural_complete) >= 3L &&
      stats::sd(sim$id_df$x[structural_complete]) > sqrt(.Machine$double.eps) &&
      stats::sd(sim$id_df$true_slope_dev[structural_complete]) >
        sqrt(.Machine$double.eps)) {
    stats::cor(
      sim$id_df$x[structural_complete],
      sim$id_df$true_slope_dev[structural_complete]
    )^2
  } else {
    NA_real_
  }

  ordered_ids <- as.character(sim$id_df$id)
  split_dat <- split(sim$dat, as.character(sim$dat$id), drop = TRUE)[ordered_ids]
  Z_list <- lapply(split_dat, function(df_i) {
    cbind(intercept = 1, slope = as.numeric(df_i$z))
  })
  R_list <- sim$R_list[ordered_ids]
  empirical_reliability <- tryCatch(
    expected_slope_reliability(G_marginal, Z_list, R_list),
    error = function(e) NA_real_
  )

  residual_g_min_eigenvalue <- min(eigen(
    G_residual,
    symmetric = TRUE,
    only.values = TRUE
  )$values)

  tibble::tibble(
    empirical_g_error = empirical_g_error,
    empirical_structural_r2 = empirical_structural_r2,
    empirical_reliability = empirical_reliability,
    residual_g_min_eigenvalue = residual_g_min_eigenvalue
  )
}

#' Resolve legacy or reliability-calibrated simulation parameters.
#'
#' Reliability-calibrated conditions distinguish the marginal intercept-slope
#' correlation targeted in the Stage-1 G matrix from the residual correlation
#' used to draw `(u0, u1)` before adding `gamma * x`. Legacy conditions contain
#' no `rho_residual` field and retain their original `rho`/`tau1` semantics.
#'
#' @param condition One-row design condition.
#' @param params Shared fixed simulation parameters.
#'
#' @return List containing `truth`, `sim_params`, and residual `tau1`.
resolve_blup_outcome_simulation_parameters <- function(condition, params) {
  truth <- as.numeric(condition$gamma_x_on_slope[[1]])
  sim_params <- params
  sim_params$gamma_x_on_slope <- truth

  calibrated <- all(c(
    "target_reliability", "structural_r2", "marginal_rho",
    "rho_residual", "tau1_residual", "calibration_tau0"
  ) %in% names(condition))

  if (calibrated) {
    sim_params$tau0 <- as.numeric(condition$calibration_tau0[[1]])
    sim_params$rho <- as.numeric(condition$rho_residual[[1]])
    tau1 <- as.numeric(condition$tau1_residual[[1]])
  } else {
    sim_params$rho <- as.numeric(condition$rho[[1]])
    tau1 <- as.numeric(condition$tau1[[1]])
  }

  list(
    calibrated = calibrated,
    truth = truth,
    sim_params = sim_params,
    tau1 = tau1
  )
}

#' Run one BLUP-outcome Monte Carlo replication.
#'
#' @details
#' This function is the core replication pipeline:
#'
#' 1. Simulate a random-intercept/random-slope dataset for one design condition.
#' 2. Fit the null mixed model used for score extraction and a direct MLM
#'    benchmark.
#' 3. Recover several subject-level slope scores: EB/BLUP, matrix prior-
#'    unweighted scores, diagonal-only corrected scores, closed-form OLS scores,
#'    and single-subject OLS slopes.
#' 4. Fit Stage 2 regressions where each recovered slope score is the outcome
#'    and the true level-2 predictor `x` is the predictor.
#' 5. In full mode, add Lai 2S-PA/PAA, Fuller, and stacked-sandwich corrected-
#'    score variants.
#' 6. Attach diagnostics and realized trial-count context to every method row.
#'
#' The target truth is `gamma_x_on_slope`, the population effect of the level-2
#' predictor on the latent random slope.
#'
#' @param condition One-row design condition tibble.
#' @param params List of fixed simulation parameters shared across conditions.
#' @param derivative_backend Derivative backend object used by the stacked
#'   sandwich estimator.
#' @param analysis_mode Character scalar. `"full"` includes Lai and stacked
#'   sandwich estimators; `"screen"` returns only the core estimators.
#'
#' @return Replication-level tibble with one row per estimator method.
run_blup_outcome_rep <- function(condition, params, derivative_backend, analysis_mode = "full") {
  resolved <- resolve_blup_outcome_simulation_parameters(condition, params)
  truth <- resolved$truth
  sim_params <- resolved$sim_params
  r_spec <- condition_to_r_spec(condition)

  sim <- simulate_dataset(
    n_id = condition$n_id[[1]],
    mean_n_trial = condition$mean_n_trial[[1]],
    params = sim_params,
    tau1 = resolved$tau1,
    sigma = condition$sigma[[1]],
    has_random_slope = TRUE,
    balanced = balance_mode_to_sim_arg(condition$balance_mode[[1]]),
    min_n_trial = condition$min_n_trial[[1]],
    highly_unbalanced_min_n_trial = condition$highly_unbalanced_min_n_trial[[1]],
    highly_unbalanced_power = condition$highly_unbalanced_power[[1]],
    r_spec = r_spec
  )

  # The null model supplies all extracted score outcomes. The direct model is a
  # one-stage benchmark for the same x-by-z interaction target.
  fit_null <- safe_lmer(y ~ 1 + z + (1 + z | id), data = sim$dat, REML = FALSE)
  fit_direct <- safe_lmer(y ~ x + z + x:z + (1 + z | id), data = sim$dat, REML = FALSE)
  fit_null_R <- NULL
  if (condition_uses_non_iid_R(condition) && requireNamespace("nlme", quietly = TRUE)) {
    fit_null_R <- safe_lme(
      fixed = y ~ z,
      random = ~1 + z | id,
      data = sim$dat,
      correlation = condition_to_nlme_correlation(condition),
      method = "ML",
      control = nlme::lmeControl(returnObject = TRUE, msMaxIter = 100L, opt = "optim")
    )
    fit_direct <- safe_lme(
      fixed = y ~ x + z + x:z,
      random = ~1 + z | id,
      data = sim$dat,
      correlation = condition_to_nlme_correlation(condition),
      method = "ML",
      control = nlme::lmeControl(returnObject = TRUE, msMaxIter = 100L, opt = "optim")
    )
  }
  if (is.null(fit_null)) {
    return(empty_blup_outcome_result(blup_outcome_methods(analysis_mode), truth) %>%
      add_empty_blup_outcome_context(sim = sim, condition = condition))
  }

  ordered_ids <- as.character(sim$id_df$id)
  split_dat <- split(sim$dat, sim$dat$id)[ordered_ids]
  use_non_iid_R <- condition_uses_non_iid_R(condition)
  stage1_score_fit <- if (use_non_iid_R) fit_null_R else fit_null

  # This table is the authoritative source for every downstream method that
  # needs Stage-1 EB means, posterior covariance, G-aware corrections, or Lai
  # measurement-model inputs. In iid conditions it is computed from the lme4
  # fit and matches ranef(..., condVar = TRUE); in non-iid conditions it is
  # computed from the nlme fit so beta, G, and R come from the same likelihood.
  stage1_components <- if (is.null(stage1_score_fit)) {
    tibble::tibble(id = ordered_ids)
  } else {
    tryCatch(
      get_stage1_eb_components(
        fit_obj = stage1_score_fit,
        data = sim$dat,
        cluster_var = "id",
        outcome_var = "y",
        within_var = "z"
      ),
      error = function(e) tibble::tibble(id = ordered_ids)
    )
  }

  closed_form_scores <- if (is.null(stage1_score_fit)) {
    tibble::tibble(id = ordered_ids)
  } else {
    tryCatch(
      get_closed_form_corrected_scores(
        fit_obj = stage1_score_fit,
        data = sim$dat,
        cluster_var = "id",
        outcome_var = "y",
        within_var = "z"
      ),
      error = function(e) tibble::tibble(id = ordered_ids)
    )
  }

  # Single-subject OLS slopes are an intentionally simple benchmark computed
  # without borrowing strength or mixed-model covariance information.
  ols_slopes <- vapply(split_dat, function(df_i) {
    fit_i <- tryCatch(stats::lm(y ~ z, data = df_i), error = function(e) NULL)
    if (is.null(fit_i) || !("z" %in% names(stats::coef(fit_i)))) {
      return(NA_real_)
    }
    unname(stats::coef(fit_i)[["z"]])
  }, numeric(1))

  # Assemble one row per subject/cluster. Missing score columns are added below
  # so each estimator can fail through complete-case logic rather than through
  # absent-column errors.
  stage2_df <- sim$id_df %>%
    dplyr::left_join(stage1_components %>% dplyr::select(-dplyr::any_of("x")), by = "id") %>%
    dplyr::left_join(closed_form_scores, by = "id") %>%
    dplyr::mutate(ols_slope = ols_slopes[as.character(.data$id)])
  stage2_df <- ensure_blup_outcome_stage2_columns(stage2_df)

  diagnostics <- make_blup_outcome_diagnostics(
    stage1_score_fit,
    stage2_df,
    sim,
    condition = condition
  )

  # Core estimators all regress a recovered random-slope outcome on the true
  # cluster-level predictor x, except `direct_mlm`, which estimates the target
  # interaction directly in the long-format mixed model.
  base_rows <- dplyr::bind_rows(
    fit_score_outcome_ols(stage2_df, outcome = "true_slope_dev", base_method = "oracle", include_hc3 = FALSE),
    fit_score_outcome_ols(stage2_df, outcome = "u1_eb", base_method = "naive_blup", include_hc3 = TRUE),
    fit_score_outcome_ols(stage2_df, outcome = "corrected_z_diag", base_method = "diag_corrected", include_hc3 = TRUE),
    fit_score_outcome_ols(stage2_df, outcome = "corrected_z", base_method = "matrix_corrected", include_hc3 = TRUE),
    fit_score_outcome_ols(stage2_df, outcome = "corrected_slope_full", base_method = "closed_form", include_hc3 = TRUE),
    fit_score_outcome_ols(stage2_df, outcome = "ols_slope", base_method = "single_subject_ols", include_hc3 = TRUE),
    fit_direct_mlm_row(fit_direct, n_id = nrow(sim$id_df))
  )

  full_rows <- if (identical(analysis_mode, "full")) {
    # Lai uses the EB measurement-model representation; the stacked sandwich
    # uses the closed-form corrected scores and propagates Stage 1 uncertainty.
    lai_rows <- dplyr::bind_rows(
      fit_lai_2spa(stage2_df, use_average = FALSE) %>%
        dplyr::mutate(method = "lai_2spa") %>%
        dplyr::select("method", dplyr::everything()),
      fit_lai_2spa(stage2_df, use_average = TRUE) %>%
        dplyr::mutate(method = "lai_2spaa") %>%
        dplyr::select("method", dplyr::everything())
    )
    
    stage2_df_fuller <- stage2_df %>% dplyr::mutate(zero = 0)
    fuller_rows <- dplyr::bind_rows(
      fit_fuller(
        stage2_df_fuller,
        outcome = "u1_eb",
        predictor_u1 = "x",
        meas22 = "zero",
        outcome_meas_var = "postvar22"
      ) %>%
        dplyr::mutate(method = "fuller_blup") %>%
        dplyr::select("method", dplyr::everything()),
      fit_fuller(
        stage2_df_fuller,
        outcome = "corrected_z_diag",
        predictor_u1 = "x",
        meas22 = "zero",
        outcome_meas_var = "corrected_z_diag_var"
      ) %>%
        dplyr::mutate(method = "fuller_diag_corrected") %>%
        dplyr::select("method", dplyr::everything()),
      fit_fuller(
        stage2_df_fuller,
        outcome = "corrected_z",
        predictor_u1 = "x",
        meas22 = "zero",
        outcome_meas_var = "corrected_z_var"
      ) %>%
        dplyr::mutate(method = "fuller_matrix_corrected") %>%
        dplyr::select("method", dplyr::everything()),
      fit_fuller(
        stage2_df_fuller,
        outcome = "corrected_slope_full",
        predictor_u1 = "x",
        meas22 = "zero",
        outcome_meas_var = "ols_var22"
      ) %>%
        dplyr::mutate(method = "fuller_closed_form") %>%
        dplyr::select("method", dplyr::everything()),
      fit_fuller_dual_alpha_stepdown(
        stage2_df_fuller,
        outcome = "u1_eb",
        predictor_u0 = NULL,
        predictor_u1 = "x",
        meas11 = NULL,
        meas12 = NULL,
        meas22 = "zero",
        outcome_meas_var = "postvar22"
      ) %>%
        dplyr::mutate(method = "fuller_alpha_blup") %>%
        dplyr::select("method", dplyr::everything()),
      fit_fuller_dual_alpha_stepdown(
        stage2_df_fuller,
        outcome = "corrected_z_diag",
        predictor_u0 = NULL,
        predictor_u1 = "x",
        meas11 = NULL,
        meas12 = NULL,
        meas22 = "zero",
        outcome_meas_var = "corrected_z_diag_var"
      ) %>%
        dplyr::mutate(method = "fuller_alpha_diag_corrected") %>%
        dplyr::select("method", dplyr::everything()),
      fit_fuller_dual_alpha_stepdown(
        stage2_df_fuller,
        outcome = "corrected_z",
        predictor_u0 = NULL,
        predictor_u1 = "x",
        meas11 = NULL,
        meas12 = NULL,
        meas22 = "zero",
        outcome_meas_var = "corrected_z_var"
      ) %>%
        dplyr::mutate(method = "fuller_alpha_matrix_corrected") %>%
        dplyr::select("method", dplyr::everything()),
      fit_fuller_dual_alpha_stepdown(
        stage2_df_fuller,
        outcome = "corrected_slope_full",
        predictor_u0 = NULL,
        predictor_u1 = "x",
        meas11 = NULL,
        meas12 = NULL,
        meas22 = "zero",
        outcome_meas_var = "ols_var22"
      ) %>%
        dplyr::mutate(method = "fuller_alpha_closed_form") %>%
        dplyr::select("method", dplyr::everything())
    )

    sandwich_rows <- tryCatch({
      stacked_R_list <- if (identical(sim$r_spec$structure, "iid")) NULL else sim$R_list
      sandwich_out <- stacked_sandwich_for_corrected_scores(
        split_dat = split_dat,
        id_df = sim$id_df,
        fit_null = fit_null,
        psi_hat = pack_psi(fit_null),
        derivative_backend = derivative_backend,
        R_list = stacked_R_list
      )
      format_stacked_sandwich_rows(
        sandwich_out = sandwich_out,
        df = nrow(sim$id_df) - 2L,
        method_prefix = "closed_form_stacked"
      ) %>%
        dplyr::mutate(status_code = dplyr::if_else(is.finite(.data$estimate) & is.finite(.data$se), 0L, NA_integer_))
    }, error = function(e) empty_stacked_rows("closed_form_stacked"))

    dplyr::bind_rows(lai_rows, fuller_rows, sandwich_rows)
  } else {
    tibble::tibble()
  }

  # Every method row receives the same replication-level diagnostics and
  # realized-trial descriptors so downstream summaries can group by method while
  # retaining information about the Stage 1 data quality.
  dplyr::bind_rows(base_rows, full_rows) %>%
    standardize_estimator_rows(truth = truth) %>%
    dplyr::bind_cols(diagnostics[rep(1L, nrow(.)), , drop = FALSE]) %>%
    dplyr::mutate(
      truth = truth,
      mean_realized_trials = sim$mean_realized_trials,
      min_realized_trials = sim$min_realized_trials,
      prop_ids_leq_2_trials = sim$prop_ids_leq_2_trials,
      prop_ids_leq_3_trials = sim$prop_ids_leq_3_trials
    )
}
