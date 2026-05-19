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
      cluster_size_x_cor = NA_real_
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
#'
#' @return `results` with diagnostic and realized-trial columns appended.
add_empty_blup_outcome_context <- function(results, sim = NULL) {
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

  results %>%
    dplyr::bind_cols(empty_blup_outcome_diagnostics()[rep(1L, nrow(results)), , drop = FALSE]) %>%
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
    "x", "true_slope_dev", "u0_eb", "u1_eb", "postvar22", "theta22", "lambda22",
    "corrected_z", "corrected_z_diag", "corrected_slope_full", "ols_slope"
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

#' Extract the direct mixed-model benchmark row.
#'
#' @details
#' The direct MLM benchmark estimates the target interaction in one model,
#' `y ~ x + z + x:z + (1 + z | id)`. The reported row uses a t critical value
#' with approximate residual degrees of freedom based on the number of clusters
#' minus the number of fixed effects.
#'
#' @param fit_direct Fitted `lme4` model or `NULL` if fitting failed.
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

  extract_lmer_stats(fit_direct, term = "x:z", use_t = TRUE, df = n_id - length(lme4::fixef(fit_direct))) %>%
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
#' @param fit_null Fitted null `lme4` model used to extract BLUP/corrected
#'   scores.
#' @param stage2_df Cluster-level Stage 2 data frame.
#' @param sim Simulated dataset list from `simulate_dataset()`.
#'
#' @return One-row tibble containing shared Stage 1 diagnostics plus
#'   BLUP-outcome-specific score diagnostics.
make_blup_outcome_diagnostics <- function(fit_null, stage2_df, sim) {
  stage1_diag <- get_stage1_diagnostics(fit_null, stage2_df, predictor_u0 = "u0_eb", predictor_u1 = "u1_eb")

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
    )
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
  truth <- as.numeric(condition$gamma_x_on_slope[[1]])
  sim_params <- params
  sim_params$gamma_x_on_slope <- truth
  sim_params$rho <- as.numeric(condition$rho[[1]])

  sim <- simulate_dataset(
    n_id = condition$n_id[[1]],
    mean_n_trial = condition$mean_n_trial[[1]],
    params = sim_params,
    tau1 = condition$tau1[[1]],
    sigma = condition$sigma[[1]],
    has_random_slope = TRUE,
    balanced = balance_mode_to_sim_arg(condition$balance_mode[[1]]),
    min_n_trial = condition$min_n_trial[[1]],
    highly_unbalanced_min_n_trial = condition$highly_unbalanced_min_n_trial[[1]],
    highly_unbalanced_power = condition$highly_unbalanced_power[[1]]
  )

  # The null model supplies all extracted score outcomes. The direct model is a
  # one-stage benchmark for the same x-by-z interaction target.
  fit_null <- safe_lmer(y ~ 1 + z + (1 + z | id), data = sim$dat, REML = FALSE)
  fit_direct <- safe_lmer(y ~ x + z + x:z + (1 + z | id), data = sim$dat, REML = FALSE)
  if (is.null(fit_null)) {
    return(empty_blup_outcome_result(blup_outcome_methods(analysis_mode), truth) %>%
      add_empty_blup_outcome_context(sim = sim))
  }

  ordered_ids <- as.character(sim$id_df$id)
  split_dat <- split(sim$dat, sim$dat$id)[ordered_ids]

  # Lai/OpenMx inputs contain EB scores plus cluster-specific loading/theta
  # definitions. Failures leave an id-only tibble so later joins still work.
  eb_inputs <- tryCatch(
    compute_lai_2spa_inputs(fit_null = fit_null, split_dat = split_dat, id_df = sim$id_df, R_list = sim$R_list),
    error = function(e) tibble::tibble(id = ordered_ids)
  )
  eb_inputs <- eb_inputs %>% dplyr::select(-dplyr::any_of("x"))

  # These three score paths represent increasingly direct ways to recover the
  # likelihood-only cluster slope: GLS-aware matrix prior-unweighting,
  # diagonal-only lme4 prior-unweighting, and direct GLS closed form. The
  # matrix path computes EB/posterior ingredients from the same R_i used to
  # generate the data so non-diagonal residual covariance is represented.
  matrix_scores <- tryCatch(
    get_gls_corrected_scores(
      fit_obj = fit_null,
      data = sim$dat,
      cluster_var = "id",
      outcome_var = "y",
      within_var = "z",
      R_list = sim$R_list
    ),
    error = function(e) tibble::tibble(id = ordered_ids)
  )
  diag_scores <- tryCatch(get_diagonal_corrected_scores(fit_null), error = function(e) tibble::tibble(id = ordered_ids))
  closed_form_scores <- get_closed_form_corrected_scores(
    fit_obj = fit_null,
    data = sim$dat,
    cluster_var = "id",
    outcome_var = "y",
    within_var = "z",
    R_list = sim$R_list
  )

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
    dplyr::left_join(eb_inputs, by = "id") %>%
    dplyr::left_join(matrix_scores %>% dplyr::select("id", dplyr::starts_with("corrected_")), by = "id") %>%
    dplyr::left_join(diag_scores, by = "id") %>%
    dplyr::left_join(closed_form_scores, by = "id") %>%
    dplyr::mutate(ols_slope = ols_slopes[as.character(.data$id)])
  stage2_df <- ensure_blup_outcome_stage2_columns(stage2_df)

  diagnostics <- make_blup_outcome_diagnostics(fit_null, stage2_df, sim)

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
        dplyr::select("method", dplyr::everything())
    )

    sandwich_rows <- tryCatch({
      sandwich_out <- stacked_sandwich_for_corrected_scores(
        split_dat = split_dat,
        id_df = sim$id_df,
        fit_null = fit_null,
        psi_hat = pack_psi(fit_null),
        derivative_backend = derivative_backend
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
