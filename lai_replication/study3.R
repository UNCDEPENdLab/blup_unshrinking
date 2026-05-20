#' Lai Simulation Study 3: disparate clustering for Y random slopes and repeated Z.
#'
#' Study 3 differs from the matched-clustering Study 1/2 setup because the
#' first-stage random effects are estimated from one repeated-measures outcome
#' (`Y`), while the downstream outcome (`Z`) is observed in a separate
#' repeated-measures process with its own low-count cluster structure. The
#' replication therefore has to fit two Stage 1 models, recover cluster-level
#' information from each, and assemble a subject-level Stage 2 data set where
#' corrected `Y` random slopes predict corrected or EB-derived `Z` outcomes.

#' Construct the operational Study 3 cluster-size profile for Y.
#'
#' @details
#' The Lai Study 3 materials do not give a single perfectly consistent
#' cluster-size definition across the paper, notebook text, and supplement
#' script. The repository design notes treat the supplement script
#' `simulation_scripts/sim3_revise.R` as the operational reference. This helper
#' encodes that profile as a repeating sequence of small clusters followed by
#' larger clusters.
#'
#' @param num_clus Integer-like scalar number of clusters/subjects.
#'
#' @return Integer vector of length `num_clus` giving the number of `Y`
#'   observations per cluster.
make_study3_cluster_sizes <- function(num_clus) {
  rep_len(
    c(rep(2:4, each = 5), rep(18:20, each = 5)),
    length.out = num_clus
  )
}

#' Simulate one Lai Study 3 data set.
#'
#' @details
#' The simulated data contain two related but distinct observation-level data
#' frames:
#'
#' - `lv1_y`: repeated `Y` observations used to estimate cluster-specific
#'   random intercepts and random slopes for predictor `x`;
#' - `lv1_z`: repeated `Z` observations used to estimate a separate noisy
#'   cluster-level outcome; and
#' - `lv2_true`: the latent subject-level random effects and latent `Z`
#'   outcome used by oracle comparisons.
#'
#' Random effects `u0` and `u1` are drawn from the condition-specific covariance
#' matrix returned by `make_covu()`. They are sorted by `u1` to mirror the
#' supplement script's deterministic relationship between random-slope rank and
#' the heteroscedastic `x` generator. The latent outcome `eta_z` is a linear
#' function of the true random intercept and slope, with residual variance
#' chosen so the total latent outcome variance remains approximately one.
#'
#' @param condition One-row Study 3 design condition containing `num_clus`,
#'   `sigma2`, `sigma_z`, `icc`, `var_u1`, `cor_u0_u1`, and `beta_zu1`.
#'
#' @return A list with `lv1_y`, `lv1_z`, and `lv2_true` tibbles.
simulate_study3 <- function(condition) {
  covu <- make_covu(condition)
  cluster_sizes <- make_study3_cluster_sizes(as.integer(condition$num_clus))
  cid <- rep(seq_len(as.integer(condition$num_clus)), cluster_sizes)

  # The supplement uses cluster-index-dependent x variance. Because the random
  # effects are sorted by slope below, this creates the intended Study 3 link
  # between cluster information quality and latent random-slope ordering.
  x <- stats::rnorm(length(cid), sd = sqrt(2.25 - 1.5 * cid / as.integer(condition$num_clus)))

  u <- MASS::mvrnorm(as.integer(condition$num_clus), mu = c(0, 0), Sigma = covu)
  u <- u[order(u[, 2]), , drop = FALSE]

  # First repeated-measures process: Y is generated from a random-intercept and
  # random-slope model and supplies the EB/corrected random-effect predictors.
  y <- fixed_params$gamma0 + fixed_params$gamma1 * x + rowSums(cbind(1, x) * u[cid, , drop = FALSE]) +
    draw_lai_level1_residuals(cluster_sizes, sigma = sqrt(condition$sigma2), condition = condition)

  # Second repeated-measures process: eta_z is the latent Stage 2 outcome
  # implied by the true random effects, then z adds sparse measurement noise.
  ev_z <- 1 - drop(t(c(fixed_params$beta_zu0, condition$beta_zu1)) %*% covu %*% c(fixed_params$beta_zu0, condition$beta_zu1))
  eta_z <- fixed_params$z_intercept + fixed_params$beta_zu0 * u[, 1] + condition$beta_zu1 * u[, 2] +
    stats::rnorm(as.integer(condition$num_clus), sd = sqrt(ev_z))
  cid_z_sizes <- rep_len(2:3, as.integer(condition$num_clus))
  cid_z <- rep(seq_len(as.integer(condition$num_clus)), cid_z_sizes)
  z <- eta_z[cid_z] + draw_lai_level1_residuals(cid_z_sizes, sigma = condition$sigma_z, condition = condition)

  list(
    lv1_y = tibble::tibble(
      cid = factor(cid),
      cid_chr = as.character(cid),
      trial_index = ave(seq_along(cid), cid, FUN = seq_along),
      x = x,
      y = y
    ),
    lv1_z = tibble::tibble(
      cid_z = factor(cid_z),
      cid_z_chr = as.character(cid_z),
      trial_index = ave(seq_along(cid_z), cid_z, FUN = seq_along),
      z = z
    ),
    lv2_true = tibble::tibble(
      id = as.character(seq_len(as.integer(condition$num_clus))),
      true_u0 = u[, 1],
      true_u1 = u[, 2],
      eta_z = eta_z
    )
  )
}

#' Run one Monte Carlo replication for Lai Study 3.
#'
#' @details
#' This is the Study 3 per-replication estimator pipeline. It fits separate
#' Stage 1 models for `Y` and `Z`, extracts EB measurement inputs and
#' likelihood-only corrected scores, merges those cluster-level summaries, and
#' evaluates the disparate-clustering methods:
#'
#' - oracle regression on the true latent `eta_z`, `u0`, and `u1`;
#' - naive dual-EB regression using EB `Z` and EB `Y` random effects;
#' - EIV variants using corrected `Y` scores and their OLS sampling covariance;
#' - ridge-stabilized dual-EB regression;
#' - corrected-score dual regression; and
#' - Lai 2S-PA/2S-PAA for disparate first-stage `Y` and `Z` measurements.
#'
#' If either Stage 1 mixed model fails, the function returns a full method set
#' with missing estimates so the runner preserves the replication row shape.
#'
#' @param condition One-row Study 3 design condition.
#'
#' @return Tibble with one row per Study 3 estimator, including estimate,
#'   standard error, confidence interval, status code, truth, study label, and
#'   Stage 1 diagnostics.
run_study3_rep <- function(condition) {
  truth <- lai_truth(condition)
  include_tempered_eiv <- condition_includes_tempered_eiv(condition)
  sim <- simulate_study3(condition)

  # Study 3 has genuinely disparate first-stage models: random intercept/slope
  # for Y, but random intercept only for the repeated Z measurements.
  sim$lv1_y <- add_lai_trial_index(sim$lv1_y, cluster_var = "cid")
  sim$lv1_z <- add_lai_trial_index(sim$lv1_z, cluster_var = "cid_z")
  fit_y <- fit_lai_stage1(y ~ x, random = ~x | cid, data = sim$lv1_y, condition = condition, cluster_var = "cid")
  fit_z <- fit_lai_stage1(z ~ 1, random = ~1 | cid_z, data = sim$lv1_z, condition = condition, cluster_var = "cid_z")

  if (is.null(fit_y) || is.null(fit_z)) {
    return(make_failed_result(condition, disparate_study_methods(include_tempered_eiv), truth))
  }

  ordered_ids <- sim$lv2_true$id
  eb_inputs_y <- get_stage1_eb_components(
    fit_obj = fit_y,
    data = sim$lv1_y,
    cluster_var = "cid",
    outcome_var = "y",
    within_var = "x"
  )
  eb_inputs_z <- get_stage1_eb_components(
    fit_obj = fit_z,
    data = sim$lv1_z,
    cluster_var = "cid_z",
    outcome_var = "z",
    within_var = NULL
  ) %>%
    select_lai_measurement_columns(n_re = 1L, prefix = "z_")

  # Corrected Y scores are likelihood-only estimates of the random intercept
  # and slope. Their OLS sampling covariance feeds the EIV estimators below.
  corrected_y <- get_closed_form_corrected_scores(
    fit_obj = fit_y,
    data = sim$lv1_y,
    cluster_var = "cid",
    outcome_var = "y",
    within_var = "x"
  ) %>%
    dplyr::select(id, corrected_intercept_full, corrected_slope_full, ols_var11, ols_var12, ols_var22)

  # Corrected Z is the likelihood-only cluster mean after removing the fitted
  # fixed intercept from the repeated-Z model.
  corrected_z <- get_closed_form_corrected_scores(
    fit_obj = fit_z,
    data = sim$lv1_z,
    cluster_var = "cid_z",
    outcome_var = "z",
    within_var = NULL
  ) %>%
    dplyr::select(id, corrected_outcome = corrected)

  # Assemble the cluster-level Stage 2 file. The `z_u0_eb` column comes from the
  # prefixed univariate EB extraction and acts as the measured downstream
  # outcome for naive/ridge EB comparisons.
  stage2_df <- sim$lv2_true %>%
    dplyr::left_join(eb_inputs_y, by = "id") %>%
    dplyr::left_join(eb_inputs_z, by = "id") %>%
    dplyr::left_join(corrected_y, by = "id") %>%
    dplyr::left_join(corrected_z, by = "id")
  stage1_diag <- get_stage1_diagnostics(fit_y, stage2_df)

  results <- dplyr::bind_rows(
    # Oracle target: latent Z outcome regressed on the true Y random effects.
    fit_observed_dual(stage2_df, outcome = "eta_z", predictor_u0 = "true_u0", predictor_u1 = "true_u1") %>%
      dplyr::filter(se_type == "naive") %>%
      dplyr::transmute(method = "oracle_dual", estimate, se, ci_low, ci_high, status_code),
    # Naive EB target: measured EB Z outcome regressed on EB Y random effects.
    finalize_ols_se_variants(fit_observed_dual(stage2_df, outcome = "z_u0_eb", predictor_u0 = "u0_eb", predictor_u1 = "u1_eb"), "naive_dual_eb"),
    # EIV target: corrected Z outcome regressed on corrected Y scores, with the
    # corrected Y score covariance subtracted from the predictor cross-products.
    fit_eiv_dual(
      stage2_df,
      outcome = "corrected_outcome",
      predictor_u0 = "corrected_intercept_full",
      predictor_u1 = "corrected_slope_full",
      meas11 = "ols_var11",
      meas12 = "ols_var12",
      meas22 = "ols_var22"
    ) %>%
      finalize_eiv_se_variants("eiv_dual_corrected"),
    fit_eiv_dual(
      stage2_df,
      outcome = "corrected_outcome",
      predictor_u0 = "corrected_intercept_full",
      predictor_u1 = "corrected_slope_full",
      meas11 = "ols_var11",
      meas12 = "ols_var12",
      meas22 = "ols_var22",
      stabilize_a_mat = TRUE
    ) %>%
      finalize_eiv_se_variants("eiv_dual_corrected_nearpd"),
    fit_eiv_dual(
      stage2_df,
      outcome = "corrected_outcome",
      predictor_u0 = "corrected_intercept_full",
      predictor_u1 = "corrected_slope_full",
      meas11 = "ols_var11",
      meas12 = "ols_var12",
      meas22 = "ols_var22",
      ridge_predictor_block = TRUE
    ) %>%
      finalize_eiv_se_variants("eiv_dual_corrected_ridge"),
    if (isTRUE(include_tempered_eiv)) {
      # Optional sensitivity path: partial measurement-error subtraction can be
      # useful when the full EIV correction is numerically too aggressive.
      fit_tempered_eiv_dual_set(
        stage2_df,
        outcome = "corrected_outcome",
        predictor_u0 = "corrected_intercept_full",
        predictor_u1 = "corrected_slope_full",
        meas11 = "ols_var11",
        meas12 = "ols_var12",
        meas22 = "ols_var22"
      )
    },
    # Ridge is a point-estimation stabilization benchmark for the dual-EB
    # predictors; it does not supply analytic standard errors.
    fit_ridge_dual(stage2_df, outcome = "z_u0_eb", predictor_u0 = "u0_eb", predictor_u1 = "u1_eb") %>%
      dplyr::mutate(method = "ridge_dual_eb") %>%
      dplyr::select(method, estimate, se, ci_low, ci_high, status_code),
    # Corrected-score regression is the simple Vig-style comparison for the
    # disparate-clustering setting.
    finalize_ols_se_variants(fit_observed_dual(stage2_df, outcome = "corrected_outcome", predictor_u0 = "corrected_intercept_full", predictor_u1 = "corrected_slope_full"), "corrected_dual"),
    # Lai's disparate OpenMx model carries both the Y and Z first-stage
    # measurement structures into the Stage 2 path model.
    fit_lai_2spa_disparate(stage2_df, use_average = FALSE) %>%
      dplyr::mutate(method = "lai_2spa") %>%
      dplyr::select(method, dplyr::everything()),
    fit_lai_2spa_disparate(stage2_df, use_average = TRUE) %>%
      dplyr::mutate(method = "lai_2spaa") %>%
      dplyr::select(method, dplyr::everything()),
    fit_fuller_dual(
      stage2_df,
      outcome = "corrected_outcome",
      predictor_u0 = "corrected_intercept_full",
      predictor_u1 = "corrected_slope_full",
      meas11 = "ols_var11",
      meas12 = "ols_var12",
      meas22 = "ols_var22"
    ) %>%
      dplyr::mutate(method = "fuller") %>%
      dplyr::select(method, dplyr::everything()),
    fit_fuller_dual_stepdown(
      stage2_df,
      outcome = "corrected_outcome",
      predictor_u0 = "corrected_intercept_full",
      predictor_u1 = "corrected_slope_full",
      meas11 = "ols_var11",
      meas12 = "ols_var12",
      meas22 = "ols_var22"
    ) %>%
      dplyr::mutate(method = "fuller_stepdown") %>%
      dplyr::select(method, dplyr::everything()),
    fit_fuller_dual_alpha_stepdown(
      stage2_df,
      outcome = "corrected_outcome",
      predictor_u0 = "corrected_intercept_full",
      predictor_u1 = "corrected_slope_full",
      meas11 = "ols_var11",
      meas12 = "ols_var12",
      meas22 = "ols_var22"
    ) %>%
      dplyr::mutate(method = "fuller_alpha_stepdown") %>%
      dplyr::select(method, dplyr::everything())
  )

  dplyr::bind_cols(results, stage1_diag[rep(1L, nrow(results)), , drop = FALSE]) %>%
    add_study_result_context(condition, truth)
}
