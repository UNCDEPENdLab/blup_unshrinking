#' Vig-Hallquist (2026) simulation Study 3: BLUP as predictor and outcome.

study3_methods <- function() {
  c(
    "oracle_dual",
    "naive_blup_on_blup",
    "closed_form_on_blup",
    "blup_on_closed_form",
    "closed_form_on_closed_form",
    "fuller_closed_form",
    "fuller_alpha_stepdown_closed_form",
    "lai_2spa", "lai_2spaa"
  )
}

simulate_study3 <- function(condition) {
  simulate_data_dual_blup(condition)
}

run_study3_rep <- function(condition) {
  truth <- as.numeric(condition$standardized_beta_target[[1]])
  reporting_scale <- as.numeric(
    condition$tau1_y[[1]] / condition$tau1_q[[1]]
  )
  sim <- simulate_study3(condition)

  fit_y <- fit_stage1(
    y ~ x,
    random = ~x | cid,
    data = sim$lv1_y,
    condition = condition,
    cluster_var = "cid"
  )
  fit_q <- fit_stage1(
    q ~ x,
    random = ~x | cid,
    data = sim$lv1_q,
    condition = condition,
    cluster_var = "cid"
  )
  if (is.null(fit_y) || is.null(fit_q)) {
    return(make_failed_result(condition, study3_methods(), truth))
  }

  y_eb <- tryCatch(
    get_stage1_eb_components(
      fit_obj = fit_y,
      data = sim$lv1_y,
      cluster_var = "cid",
      outcome_var = "y",
      within_var = "x"
    ),
    error = function(e) NULL
  )
  q_eb <- tryCatch(
    get_stage1_eb_components(
      fit_obj = fit_q,
      data = sim$lv1_q,
      cluster_var = "cid",
      outcome_var = "q",
      within_var = "x"
    ) %>%
      select_lai_measurement_columns(n_re = 2L, prefix = "q_"),
    error = function(e) NULL
  )
  y_corrected <- tryCatch(
    get_closed_form_corrected_scores(
      fit_obj = fit_y,
      data = sim$lv1_y,
      cluster_var = "cid",
      outcome_var = "y",
      within_var = "x"
    ) %>%
      dplyr::select(
        id,
        corrected_intercept_full,
        corrected_slope_full,
        ols_var11,
        ols_var12,
        ols_var22
      ),
    error = function(e) NULL
  )
  q_corrected <- tryCatch(
    get_closed_form_corrected_scores(
      fit_obj = fit_q,
      data = sim$lv1_q,
      cluster_var = "cid",
      outcome_var = "q",
      within_var = "x"
    ) %>%
      dplyr::transmute(
        id,
        q_corrected_intercept = corrected_intercept_full,
        q_corrected_slope = corrected_slope_full,
        q_ols_var11 = ols_var11,
        q_ols_var12 = ols_var12,
        q_ols_var22 = ols_var22
      ),
    error = function(e) NULL
  )
  if (is.null(y_eb) || is.null(q_eb) ||
      is.null(y_corrected) || is.null(q_corrected)) {
    return(make_failed_result(condition, study3_methods(), truth))
  }

  stage2_df <- sim$lv2_true %>%
    dplyr::left_join(y_eb, by = "id") %>%
    dplyr::left_join(q_eb, by = "id") %>%
    dplyr::left_join(y_corrected, by = "id") %>%
    dplyr::left_join(q_corrected, by = "id")
  stage1_diag <- combine_dual_stage1_diagnostics(
    get_stage1_diagnostics(
      fit_y,
      stage2_df,
      predictor_u0 = "u0_eb",
      predictor_u1 = "u1_eb"
    ),
    get_stage1_diagnostics(
      fit_q,
      stage2_df,
      predictor_u0 = "q_u0_eb",
      predictor_u1 = "q_u1_eb"
    )
  )

  results <- dplyr::bind_rows(
    fit_observed_dual(
      stage2_df,
      outcome = "true_q1",
      predictor_u0 = "true_y0",
      predictor_u1 = "true_y1",
      reporting_scale = reporting_scale
    ) %>%
      dplyr::filter(se_type == "naive") %>%
      dplyr::transmute(
        method = "oracle_dual",
        estimate, se, ci_low, ci_high, status_code
      ),
    finalize_ols_se_variants(
      fit_observed_dual(
        stage2_df,
        outcome = "q_u1_eb",
        predictor_u0 = "u0_eb",
        predictor_u1 = "u1_eb",
        reporting_scale = reporting_scale
      ),
      "naive_blup_on_blup"
    ),
    finalize_ols_se_variants(
      fit_observed_dual(
        stage2_df,
        outcome = "q_corrected_slope",
        predictor_u0 = "u0_eb",
        predictor_u1 = "u1_eb",
        reporting_scale = reporting_scale
      ),
      "closed_form_on_blup"
    ),
    finalize_ols_se_variants(
      fit_observed_dual(
        stage2_df,
        outcome = "q_u1_eb",
        predictor_u0 = "corrected_intercept_full",
        predictor_u1 = "corrected_slope_full",
        reporting_scale = reporting_scale
      ),
      "blup_on_closed_form"
    ),
    finalize_ols_se_variants(
      fit_observed_dual(
        stage2_df,
        outcome = "q_corrected_slope",
        predictor_u0 = "corrected_intercept_full",
        predictor_u1 = "corrected_slope_full",
        reporting_scale = reporting_scale
      ),
      "closed_form_on_closed_form"
    ),
    fit_fuller_dual(
      stage2_df,
      outcome = "q_corrected_slope",
      predictor_u0 = "corrected_intercept_full",
      predictor_u1 = "corrected_slope_full",
      meas11 = "ols_var11",
      meas12 = "ols_var12",
      meas22 = "ols_var22",
      outcome_meas_var = "q_ols_var22"
    ) %>%
      rescale_fuller_to_population_sd(reporting_scale) %>%
      dplyr::mutate(method = "fuller_closed_form") %>%
      dplyr::select(method, dplyr::everything()),
    fit_fuller_dual_alpha_stepdown(
      stage2_df,
      outcome = "q_corrected_slope",
      predictor_u0 = "corrected_intercept_full",
      predictor_u1 = "corrected_slope_full",
      meas11 = "ols_var11",
      meas12 = "ols_var12",
      meas22 = "ols_var22",
      outcome_meas_var = "q_ols_var22"
    ) %>%
      rescale_fuller_to_population_sd(reporting_scale) %>%
      dplyr::mutate(method = "fuller_alpha_stepdown_closed_form") %>%
      dplyr::select(method, dplyr::everything()),
    fit_lai_2spa_dual_process(
      stage2_df,
      use_average = FALSE,
      theta0_start = condition$theta0[[1]],
      theta1_start = condition$theta1[[1]],
      reporting_scale = reporting_scale
    ) %>%
      dplyr::mutate(method = "lai_2spa") %>%
      dplyr::select(method, dplyr::everything()),
    fit_lai_2spa_dual_process(
      stage2_df,
      use_average = TRUE,
      theta0_start = condition$theta0[[1]],
      theta1_start = condition$theta1[[1]],
      reporting_scale = reporting_scale
    ) %>%
      dplyr::mutate(method = "lai_2spaa") %>%
      dplyr::select(method, dplyr::everything())
  )

  dplyr::bind_cols(
    results,
    stage1_diag[rep(1L, nrow(results)), , drop = FALSE]
  ) %>%
    add_study_result_context(condition, truth)
}
