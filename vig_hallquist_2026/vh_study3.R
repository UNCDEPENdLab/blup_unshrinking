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
    "lai_2spa",
    "sem"
  )
}

#' Classify Study-3 rows for point and interval performance summaries.
#'
#' `analysis_eligible` is deliberately a point-estimate criterion: it controls
#' bias and RMSE summaries. `interval_eligible` additionally requires a finite,
#' positive standard error and therefore excludes results without a usable
#' interval estimate from variance and coverage summaries without discarding
#' their point estimates.
add_study3_analysis_eligibility <- function(results) {
  get_column <- function(name, default) {
    if (name %in% names(results)) results[[name]] else rep(default, nrow(results))
  }
  set_reason <- function(reason, when, value) {
    index <- is.na(reason) & !is.na(when) & when
    if (length(value) == 1L) {
      reason[index] <- value
    } else {
      reason[index] <- value[index]
    }
    reason
  }

  method <- as.character(get_column("method", NA_character_))
  status_code <- suppressWarnings(as.integer(get_column("status_code", NA_integer_)))
  estimate <- suppressWarnings(as.numeric(get_column("estimate", NA_real_)))
  se <- suppressWarnings(as.numeric(get_column("se", NA_real_)))
  dual_eligible <- as.logical(get_column("analysis_eligible", NA))
  dual_reason <- as.character(get_column("analysis_exclusion_reason", NA_character_))
  fuller_guard_pass <- as.logical(get_column("fuller_auto_guard_pass", NA))
  fuller_guard_reason <- as.character(get_column("fuller_auto_guard_reason", NA_character_))
  mx_issue_class <- as.character(get_column("mx_issue_class", NA_character_))
  mx_info_definite <- as.logical(get_column("mx_info_definite", NA))
  mx_condition_number <- suppressWarnings(as.numeric(get_column("mx_condition_number", NA_real_)))
  mplus_critical_warning <- as.logical(get_column("mplus_critical_warning", NA))
  mplus_target_parameter_count <- suppressWarnings(as.integer(get_column("mplus_target_parameter_count", NA_integer_)))

  point_reason <- rep(NA_character_, nrow(results))
  point_reason <- set_reason(point_reason, is.na(status_code) | is.na(estimate), "estimation_unavailable")
  point_reason <- set_reason(point_reason, !is.na(status_code) & status_code != 0L, "estimation_status_nonzero")
  point_reason <- set_reason(point_reason, !is.finite(estimate), "nonfinite_estimate")

  dual_ols_methods <- c(
    "oracle_dual", "naive_blup_on_blup", "closed_form_on_blup",
    "blup_on_closed_form", "closed_form_on_closed_form"
  )
  dual_bad <- method %in% dual_ols_methods & !is.na(dual_eligible) & !dual_eligible
  point_reason <- set_reason(
    point_reason,
    dual_bad,
    ifelse(is.na(dual_reason), "stage2_design_ineligible", dual_reason)
  )

  alpha_fuller <- method == "fuller_alpha_stepdown_closed_form"
  point_reason <- set_reason(
    point_reason,
    alpha_fuller & (is.na(fuller_guard_pass) | !fuller_guard_pass),
    ifelse(
      is.na(fuller_guard_reason) | fuller_guard_reason == "",
      "fuller_guard_failed",
      paste0("fuller_guard_", fuller_guard_reason)
    )
  )

  lai <- method == "lai_2spa"
  point_reason <- set_reason(point_reason, lai & !is.na(mx_issue_class) & mx_issue_class != "ok", "openmx_issue")
  point_reason <- set_reason(point_reason, lai & !is.na(mx_info_definite) & !mx_info_definite, "openmx_information_not_definite")
  point_reason <- set_reason(
    point_reason,
    lai & is.finite(mx_condition_number) & mx_condition_number > 1e12,
    "openmx_condition_number_excessive"
  )

  sem <- method == "sem"
  point_reason <- set_reason(
    point_reason,
    sem & !is.na(mplus_critical_warning) & mplus_critical_warning,
    "mplus_critical_warning"
  )
  point_reason <- set_reason(
    point_reason,
    sem & !is.na(mplus_target_parameter_count) & mplus_target_parameter_count != 1L,
    "mplus_target_parameter_not_unique"
  )

  interval_reason <- point_reason
  interval_reason <- set_reason(interval_reason, !is.finite(se) | se <= 0, "invalid_standard_error")

  dplyr::mutate(
    results,
    analysis_eligible = is.na(point_reason),
    analysis_exclusion_reason = point_reason,
    interval_eligible = is.na(interval_reason),
    interval_exclusion_reason = interval_reason
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
    return(
      make_failed_result(condition, study3_methods(), truth) %>%
        add_study3_analysis_eligibility()
    )
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
    return(
      make_failed_result(condition, study3_methods(), truth) %>%
        add_study3_analysis_eligibility()
    )
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
      dplyr::mutate(method = "oracle_dual") %>%
      dplyr::select(method, -se_type),
    fit_observed_dual(
      stage2_df,
      outcome = "q_u1_eb",
      predictor_u0 = "u0_eb",
      predictor_u1 = "u1_eb",
      reporting_scale = reporting_scale
    ) %>%
      dplyr::filter(se_type == "naive") %>%
      dplyr::mutate(method = "naive_blup_on_blup") %>%
      dplyr::select(method, -se_type),
    fit_observed_dual(
      stage2_df,
      outcome = "q_corrected_slope",
      predictor_u0 = "u0_eb",
      predictor_u1 = "u1_eb",
      reporting_scale = reporting_scale
    ) %>%
      dplyr::filter(se_type == "naive") %>%
      dplyr::mutate(method = "closed_form_on_blup") %>%
      dplyr::select(method, -se_type),
    fit_observed_dual(
      stage2_df,
      outcome = "q_u1_eb",
      predictor_u0 = "corrected_intercept_full",
      predictor_u1 = "corrected_slope_full",
      reporting_scale = reporting_scale
    ) %>%
      dplyr::filter(se_type == "naive") %>%
      dplyr::mutate(method = "blup_on_closed_form") %>%
      dplyr::select(method, -se_type),
    fit_observed_dual(
      stage2_df,
      outcome = "q_corrected_slope",
      predictor_u0 = "corrected_intercept_full",
      predictor_u1 = "corrected_slope_full",
      reporting_scale = reporting_scale
    ) %>%
      dplyr::filter(se_type == "naive") %>%
      dplyr::mutate(method = "closed_form_on_closed_form") %>%
      dplyr::select(method, -se_type),
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
    fit_mplus_dual_process(
      proc1_data = sim$lv1_y,
      proc2_data = sim$lv1_q,
      outcome1_var = "y", 
      outcome2_var = "q",
      cluster_id = "cid",
      time_index_var = "trial_index",
      time_value_var = "x",
      reporting_scale = reporting_scale
    )  %>%
      dplyr::mutate(method = "sem") %>%
      dplyr::select(method, dplyr::everything())
  )

  results <- add_study3_analysis_eligibility(results)

  dplyr::bind_cols(
    results,
    stage1_diag[rep(1L, nrow(results)), , drop = FALSE]
  ) %>%
    add_study_result_context(condition, truth)
}
