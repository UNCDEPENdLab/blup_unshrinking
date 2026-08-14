#' Vig-Hallquist (2026) simulation study 2: BLUP as predictor
#' 
study2_methods <- function() {
  c(
    "oracle_dual",
    "naive_dual_blup",
    "closed_form_dual",
    "fuller_closed_form",
    "fuller_alpha_stepdown_closed_form",
    "lai_2spa",
    "msem"
  )
}

add_study2_method_roles <- function(results) {
  dplyr::mutate(
    results,
    method_role = dplyr::case_when(
      method == "oracle_dual" ~ "benchmark",
      grepl("_slope", method) ~ "slope_only_diagnostic",
      TRUE ~ "primary_dual"
    )
  )
}

#' Classify Study-2 rows for primary performance summaries.
#'
#' The raw estimate is never altered here. `analysis_eligible` records whether
#' prespecified method-specific diagnostics support including the row in bias,
#' variance, RMSE, or coverage summaries, while `analysis_exclusion_reason`
#' retains the first applicable reason for exclusion.
add_study2_analysis_eligibility <- function(results) {
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
  # fuller_guard_pass <- as.logical(get_column("fuller_auto_guard_pass", NA))
  # fuller_guard_reason <- as.character(get_column("fuller_auto_guard_reason", NA_character_))
  mx_issue_class <- as.character(get_column("mx_issue_class", NA_character_))
  mx_info_definite <- as.logical(get_column("mx_info_definite", NA))
  mx_condition_number <- suppressWarnings(as.numeric(get_column("mx_condition_number", NA_real_)))
  mplus_critical_warning <- as.logical(get_column("mplus_critical_warning", NA))
  mplus_target_parameter_count <- suppressWarnings(as.integer(get_column("mplus_target_parameter_count", NA_integer_)))

  reason <- rep(NA_character_, nrow(results))
  reason <- set_reason(reason, is.na(status_code) | is.na(estimate), "estimation_unavailable")
  reason <- set_reason(reason, !is.na(status_code) & status_code != 0L, "estimation_status_nonzero")
  reason <- set_reason(reason, !is.finite(estimate), "nonfinite_estimate")
  reason <- set_reason(reason, !is.finite(se) | se <= 0, "invalid_standard_error")

  dual_ols_methods <- c("oracle_dual", "naive_dual_blup", "closed_form_dual")
  dual_bad <- method %in% dual_ols_methods & !is.na(dual_eligible) & !dual_eligible
  reason <- set_reason(
    reason,
    dual_bad,
    ifelse(is.na(dual_reason), "stage2_design_ineligible", dual_reason)
  )

  # ZV: the auto guard was from the previous implementation of the stepdown so this incorrectly filters
  # all rows with alpha_stepdown as the method
  # alpha_fuller <- method == "fuller_alpha_stepdown_closed_form"
  # reason <- set_reason(
  #   reason,
  #   alpha_fuller & (is.na(fuller_guard_pass) | !fuller_guard_pass),
  #   ifelse(
  #     is.na(fuller_guard_reason) | fuller_guard_reason == "",
  #     "fuller_guard_failed",
  #     paste0("fuller_guard_", fuller_guard_reason)
  #   )
  # )

  lai <- method == "lai_2spa"
  reason <- set_reason(reason, lai & !is.na(mx_issue_class) & mx_issue_class != "ok", "openmx_issue")
  reason <- set_reason(reason, lai & !is.na(mx_info_definite) & !mx_info_definite, "openmx_information_not_definite")
  # OpenMx reports this diagnostic only for some fits; when available, 1e12 is
  # a deliberately conservative predeclared condition-number cap.
  reason <- set_reason(
    reason,
    lai & is.finite(mx_condition_number) & mx_condition_number > 1e12,
    "openmx_condition_number_excessive"
  )

  msem <- method == "msem"
  reason <- set_reason(
    reason,
    msem & !is.na(mplus_critical_warning) & mplus_critical_warning,
    "mplus_critical_warning"
  )
  reason <- set_reason(
    reason,
    msem & !is.na(mplus_target_parameter_count) & mplus_target_parameter_count != 1L,
    "mplus_target_parameter_not_unique"
  )

  dplyr::mutate(
    results,
    analysis_eligible = is.na(reason),
    analysis_exclusion_reason = reason
  )
}

simulate_study2 <- function(condition) {
  simulate_data_blup_as_predictor(condition)
}

run_study2_rep <- function(condition) {
  truth <- as.numeric(condition$standardized_beta_target[[1]])
  latent_slope_sd <- as.numeric(condition$tau1[[1]])
  sim <- simulate_study2(condition)
  fit_y <- fit_stage1(
    y ~ x,
    random = ~x | cid,
    data = sim$lv1,
    condition = condition,
    cluster_var = "cid"
  )
  if (is.null(fit_y)) {
    return(
      make_failed_result(condition, study2_methods(), truth) %>%
        dplyr::mutate(
          analysis_eligible = FALSE,
          analysis_exclusion_reason = "estimation_unavailable"
        ) %>%
        add_study2_method_roles()
    )
  }

  ordered_ids <- sim$lv2_true$id
  stage1_y <- tryCatch(
    get_stage1_eb_components(
      fit_obj = fit_y,
      data = sim$lv1,
      cluster_var = "cid",
      outcome_var = "y",
      within_var = "x"
    ),
    error = function(e) NULL
  )
  corrected_y <- tryCatch(
    get_closed_form_corrected_scores(
      fit_obj = fit_y,
      data = sim$lv1,
      cluster_var = "cid",
      outcome_var = "y",
      within_var = "x"
    ),
    error = function(e) NULL
  )
  if (is.null(stage1_y) || is.null(corrected_y)) {
    return(
      make_failed_result(condition, study2_methods(), truth) %>%
        dplyr::mutate(
          analysis_eligible = FALSE,
          analysis_exclusion_reason = "estimation_unavailable"
        ) %>%
        add_study2_method_roles()
    )
  }

  centered_dat <- dplyr::mutate(
    sim$lv1,
    yc = y - ave(y, cid),
    xc = x - ave(x, cid)
  )
  fit_centered <- fit_stage1(
    yc ~ 0 + xc,
    random = ~0 + xc | cid,
    data = centered_dat,
    condition = condition,
    cluster_var = "cid"
  )
  centered_u1 <- extract_centered_slope_eb(
    fit_centered,
    ordered_ids = ordered_ids
  )

  stage2_df <- sim$lv2_true %>%
    dplyr::left_join(stage1_y, by = "id") %>%
    dplyr::left_join(
      corrected_y %>%
        dplyr::select(
          id,
          corrected_intercept_full,
          corrected_slope_full,
          ols_var11,
          ols_var12,
          ols_var22
        ),
      by = "id"
    ) %>%
    dplyr::left_join(centered_u1, by = "id")
  stage1_diag <- get_stage1_diagnostics(fit_y, stage2_df)

  results <- dplyr::bind_rows(
    fit_observed_dual(
      stage2_df,
      outcome = "z",
      predictor_u0 = "true_u0",
      predictor_u1 = "true_u1",
      reporting_scale = latent_slope_sd
    ) %>%
      dplyr::filter(se_type == "naive") %>%
      dplyr::transmute(
        method = "oracle_dual",
        estimate, se, ci_low, ci_high, status_code
      ),
    fit_observed_dual(
      stage2_df,
      outcome = "z",
      predictor_u0 = "u0_eb",
      predictor_u1 = "u1_eb",
      reporting_scale = latent_slope_sd
    ) %>%
      dplyr::filter(se_type == "naive") %>%
      dplyr::transmute(
        method = "naive_dual_blup",
        estimate, se, ci_low, ci_high, status_code
      ),
    fit_observed_dual(
      stage2_df,
      outcome = "z",
      predictor_u0 = "corrected_intercept_full",
      predictor_u1 = "corrected_slope_full",
      reporting_scale = latent_slope_sd
    ) %>%
      dplyr::filter(se_type == "naive") %>%
      dplyr::transmute(
        method = "closed_form_dual",
        estimate, se, ci_low, ci_high, status_code
      ),
    fit_fuller_dual(
      stage2_df,
      outcome = "z",
      predictor_u0 = "corrected_intercept_full",
      predictor_u1 = "corrected_slope_full",
      meas11 = "ols_var11",
      meas12 = "ols_var12",
      meas22 = "ols_var22"
    ) %>%
      rescale_fuller_to_population_sd(latent_slope_sd) %>%
      dplyr::mutate(method = "fuller_closed_form") %>%
      dplyr::select(method, dplyr::everything()),
    fit_fuller_dual_alpha_stepdown(
      stage2_df,
      outcome = "z",
      predictor_u0 = "corrected_intercept_full",
      predictor_u1 = "corrected_slope_full",
      meas11 = "ols_var11",
      meas12 = "ols_var12",
      meas22 = "ols_var22"
    ) %>%
      rescale_fuller_to_population_sd(latent_slope_sd) %>%
      dplyr::mutate(method = "fuller_alpha_stepdown_closed_form") %>%
      dplyr::select(method, dplyr::everything()),
    fit_lai_2spa_observed_outcome(
      stage2_df,
      use_average = FALSE,
      u0_start = condition$beta1z[[1]],
      reporting_scale = latent_slope_sd
    ) %>%
      dplyr::mutate(method = "lai_2spa") %>%
      dplyr::select(method, dplyr::everything()),
    fit_mplus_blup_predictor(
      level1_data = sim$lv1,
      level2_data = sim$lv2_true,
      outcome_variable = "y",
      within_component = "x",
      between_component = "z",
      cluster_id = "cid",
      reporting_scale = latent_slope_sd
    ) %>%
      dplyr::mutate(method = "msem") %>%
      dplyr::select(method, dplyr::everything())
  )

  results <- results %>%
    add_stage1_estimates(
      fit_obj = fit_y,
      data = sim$lv1,
      cluster_var = "cid",
      within_var = "x"
    ) %>%
    add_study2_analysis_eligibility() %>%
    add_study2_method_roles()

  dplyr::bind_cols(
    results,
    stage1_diag[rep(1L, nrow(results)), , drop = FALSE]
  ) %>%
    add_study_result_context(condition, truth)
}
