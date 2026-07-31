# Current estimator bundle applied to the Lai & Liu Study 1 DGM.

lai_study1_vh_methods <- function() {
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

lai_study1_vh_failure <- function(method, issue = "estimation_unavailable", detail = NA_character_) {
  tibble::tibble(
    method = method,
    estimate = NA_real_,
    se = NA_real_,
    ci_low = NA_real_,
    ci_high = NA_real_,
    status_code = NA_integer_,
    mx_issue_class = issue,
    mx_issue_detail = detail,
    analysis_eligible = FALSE,
    analysis_exclusion_reason = issue
  )
}

lai_study1_vh_as_method <- function(result, method) {
  result <- tibble::as_tibble(result)
  if (nrow(result) != 1L) {
    stop("Each Lai Study 1 VH estimator must return exactly one row.")
  }
  dplyr::mutate(result, method = method, .before = 1L)
}

lai_study1_vh_safe_method <- function(method, expression) {
  tryCatch(
    lai_study1_vh_as_method(expression, method),
    error = function(e) lai_study1_vh_failure(method, "estimation_error", conditionMessage(e))
  )
}

lai_study1_vh_fit_stage1 <- function(data) {
  safe_lmer(y ~ x + (x | cid), data = data)
}

#' Return an NA-stable template for replication-specific slope diagnostics.
#'
#' These quantities are all on the random-slope's raw units.  In particular,
#' `stage1_mean_lambda22` is retained as Lai's loading diagnostic, but is not
#' called reliability: with correlated random intercepts and slopes, the full
#' posterior covariance is required for a slope-reliability interpretation.
empty_lai_study1_vh_slope_diagnostics <- function() {
  tibble::tibble(
    stage1_realized_true_slope_sd = NA_real_,
    stage1_eb_slope_sd = NA_real_,
    stage1_corrected_slope_sd = NA_real_,
    stage1_fitted_slope_sd = NA_real_,
    stage1_fitted_residual_sd = NA_real_,
    stage1_mean_posterior_slope_variance = NA_real_,
    stage1_mean_lambda22 = NA_real_,
    stage1_fitted_posterior_slope_reliability = NA_real_,
    stage1_eb_to_population_slope_sd = NA_real_,
    stage1_corrected_to_population_slope_sd = NA_real_,
    stage1_fitted_to_population_slope_sd = NA_real_
  )
}

lai_study1_vh_sample_sd <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) < 2L) NA_real_ else stats::sd(x)
}

lai_study1_vh_finite_mean <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else mean(x)
}

#' Extract replication-specific slope-SD and posterior-reliability diagnostics.
#'
#' The fitted reliability is `mean_i(1 - V_post,i[2,2] / G_hat[2,2])`, using
#' exactly the posterior variances used to form the BLUPs.  It is therefore a
#' reliability diagnostic for the fitted Stage-1 model rather than the
#' `lambda22` entry of Lai's bivariate measurement matrix.
lai_study1_vh_slope_diagnostics <- function(fit_y, stage1_data, stage2_df,
                                             condition) {
  out <- empty_lai_study1_vh_slope_diagnostics()
  population_slope_sd <- sqrt(as.numeric(condition$var_u1[[1]]))
  if (!is.finite(population_slope_sd) || population_slope_sd <= 0) {
    return(out)
  }

  realized_true_slope_sd <- lai_study1_vh_sample_sd(stage2_df$true_u1)
  eb_slope_sd <- lai_study1_vh_sample_sd(stage2_df$u1_eb)
  corrected_slope_sd <- lai_study1_vh_sample_sd(stage2_df$corrected_slope_full)
  mean_postvar22 <- lai_study1_vh_finite_mean(stage2_df$postvar22)
  mean_lambda22 <- lai_study1_vh_finite_mean(stage2_df$lambda22)

  fitted <- tryCatch(
    extract_stage1_components(
      fit_obj = fit_y,
      data = stage1_data,
      cluster_var = "cid",
      within_var = "x"
    ),
    error = function(e) NULL
  )
  if (is.null(fitted) || nrow(fitted$G_hat) < 2L || ncol(fitted$G_hat) < 2L) {
    return(dplyr::mutate(
      out,
      stage1_realized_true_slope_sd = realized_true_slope_sd,
      stage1_eb_slope_sd = eb_slope_sd,
      stage1_corrected_slope_sd = corrected_slope_sd,
      stage1_mean_posterior_slope_variance = mean_postvar22,
      stage1_mean_lambda22 = mean_lambda22,
      stage1_eb_to_population_slope_sd = eb_slope_sd / population_slope_sd,
      stage1_corrected_to_population_slope_sd = corrected_slope_sd / population_slope_sd
    ))
  }

  fitted_slope_variance <- as.numeric(fitted$G_hat[2L, 2L])
  fitted_slope_sd <- if (is.finite(fitted_slope_variance) && fitted_slope_variance > 0) {
    sqrt(fitted_slope_variance)
  } else {
    NA_real_
  }
  residual_variances <- unlist(lapply(fitted$R_list, diag), use.names = FALSE)
  fitted_residual_sd <- if (length(residual_variances) > 0L &&
      all(is.finite(residual_variances)) && mean(residual_variances) >= 0) {
    sqrt(mean(residual_variances))
  } else {
    NA_real_
  }
  fitted_reliability <- if (is.finite(fitted_slope_variance) && fitted_slope_variance > 0 &&
      is.finite(mean_postvar22)) {
    1 - mean_postvar22 / fitted_slope_variance
  } else {
    NA_real_
  }

  dplyr::mutate(
    out,
    stage1_realized_true_slope_sd = realized_true_slope_sd,
    stage1_eb_slope_sd = eb_slope_sd,
    stage1_corrected_slope_sd = corrected_slope_sd,
    stage1_fitted_slope_sd = fitted_slope_sd,
    stage1_fitted_residual_sd = fitted_residual_sd,
    stage1_mean_posterior_slope_variance = mean_postvar22,
    stage1_mean_lambda22 = mean_lambda22,
    stage1_fitted_posterior_slope_reliability = fitted_reliability,
    stage1_eb_to_population_slope_sd = eb_slope_sd / population_slope_sd,
    stage1_corrected_to_population_slope_sd = corrected_slope_sd / population_slope_sd,
    stage1_fitted_to_population_slope_sd = fitted_slope_sd / population_slope_sd
  )
}

#' Add historical-Lai scale diagnostics before creating reporting views.
#'
#' Lai's naïve dual-EB model standardized its slope path by the observed EB
#' slope SD estimated under raw-data ML.  This equals the centered root mean
#' square (denominator `J`), whereas `stats::sd()` uses `J - 1`; retain both so
#' the historical multiplier and the familiar sample-SD diagnostic are clear.
#' Lai's 2S-PA and MSEM paths instead used each fitted model's latent slope SD.
add_lai_study1_vh_original_scale_diagnostics <- function(raw_results, stage2_df) {
  get_column <- function(name, default = NA_real_) {
    if (name %in% names(raw_results)) raw_results[[name]] else rep(default, nrow(raw_results))
  }
  eb_values <- as.numeric(stage2_df$u1_eb)
  eb_values <- eb_values[is.finite(eb_values)]
  naive_eb_mle_sd <- if (length(eb_values) >= 1L) {
    sqrt(mean((eb_values - mean(eb_values))^2))
  } else {
    NA_real_
  }
  naive_eb_sample_sd <- lai_study1_vh_sample_sd(eb_values)
  lai_2spa_slope_sd <- suppressWarnings(as.numeric(get_column("lai_fitted_latent_slope_sd")))
  msem_slope_sd <- suppressWarnings(as.numeric(get_column("mplus_fitted_latent_slope_sd")))

  dplyr::mutate(
    raw_results,
    lai_original_naive_eb_slope_sd = naive_eb_mle_sd,
    lai_original_naive_eb_slope_sample_sd = naive_eb_sample_sd,
    lai_original_method_multiplier = dplyr::case_when(
      method == "naive_dual_blup" ~ naive_eb_mle_sd,
      method == "lai_2spa" ~ lai_2spa_slope_sd,
      method == "msem" ~ msem_slope_sd,
      TRUE ~ NA_real_
    ),
    lai_original_multiplier_source = dplyr::case_when(
      method == "naive_dual_blup" ~ "observed_eb_mle_slope_sd",
      method == "lai_2spa" ~ "fitted_2spa_latent_slope_sd",
      method == "msem" ~ "fitted_msem_latent_slope_sd",
      TRUE ~ NA_character_
    ),
    lai_original_reporting_available = is.finite(lai_original_method_multiplier) &
      lai_original_method_multiplier > 0
  )
}

#' Add raw, common-latent-SD, and historical-Lai reporting views.
#'
#' Estimators are always fit and extracted on the raw coefficient scale.  This
#' function is the only place where reporting units are changed.  The `raw`
#' and `latent_sd` views use a common multiplier for every method; the
#' `lai_original_standardized` view is emitted only for methods with a direct
#' historical counterpart and preserves Lai's method-specific multiplier.
add_lai_study1_vh_reporting_scales <- function(raw_results, condition) {
  latent_sd <- sqrt(as.numeric(condition$var_u1[[1]]))
  raw_truth <- as.numeric(condition$beta_zu1[[1]])
  if (!is.finite(latent_sd) || latent_sd <= 0 || !is.finite(raw_truth)) {
    stop("Condition has an invalid Study 1 reporting scale or truth.")
  }

  raw_view <- dplyr::mutate(
    raw_results,
    reporting_scale = "raw",
    reporting_multiplier = 1,
    truth = raw_truth
  )
  latent_view <- dplyr::mutate(
    raw_results,
    reporting_scale = "latent_sd",
    reporting_multiplier = latent_sd,
    estimate = estimate * latent_sd,
    se = se * latent_sd,
    ci_low = ci_low * latent_sd,
    ci_high = ci_high * latent_sd,
    truth = raw_truth * latent_sd
  )
  lai_original_view <- raw_results %>%
    dplyr::filter(lai_original_reporting_available) %>%
    dplyr::mutate(
      reporting_scale = "lai_original_standardized",
      reporting_multiplier = lai_original_method_multiplier,
      estimate = estimate * reporting_multiplier,
      se = se * reporting_multiplier,
      ci_low = ci_low * reporting_multiplier,
      ci_high = ci_high * reporting_multiplier,
      # The original supplement used this common population benchmark even
      # though the naïve EB estimate used an observed-score multiplier.
      truth = raw_truth * latent_sd
    )
  dplyr::bind_rows(raw_view, latent_view, lai_original_view)
}

#' Record the historical Lai Study 1 convergence/retention convention.
#'
#' The original `sim1.R` records each method's optimizer `code` and retains a
#' replication in its converged summaries exactly when `code == 0`.  It does
#' not additionally screen standard errors, OpenMx information diagnostics,
#' Mplus warnings, or OLS design collinearity.  These fields intentionally
#' preserve that historical rule, even when it differs from VH eligibility.
add_lai_study1_original_eligibility <- function(results) {
  status_code <- suppressWarnings(as.integer(results$status_code))
  reason <- dplyr::case_when(
    is.na(status_code) ~ "status_code_missing",
    status_code != 0L ~ "status_code_nonzero",
    TRUE ~ NA_character_
  )
  dplyr::mutate(
    results,
    lai_original_eligible = is.na(reason),
    lai_original_exclusion_reason = reason
  )
}

#' Apply the exact VH Study 2 primary eligibility screen to the refreshed DGM.
#'
#' This is deliberately parallel to `add_study2_analysis_eligibility()`.  The
#' pre-existing `analysis_eligible` fields emitted by OLS are used as that
#' method's design diagnostic before this function overwrites the public VH
#' eligibility fields.
add_lai_study1_vh_analysis_eligibility <- function(results) {
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

  results <- add_lai_study1_original_eligibility(results)
  method <- as.character(get_column("method", NA_character_))
  status_code <- suppressWarnings(as.integer(get_column("status_code", NA_integer_)))
  estimate <- suppressWarnings(as.numeric(get_column("estimate", NA_real_)))
  se <- suppressWarnings(as.numeric(get_column("se", NA_real_)))
  dual_eligible <- as.logical(get_column("analysis_eligible", NA))
  dual_reason <- as.character(get_column("analysis_exclusion_reason", NA_character_))
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

  lai <- method == "lai_2spa"
  reason <- set_reason(reason, lai & !is.na(mx_issue_class) & mx_issue_class != "ok", "openmx_issue")
  reason <- set_reason(reason, lai & !is.na(mx_info_definite) & !mx_info_definite, "openmx_information_not_definite")
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
    vh_analysis_eligible = is.na(reason),
    vh_analysis_exclusion_reason = reason,
    # Maintain the existing public names as the VH-primary analysis rule.
    analysis_eligible = vh_analysis_eligible,
    analysis_exclusion_reason = vh_analysis_exclusion_reason,
    eligibility_comparison = dplyr::case_when(
      vh_analysis_eligible & lai_original_eligible ~ "both",
      lai_original_eligible & !vh_analysis_eligible ~ "lai_original_only",
      vh_analysis_eligible & !lai_original_eligible ~ "vh_only",
      TRUE ~ "neither"
    ),
    interval_eligible = vh_analysis_eligible & is.finite(se) & se > 0 &
      is.finite(ci_low) & is.finite(ci_high) & ci_low <= ci_high
  )
}

#' Fit the current estimator bundle on one Lai Study 1 replication.
#'
#' All method calls use reporting scale one (or raw Fuller output).  The caller
#' subsequently applies `add_lai_study1_vh_reporting_scales()`.
fit_lai_study1_vh_estimators <- function(condition, sim, methods = lai_study1_vh_methods()) {
  methods <- unique(as.character(methods))
  unknown <- setdiff(methods, lai_study1_vh_methods())
  if (length(unknown) > 0L) {
    stop("Unknown Lai Study 1 VH method(s): ", paste(unknown, collapse = ", "))
  }

  fit_y <- lai_study1_vh_fit_stage1(sim$lv1)
  if (is.null(fit_y)) {
    return(
      dplyr::bind_cols(
        dplyr::bind_rows(lapply(methods, lai_study1_vh_failure)),
        empty_stage1_diagnostics(),
        empty_lai_study1_vh_slope_diagnostics()
      ) %>%
        add_lai_study1_vh_analysis_eligibility() %>%
        add_lai_study1_vh_original_scale_diagnostics(tibble::tibble(u1_eb = numeric()))
    )
  }

  stage1_y <- tryCatch(
    get_stage1_eb_components(fit_y, sim$lv1, "cid", "y", "x"),
    error = function(e) NULL
  )
  corrected_y <- tryCatch(
    get_closed_form_corrected_scores(fit_y, sim$lv1, "cid", "y", "x"),
    error = function(e) NULL
  )
  if (is.null(stage1_y) || is.null(corrected_y)) {
    return(
      dplyr::bind_cols(
        dplyr::bind_rows(lapply(methods, lai_study1_vh_failure)),
        empty_stage1_diagnostics(),
        empty_lai_study1_vh_slope_diagnostics()
      ) %>%
        add_lai_study1_vh_analysis_eligibility() %>%
        add_lai_study1_vh_original_scale_diagnostics(tibble::tibble(u1_eb = numeric()))
    )
  }

  stage2_df <- sim$lv2_true %>%
    dplyr::left_join(stage1_y, by = "id") %>%
    dplyr::left_join(
      dplyr::select(
        corrected_y, id, corrected_intercept_full, corrected_slope_full,
        ols_var11, ols_var12, ols_var22
      ),
      by = "id"
    )
  stage1_diag <- get_stage1_diagnostics(fit_y, stage2_df)
  slope_diag <- lai_study1_vh_slope_diagnostics(
    fit_y = fit_y,
    stage1_data = sim$lv1,
    stage2_df = stage2_df,
    condition = condition
  )

  requested <- function(name) name %in% methods
  fitted <- list()

  if (requested("oracle_dual")) {
    fitted[["oracle_dual"]] <- lai_study1_vh_safe_method("oracle_dual", {
      fit_observed_dual(
        stage2_df, "z", "true_u0", "true_u1", reporting_scale = 1
      ) %>% dplyr::filter(se_type == "naive") %>% dplyr::select(-se_type)
    })
  }
  if (requested("naive_dual_blup")) {
    fitted[["naive_dual_blup"]] <- lai_study1_vh_safe_method("naive_dual_blup", {
      fit_observed_dual(stage2_df, "z", "u0_eb", "u1_eb", reporting_scale = 1) %>%
        dplyr::filter(se_type == "naive") %>% dplyr::select(-se_type)
    })
  }
  if (requested("closed_form_dual")) {
    fitted[["closed_form_dual"]] <- lai_study1_vh_safe_method("closed_form_dual", {
      fit_observed_dual(
        stage2_df, "z", "corrected_intercept_full", "corrected_slope_full",
        reporting_scale = 1
      ) %>% dplyr::filter(se_type == "naive") %>% dplyr::select(-se_type)
    })
  }
  fuller_args <- list(
    stage2_df = stage2_df, outcome = "z",
    predictor_u0 = "corrected_intercept_full",
    predictor_u1 = "corrected_slope_full",
    meas11 = "ols_var11", meas12 = "ols_var12", meas22 = "ols_var22"
  )
  if (requested("fuller_closed_form")) {
    fitted[["fuller_closed_form"]] <- lai_study1_vh_safe_method("fuller_closed_form", {
      do.call(fit_fuller_dual, fuller_args)
    })
  }
  if (requested("fuller_alpha_stepdown_closed_form")) {
    fitted[["fuller_alpha_stepdown_closed_form"]] <- lai_study1_vh_safe_method("fuller_alpha_stepdown_closed_form", {
      do.call(fit_fuller_dual_alpha_stepdown, fuller_args)
    })
  }
  if (requested("lai_2spa")) {
    fitted[["lai_2spa"]] <- lai_study1_vh_safe_method("lai_2spa", {
      fit_lai_2spa_observed_outcome(
        stage2_df,
        use_average = FALSE,
        # Match the original Study 1 and VH Study 2 convention: this is an
        # OpenMx start value for the free u0 -> z path, not a fixed parameter.
        u0_start = lai_study1_vh_fixed_params$beta_zu0,
        reporting_scale = 1
      )
    })
  }
  if (requested("msem")) {
    fitted[["msem"]] <- lai_study1_vh_safe_method("msem", {
      fit_mplus_blup_predictor(
        level1_data = sim$lv1,
        level2_data = sim$lv2_true,
        outcome_variable = "y",
        within_component = "x",
        between_component = "z",
        cluster_id = "cid",
        reporting_scale = 1
      )
    })
  }

  raw_results <- dplyr::bind_rows(fitted[methods]) %>%
    add_lai_study1_vh_analysis_eligibility() %>%
    add_lai_study1_vh_original_scale_diagnostics(stage2_df)
  dplyr::bind_cols(
    raw_results,
    stage1_diag[rep(1L, nrow(raw_results)), , drop = FALSE],
    slope_diag[rep(1L, nrow(raw_results)), , drop = FALSE]
  )
}
