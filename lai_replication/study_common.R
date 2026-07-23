#' Shared Lai study machinery that is not part of the general score helpers.

lai_truth <- function(condition) {
  condition$beta_zu1 * sqrt(condition$var_u1)
}

make_covu <- function(condition) {
  matrix(
    c(
      condition$icc,
      condition$cor_u0_u1 * sqrt(condition$icc * condition$var_u1),
      condition$cor_u0_u1 * sqrt(condition$icc * condition$var_u1),
      condition$var_u1
    ),
    nrow = 2L,
    byrow = TRUE
  )
}

tempered_eiv_methods <- function() {
  unlist(lapply(
    c(
      "tempered_eiv_dual_corrected_l25",
      "tempered_eiv_dual_corrected_l50",
      "tempered_eiv_dual_corrected_l75"
    ),
    eiv_se_variant_methods
  ), use.names = FALSE)
}

eiv_se_variant_methods <- function(base_method) {
  c(base_method, paste0(base_method, "_hc0"), paste0(base_method, "_hc3"))
}

condition_includes_tempered_eiv <- function(condition) {
  if (!("include_tempered_eiv" %in% names(condition))) {
    return(FALSE)
  }
  value <- condition$include_tempered_eiv[[1]]
  if (is.logical(value)) {
    return(isTRUE(value))
  }
  if (is.numeric(value)) {
    return(is.finite(value) && value != 0)
  }
  tolower(trimws(as.character(value))) %in% c("1", "true", "t", "yes", "y")
}

lai_condition_to_r_spec <- function(condition) {
  r_structure <- if ("r_structure" %in% names(condition)) {
    as.character(condition$r_structure[[1]])
  } else if ("residual_structure" %in% names(condition)) {
    as.character(condition$residual_structure[[1]])
  } else {
    "iid"
  }

  switch(
    tolower(r_structure),
    iid = list(structure = "iid"),
    ar1 = {
      rho <- if ("r_rho" %in% names(condition)) condition$r_rho[[1]] else condition$residual_rho[[1]]
      list(structure = "ar1", rho = as.numeric(rho))
    },
    toeplitz = {
      correlations <- if ("r_correlations" %in% names(condition)) condition$r_correlations[[1]] else condition$residual_correlations[[1]]
      list(structure = "toeplitz", correlations = as.numeric(correlations))
    },
    stop("Unsupported residual structure for Lai replication: ", r_structure)
  )
}

lai_condition_uses_non_iid_R <- function(condition) {
  !identical(normalize_r_spec(lai_condition_to_r_spec(condition))$structure, "iid")
}

lai_condition_to_nlme_correlation <- function(condition, cluster_var, index_var = "trial_index") {
  r_spec <- normalize_r_spec(lai_condition_to_r_spec(condition))
  switch(
    r_spec$structure,
    iid = NULL,
    ar1 = {
      if (!requireNamespace("nlme", quietly = TRUE)) {
        stop("The `nlme` package is required for AR(1) Stage-1 residual covariance fits.")
      }
      nlme::corAR1(form = stats::as.formula(sprintf("~%s | %s", index_var, cluster_var)))
    },
    stop("No nlme residual-correlation adapter is implemented for `", r_spec$structure, "`.")
  )
}

add_lai_trial_index <- function(data, cluster_var, index_var = "trial_index") {
  cluster_ids <- as.character(data[[cluster_var]])
  data[[index_var]] <- ave(seq_len(nrow(data)), cluster_ids, FUN = seq_along)
  data
}

draw_lai_level1_residuals <- function(cluster_sizes, sigma, condition) {
  r_spec <- lai_condition_to_r_spec(condition)
  unlist(lapply(cluster_sizes, function(n_i) {
    draw_residuals_from_R(make_R_matrix(n_i, sigma = sigma, r_spec = r_spec))
  }), use.names = FALSE)
}

fit_lai_stage1 <- function(fixed, random, data, condition, cluster_var) {
  if (lai_condition_uses_non_iid_R(condition)) {
    if (!requireNamespace("nlme", quietly = TRUE)) {
      return(NULL)
    }
    return(safe_lme(
      fixed = fixed,
      random = random,
      data = data,
      correlation = lai_condition_to_nlme_correlation(condition, cluster_var = cluster_var),
      method = "REML",
      control = nlme::lmeControl(returnObject = TRUE, msMaxIter = 100L, opt = "optim")
    ))
  }

  fixed_txt <- paste(deparse(fixed), collapse = " ")
  random_txt <- sub("^~", "", paste(deparse(random), collapse = " "))
  safe_lmer(stats::as.formula(paste0(fixed_txt, " + (", random_txt, ")")), data = data)
}

extract_centered_slope_eb <- function(fit_obj, ordered_ids) {
  if (is.null(fit_obj)) {
    return(tibble::tibble(id = ordered_ids, centered_u1_eb = NA_real_))
  }

  re_df <- tryCatch({
    if (inherits(fit_obj, "lme")) {
      nlme::ranef(fit_obj)
    } else {
      lme4::ranef(fit_obj)[[1]]
    }
  }, error = function(e) NULL)

  if (is.null(re_df)) {
    return(tibble::tibble(id = ordered_ids, centered_u1_eb = NA_real_))
  }

  tibble::tibble(id = rownames(re_df), centered_u1_eb = as.numeric(re_df[[ncol(re_df)]]))
}

matched_study_methods <- function(include_tempered_eiv = FALSE) {
  methods <- c(
    "oracle_dual",
    "naive_slope_only", "naive_slope_only_hc3",
    "centered_slope_only", "centered_slope_only_hc3",
    "naive_dual_eb", "naive_dual_eb_hc3",
    eiv_se_variant_methods("eiv_dual_corrected"),
    eiv_se_variant_methods("eiv_dual_corrected_nearpd"),
    eiv_se_variant_methods("eiv_dual_corrected_ridge"),
    "ridge_dual_eb",
    "corrected_slope_only", "corrected_slope_only_hc3",
    "corrected_dual", "corrected_dual_hc3",
    "lai_2spa", "lai_2spaa", "fuller", "fuller_stepdown", "fuller_alpha_stepdown"
  )
  if (isTRUE(include_tempered_eiv)) {
    methods <- append(methods, tempered_eiv_methods(), after = match("eiv_dual_corrected_ridge", methods))
  }
  methods
}

disparate_study_methods <- function(include_tempered_eiv = FALSE) {
  methods <- c(
    "oracle_dual",
    "naive_dual_eb", "naive_dual_eb_hc3",
    eiv_se_variant_methods("eiv_dual_corrected"),
    eiv_se_variant_methods("eiv_dual_corrected_nearpd"),
    eiv_se_variant_methods("eiv_dual_corrected_ridge"),
    "ridge_dual_eb",
    "corrected_dual", "corrected_dual_hc3",
    "lai_2spa", "lai_2spaa", "fuller", "fuller_stepdown", "fuller_alpha_stepdown"
  )
  if (isTRUE(include_tempered_eiv)) {
    methods <- append(methods, tempered_eiv_methods(), after = match("eiv_dual_corrected_ridge", methods))
  }
  methods
}

study_methods_for_condition <- function(condition) {
  include_tempered_eiv <- condition_includes_tempered_eiv(condition)
  study_key <- as.character(condition$study[[1]])
  switch(
    study_key,
    study1 = matched_study_methods(include_tempered_eiv),
    study2 = matched_study_methods(include_tempered_eiv),
    study3 = disparate_study_methods(include_tempered_eiv),
    stop("Unsupported Lai study key: ", study_key)
  )
}

#' Classify matched-outcome replication results for performance summaries.
#'
#' This pathway is used by Lai Studies 1, 2, and 4.  Dual-predictor OLS fits
#' retain their raw estimates, but a rank-deficient or near-collinear Stage-2
#' design is not eligible for the primary bias and RMSE summaries.  Interval
#' eligibility is stricter because coverage additionally needs a valid standard
#' error and confidence interval.
add_matched_outcome_analysis_eligibility <- function(results) {
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
  ci_low <- suppressWarnings(as.numeric(get_column("ci_low", NA_real_)))
  ci_high <- suppressWarnings(as.numeric(get_column("ci_high", NA_real_)))
  dual_eligible <- as.logical(get_column("analysis_eligible", NA))
  dual_reason <- as.character(get_column("analysis_exclusion_reason", NA_character_))
  fuller_guard_pass <- as.logical(get_column("fuller_auto_guard_pass", NA))
  fuller_guard_reason <- as.character(get_column("fuller_auto_guard_reason", NA_character_))
  mx_issue_class <- as.character(get_column("mx_issue_class", NA_character_))
  mx_info_definite <- as.logical(get_column("mx_info_definite", NA))
  mx_condition_number <- suppressWarnings(as.numeric(get_column("mx_condition_number", NA_real_)))

  point_reason <- rep(NA_character_, nrow(results))
  point_reason <- set_reason(point_reason, is.na(status_code) | is.na(estimate), "estimation_unavailable")
  point_reason <- set_reason(point_reason, !is.na(status_code) & status_code != 0L, "estimation_status_nonzero")
  point_reason <- set_reason(point_reason, !is.finite(estimate), "nonfinite_estimate")

  dual_ols_methods <- c(
    "oracle_dual", "naive_dual_eb", "naive_dual_eb_hc3",
    "corrected_dual", "corrected_dual_hc3"
  )
  dual_bad <- method %in% dual_ols_methods & !is.na(dual_eligible) & !dual_eligible
  point_reason <- set_reason(
    point_reason,
    dual_bad,
    ifelse(is.na(dual_reason), "stage2_design_ineligible", dual_reason)
  )

  alpha_fuller <- method == "fuller_alpha_stepdown"
  point_reason <- set_reason(
    point_reason,
    alpha_fuller & (is.na(fuller_guard_pass) | !fuller_guard_pass),
    ifelse(
      is.na(fuller_guard_reason) | fuller_guard_reason == "",
      "fuller_guard_failed",
      paste0("fuller_guard_", fuller_guard_reason)
    )
  )

  lai <- method %in% c("lai_2spa", "lai_2spaa")
  point_reason <- set_reason(point_reason, lai & !is.na(mx_issue_class) & mx_issue_class != "ok", "openmx_issue")
  point_reason <- set_reason(point_reason, lai & !is.na(mx_info_definite) & !mx_info_definite, "openmx_information_not_definite")
  point_reason <- set_reason(
    point_reason,
    lai & is.finite(mx_condition_number) & mx_condition_number > 1e12,
    "openmx_condition_number_excessive"
  )

  interval_reason <- point_reason
  interval_reason <- set_reason(interval_reason, !is.finite(se) | se <= 0, "invalid_standard_error")
  interval_reason <- set_reason(
    interval_reason,
    !is.finite(ci_low) | !is.finite(ci_high) | ci_low > ci_high,
    "invalid_confidence_interval"
  )

  dplyr::mutate(
    results,
    analysis_eligible = is.na(point_reason),
    analysis_exclusion_reason = point_reason,
    interval_eligible = is.na(interval_reason),
    interval_exclusion_reason = interval_reason
  )
}

make_failed_result <- function(condition, methods, truth) {
  dplyr::bind_cols(
    tibble::tibble(
      study = condition$study,
      method = methods,
      estimate = NA_real_,
      se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      status_code = NA_integer_,
      truth = truth
    ),
    empty_stage1_diagnostics()[rep(1L, length(methods)), , drop = FALSE]
  )
}

add_study_result_context <- function(results, condition, truth) {
  dplyr::mutate(results, study = condition$study, truth = truth)
}

fit_tempered_eiv_dual_set <- function(stage2_df,
                                      outcome,
                                      predictor_u0,
                                      predictor_u1,
                                      meas11,
                                      meas12,
                                      meas22,
                                      outcome_meas_var = NULL,
                                      weights = c(0.25, 0.50, 0.75)) {
  purrr::map_dfr(weights, function(weight) {
    suffix <- sprintf("l%02d", as.integer(round(weight * 100)))
    fit_eiv_dual(
      stage2_df,
      outcome = outcome,
      predictor_u0 = predictor_u0,
      predictor_u1 = predictor_u1,
      meas11 = meas11,
      meas12 = meas12,
      meas22 = meas22,
      outcome_meas_var = outcome_meas_var,
      measurement_weight = weight
    ) %>%
      finalize_eiv_se_variants(paste0("tempered_eiv_dual_corrected_", suffix))
  })
}

run_matched_outcome_rep <- function(condition, sim) {
  truth <- lai_truth(condition)
  include_tempered_eiv <- condition_includes_tempered_eiv(condition)
  sim$lv1 <- add_lai_trial_index(sim$lv1, cluster_var = "cid")
  fit_y <- fit_lai_stage1(y ~ x, random = ~x | cid, data = sim$lv1, condition = condition, cluster_var = "cid")
  if (is.null(fit_y)) {
    return(
      make_failed_result(condition, matched_study_methods(include_tempered_eiv), truth) %>%
        add_matched_outcome_analysis_eligibility()
    )
  }

  ordered_ids <- sim$lv2_true$id
  stage1_y <- get_stage1_eb_components(
    fit_obj = fit_y,
    data = sim$lv1,
    cluster_var = "cid",
    outcome_var = "y",
    within_var = "x"
  )
  corrected_y <- get_closed_form_corrected_scores(
    fit_obj = fit_y,
    data = sim$lv1,
    cluster_var = "cid",
    outcome_var = "y",
    within_var = "x"
  )

  yc <- sim$lv1$y - ave(sim$lv1$y, sim$lv1$cid)
  xc <- sim$lv1$x - ave(sim$lv1$x, sim$lv1$cid)
  centered_dat <- dplyr::mutate(sim$lv1, yc = yc, xc = xc)
  fit_centered <- fit_lai_stage1(yc ~ 0 + xc, random = ~0 + xc | cid, data = centered_dat, condition = condition, cluster_var = "cid")
  centered_u1 <- extract_centered_slope_eb(fit_centered, ordered_ids = ordered_ids)

  stage2_df <- sim$lv2_true %>%
    dplyr::left_join(stage1_y, by = "id") %>%
    dplyr::left_join(corrected_y %>% dplyr::select(id, corrected_intercept_full, corrected_slope_full, ols_var11, ols_var12, ols_var22), by = "id") %>%
    dplyr::left_join(centered_u1, by = "id")
  stage1_diag <- get_stage1_diagnostics(fit_y, stage2_df)

  results <- dplyr::bind_rows(
    fit_observed_dual(stage2_df, outcome = "z", predictor_u0 = "true_u0", predictor_u1 = "true_u1") %>%
      dplyr::filter(se_type == "naive") %>%
      dplyr::mutate(method = "oracle_dual") %>%
      dplyr::select(method, -se_type),
    finalize_ols_se_variants(fit_observed_single(stage2_df, outcome = "z", predictor = "u1_eb"), "naive_slope_only"),
    finalize_ols_se_variants(fit_observed_single(stage2_df, outcome = "z", predictor = "centered_u1_eb"), "centered_slope_only"),
    finalize_ols_se_variants(fit_observed_dual(stage2_df, outcome = "z", predictor_u0 = "u0_eb", predictor_u1 = "u1_eb"), "naive_dual_eb"),
    # IMPORTANT: The Errors-in-Variables (EIV) approach is mathematically formulated to correct
    # the variance of unshrunken estimates by subtracting their sampling error variance. 
    # Because BLUPs are already shrunken toward the population mean (meaning Var(BLUP) < Var(True)),
    # subtracting measurement error from their variance artificially deflates it further, inflating
    # the Stage-2 coefficients to incorrect levels. Therefore, we use the unshrunken OLS likelihood 
    # scores (corrected_intercept_full, corrected_slope_full) and their corresponding sampling 
    # variances (ols_var11, ols_var12, ols_var22) rather than BLUPs (u0_eb, u1_eb) and posterior variances.
    fit_eiv_dual(
      stage2_df,
      outcome = "z",
      predictor_u0 = "corrected_intercept_full",
      predictor_u1 = "corrected_slope_full",
      meas11 = "ols_var11",
      meas12 = "ols_var12",
      meas22 = "ols_var22"
    ) %>%
      finalize_eiv_se_variants("eiv_dual_corrected"),
    fit_eiv_dual(
      stage2_df,
      outcome = "z",
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
      outcome = "z",
      predictor_u0 = "corrected_intercept_full",
      predictor_u1 = "corrected_slope_full",
      meas11 = "ols_var11",
      meas12 = "ols_var12",
      meas22 = "ols_var22",
      ridge_predictor_block = TRUE
    ) %>%
      finalize_eiv_se_variants("eiv_dual_corrected_ridge"),
    if (isTRUE(include_tempered_eiv)) {
      fit_tempered_eiv_dual_set(
        stage2_df,
        outcome = "z",
        predictor_u0 = "corrected_intercept_full",
        predictor_u1 = "corrected_slope_full",
        meas11 = "ols_var11",
        meas12 = "ols_var12",
        meas22 = "ols_var22"
      )
    },
    fit_ridge_dual(stage2_df, outcome = "z", predictor_u0 = "u0_eb", predictor_u1 = "u1_eb") %>%
      dplyr::mutate(method = "ridge_dual_eb") %>%
      dplyr::select(method, dplyr::everything()),
    finalize_ols_se_variants(fit_observed_single(stage2_df, outcome = "z", predictor = "corrected_slope_full"), "corrected_slope_only"),
    finalize_ols_se_variants(fit_observed_dual(stage2_df, outcome = "z", predictor_u0 = "corrected_intercept_full", predictor_u1 = "corrected_slope_full"), "corrected_dual"),
    fit_lai_2spa_observed_outcome(stage2_df, use_average = FALSE) %>%
      dplyr::mutate(method = "lai_2spa") %>%
      dplyr::select(method, dplyr::everything()),
    fit_lai_2spa_observed_outcome(stage2_df, use_average = TRUE) %>%
      dplyr::mutate(method = "lai_2spaa") %>%
      dplyr::select(method, dplyr::everything()),
    fit_fuller_dual(
      stage2_df,
      outcome = "z",
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
      outcome = "z",
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
      outcome = "z",
      predictor_u0 = "corrected_intercept_full",
      predictor_u1 = "corrected_slope_full",
      meas11 = "ols_var11",
      meas12 = "ols_var12",
      meas22 = "ols_var22"
    ) %>%
      dplyr::mutate(method = "fuller_alpha_stepdown") %>%
      dplyr::select(method, dplyr::everything())
  )

  results <- add_matched_outcome_analysis_eligibility(results)

  dplyr::bind_cols(results, stage1_diag[rep(1L, nrow(results)), , drop = FALSE]) %>%
    add_study_result_context(condition, truth)
}
