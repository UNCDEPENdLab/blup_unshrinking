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
  fit_y <- safe_lmer(y ~ x + (x | cid), data = sim$lv1)
  if (is.null(fit_y)) {
    return(make_failed_result(condition, matched_study_methods(include_tempered_eiv), truth))
  }

  split_y <- split(sim$lv1, sim$lv1$cid_chr)
  ordered_ids <- sim$lv2_true$id
  eb_inputs <- compute_bivariate_eb_inputs(fit_y, split_y, ordered_ids, within_var = "x")
  corrected_y <- get_closed_form_corrected_scores(
    fit_obj = fit_y,
    data = sim$lv1,
    cluster_var = "cid",
    outcome_var = "y",
    within_var = "x"
  )

  yc <- sim$lv1$y - ave(sim$lv1$y, sim$lv1$cid)
  xc <- sim$lv1$x - ave(sim$lv1$x, sim$lv1$cid)
  fit_centered <- safe_lmer(yc ~ 0 + xc + (0 + xc | cid), data = dplyr::mutate(sim$lv1, yc = yc, xc = xc))
  centered_u1 <- if (is.null(fit_centered)) {
    tibble::tibble(id = ordered_ids, centered_u1_eb = NA_real_)
  } else {
    centered_re <- lme4::ranef(fit_centered)[[1]]
    tibble::tibble(id = rownames(centered_re), centered_u1_eb = centered_re[[1]])
  }

  stage2_df <- sim$lv2_true %>%
    dplyr::left_join(eb_inputs, by = "id") %>%
    dplyr::left_join(corrected_y %>% dplyr::select(id, corrected_intercept_full, corrected_slope_full, ols_var11, ols_var12, ols_var22), by = "id") %>%
    dplyr::left_join(centered_u1, by = "id")
  stage1_diag <- get_stage1_diagnostics(fit_y, stage2_df)

  results <- dplyr::bind_rows(
    fit_observed_dual(stage2_df, outcome = "z", predictor_u0 = "true_u0", predictor_u1 = "true_u1") %>%
      dplyr::filter(se_type == "naive") %>%
      dplyr::transmute(method = "oracle_dual", estimate, se, ci_low, ci_high, status_code),
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
      dplyr::select(method, estimate, se, ci_low, ci_high, status_code),
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

  dplyr::bind_cols(results, stage1_diag[rep(1L, nrow(results)), , drop = FALSE]) %>%
    add_study_result_context(condition, truth)
}
