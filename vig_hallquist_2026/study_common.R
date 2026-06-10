#' Shared study machinery that is not part of the general score helpers.

vig_hallquist_truth <- function(condition) {
  condition$beta_zu1 * sqrt(condition$var_u1)
}

# TODO: do we have to use ICC here?
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

# tempered_eiv_methods <- function() {
#   unlist(lapply(
#     c(
#       "tempered_eiv_dual_corrected_l25",
#       "tempered_eiv_dual_corrected_l50",
#       "tempered_eiv_dual_corrected_l75"
#     ),
#     eiv_se_variant_methods
#   ), use.names = FALSE)
# }

# eiv_se_variant_methods <- function(base_method) {
#   c(base_method, paste0(base_method, "_hc0"), paste0(base_method, "_hc3"))
# }

# condition_includes_tempered_eiv <- function(condition) {
#   if (!("include_tempered_eiv" %in% names(condition))) {
#     return(FALSE)
#   }
#   value <- condition$include_tempered_eiv[[1]]
#   if (is.logical(value)) {
#     return(isTRUE(value))
#   }
#   if (is.numeric(value)) {
#     return(is.finite(value) && value != 0)
#   }
#   tolower(trimws(as.character(value))) %in% c("1", "true", "t", "yes", "y")
# }

# lai_condition_to_r_spec <- function(condition) {
#   r_structure <- if ("r_structure" %in% names(condition)) {
#     as.character(condition$r_structure[[1]])
#   } else if ("residual_structure" %in% names(condition)) {
#     as.character(condition$residual_structure[[1]])
#   } else {
#     "iid"
#   }

#   switch(
#     tolower(r_structure),
#     iid = list(structure = "iid"),
#     ar1 = {
#       rho <- if ("r_rho" %in% names(condition)) condition$r_rho[[1]] else condition$residual_rho[[1]]
#       list(structure = "ar1", rho = as.numeric(rho))
#     },
#     toeplitz = {
#       correlations <- if ("r_correlations" %in% names(condition)) condition$r_correlations[[1]] else condition$residual_correlations[[1]]
#       list(structure = "toeplitz", correlations = as.numeric(correlations))
#     },
#     stop("Unsupported residual structure for Lai replication: ", r_structure)
#   )
# }

# lai_condition_uses_non_iid_R <- function(condition) {
#   !identical(normalize_r_spec(lai_condition_to_r_spec(condition))$structure, "iid")
# }

# lai_condition_to_nlme_correlation <- function(condition, cluster_var, index_var = "trial_index") {
#   r_spec <- normalize_r_spec(lai_condition_to_r_spec(condition))
#   switch(
#     r_spec$structure,
#     iid = NULL,
#     ar1 = {
#       if (!requireNamespace("nlme", quietly = TRUE)) {
#         stop("The `nlme` package is required for AR(1) Stage-1 residual covariance fits.")
#       }
#       nlme::corAR1(form = stats::as.formula(sprintf("~%s | %s", index_var, cluster_var)))
#     },
#     stop("No nlme residual-correlation adapter is implemented for `", r_spec$structure, "`.")
#   )
# }

add_trial_index <- function(data, cluster_var, index_var = "trial_index") {
  cluster_ids <- as.character(data[[cluster_var]])
  data[[index_var]] <- ave(seq_len(nrow(data)), cluster_ids, FUN = seq_along)
  data
}

draw_level1_residuals <- function(cluster_sizes, sigma, condition) {
  r_spec <- list(structure = "iid")
  unlist(lapply(cluster_sizes, function(n_i) {
    draw_residuals_from_R(make_R_matrix(n_i, sigma = sigma, r_spec = r_spec))
  }), use.names = FALSE)
}

fit_stage1 <- function(fixed, random, data, condition, cluster_var) {
  # if (lai_condition_uses_non_iid_R(condition)) {
  #   if (!requireNamespace("nlme", quietly = TRUE)) {
  #     return(NULL)
  #   }
  #   return(safe_lme(
  #     fixed = fixed,
  #     random = random,
  #     data = data,
  #     correlation = lai_condition_to_nlme_correlation(condition, cluster_var = cluster_var),
  #     method = "REML",
  #     control = nlme::lmeControl(returnObject = TRUE, msMaxIter = 100L, opt = "optim")
  #   ))
  # }

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



# matched_study_methods <- function() {
#   # TODO: update these methods
#   c(
#     "oracle_dual",
#     "naive_slope_only", "naive_slope_only_hc3",
#     "centered_slope_only", "centered_slope_only_hc3",
#     "naive_dual_eb", "naive_dual_eb_hc3",
#     "ridge_dual_eb",
#     "corrected_slope_only", "corrected_slope_only_hc3",
#     "corrected_dual", "corrected_dual_hc3",
#     "lai_2spa", "lai_2spaa", "fuller", "fuller_stepdown", "fuller_alpha_stepdown"
#   )
# }

# disparate_study_methods <- function() {
#   # TODO: update these methods
#   c(
#     "oracle_dual",
#     "naive_dual_eb", "naive_dual_eb_hc3",
#     "ridge_dual_eb",
#     "corrected_dual", "corrected_dual_hc3",
#     "lai_2spa", "lai_2spaa", "fuller", "fuller_stepdown", "fuller_alpha_stepdown"
#   )
# }

study_methods_for_condition <- function(condition) {
  study_key <- as.character(condition$study[[1]])
  switch(
    study_key,
    study1 = study1_methods(),
    study2 = study2_methods(),
    study3 = study3_methods(),
    study4 = study4_methods(),
    stop("Unsupported study key: ", study_key)
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

# fit_tempered_eiv_dual_set <- function(stage2_df,
#                                       outcome,
#                                       predictor_u0,
#                                       predictor_u1,
#                                       meas11,
#                                       meas12,
#                                       meas22,
#                                       outcome_meas_var = NULL,
#                                       weights = c(0.25, 0.50, 0.75)) {
#   purrr::map_dfr(weights, function(weight) {
#     suffix <- sprintf("l%02d", as.integer(round(weight * 100)))
#     fit_eiv_dual(
#       stage2_df,
#       outcome = outcome,
#       predictor_u0 = predictor_u0,
#       predictor_u1 = predictor_u1,
#       meas11 = meas11,
#       meas12 = meas12,
#       meas22 = meas22,
#       outcome_meas_var = outcome_meas_var,
#       measurement_weight = weight
#     ) %>%
#       finalize_eiv_se_variants(paste0("tempered_eiv_dual_corrected_", suffix))
#   })
# }

run_matched_outcome_rep <- function(condition, sim) {
  truth <- vig_hallquist_truth(condition)
  sim$lv1 <- add_trial_index(sim$lv1, cluster_var = "cid")
  fit_y <- fit_stage1(y ~ x, random = ~x | cid, data = sim$lv1, condition = condition, cluster_var = "cid")
  if (is.null(fit_y)) {
    return(make_failed_result(condition, matched_study_methods(), truth))
  }

  ordered_ids <- sim$lv2_true$id
  # TODO: how are scoring failures handled in stage 2
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
  fit_centered <- fit_stage1(yc ~ 0 + xc, random = ~0 + xc | cid, data = centered_dat, condition = condition, cluster_var = "cid")
  centered_u1 <- extract_centered_slope_eb(fit_centered, ordered_ids = ordered_ids)

  stage2_df <- sim$lv2_true %>%
    dplyr::left_join(stage1_y, by = "id") %>%
    dplyr::left_join(corrected_y %>% dplyr::select(id, corrected_intercept_full, corrected_slope_full, ols_var11, ols_var12, ols_var22), by = "id") %>%
    dplyr::left_join(centered_u1, by = "id")
  stage1_diag <- get_stage1_diagnostics(fit_y, stage2_df)

  # TODO: update these methods
  results <- dplyr::bind_rows(
    fit_observed_dual(stage2_df, outcome = "z", predictor_u0 = "true_u0", predictor_u1 = "true_u1") %>%
      dplyr::filter(se_type == "naive") %>%
      dplyr::transmute(method = "oracle_dual", estimate, se, ci_low, ci_high, status_code),
    finalize_ols_se_variants(fit_observed_single(stage2_df, outcome = "z", predictor = "u1_eb"), "naive_slope_only"),
    finalize_ols_se_variants(fit_observed_single(stage2_df, outcome = "z", predictor = "centered_u1_eb"), "centered_slope_only"),
    finalize_ols_se_variants(fit_observed_dual(stage2_df, outcome = "z", predictor_u0 = "u0_eb", predictor_u1 = "u1_eb"), "naive_dual_eb"),
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
