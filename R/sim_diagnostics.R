#' Diagnostics for first-stage mixed-model fits used in stage-2 score analyses.

#' Return an empty stage-1 diagnostics template.
#'
#' @details
#' When upstream model fitting fails, this helper returns a row filled with
#' missing values. This keeps downstream calls to `dplyr::bind_rows()` stable
#' by ensuring the expected columns are always present.
#'
#' @return A `tibble` with a single row of NA values for diagnostic columns.
empty_stage1_diagnostics <- function() {
  tibble::tibble(
    stage1_singular_problem = NA,
    stage1_problem_detail = NA_character_,
    stage1_lmer_singular = NA,
    stage1_re_corr = NA_real_,
    stage1_eb_corr = NA_real_,
    stage1_design_kappa = NA_real_
  )
}


#' Extract diagnostics from a first-stage lme4 fit.
#'
#' @details
#' This function extracts random-effect correlations and checks for boundary
#' or near-singular estimates that could cause numerical issues in stage-2
#' score estimators. It flags fits as problematic if `lme4::isSingular` is true,
#' or if the estimated random-effect correlation or empirical Bayes correlation
#' is perfectly colinear (absolute value > 0.999).
#'
#' @param fit_obj The `lmerMod` fitted object from `lme4`. Can be `NULL` if
#' the fit failed.
#' @param stage2_df A data frame containing the cluster-level EB scores or
#' predictors used in the stage-2 model.
#' @param predictor_u0 The column name in `stage2_df` for the intercept-like score.
#' @param predictor_u1 The column name in `stage2_df` for the slope-like score.
#'
#' @return A `tibble` with 1 row containing summary diagnostics for the fit,
#' including singularity flags and correlation bounds.
get_stage1_diagnostics <- function(fit_obj, stage2_df, predictor_u0 = "u0_eb", predictor_u1 = "u1_eb") {
  if (is.null(fit_obj)) {
    return(empty_stage1_diagnostics())
  }

  vcov_re <- tryCatch(as.matrix(lme4::VarCorr(fit_obj)[[1]]), error = function(e) matrix(NA_real_, nrow = 2L, ncol = 2L))
  re_corr <- tryCatch(stats::cov2cor(vcov_re)[1, 2], error = function(e) NA_real_)
  lmer_singular <- tryCatch(lme4::isSingular(fit_obj, tol = 1e-5), error = function(e) NA)

  # Filter to complete cases to ensure the empirical Bayes correlation and
  # design kappa calculations match the exact data available for stage-2 analysis.
  dat <- stage2_df[, c(predictor_u0, predictor_u1), drop = FALSE]
  dat <- dat[stats::complete.cases(dat), , drop = FALSE]

  eb_corr <- if (nrow(dat) >= 2L &&
    is.finite(stats::sd(dat[[predictor_u0]])) &&
    is.finite(stats::sd(dat[[predictor_u1]])) &&
    stats::sd(dat[[predictor_u0]]) > sqrt(.Machine$double.eps) &&
    stats::sd(dat[[predictor_u1]]) > sqrt(.Machine$double.eps)) {
    stats::cor(dat[[predictor_u0]], dat[[predictor_u1]])
  } else {
    NA_real_
  }

  # Compute the condition number (kappa) of the design matrix. A kappa > 1e6
  # suggests extreme multicollinearity that can destabilize stage-2 estimators.
  design_kappa <- if (nrow(dat) >= 3L) {
    tryCatch(kappa(cbind(1, as.matrix(dat))), error = function(e) NA_real_)
  } else {
    NA_real_
  }

  # Check for near-perfect correlations which often indicate structural
  # non-identifiability or zero variance components in the simulation.
  boundary_re_corr <- is.finite(re_corr) && abs(re_corr) > 0.999
  boundary_eb_corr <- is.finite(eb_corr) && abs(eb_corr) > 0.999
  singular_problem <- isTRUE(lmer_singular) || boundary_re_corr || boundary_eb_corr

  detail_parts <- c(
    if (isTRUE(lmer_singular)) "lmer_is_singular" else NA_character_,
    if (boundary_re_corr) sprintf("random_effect_corr=%0.6f", re_corr) else NA_character_,
    if (boundary_eb_corr) sprintf("eb_corr=%0.6f", eb_corr) else NA_character_,
    if (is.finite(design_kappa) && design_kappa > 1e6) sprintf("design_kappa=%0.3e", design_kappa) else NA_character_
  )

  tibble::tibble(
    stage1_singular_problem = singular_problem,
    stage1_problem_detail = if (isTRUE(singular_problem)) compact_message(detail_parts) else "ok",
    stage1_lmer_singular = lmer_singular,
    stage1_re_corr = re_corr,
    stage1_eb_corr = eb_corr,
    stage1_design_kappa = design_kappa
  )
}

