#' Stage-2 estimators for observed EB/corrected-score regressions.
#'
#' The functions in this file fit second-stage regressions where the predictor
#' of substantive interest is a cluster-level random slope estimate. To keep
#' simulation summaries directly comparable across estimator variants, the
#' reported slope is scaled by the empirical standard deviation of the slope
#' predictor. A one-unit estimate therefore corresponds to a one-SD change in
#' the stage-2 slope predictor.

#' Construct the empty OLS result template.
#'
#' @details
#' OLS-based estimators return both model-based (`naive`) and HC3 robust
#' standard-error variants. This helper centralizes the NA-filled return shape
#' used when the model cannot be estimated, keeping downstream row-binding
#' stable.
#'
#' @return
#' A tibble with one row for each `se_type` (`naive`, `hc3`) and columns
#' `estimate`, `se`, `ci_low`, `ci_high`, and `status_code`, all filled with
#' missing values.
empty_lm_variants <- function() {
  dplyr::bind_rows(
    naive = tibble::tibble(
      estimate = NA_real_,
      se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      status_code = NA_integer_
    ),
    hc3 = tibble::tibble(
      estimate = NA_real_,
      se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      status_code = NA_integer_
    ),
    .id = "se_type"
  )
}

#' Fit a single-predictor observed-score stage-2 OLS regression.
#'
#' @details
#' Fits `outcome ~ predictor` after complete-case filtering. The coefficient
#' for `predictor` and its standard errors are multiplied by the observed
#' standard deviation of `predictor`, so the reported estimand is the expected
#' change in `outcome` per one observed SD of the predictor.
#'
#' Both conventional model-based and HC3 robust standard errors are returned.
#' Failures are represented by the `empty_lm_variants()` template rather than
#' by thrown errors, because these functions are used inside simulation loops.
#'
#' @param stage2_df Data frame containing the outcome and stage-2 predictor.
#' @param outcome Character scalar naming the outcome column.
#' @param predictor Character scalar naming the single predictor column.
#'
#' @return
#' A two-row tibble keyed by `se_type`, with scaled `estimate`, `se`,
#' Wald-normal confidence limits, and `status_code` (`0L` for successfully
#' estimated rows).
fit_observed_single <- function(stage2_df, outcome, predictor) {
  empty_out <- empty_lm_variants()

  # Restrict to the variables used by this estimator before applying the
  # complete-case rule; missingness in unrelated simulation columns is ignored.
  dat <- stage2_df[, c(outcome, predictor), drop = FALSE]
  dat <- dat[stats::complete.cases(dat), , drop = FALSE]

  # Need at least three rows for an intercept, slope, and residual degree of
  # freedom, and a non-degenerate predictor to make the scaled slope meaningful.
  if (nrow(dat) < 3L || !is.finite(stats::sd(dat[[predictor]])) ||
    stats::sd(dat[[predictor]]) <= sqrt(.Machine$double.eps)) {
    return(empty_out)
  }

  fit <- tryCatch(stats::lm(stats::as.formula(paste(outcome, "~", predictor)), data = dat), error = function(e) NULL)
  if (is.null(fit)) {
    return(empty_out)
  }

  coef_tab <- summary(fit)$coefficients
  if (!(predictor %in% rownames(coef_tab))) {
    return(empty_out)
  }

  # Report the effect per observed SD of the stage-2 predictor. This keeps
  # observed EB, corrected-score, and true-score estimators on the same scale.
  scale_u1 <- stats::sd(dat[[predictor]])
  est <- unname(coef_tab[predictor, "Estimate"]) * scale_u1
  se_naive <- unname(coef_tab[predictor, "Std. Error"]) * scale_u1

  # HC3 is the leverage-adjusted sandwich variance used for the robust OLS
  # variant. If sandwich fails, only the HC3 row is marked unavailable.
  vcov_hc3 <- tryCatch(sandwich::vcovHC(fit, type = "HC3"), error = function(e) NULL)
  se_hc3 <- if (!is.null(vcov_hc3) && predictor %in% rownames(vcov_hc3)) {
    sqrt(unname(vcov_hc3[predictor, predictor])) * scale_u1
  } else {
    NA_real_
  }

  dplyr::bind_rows(
    naive = tibble::tibble(
      estimate = est,
      se = se_naive,
      ci_low = est - stats::qnorm(0.975) * se_naive,
      ci_high = est + stats::qnorm(0.975) * se_naive,
      status_code = 0L
    ),
    hc3 = tibble::tibble(
      estimate = est,
      se = se_hc3,
      ci_low = est - stats::qnorm(0.975) * se_hc3,
      ci_high = est + stats::qnorm(0.975) * se_hc3,
      status_code = ifelse(is.finite(se_hc3), 0L, NA_integer_)
    ),
    .id = "se_type"
  )
}

#' Fit a dual-predictor observed-score stage-2 OLS regression.
#'
#' @details
#' Fits `outcome ~ predictor_u0 + predictor_u1` after complete-case filtering
#' and reports the coefficient for `predictor_u1`. The intercept-like random
#' effect proxy (`predictor_u0`) is included as an adjustment variable, while
#' the slope-like random effect proxy (`predictor_u1`) is the target estimand.
#'
#' The returned estimate and standard errors are scaled by the observed
#' standard deviation of `predictor_u1`. This mirrors `fit_observed_single()`
#' and supports direct comparisons across simulation conditions.
#'
#' @param stage2_df Data frame containing the outcome and both predictors.
#' @param outcome Character scalar naming the outcome column.
#' @param predictor_u0 Character scalar naming the intercept-like predictor.
#' @param predictor_u1 Character scalar naming the slope-like predictor whose
#' coefficient is reported.
#'
#' @return
#' A two-row tibble keyed by `se_type`, with scaled `estimate`, `se`,
#' Wald-normal confidence limits, and `status_code` (`0L` for successfully
#' estimated rows).
fit_observed_dual <- function(stage2_df, outcome, predictor_u0, predictor_u1) {
  empty_out <- empty_lm_variants()

  # Complete-case filtering is estimator-specific: rows only need to be
  # observed on the outcome and the two predictors used here.
  dat <- stage2_df[, c(outcome, predictor_u0, predictor_u1), drop = FALSE]
  dat <- dat[stats::complete.cases(dat), , drop = FALSE]

  # Four rows are the minimum for a two-predictor model with residual degrees
  # of freedom. Degenerate slope predictors cannot be converted to one-SD units.
  if (nrow(dat) < 4L || !is.finite(stats::sd(dat[[predictor_u1]])) ||
    stats::sd(dat[[predictor_u1]]) <= sqrt(.Machine$double.eps)) {
    return(empty_out)
  }

  fit <- tryCatch(
    stats::lm(stats::as.formula(paste(outcome, "~", predictor_u0, "+", predictor_u1)), data = dat),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(empty_out)
  }

  coef_tab <- summary(fit)$coefficients
  if (!(predictor_u1 %in% rownames(coef_tab))) {
    return(empty_out)
  }

  # Only the slope-like predictor is reported; the u0 proxy is an adjustment
  # covariate used to reduce confounding between random intercepts and slopes.
  scale_u1 <- stats::sd(dat[[predictor_u1]])
  est <- unname(coef_tab[predictor_u1, "Estimate"]) * scale_u1
  se_naive <- unname(coef_tab[predictor_u1, "Std. Error"]) * scale_u1

  # Compute the HC3 standard error for the target coefficient on the same
  # one-SD scale as the conventional OLS standard error.
  vcov_hc3 <- tryCatch(sandwich::vcovHC(fit, type = "HC3"), error = function(e) NULL)
  se_hc3 <- if (!is.null(vcov_hc3) && predictor_u1 %in% rownames(vcov_hc3)) {
    sqrt(unname(vcov_hc3[predictor_u1, predictor_u1])) * scale_u1
  } else {
    NA_real_
  }

  dplyr::bind_rows(
    naive = tibble::tibble(
      estimate = est,
      se = se_naive,
      ci_low = est - stats::qnorm(0.975) * se_naive,
      ci_high = est + stats::qnorm(0.975) * se_naive,
      status_code = 0L
    ),
    hc3 = tibble::tibble(
      estimate = est,
      se = se_hc3,
      ci_low = est - stats::qnorm(0.975) * se_hc3,
      ci_high = est + stats::qnorm(0.975) * se_hc3,
      status_code = ifelse(is.finite(se_hc3), 0L, NA_integer_)
    ),
    .id = "se_type"
  )
}

#' Convert OLS standard-error variants to method-labelled estimator rows.
#'
#' @details
#' `fit_observed_single()` and `fit_observed_dual()` return rows keyed by
#' `se_type` so the model fit and the robust variance calculation can stay
#' together. Simulation summaries instead group by `method`, so this helper
#' maps the `naive` row to `base_method` and the `hc3` row to
#' `paste0(base_method, "_hc3")`.
#'
#' @param fit_tbl Tibble returned by an OLS stage-2 fitting helper.
#' @param base_method Character scalar naming the non-HC3 method.
#'
#' @return
#' A tibble with `method`, `estimate`, `se`, `ci_low`, `ci_high`, and
#' `status_code`.
finalize_ols_se_variants <- function(fit_tbl, base_method) {
  fit_tbl <- dplyr::mutate(
    fit_tbl,
    method = ifelse(se_type == "hc3", paste0(base_method, "_hc3"), base_method)
  )
  dplyr::select(fit_tbl, method, estimate, se, ci_low, ci_high, status_code)
}

#' Fit a ridge-regularized dual-predictor stage-2 regression.
#'
#' @details
#' Fits a Gaussian ridge regression with `glmnet::cv.glmnet()` using
#' `predictor_u0` and `predictor_u1`, then reports the cross-validated
#' coefficient for `predictor_u1` scaled by the observed SD of `predictor_u1`.
#' Ridge is used as a point-estimation fallback for collinear stage-2 predictors;
#' no analytic standard error or confidence interval is computed here.
#'
#' The number of folds is capped at 10 and floored at 3, with the sample size
#' used when fewer than 10 complete cases are available.
#'
#' @param stage2_df Data frame containing the outcome and both predictors.
#' @param outcome Character scalar naming the outcome column.
#' @param predictor_u0 Character scalar naming the intercept-like predictor.
#' @param predictor_u1 Character scalar naming the slope-like predictor whose
#' coefficient is reported.
#' @param lambda_rule Character scalar naming the `cv.glmnet` lambda selector,
#' typically `"lambda.1se"` or `"lambda.min"`.
#'
#' @return
#' A one-row tibble with scaled `estimate`, missing `se` and confidence limits,
#' `status_code`, and the selected `ridge_lambda` plus `ridge_lambda_rule`.
fit_ridge_dual <- function(stage2_df, outcome, predictor_u0, predictor_u1, lambda_rule = "lambda.1se") {
  # Ridge uses the same complete-case and non-degeneracy checks as the dual OLS
  # estimator so comparisons are driven by estimator behavior, not row sets.
  dat <- stage2_df[, c(outcome, predictor_u0, predictor_u1), drop = FALSE]
  dat <- dat[stats::complete.cases(dat), , drop = FALSE]

  if (nrow(dat) < 4L || !is.finite(stats::sd(dat[[predictor_u1]])) ||
    stats::sd(dat[[predictor_u1]]) <= sqrt(.Machine$double.eps)) {
    return(tibble::tibble(
      estimate = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      status_code = NA_integer_, ridge_lambda = NA_real_, ridge_lambda_rule = lambda_rule
    ))
  }

  x_mat <- as.matrix(dat[, c(predictor_u0, predictor_u1), drop = FALSE])
  y_vec <- dat[[outcome]]

  # Cross-validation needs at least three folds; with small simulation samples,
  # use leave-one-ish folds up to the usual 10-fold cap.
  nfolds <- max(3L, min(10L, nrow(dat)))
  cv_fit <- tryCatch(
    glmnet::cv.glmnet(
      x = x_mat,
      y = y_vec,
      family = "gaussian",
      alpha = 0,
      nfolds = nfolds,
      standardize = TRUE,
      intercept = TRUE,
      type.measure = "mse"
    ),
    error = function(e) NULL
  )

  if (is.null(cv_fit)) {
    return(tibble::tibble(
      estimate = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      status_code = NA_integer_, ridge_lambda = NA_real_, ridge_lambda_rule = lambda_rule
    ))
  }

  lambda_value <- cv_fit[[lambda_rule]]

  # Extracting by the lambda rule rather than by the numeric lambda preserves
  # glmnet's handling of named selectors such as "lambda.1se".
  coef_mat <- tryCatch(as.matrix(stats::coef(cv_fit, s = lambda_rule)), error = function(e) NULL)
  if (is.null(coef_mat) || !(predictor_u1 %in% rownames(coef_mat))) {
    return(tibble::tibble(
      estimate = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      status_code = NA_integer_, ridge_lambda = lambda_value, ridge_lambda_rule = lambda_rule
    ))
  }

  scale_u1 <- stats::sd(dat[[predictor_u1]])
  tibble::tibble(
    estimate = unname(coef_mat[predictor_u1, 1]) * scale_u1,
    se = NA_real_,
    ci_low = NA_real_,
    ci_high = NA_real_,
    status_code = 0L,
    ridge_lambda = lambda_value,
    ridge_lambda_rule = lambda_rule
  )
}

#' Fit a dual-predictor errors-in-variables corrected-score estimator.
#'
#' @details
#' This estimator treats the observed stage-2 predictors as noisy measurements
#' of latent random-effect scores. For each complete case it subtracts the
#' supplied measurement-error covariance matrix from `x_i x_i'` before solving
#' the corrected normal equations. The reported target is the `predictor_u1`
#' coefficient multiplied by an estimated latent SD for `predictor_u1`.
#'
#' The estimating equation uses predictors ordered as
#' `(intercept, predictor_u0, predictor_u1)`. The measurement-error variances
#' and covariance (`meas11`, `meas12`, `meas22`) fill the predictor block only;
#' the intercept is treated as measured without error. `measurement_weight`
#' scales the supplied covariance terms, enabling sensitivity checks that
#' attenuate or amplify the correction.
#'
#' A sandwich variance is computed from the empirical estimating-function
#' residuals. If `stabilize_a_mat` is `TRUE`, the full corrected cross-product
#' matrix is projected to positive definite before solving. If
#' `ridge_predictor_block` is `TRUE`, only the two-predictor block is ridged up
#' to `ridge_min_eigen`.
#'
#' @param stage2_df Data frame containing outcome, predictors, and measurement
#' error covariance columns.
#' @param outcome Character scalar naming the outcome column.
#' @param predictor_u0 Character scalar naming the intercept-like observed
#' predictor.
#' @param predictor_u1 Character scalar naming the slope-like observed
#' predictor whose coefficient is reported.
#' @param meas11 Character scalar naming the measurement-error variance column
#' for `predictor_u0`.
#' @param meas12 Character scalar naming the measurement-error covariance column
#' between `predictor_u0` and `predictor_u1`.
#' @param meas22 Character scalar naming the measurement-error variance column
#' for `predictor_u1`.
#' @param outcome_meas_var Optional character scalar naming an outcome
#' measurement-error variance column. The current estimator does not use this
#' variance in the estimating equations, but includes it in complete-case
#' filtering when supplied so corrected-outcome workflows use a consistent row
#' set.
#' @param stabilize_a_mat Logical; project the corrected normal-equation matrix
#' to positive definite before solving.
#' @param min_eigen Numeric lower bound used by positive-definite projections.
#' @param ridge_predictor_block Logical; add ridge only to the corrected
#' two-predictor cross-product block.
#' @param ridge_min_eigen Numeric minimum eigenvalue targeted by the predictor
#' block ridge.
#' @param measurement_weight Numeric multiplier applied to the supplied
#' predictor measurement-error covariance terms.
#'
#' @return
#' A one-row tibble with scaled `estimate`, sandwich `se`, Wald-normal
#' confidence limits, and `status_code`. Status `0L` indicates success, `1L`
#' indicates that the corrected normal equations could not be solved, and `2L`
#' indicates a non-finite scaled estimate or standard error.
fit_eiv_dual <- function(stage2_df,
                         outcome,
                         predictor_u0,
                         predictor_u1,
                         meas11,
                         meas12,
                         meas22,
                         outcome_meas_var = NULL,
                         stabilize_a_mat = FALSE,
                         min_eigen = 1e-6,
                         ridge_predictor_block = FALSE,
                         ridge_min_eigen = 1e-4,
                         measurement_weight = 1) {
  out_fail <- tibble::tibble(
    estimate = NA_real_,
    se = NA_real_,
    ci_low = NA_real_,
    ci_high = NA_real_,
    status_code = NA_integer_
  )

  # Include the optional outcome measurement-variance column only to enforce a
  # consistent complete-case sample when upstream simulations create it.
  cols_needed <- c(outcome, predictor_u0, predictor_u1, meas11, meas12, meas22)
  if (!is.null(outcome_meas_var)) {
    cols_needed <- c(cols_needed, outcome_meas_var)
  }

  dat <- stage2_df[, cols_needed, drop = FALSE]
  dat <- dat[stats::complete.cases(dat), , drop = FALSE]

  # EIV needs enough rows to estimate a three-parameter equation and its
  # sandwich variance with some stability. The target predictor must also vary.
  if (nrow(dat) < 8L || !is.finite(stats::sd(dat[[predictor_u1]])) ||
    stats::sd(dat[[predictor_u1]]) <= sqrt(.Machine$double.eps)) {
    return(out_fail)
  }

  w_mat <- as.matrix(dat[, c(predictor_u0, predictor_u1), drop = FALSE])
  y_vec <- dat[[outcome]]

  # Variance components should not reduce observed second moments because of
  # negative numerical estimates; clamp variances at zero but retain covariance.
  s11 <- pmax(dat[[meas11]], 0)
  s12 <- dat[[meas12]]
  s22 <- pmax(dat[[meas22]], 0)

  # Build the per-row measurement-error covariance array for
  # (intercept, predictor_u0, predictor_u1). The intercept row/column remains 0.
  x_mat <- cbind(1, w_mat)

  # Corrected normal equations:
  #   sum_i (x_i x_i' - S_i) beta = sum_i x_i y_i
  # where S_i is the supplied measurement-error covariance for the predictors.
  a_mat <- crossprod(x_mat)
  a_mat[2, 2] <- a_mat[2, 2] - measurement_weight * sum(s11)
  a_mat[2, 3] <- a_mat[2, 3] - measurement_weight * sum(s12)
  a_mat[3, 2] <- a_mat[3, 2] - measurement_weight * sum(s12)
  a_mat[3, 3] <- a_mat[3, 3] - measurement_weight * sum(s22)

  b_vec <- as.vector(crossprod(x_mat, y_vec))

  a_mat_use <- a_mat
  if (isTRUE(stabilize_a_mat)) {
    a_mat_use <- project_to_pd(a_mat_use, min_eigen = min_eigen)
  }
  if (isTRUE(ridge_predictor_block)) {
    # Ridge only the observed predictor block, leaving the intercept coupling
    # intact while preventing near-singular corrected predictor covariance.
    pred_block <- (a_mat_use[2:3, 2:3, drop = FALSE] + t(a_mat_use[2:3, 2:3, drop = FALSE])) / 2
    eig_min <- min(eigen(pred_block, symmetric = TRUE, only.values = TRUE)$values)
    ridge_lambda <- max(0, ridge_min_eigen - eig_min)
    a_mat_use[2:3, 2:3] <- pred_block + diag(ridge_lambda, nrow = 2L)
  }

  beta_hat <- tryCatch(solve(a_mat_use, b_vec), error = function(e) NULL)
  if (is.null(beta_hat) || any(!is.finite(beta_hat))) {
    return(dplyr::mutate(out_fail, status_code = 1L))
  }

  # If stabilization changed the estimating matrix, include the average matrix
  # perturbation in each influence contribution so the sandwich variance matches
  # the equation actually solved.
  a_adjustment <- (a_mat_use - a_mat) / nrow(dat)
  
  term1 <- x_mat * y_vec
  term2 <- x_mat * as.vector(x_mat %*% beta_hat)
  
  term3 <- matrix(0, nrow = nrow(dat), ncol = 3L)
  term3[, 2] <- measurement_weight * (s11 * beta_hat[2] + s12 * beta_hat[3])
  term3[, 3] <- measurement_weight * (s12 * beta_hat[2] + s22 * beta_hat[3])
  
  adj_beta <- as.vector(a_adjustment %*% beta_hat)
  term4 <- matrix(adj_beta, nrow = nrow(dat), ncol = 3L, byrow = TRUE)
  
  psi_mat <- term1 - term2 + term3 - term4

  # Empirical sandwich variance for beta_hat. The generalized inverse fallback
  # preserves a finite variance estimate in borderline stabilized cases.
  bread_inv <- tryCatch(solve(a_mat_use), error = function(e) MASS::ginv(a_mat_use))
  meat <- crossprod(psi_mat)
  vcov_beta <- bread_inv %*% meat %*% bread_inv

  # Estimate the latent predictor covariance by subtracting the aggregate
  # measurement-error covariance from the centered observed cross-product.
  w_centered <- scale(w_mat, center = TRUE, scale = FALSE)
  sum_s <- measurement_weight * matrix(c(sum(s11), sum(s12), sum(s12), sum(s22)), nrow = 2L, byrow = TRUE)
  sigma_x_hat <- (crossprod(w_centered) - sum_s) / max(1, nrow(dat) - 1L)
  sigma_x_hat <- project_to_pd(sigma_x_hat, min_eigen = min_eigen)
  scale_u1 <- sqrt(sigma_x_hat[2, 2])

  # Convert the raw latent-slope coefficient and its sandwich SE to the same
  # one-SD target scale used by the observed-score estimators.
  est <- unname(beta_hat[[3]]) * scale_u1
  se_beta1 <- if (nrow(vcov_beta) >= 3L) sqrt(unname(vcov_beta[3, 3])) else NA_real_
  se <- se_beta1 * scale_u1

  if (!is.finite(est) || !is.finite(se)) {
    return(dplyr::mutate(out_fail, status_code = 2L))
  }

  tibble::tibble(
    estimate = est,
    se = se,
    ci_low = est - stats::qnorm(0.975) * se,
    ci_high = est + stats::qnorm(0.975) * se,
    status_code = 0L
  )
}

#' Format stacked-sandwich covariance variants as estimator rows.
#'
#' @details
#' `stacked_sandwich_for_corrected_scores()` returns a fitted stage-2
#' coefficient vector and several covariance matrix variants. Simulation
#' summaries expect one row per method, so this helper converts HC0-HC3
#' covariance entries into method-labelled rows with Wald confidence limits.
#'
#' @param sandwich_out List returned by `stacked_sandwich_for_corrected_scores()`.
#' @param df Degrees of freedom for the `t` critical value.
#' @param alpha_names Character vector naming the entries of
#' `sandwich_out$alpha_hat`. Defaults to `c("(Intercept)", "x")`.
#' @param term Character scalar naming the coefficient to report.
#' @param method_prefix Prefix used to build method names. Variants append
#' `"_hc0"` through `"_hc3"`.
#'
#' @return A tibble with `method`, `estimate`, `se`, `ci_low`, and `ci_high`.
format_stacked_sandwich_rows <- function(sandwich_out,
                                         df,
                                         alpha_names = c("(Intercept)", "x"),
                                         term = "x",
                                         method_prefix = "corrected_full_stacked") {
  variants <- c("hc0", "hc1", "hc2", "hc3")
  x_idx <- which(alpha_names == term)

  if (length(x_idx) != 1L) {
    return(purrr::map_dfr(variants, function(variant) {
      tibble::tibble(
        method = paste0(method_prefix, "_", variant),
        estimate = NA_real_,
        se = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_
      )
    }))
  }

  crit <- stats::qt(0.975, df)
  purrr::map_dfr(variants, function(variant) {
    vcov_name <- paste0("vcov_", variant)
    vcov_mat <- sandwich_out[[vcov_name]]
    se <- if (!is.null(vcov_mat) && nrow(vcov_mat) >= x_idx && ncol(vcov_mat) >= x_idx) {
      sqrt(unname(vcov_mat[x_idx, x_idx]))
    } else {
      NA_real_
    }
    est <- unname(sandwich_out$alpha_hat[x_idx])

    tibble::tibble(
      method = paste0(method_prefix, "_", variant),
      estimate = est,
      se = se,
      ci_low = est - crit * se,
      ci_high = est + crit * se
    )
  })
}
