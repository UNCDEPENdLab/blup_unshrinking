#' Fit BLUP-Adjusted Dual
#'
#' @description
#' NOTE: This function has been archived because it demonstrates a theoretical flaw
#' in attempting to apply Errors-in-Variables (EIV) corrections directly to BLUPs.
#' 
#' **Why this approach fails (Attenuation Bias):**
#' If you start with shrunken Empirical Bayes estimates (BLUPs), the variance of the
#' BLUP is strictly smaller than the true variance. To "un-shrink" the variance, this 
#' function correctly *adds* the posterior variance to the predictor cross-product matrix 
#' (the denominator). However, the numerator (the covariance between the BLUP and the outcome)
#' remains shrunken. Dividing a shrunken numerator by an un-shrunken denominator 
#' results in a slope estimator that is severely attenuated (biased toward zero).
#' 
#' **What is required to correctly use BLUPs:**
#' To properly avoid attenuation bias when using BLUPs, one must "un-shrink" BOTH the 
#' denominator (variance) and the numerator (covariance). This requires extracting the 
#' cluster-specific reliability matrices ($\Lambda_i$) and pre/post-multiplying the 
#' sufficient statistics by $\Lambda_i^{-1}$. This is exactly what Lai's 2S-PA method does.
#' 
#' Alternatively, and far more simply, one can just use the unshrunken Maximum Likelihood 
#' (OLS) scores and apply conventional EIV (subtracting the sampling variance), which is 
#' mathematically isomorphic to Lai's 2S-PA but computationally much faster.
#'
#' @param stage2_df Data frame containing outcome, predictors (BLUPs), and posterior
#' variance columns.
#' @param outcome Character scalar naming the outcome column.
#' @param predictor_u0 Character scalar naming the intercept-like BLUP predictor.
#' @param predictor_u1 Character scalar naming the slope-like BLUP predictor whose
#' coefficient is reported.
#' @param meas11 Character scalar naming the posterior variance column for `predictor_u0`.
#' @param meas12 Character scalar naming the posterior covariance column.
#' @param meas22 Character scalar naming the posterior variance column for `predictor_u1`.
#' @param outcome_meas_var Optional character scalar naming an outcome measurement-error
#' variance column (used only for consistent complete-case filtering).
#' @param stabilize_a_mat Logical; project the adjusted normal-equation matrix
#' to positive definite before solving.
#' @param min_eigen Numeric lower bound used by positive-definite projections.
#' @param ridge_predictor_block Logical; add ridge only to the corrected predictor block.
#' @param ridge_min_eigen Numeric minimum eigenvalue targeted by the predictor block ridge.
#' @param measurement_weight Numeric multiplier applied to the supplied posterior
#' covariance terms.
#'
#' @return
#' A one-row tibble with scaled `estimate`, sandwich `se`, Wald-normal
#' confidence limits, and `status_code`.
fit_blup_adjusted_dual <- function(stage2_df,
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

  cols_needed <- c(outcome, predictor_u0, predictor_u1, meas11, meas12, meas22)
  if (!is.null(outcome_meas_var)) {
    cols_needed <- c(cols_needed, outcome_meas_var)
  }

  dat <- stage2_df[, cols_needed, drop = FALSE]
  dat <- dat[stats::complete.cases(dat), , drop = FALSE]

  if (nrow(dat) < 8L || !is.finite(stats::sd(dat[[predictor_u1]])) ||
    stats::sd(dat[[predictor_u1]]) <= sqrt(.Machine$double.eps)) {
    return(out_fail)
  }

  w_mat <- as.matrix(dat[, c(predictor_u0, predictor_u1), drop = FALSE])
  y_vec <- dat[[outcome]]

  s11 <- pmax(dat[[meas11]], 0)
  s12 <- dat[[meas12]]
  s22 <- pmax(dat[[meas22]], 0)

  x_mat <- cbind(1, w_mat)

  # Adjusted normal equations for BLUPs:
  #   sum_i (x_i x_i' + S_i) beta = sum_i x_i y_i
  # where S_i is the posterior covariance matrix for the BLUPs.
  a_mat <- crossprod(x_mat)
  a_mat[2, 2] <- a_mat[2, 2] + measurement_weight * sum(s11)
  a_mat[2, 3] <- a_mat[2, 3] + measurement_weight * sum(s12)
  a_mat[3, 2] <- a_mat[3, 2] + measurement_weight * sum(s12)
  a_mat[3, 3] <- a_mat[3, 3] + measurement_weight * sum(s22)

  b_vec <- as.vector(crossprod(x_mat, y_vec))

  a_mat_use <- a_mat
  if (isTRUE(stabilize_a_mat)) {
    a_mat_use <- project_to_pd(a_mat_use, min_eigen = min_eigen)
  }
  if (isTRUE(ridge_predictor_block)) {
    pred_block <- (a_mat_use[2:3, 2:3, drop = FALSE] + t(a_mat_use[2:3, 2:3, drop = FALSE])) / 2
    eig_min <- min(eigen(pred_block, symmetric = TRUE, only.values = TRUE)$values)
    ridge_lambda <- max(0, ridge_min_eigen - eig_min)
    a_mat_use[2:3, 2:3] <- pred_block + diag(ridge_lambda, nrow = 2L)
  }

  beta_hat <- tryCatch(solve(a_mat_use, b_vec), error = function(e) NULL)
  if (is.null(beta_hat) || any(!is.finite(beta_hat))) {
    return(dplyr::mutate(out_fail, status_code = 1L))
  }

  a_adjustment <- (a_mat_use - a_mat) / nrow(dat)
  
  term1 <- x_mat * y_vec
  term2 <- x_mat * as.vector(x_mat %*% beta_hat)
  
  term3 <- matrix(0, nrow = nrow(dat), ncol = 3L)
  term3[, 2] <- measurement_weight * (s11 * beta_hat[2] + s12 * beta_hat[3])
  term3[, 3] <- measurement_weight * (s12 * beta_hat[2] + s22 * beta_hat[3])
  
  adj_beta <- as.vector(a_adjustment %*% beta_hat)
  term4 <- matrix(adj_beta, nrow = nrow(dat), ncol = 3L, byrow = TRUE)
  
  # For BLUPs, the estimating equation was X_i y_i - (X_i X_i' + S_i + A_adj) beta
  psi_mat <- term1 - term2 - term3 - term4

  bread_inv <- tryCatch(solve(a_mat_use), error = function(e) MASS::ginv(a_mat_use))
  meat <- crossprod(psi_mat)
  vcov_beta <- bread_inv %*% meat %*% bread_inv

  # Estimate the latent predictor covariance by adding the aggregate
  # posterior covariance to the centered observed cross-product.
  w_centered <- scale(w_mat, center = TRUE, scale = FALSE)
  sum_s <- measurement_weight * matrix(c(sum(s11), sum(s12), sum(s12), sum(s22)), nrow = 2L, byrow = TRUE)
  sigma_x_hat <- (crossprod(w_centered) + sum_s) / max(1, nrow(dat) - 1L)
  sigma_x_hat <- project_to_pd(sigma_x_hat, min_eigen = min_eigen)
  scale_u1 <- sqrt(sigma_x_hat[2, 2])

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
