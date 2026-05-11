#' Stacked sandwich helpers for corrected-score two-stage estimators.
#'
#' These helpers implement the reusable estimation pieces behind the
#' random-slope sandwich coverage simulation: mixed-model parameter packing,
#' cluster-level likelihood evaluation, likelihood-only corrected scores,
#' stage-1 score/bread extraction, and full stacked M-estimation covariance
#' assembly.

#' Pack an `lme4` fit into an unconstrained standard-deviation parameter vector.
#'
#' @details
#' The finite-difference backend works on an unconstrained scale:
#' fixed effects are left unchanged, random-effect standard deviations and the
#' residual standard deviation are log-transformed, and the random-effect
#' correlation is Fisher-z transformed.
#'
#' @param fit Fitted `lme4` model with fixed effects `(Intercept)` and `z`, and
#' a random intercept/slope covariance matrix named `id`.
#'
#' @return Numeric vector `c(beta0, beta1, log_sd0, log_sd1, atanh_rho,
#' log_sigma)`.
pack_psi <- function(fit) {
  beta_hat <- unname(lme4::fixef(fit)[c("(Intercept)", "z")])
  g_hat <- as.matrix(lme4::VarCorr(fit)$id)
  sd0 <- sqrt(g_hat[1, 1])
  sd1 <- sqrt(g_hat[2, 2])
  rho <- g_hat[1, 2] / (sd0 * sd1)
  rho <- max(min(rho, 0.999), -0.999)
  c(
    beta_hat[1],
    beta_hat[2],
    log(sd0),
    log(sd1),
    atanh(rho),
    log(stats::sigma(fit))
  )
}

#' Pack an `lme4` fit into a variance-scale parameter vector.
#'
#' @details
#' Analytical, `merDeriv`, and TMB paths operate on the variance scale rather
#' than the unconstrained finite-difference scale.
#'
#' @param fit Fitted `lme4` model with fixed effects `(Intercept)` and `z`, and
#' a random intercept/slope covariance matrix named `id`.
#'
#' @return Numeric vector `c(beta0, beta1, var0, cov01, var1, sigma2)`.
pack_psi_var <- function(fit) {
  beta_hat <- unname(lme4::fixef(fit)[c("(Intercept)", "z")])
  g_hat <- unclass(lme4::VarCorr(fit)[["id"]])
  c(
    beta_hat[1],
    beta_hat[2],
    g_hat[1, 1],
    g_hat[2, 1],
    g_hat[2, 2],
    stats::sigma(fit)^2
  )
}

#' Unpack an unconstrained parameter vector.
#'
#' @param psi_vec Numeric vector from `pack_psi()`.
#'
#' @return A list with fixed effects `beta`, random-effect covariance `G`, and
#' residual standard deviation `sigma`.
unpack_psi <- function(psi_vec) {
  beta <- psi_vec[1:2]
  sd0 <- exp(psi_vec[3])
  sd1 <- exp(psi_vec[4])
  rho <- tanh(psi_vec[5])
  sigma <- exp(psi_vec[6])
  g_mat <- matrix(
    c(sd0^2, rho * sd0 * sd1, rho * sd0 * sd1, sd1^2),
    nrow = 2
  )
  list(beta = beta, G = g_mat, sigma = sigma)
}

#' Unpack a variance-scale parameter vector.
#'
#' @param psi_vec Numeric vector from `pack_psi_var()`.
#'
#' @return A list with fixed effects `beta`, random-effect covariance `G`, and
#' residual standard deviation `sigma`.
unpack_psi_var <- function(psi_vec) {
  psi_vec <- as.numeric(psi_vec)
  beta <- psi_vec[1:2]
  g_mat <- matrix(
    c(psi_vec[3], psi_vec[4], psi_vec[4], psi_vec[5]),
    nrow = 2
  )
  sigma2 <- psi_vec[6]
  list(beta = beta, G = g_mat, sigma = sqrt(sigma2))
}

#' Precompute cluster-level matrices for Gaussian random-slope calculations.
#'
#' @details
#' The sandwich estimator repeatedly evaluates cluster likelihoods and
#' corrected scores. This helper caches design matrices and cross-products so
#' derivative backends do not rebuild them for every parameter perturbation.
#'
#' @param split_dat Named list of cluster data frames with columns `y` and `z`.
#'
#' @return A named list of cluster objects containing `y`, `X`, `Z`,
#' cross-products, cluster size, cluster OLS coefficients, and the derivative
#' of the slope score with respect to fixed effects.
prepare_cluster_objects <- function(split_dat) {
  purrr::map(split_dat, function(df_i) {
    x_mat <- stats::model.matrix(~z, data = df_i)
    y_vec <- df_i$y

    ZtZ <- crossprod(x_mat)
    ZtX <- crossprod(x_mat, x_mat)
    ZtY <- crossprod(x_mat, y_vec)
    list(
      y = y_vec,
      X = x_mat,
      Z = x_mat,
      ZtZ = ZtZ,
      ZtX = ZtX,
      ZtY = ZtY,
      n = nrow(x_mat),
      ols_coef = drop(solve(ZtZ, ZtY)),
      dh_beta = -drop(solve(ZtZ, ZtX)[2, ])
    )
  })
}

#' Evaluate one cluster's marginal Gaussian log-likelihood.
#'
#' @param cluster_obj Cluster object from `prepare_cluster_objects()`.
#' @param psi_vec Numeric parameter vector.
#' @param parameterization Either `"sd"` for `pack_psi()` scale or `"var"` for
#' `pack_psi_var()` scale.
#'
#' @return Scalar log-likelihood. Returns `-Inf` if the marginal covariance is
#' not positive definite.
cluster_loglik_precomputed <- function(cluster_obj, psi_vec, parameterization = c("sd", "var")) {
  parameterization <- match.arg(parameterization)
  pars <- if (identical(parameterization, "sd")) unpack_psi(psi_vec) else unpack_psi_var(psi_vec)
  resid_vec <- cluster_obj$y - drop(cluster_obj$X %*% pars$beta)
  vy <- cluster_obj$Z %*% pars$G %*% t(cluster_obj$Z) + (pars$sigma^2) * diag(cluster_obj$n)
  chol_vy <- tryCatch(chol(vy), error = function(e) NULL)
  if (is.null(chol_vy)) {
    return(-Inf)
  }
  log_det <- 2 * sum(log(diag(chol_vy)))
  quad <- sum(backsolve(chol_vy, resid_vec, transpose = TRUE)^2)
  -0.5 * (cluster_obj$n * log(2 * pi) + log_det + quad)
}

#' Evaluate the total marginal Gaussian log-likelihood.
#'
#' @param cluster_objects List from `prepare_cluster_objects()`.
#' @param psi_vec Numeric parameter vector.
#' @param parameterization Either `"sd"` or `"var"`.
#'
#' @return Scalar summed log-likelihood across clusters.
total_loglik_precomputed <- function(cluster_objects, psi_vec, parameterization = c("sd", "var")) {
  parameterization <- match.arg(parameterization)
  sum(vapply(
    cluster_objects,
    cluster_loglik_precomputed,
    numeric(1),
    psi_vec = psi_vec,
    parameterization = parameterization
  ))
}

#' Compute a likelihood-only corrected random-slope score.
#'
#' @details
#' For this Gaussian random-slope simulation, the likelihood-only cluster slope
#' is the within-cluster OLS slope after subtracting the fixed-effect slope.
#'
#' @param cluster_obj Cluster object from `prepare_cluster_objects()`.
#' @param psi_vec Numeric parameter vector.
#' @param parameterization Either `"sd"` or `"var"`.
#'
#' @return Numeric scalar corrected random-slope score.
corrected_slope_from_precomputed <- function(cluster_obj, psi_vec, parameterization = c("sd", "var")) {
  parameterization <- match.arg(parameterization)
  pars <- if (identical(parameterization, "sd")) unpack_psi(psi_vec) else unpack_psi_var(psi_vec)
  unname(cluster_obj$ols_coef[2] - pars$beta[2])
}

#' Extract stage-1 score and bread components for the stacked sandwich.
#'
#' @details
#' This dispatcher routes to the analytical, `merDeriv`, TMB, or finite
#' difference path while returning a common structure consumed by
#' `stacked_sandwich_for_corrected_scores()`.
#'
#' @param fit_null Fitted first-stage `lme4` model.
#' @param split_dat Named list of cluster data frames.
#' @param cluster_objects List from `prepare_cluster_objects()`.
#' @param psi_hat Parameter vector on the finite-difference scale.
#' @param derivative_backend Backend object from `make_derivative_backend()`.
#'
#' @return A list containing `s1_mat`, `a11_hat`, `psi_stage1`,
#' `corrected_score_fn`, and `psi_parameterization`.
get_stage1_sandwich_inputs <- function(fit_null, split_dat, cluster_objects, psi_hat, derivative_backend) {
  if (identical(derivative_backend$name, "analytical")) {
    psi_stage1 <- pack_psi_var(fit_null)
    res <- get_exact_stage1_matrices(cluster_objects, psi_stage1)
    return(list(
      s1_mat = res$s1_mat,
      a11_hat = res$a11_hat,
      psi_stage1 = psi_stage1,
      corrected_score_fn = corrected_slope_from_precomputed,
      psi_parameterization = "var"
    ))
  }

  if (identical(derivative_backend$name, "merDeriv")) {
    if (!requireNamespace("merDeriv", quietly = TRUE)) {
      stop("The `merDeriv` package is required for the `merDeriv` backend.")
    }

    s1_mat <- sandwich::estfun(fit_null, level = 2, ranpar = "var")
    bread_mat <- sandwich::bread(fit_null, full = TRUE, information = "observed", ranpar = "var")
    a11_hat <- solve(bread_mat)
    psi_stage1 <- pack_psi_var(fit_null)

    return(list(
      s1_mat = unname(as.matrix(s1_mat)),
      a11_hat = unname(as.matrix(a11_hat)),
      psi_stage1 = psi_stage1,
      corrected_score_fn = corrected_slope_from_precomputed,
      psi_parameterization = "var"
    ))
  }

  if (identical(derivative_backend$name, "tmb")) {
    if (!requireNamespace("TMB", quietly = TRUE)) {
      stop("The `TMB` package is required for the `tmb` backend.")
    }

    tmb_obj <- make_tmb_stage1_object(
      cluster_objects = cluster_objects,
      fit_null = fit_null
    )

    h_nll <- get_tmb_stage1_hessian(tmb_obj)
    psi_stage1 <- pack_psi_var(fit_null)

    if (requireNamespace("merDeriv", quietly = TRUE)) {
      s1_mat <- sandwich::estfun(fit_null, level = 2, ranpar = "var")
      s1_mat <- unname(as.matrix(s1_mat))
    } else {
      s1_mat <- t(vapply(split_dat, function(df_i) {
        cluster_obj <- TMB::MakeADFun(
          data = make_tmb_stage1_data(prepare_cluster_objects(list(df_i))),
          parameters = list(
            beta0 = psi_stage1[[1]],
            beta1 = psi_stage1[[2]],
            var0 = psi_stage1[[3]],
            cov01 = psi_stage1[[4]],
            var1 = psi_stage1[[5]],
            sigma2 = psi_stage1[[6]]
          ),
          DLL = "tmb_stage1_gaussian",
          silent = TRUE
        )

        -cluster_obj$gr(cluster_obj$par)
      }, numeric(length(psi_stage1))))
    }

    return(list(
      s1_mat = s1_mat,
      a11_hat = h_nll / nrow(s1_mat),
      psi_stage1 = psi_stage1,
      corrected_score_fn = corrected_slope_from_precomputed,
      psi_parameterization = "var"
    ))
  }

  s1_mat <- t(vapply(cluster_objects, function(cluster_obj) {
    derivative_backend$gradient(
      fn = function(p) cluster_loglik_precomputed(cluster_obj, p, parameterization = "sd"),
      x = psi_hat
    )
  }, numeric(length(psi_hat))))

  h_total <- derivative_backend$hessian(
    fn = function(p) total_loglik_precomputed(cluster_objects, p, parameterization = "sd"),
    x = psi_hat
  )

  list(
    s1_mat = s1_mat,
    a11_hat = -h_total / nrow(s1_mat),
    psi_stage1 = psi_hat,
    corrected_score_fn = corrected_slope_from_precomputed,
    psi_parameterization = "sd"
  )
}

#' Compute the full stacked sandwich covariance for corrected scores.
#'
#' @details
#' This function stacks stage-1 estimating equations for the first-stage
#' Gaussian mixed model with stage-2 OLS estimating equations for corrected
#' random-slope scores. It returns HC0, HC1, HC2, and HC3 variants for the
#' stage-2 coefficient block.
#'
#' @param split_dat Named list of cluster data frames.
#' @param id_df Level-2 data frame containing `id` and `x`.
#' @param fit_null Fitted first-stage `lme4` model.
#' @param psi_hat Parameter vector from `pack_psi(fit_null)`.
#' @param derivative_backend Backend object from `make_derivative_backend()`.
#'
#' @return A list with `alpha_hat`, `corrected_scores`, and stage-2 covariance
#' matrices `vcov_hc0`, `vcov_hc1`, `vcov_hc2`, and `vcov_hc3`.
stacked_sandwich_for_corrected_scores <- function(split_dat, id_df, fit_null, psi_hat, derivative_backend) {
  n_id <- nrow(id_df)
  ordered_ids <- as.character(id_df$id)
  split_dat <- split_dat[ordered_ids]

  if (any(vapply(split_dat, is.null, logical(1)))) {
    stop("`split_dat` does not contain one cluster per subject in `id_df`.")
  }

  cluster_objects <- prepare_cluster_objects(split_dat)
  corrected_scores <- vapply(
    cluster_objects,
    corrected_slope_from_precomputed,
    numeric(1),
    psi_vec = psi_hat,
    parameterization = "sd"
  )
  stage1_parts <- get_stage1_sandwich_inputs(
    fit_null = fit_null,
    split_dat = split_dat,
    cluster_objects = cluster_objects,
    psi_hat = psi_hat,
    derivative_backend = derivative_backend
  )
  corrected_fn <- stage1_parts$corrected_score_fn

  if (identical(derivative_backend$name, "merDeriv")) {
    corrected_scores <- vapply(
      cluster_objects,
      corrected_fn,
      numeric(1),
      psi_vec = stage1_parts$psi_stage1,
      parameterization = stage1_parts$psi_parameterization
    )
  }

  x_stage2 <- cbind(1, id_df$x)
  alpha_hat <- drop(solve(crossprod(x_stage2), crossprod(x_stage2, corrected_scores)))
  resid_stage2 <- corrected_scores - drop(x_stage2 %*% alpha_hat)

  s1_mat <- stage1_parts$s1_mat
  s2_mat <- resid_stage2 * x_stage2
  g_mat <- cbind(s1_mat, s2_mat)
  b_hat <- crossprod(g_mat) / n_id

  a11_hat <- stage1_parts$a11_hat
  dh_mat <- t(vapply(cluster_objects, function(cluster_obj) {
    c(cluster_obj$dh_beta, rep(0, length(stage1_parts$psi_stage1) - 2L))
  }, numeric(length(stage1_parts$psi_stage1))))

  a21_hat <- -crossprod(x_stage2, dh_mat) / n_id
  a22_hat <- crossprod(x_stage2) / n_id
  a_hat <- rbind(
    cbind(a11_hat, matrix(0, nrow = nrow(a11_hat), ncol = ncol(a22_hat))),
    cbind(a21_hat, a22_hat)
  )

  vcov_hc0 <- tryCatch(
    solve(a_hat, b_hat) %*% t(solve(a_hat)) / n_id,
    error = function(e) matrix(NA_real_, nrow(a_hat), ncol(a_hat))
  )
  p_eta <- ncol(a_hat)
  hc1_scale <- n_id / max(n_id - p_eta, 1)
  vcov_hc1 <- vcov_hc0 * hc1_scale

  H_stage2 <- x_stage2 %*% solve(crossprod(x_stage2)) %*% t(x_stage2)
  h_ii <- pmin(diag(H_stage2), 0.999)

  s2_mat_hc2 <- (resid_stage2 / sqrt(1 - h_ii)) * x_stage2
  g_mat_hc2 <- cbind(s1_mat, s2_mat_hc2)
  b_hat_hc2 <- crossprod(g_mat_hc2) / n_id
  vcov_hc2 <- tryCatch(
    solve(a_hat, b_hat_hc2) %*% t(solve(a_hat)) / n_id,
    error = function(e) matrix(NA_real_, nrow(a_hat), ncol(a_hat))
  )

  s2_mat_hc3 <- (resid_stage2 / (1 - h_ii)) * x_stage2
  g_mat_hc3 <- cbind(s1_mat, s2_mat_hc3)
  b_hat_hc3 <- crossprod(g_mat_hc3) / n_id
  vcov_hc3 <- tryCatch(
    solve(a_hat, b_hat_hc3) %*% t(solve(a_hat)) / n_id,
    error = function(e) matrix(NA_real_, nrow(a_hat), ncol(a_hat))
  )

  alpha_slice <- (length(psi_hat) + 1L):(length(psi_hat) + 2L)

  list(
    alpha_hat = alpha_hat,
    corrected_scores = corrected_scores,
    vcov_hc0 = vcov_hc0[alpha_slice, alpha_slice, drop = FALSE],
    vcov_hc1 = vcov_hc1[alpha_slice, alpha_slice, drop = FALSE],
    vcov_hc2 = vcov_hc2[alpha_slice, alpha_slice, drop = FALSE],
    vcov_hc3 = vcov_hc3[alpha_slice, alpha_slice, drop = FALSE]
  )
}
