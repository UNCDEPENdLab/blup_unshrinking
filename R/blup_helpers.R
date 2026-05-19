#' BLUP and corrected-score helpers.
#'
#' These helpers extract empirical Bayes random-effect predictions from
#' `lme4` models and construct corrected cluster scores used in downstream
#' stage-2 regressions. The main correction implemented here removes the
#' Gaussian prior contribution from the posterior BLUP mean using the posterior
#' and prior precision matrices.

#' Return a fallback value when an object is `NULL`.
#'
#' @details
#' This small infix helper is used for optional arguments where `NULL` means
#' "use the default discovered from the fitted model." Unlike `%in%`-style
#' vectorized helpers, it tests only whether the left-hand side object itself is
#' `NULL`.
#'
#' @param x Candidate value.
#' @param y Fallback value used when `x` is `NULL`.
#'
#' @return
#' `x` when it is not `NULL`; otherwise `y`.
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' Convert random-effect names into stable column-name fragments.
#'
#' @details
#' `lme4` random-effect names can contain punctuation such as `(Intercept)` or
#' operators from model formulas. This helper maps those names to simple
#' alphanumeric-plus-underscore fragments so score columns can be named
#' predictably, for example `blup_intercept` and `corrected_slope`.
#'
#' @param x Character vector of raw random-effect names, usually from
#' `colnames(lme4::ranef(fit)[[group]])`.
#'
#' @return
#' A character vector with `(Intercept)` mapped to `intercept`, non-alphanumeric
#' runs replaced by underscores, leading/trailing underscores stripped, and any
#' empty result replaced by `"re"`.
sanitize_re_name <- function(x) {
  # Preserve the conventional intercept label before removing punctuation.
  x <- gsub("^\\(Intercept\\)$", "intercept", x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  ifelse(nzchar(x), x, "re")
}

#' Unweight random effects from the prior using matrix inversion.
#'
#' @details
#' In the Gaussian mixed-model setup, each empirical Bayes random-effect
#' prediction is a posterior mean that combines the cluster likelihood with the
#' population random-effect prior. If
#' `V_post` is the conditional posterior covariance, `V_prior` is the
#' population random-effect covariance, `m_post` is the BLUP/posterior mean, and
#' `m_prior` is the prior mean, this helper solves
#'
#' `m_like = (V_post^{-1} - V_prior^{-1})^{-1}
#'           (V_post^{-1} m_post - V_prior^{-1} m_prior)`.
#'
#' The result is a likelihood-only or "unweighted" score estimate under this
#' normal-normal algebra. Singular or numerically invalid systems return an
#' all-`NA_real_` vector of the same length as `post_mean`, allowing simulation
#' loops to continue without special error handling.
#'
#' @param post_mean Numeric vector of posterior means, typically one row of
#' `lme4::ranef(fit, condVar = TRUE)[[group]]`.
#' @param post_vcov Numeric posterior covariance matrix for the same cluster.
#' @param prior_mean Numeric prior mean vector for the random effects, usually
#' zeros in `lme4` fits.
#' @param prior_vcov Numeric prior covariance matrix from
#' `lme4::VarCorr(fit)[[group]]`.
#' @param return_var Logical scalar indicating whether to return the variance of
#' the unweighted score. If `TRUE`, the function returns a list with `scores
#' = corrected score` and `var = (V_post^{-1} - V_prior^{-1})^{-1}`. If `FALSE`,
#' only the score is returned.
#'
#' @return
#' Numeric vector of corrected likelihood-only random-effect scores, or an
#' `NA_real_` vector if the precision subtraction or solve fails.
unweight_random_effects <- function(post_mean, post_vcov, prior_mean, prior_vcov, return_var = FALSE) {
  out <- tryCatch({
    # Subtract the prior precision contribution from the posterior precision
    # and solve for the score implied by the cluster likelihood alone.
    solve(
      solve(post_vcov) - solve(prior_vcov),
      solve(post_vcov, post_mean) - solve(prior_vcov, prior_mean)
    )
  }, error = function(e) {
    rep(NA_real_, length(post_mean))
  })
  scores <- as.numeric(out)

  if (return_var) {
    out_var <- tryCatch({
      solve(solve(post_vcov) - solve(prior_vcov))
    }, error = function(e) {
      matrix(NA_real_, nrow = length(post_mean), ncol = length(post_mean))
    })
    return(list(scores = scores, vars = diag(out_var)))
  } else {
    return(scores)
  }
}

#' Extract EB/BLUPs and Vig-style corrected scores from an lme4 fit.
#'
#' @details
#' Calls `lme4::ranef(fit_null, condVar = TRUE)` to obtain cluster-level BLUPs
#' and their conditional posterior covariance matrices, then applies
#' `unweight_random_effects()` cluster by cluster. The returned tibble uses
#' model-derived columns of the form `blup_<effect>` and `corrected_<effect>`,
#' where `<effect>` is sanitized by `sanitize_re_name()`. For example,
#' `(Intercept)` becomes `blup_intercept` / `corrected_intercept`, while a
#' random slope for `z` becomes `blup_z` / `corrected_z`.
#'
#' @param fit_null Fitted `lme4` mixed model with random effects available via
#' `ranef(..., condVar = TRUE)`.
#' @param group Optional character scalar naming the grouping factor to extract.
#' If `NULL`, the first grouping factor in the `ranef()` list is used.
#' @return
#' A tibble with one row per group level. The `id` column contains group-level
#' row names from `ranef()`, followed by BLUP and corrected-score columns.
get_corrected_scores <- function(fit_null, group = NULL) {
  re_list <- lme4::ranef(fit_null, condVar = TRUE)
  group <- group %||% names(re_list)[[1]]
  re_df <- re_list[[group]]

  # `postVar` is a p x p x n array of conditional posterior covariance
  # matrices aligned with the rows of the random-effects data frame.
  post_var_arr <- attr(re_df, "postVar")
  prior_vcov <- as.matrix(lme4::VarCorr(fit_null)[[group]])
  n_re <- ncol(re_df)
  prior_mean <- rep(0, n_re)
  re_names <- sanitize_re_name(colnames(re_df))

  purrr::map_dfr(seq_len(nrow(re_df)), function(i) {
    post_mean <- as.numeric(re_df[i, ])
    post_vcov <- post_var_arr[, , i]

    # For one-dimensional random-effect structures, subsetting can drop the
    # posterior covariance to a scalar; restore matrix shape for the algebra.
    if (!is.matrix(post_vcov)) {
      post_vcov <- matrix(post_vcov, nrow = 1L, ncol = 1L)
    }
    corrected <- unweight_random_effects(
      post_mean = post_mean,
      post_vcov = post_vcov,
      prior_mean = prior_mean,
      prior_vcov = prior_vcov,
      return_var = TRUE
    )

    out <- tibble::tibble(id = rownames(re_df)[[i]])

    # Always emit formula-derived names so this helper also works for models
    # whose random effects are not simply intercept/slope.
    for (j in seq_len(n_re)) {
      out[[paste0("blup_", re_names[[j]])]] <- post_mean[[j]]
      out[[paste0("corrected_", re_names[[j]])]] <- corrected$scores[[j]]
      out[[paste0("corrected_", re_names[[j]], "_var")]] <- corrected$vars[[j]]
    }

    out
  })
}

#' Extract diagonal-only corrected random-effect scores from an lme4 fit.
#'
#' @details
#' This helper is intentionally separate from `get_corrected_scores()` because
#' the diagonal-only correction is a diagnostic comparison, not the full
#' corrected-score estimator. It applies the scalar prior-unweighting formula to
#' each random-effect component independently and ignores posterior/prior
#' covariance among random effects.
#'
#' For a component `j`, the returned value is
#' `(V_post[j,j]^{-1} - V_prior[j,j]^{-1})^{-1}
#'  V_post[j,j]^{-1} m_post[j]`.
#'
#' @param fit_null Fitted `lme4` mixed model with random effects available via
#' `ranef(..., condVar = TRUE)`.
#' @param group Optional character scalar naming the grouping factor to extract.
#' If `NULL`, the first grouping factor in the `ranef()` list is used.
#'
#' @return
#' A tibble with one row per group level, an `id` column, and columns named
#' `corrected_<effect>_diag` for each random-effect component.
get_diagonal_corrected_scores <- function(fit_null, group = NULL) {
  re_list <- lme4::ranef(fit_null, condVar = TRUE)
  group <- group %||% names(re_list)[[1]]
  re_df <- re_list[[group]]
  post_var_arr <- attr(re_df, "postVar")
  prior_vcov <- as.matrix(lme4::VarCorr(fit_null)[[group]])
  n_re <- ncol(re_df)
  re_names <- sanitize_re_name(colnames(re_df))

  purrr::map_dfr(seq_len(nrow(re_df)), function(i) {
    post_mean <- as.numeric(re_df[i, ])
    post_vcov <- post_var_arr[, , i]
    if (!is.matrix(post_vcov)) {
      post_vcov <- matrix(post_vcov, nrow = 1L, ncol = 1L)
    }

    out <- tibble::tibble(id = rownames(re_df)[[i]])
    for (j in seq_len(n_re)) {
      diag_corrected <- tryCatch({
        denom <- (1 / post_vcov[j, j]) - (1 / prior_vcov[j, j])
        if (!is.finite(denom) || denom <= sqrt(.Machine$double.eps)) {
          NA_real_
        } else {
          (1 / denom) * (post_mean[[j]] / post_vcov[j, j])
        }
      }, error = function(e) {
        NA_real_
      })
      diag_var <- tryCatch({
        denom <- (1 / post_vcov[j, j]) - (1 / prior_vcov[j, j])
        if (!is.finite(denom) || denom <= sqrt(.Machine$double.eps)) {
          NA_real_
        } else {
          1 / denom
        }
      }, error = function(e) {
        NA_real_
      })

      out[[paste0("corrected_", re_names[[j]], "_diag")]] <- diag_corrected
      out[[paste0("corrected_", re_names[[j]], "_diag_var")]] <- diag_var
    }

    out
  })
}

#' Extract GLS-aware EB/BLUPs, posterior covariance, and MLE scores.
#'
#' @details
#' This helper is the residual-covariance-aware counterpart to
#' `get_corrected_scores()`. It does not use `lme4::ranef(..., condVar = TRUE)`
#' for the cluster predictions, because those conditional modes and posterior
#' covariance matrices assume the residual covariance structure fitted by
#' `lme4`. Instead, it rebuilds the Gaussian score ingredients directly from
#' plug-in first-stage parameters and each supplied cluster-level `R_i`.
#'
#' For cluster `i`, let `r_i = y_i - X_i beta_hat`, `Z_i` be the random-effect
#' design, `G` be the fitted random-effect covariance, `R_i` be the residual
#' covariance, and `Sigma_i = Z_i G Z_i' + R_i` be the marginal covariance.
#' The EB/BLUP scoring matrix in Vig and Hallquist (V&H) equation 29 is
#'
#' `A_Ei = G Z_i' Sigma_i^{-1}`,
#'
#' giving the EB/BLUP prediction
#'
#' `b_Ei = A_Ei r_i = G Z_i' Sigma_i^{-1} r_i`.
#'
#' The conditional posterior covariance of `b_i | y_i` is
#'
#' `V_post_i = (G^{-1} + Z_i' R_i^{-1} Z_i)^{-1}`.
#'
#' The MLE/Bartlett scoring matrix in V&H equation 30 is
#'
#' `A_Mi = (Z_i' R_i^{-1} Z_i)^{-1} Z_i' R_i^{-1}`,
#'
#' giving the likelihood-only score
#'
#' `b_Mi = A_Mi r_i`.
#'
#' The V&H equation 33 gives the same MLE score as a conversion from
#' the EB/BLUP score:
#'
#' `b_Mi = ((G Z_i' R_i^{-1} Z_i)^{-1} + I) b_Ei`.
#'
#' This implementation uses the equivalent posterior-minus-prior precision
#' identity already used by `unweight_random_effects()`:
#'
#' `b_Mi = (V_post_i^{-1} - G^{-1})^{-1} V_post_i^{-1} b_Ei`.
#'
#' The returned `corrected_*` columns are these matrix-converted MLE/Bartlett
#' scores. The returned `mle_*` columns are the direct GLS scores from equation
#' 30, included for diagnostics; the two should agree up to numerical
#' tolerance when `G` and `R_i` are nonsingular. The `corrected_*_var` columns
#' are the diagonal elements of `(Z_i' R_i^{-1} Z_i)^{-1}`, matching the
#' conditional variance expression in V&H equation 38.
#'
#' @param fit_obj Fitted `lme4` model supplying plug-in `beta_hat` and `G_hat`.
#' @param data Original long-format data frame.
#' @param cluster_var Character scalar naming the cluster/grouping column.
#' @param outcome_var Character scalar naming the outcome column.
#' @param within_var Optional character scalar naming the random-slope
#'   predictor. If `NULL`, an intercept-only random-effect design is used.
#' @param R_list Optional named list of cluster-level residual covariance
#'   matrices. Names should match `cluster_var` values. If unnamed, matrices are
#'   matched to clusters in data order. If `NULL`, iid covariance
#'   `sigma(fit_obj)^2 I` is used.
#' @param group Optional grouping-factor name for extracting `G_hat` and
#'   random-effect names from `fit_obj`. Defaults to the first grouping factor.
#'
#' @return
#' A tibble with one row per cluster. It includes EB/BLUP columns
#' `blup_<effect>`, posterior covariance columns `postvar*`, direct
#' MLE/Bartlett columns `mle_<effect>`, matrix-converted score columns
#' `corrected_<effect>`, and measurement-variance columns
#' `corrected_<effect>_var`.
get_gls_corrected_scores <- function(fit_obj, data, cluster_var, outcome_var, within_var = NULL,
                                     R_list = NULL, group = NULL) {
  cluster_ids <- unique(as.character(data[[cluster_var]]))
  split_dat <- split(data, as.character(data[[cluster_var]]), drop = TRUE)[cluster_ids]
  beta_hat <- lme4::fixef(fit_obj)
  group <- group %||% names(lme4::ranef(fit_obj))[[1]]
  g_hat <- as.matrix(lme4::VarCorr(fit_obj)[[group]])
  sigma2_hat <- stats::sigma(fit_obj)^2
  n_re <- if (is.null(within_var)) 1L else 2L

  re_names_raw <- tryCatch(colnames(lme4::ranef(fit_obj)[[group]]), error = function(e) NULL)
  if (is.null(re_names_raw) || length(re_names_raw) != n_re) {
    re_names_raw <- if (is.null(within_var)) "(Intercept)" else c("(Intercept)", within_var)
  }
  re_names <- sanitize_re_name(re_names_raw)

  if (!is.null(R_list) && !is.list(R_list)) {
    stop("`R_list` must be NULL or a list of cluster-level residual covariance matrices.")
  }
  if (!is.null(R_list) && is.null(names(R_list))) {
    if (length(R_list) != length(cluster_ids)) {
      stop("Unnamed `R_list` must have one matrix per cluster in `data` order.")
    }
    names(R_list) <- cluster_ids
  }

  if (nrow(g_hat) != n_re || ncol(g_hat) != n_re) {
    stop("The fitted random-effect covariance dimension does not match `within_var`.")
  }

  prior_mean <- rep(0, n_re)
  g_inv <- tryCatch(solve(g_hat), error = function(e) NULL)

  purrr::map_dfr(cluster_ids, function(cluster_id) {
    df_i <- split_dat[[cluster_id]]

    if (is.null(within_var)) {
      z_mat <- matrix(1, nrow = nrow(df_i), ncol = 1L)
      x_mat <- z_mat
      beta_vec <- beta_hat[[1]]
    } else {
      z_vec <- df_i[[within_var]]
      z_mat <- cbind(1, z_vec)
      x_mat <- z_mat
      beta_vec <- c(beta_hat[[1]], beta_hat[[within_var]])
    }

    resid_i <- df_i[[outcome_var]] - as.numeric(x_mat %*% beta_vec)
    R_i <- if (is.null(R_list)) {
      sigma2_hat * diag(nrow(df_i))
    } else {
      as.matrix(R_list[[cluster_id]])
    }

    pieces <- tryCatch({
      if (is.null(g_inv)) {
        stop("G is singular.")
      }
      if (!is.matrix(R_i) || nrow(R_i) != nrow(df_i) || ncol(R_i) != nrow(df_i)) {
        stop("R_i has incompatible dimensions.")
      }

      # EB/BLUP from V&H equations 29 and 31:
      # b_Ei = G Z_i' Sigma_i^-1 r_i, Sigma_i = Z_i G Z_i' + R_i.
      sigma_y_i <- z_mat %*% g_hat %*% t(z_mat) + R_i
      sigma_inv_resid <- solve(sigma_y_i, resid_i)
      eb <- as.numeric(g_hat %*% crossprod(z_mat, sigma_inv_resid))

      # Posterior covariance and MLE/Bartlett information. Avoid explicitly
      # forming R_i^-1; solve(R_i, Z_i) and solve(R_i, r_i) provide the needed
      # products for Z_i' R_i^-1 Z_i and Z_i' R_i^-1 r_i.
      R_inv_Z <- solve(R_i, z_mat)
      R_inv_resid <- solve(R_i, resid_i)
      info_like <- crossprod(z_mat, R_inv_Z)
      post_vcov <- solve(g_inv + info_like)

      # Direct MLE/Bartlett score from V&H equation 30.
      mle_vcov <- solve(info_like)
      mle <- as.numeric(mle_vcov %*% crossprod(z_mat, R_inv_resid))

      # Matrix conversion from EB/BLUP to MLE/Bartlett. This is the precision subtraction form of V&H equation 33.
      corrected <- unweight_random_effects(
        post_mean = eb,
        post_vcov = post_vcov,
        prior_mean = prior_mean,
        prior_vcov = g_hat,
        return_var = TRUE
      )

      list(
        eb = eb,
        post_vcov = post_vcov,
        mle = mle,
        mle_vcov = mle_vcov,
        corrected = corrected$scores,
        corrected_vars = corrected$vars
      )
    }, error = function(e) {
      list(
        eb = rep(NA_real_, n_re),
        post_vcov = matrix(NA_real_, nrow = n_re, ncol = n_re),
        mle = rep(NA_real_, n_re),
        mle_vcov = matrix(NA_real_, nrow = n_re, ncol = n_re),
        corrected = rep(NA_real_, n_re),
        corrected_vars = rep(NA_real_, n_re)
      )
    })

    out <- tibble::tibble(id = cluster_id)
    for (j in seq_len(n_re)) {
      out[[paste0("blup_", re_names[[j]])]] <- pieces$eb[[j]]
      out[[paste0("mle_", re_names[[j]])]] <- pieces$mle[[j]]
      out[[paste0("corrected_", re_names[[j]])]] <- pieces$corrected[[j]]
      out[[paste0("corrected_", re_names[[j]], "_var")]] <- pieces$corrected_vars[[j]]
      out[[paste0("mle_", re_names[[j]], "_var")]] <- pieces$mle_vcov[j, j]
    }

    if (n_re == 1L) {
      out$postvar11 <- pieces$post_vcov[1, 1]
      out$mle_var11 <- pieces$mle_vcov[1, 1]
    } else {
      out$postvar11 <- pieces$post_vcov[1, 1]
      out$postvar12 <- pieces$post_vcov[1, 2]
      out$postvar22 <- pieces$post_vcov[2, 2]
      out$mle_var11 <- pieces$mle_vcov[1, 1]
      out$mle_var12 <- pieces$mle_vcov[1, 2]
      out$mle_var22 <- pieces$mle_vcov[2, 2]
    }

    out
  })
}

#' Compute likelihood-only cluster scores directly from residualized cluster GLS.
#'
#' @details
#' This helper computes likelihood-only random-effect scores without using
#' posterior BLUP covariance matrices. For cluster `i`, let `r_i = y_i - X_i
#' beta` be the fixed-effect residual, `Z_i` be the random-effect design, and
#' `R_i` be the level-1 residual covariance matrix. The MLE/Bartlett score is
#' the GLS coefficient from regressing `r_i` on `Z_i`:
#'
#' `b_Mi = (Z_i' R_i^{-1} Z_i)^{-1} Z_i' R_i^{-1} (y_i - X_i beta)`.
#'
#' Its conditional sampling covariance is `(Z_i' R_i^{-1} Z_i)^{-1}`. When
#' `R_list` is `NULL`, `R_i` defaults to `sigma^2 I`, the `sigma^2` factor
#' cancels from the coefficient equation, and the score reduces to the original
#' OLS closed form `(Z_i'Z_i)^{-1} Z_i'r_i`; the covariance becomes
#' `sigma^2 (Z_i'Z_i)^{-1}`.
#'
#' The iid path is performance-specialized because it is the common simulation
#' case. It batches clusters by identical `Z_i`, solves `(Z_i'Z_i)^{-1}` once
#' per unique design, and multiplies that projection into all residual vectors
#' for the design group at once. The explicit-`R_list` path keeps the direct GLS
#' formula cluster by cluster because arbitrary non-diagonal `R_i` generally
#' differs across clusters.
#'
#' For a random-intercept-only model, the returned score is the GLS intercept of
#' the fixed-effect residuals. For a random-intercept/random-slope model, the
#' returned scores are the GLS intercept and within-cluster slope from residuals
#' on `(1, within_var)` inside each cluster.
#'
#' The function is useful for balanced comparisons where the likelihood-only
#' cluster score has a direct closed form and should not depend on the
#' prior-unweighting matrix algebra used by `get_corrected_scores()`. The
#' returned `ols_var*` columns are retained as compatibility aliases for the
#' measurement-error estimators; under non-diagonal `R_i` they contain GLS/MLE
#' score variances.
#'
#' @param fit_obj Fitted `lme4` model whose fixed effects are used to residualize
#' the outcome.
#' @param data Original long-format data frame used for the mixed model.
#' @param cluster_var Character scalar naming the cluster/grouping column.
#' @param outcome_var Character scalar naming the outcome column.
#' @param within_var Optional character scalar naming the within-cluster
#' predictor. If `NULL`, an intercept-only cluster score is computed.
#' @param R_list Optional named list of cluster-level residual covariance
#' matrices. Names should match `cluster_var` values. If unnamed, matrices are
#' matched to clusters in the order they appear in `data`. If `NULL`, iid
#' residual covariance `sigma(fit_obj)^2 I` is used for each cluster.
#'
#' @return
#' A tibble with one row per cluster and an `id` column. Intercept-only calls
#' return `corrected` and `corrected_intercept_full`; calls with `within_var`
#' return `corrected_intercept_full` and `corrected_slope_full`. Singular or
#' numerically invalid cluster systems return `NA_real_` values for that
#' cluster's corrected scores and variances.
get_closed_form_corrected_scores <- function(fit_obj, data, cluster_var, outcome_var, within_var = NULL,
                                             R_list = NULL) {
  cluster_ids <- unique(as.character(data[[cluster_var]]))
  split_dat <- split(data, as.character(data[[cluster_var]]), drop = TRUE)[cluster_ids]
  beta_hat <- lme4::fixef(fit_obj)
  sigma2_hat <- stats::sigma(fit_obj)^2

  if (!is.null(R_list) && !is.list(R_list)) {
    stop("`R_list` must be NULL or a list of cluster-level residual covariance matrices.")
  }
  if (!is.null(R_list) && is.null(names(R_list))) {
    if (length(R_list) != length(cluster_ids)) {
      stop("Unnamed `R_list` must have one matrix per cluster in `data` order.")
    }
    names(R_list) <- cluster_ids
  }

  if (is.null(R_list)) {
    # Fast iid branch: R_i = sigma^2 I. The scalar sigma^2 does not affect the
    # coefficient estimate, so the score is OLS on fixed-effect residuals. We
    # still use sigma^2 below for the score covariance.
    id_vec <- as.character(data[[cluster_var]])

    if (is.null(within_var)) {
      # Random intercept only:
      #   Z_i = 1_i
      #   r_i = y_i - beta_0
      resid_all <- data[[outcome_var]] - beta_hat[[1]]
      resid_by_id <- split(resid_all, id_vec, drop = TRUE)[cluster_ids]
      z_by_id <- lapply(resid_by_id, function(resid_i) matrix(1, nrow = length(resid_i), ncol = 1L))
      n_re <- 1L
    } else {
      # Random intercept and slope:
      #   Z_i = X_i = [1, z_i]
      #   r_i = y_i - beta_0 - beta_z z_i
      # This assumes the fixed-effect and random-effect designs match, as in
      # the current random-slope simulations.
      z_all <- data[[within_var]]
      resid_all <- data[[outcome_var]] - beta_hat[[1]] - beta_hat[[within_var]] * z_all
      resid_by_id <- split(resid_all, id_vec, drop = TRUE)[cluster_ids]
      z_vec_by_id <- split(z_all, id_vec, drop = TRUE)[cluster_ids]
      z_by_id <- lapply(z_vec_by_id, function(z_vec) cbind(1, z_vec))
      n_re <- 2L
    }

    # Bucket clusters by exact Z design. Balanced simulations usually have one
    # bucket; unbalanced simulations often have one bucket per realized trial
    # count because z is generated deterministically from the count.
    design_buckets <- new.env(parent = emptyenv())
    for (i in seq_along(cluster_ids)) {
      z_mat <- z_by_id[[i]]
      z_values <- as.numeric(z_mat)
      # First partition by dimensions to avoid comparing vectors of different
      # lengths, then compare the full numeric design for exact reuse.
      bucket_key <- paste(nrow(z_mat), ncol(z_mat), sep = "|")
      bucket <- if (exists(bucket_key, envir = design_buckets, inherits = FALSE)) {
        get(bucket_key, envir = design_buckets, inherits = FALSE)
      } else {
        list()
      }

      matched <- FALSE
      if (length(bucket) > 0L) {
        for (j in seq_along(bucket)) {
          if (identical(z_values, bucket[[j]]$z_values)) {
            bucket[[j]]$indices <- c(bucket[[j]]$indices, i)
            matched <- TRUE
            break
          }
        }
      }
      if (!matched) {
        bucket[[length(bucket) + 1L]] <- list(
          z_values = z_values,
          z_mat = z_mat,
          indices = i
        )
      }
      assign(bucket_key, bucket, envir = design_buckets)
    }

    corrected <- matrix(NA_real_, nrow = n_re, ncol = length(cluster_ids))
    vcov_arr <- array(NA_real_, dim = c(n_re, n_re, length(cluster_ids)))

    # For each unique design:
    #   P_Z = (Z'Z)^-1 Z'
    #   B_group = P_Z R_group
    # where R_group is the matrix of residual vectors for all clusters sharing Z.
    for (bucket_key in ls(design_buckets, all.names = TRUE)) {
      bucket <- get(bucket_key, envir = design_buckets, inherits = FALSE)
      for (entry in bucket) {
        idx <- entry$indices
        z_mat <- entry$z_mat
        ztz_inv <- tryCatch(solve(crossprod(z_mat)), error = function(e) NULL)
        if (is.null(ztz_inv)) {
          next
        }

        resid_mat <- do.call(cbind, resid_by_id[idx])
        corrected[, idx] <- ztz_inv %*% crossprod(z_mat, resid_mat)
        vcov_i <- ztz_inv * sigma2_hat
        # Every cluster in this design group has the same iid score covariance.
        for (cluster_pos in idx) {
          vcov_arr[, , cluster_pos] <- vcov_i
        }
      }
    }

    out <- tibble::tibble(id = cluster_ids)
    if (n_re == 1L) {
      out$corrected <- corrected[1, ]
      out$corrected_intercept_full <- corrected[1, ]
      out$gls_var11 <- vcov_arr[1, 1, ]
      out$ols_var11 <- vcov_arr[1, 1, ]
    } else {
      out$corrected_intercept_full <- corrected[1, ]
      out$corrected_slope_full <- corrected[2, ]
      out$gls_var11 <- vcov_arr[1, 1, ]
      out$gls_var12 <- vcov_arr[1, 2, ]
      out$gls_var22 <- vcov_arr[2, 2, ]
      out$ols_var11 <- vcov_arr[1, 1, ]
      out$ols_var12 <- vcov_arr[1, 2, ]
      out$ols_var22 <- vcov_arr[2, 2, ]
    }
    return(out)
  }

  # General GLS branch for explicit R_i. This is the path used for AR(1),
  # Toeplitz, or any other non-diagonal residual covariance supplied by the
  # simulator. It intentionally mirrors the formula in the roxygen block.
  purrr::map_dfr(cluster_ids, function(cluster_id) {
    df_i <- split_dat[[cluster_id]]

    if (is.null(within_var)) {
      # Random-intercept-only case: the cluster score is the intercept from a
      # regression of fixed-effect residuals on a column of ones.
      z_mat <- matrix(1, nrow = nrow(df_i), ncol = 1L)
      x_mat <- z_mat
      beta_vec <- beta_hat[[1]]
    } else {
      # Random intercept and slope case: use the same fixed-effect design for
      # residualization, then solve the within-cluster random-effect design.
      z_vec <- df_i[[within_var]]
      z_mat <- cbind(1, z_vec)
      x_mat <- z_mat
      beta_vec <- c(beta_hat[[1]], beta_hat[[within_var]])
    }

    # Residualize using fixed effects only, then estimate the cluster-specific
    # likelihood-only coefficients from the random-effect design. With the
    # default iid R_i this is exactly the original OLS closed form; with
    # non-diagonal R_i it is the MLE/Bartlett GLS score.
    resid_i <- df_i[[outcome_var]] - as.numeric(x_mat %*% beta_vec)

    R_i <- as.matrix(R_list[[cluster_id]])

    gls <- tryCatch({
      if (!is.matrix(R_i) || nrow(R_i) != nrow(df_i) || ncol(R_i) != nrow(df_i)) {
        stop("R_i has incompatible dimensions.")
      }
      # Avoid forming R_i^-1 explicitly. solve(R_i, Z_i) and solve(R_i, r_i)
      # supply R_i^-1 Z_i and R_i^-1 r_i for the GLS normal equations.
      R_inv_Z <- solve(R_i, z_mat)
      R_inv_resid <- solve(R_i, resid_i)
      z_Rinv_z <- crossprod(z_mat, R_inv_Z)
      gls_vcov <- solve(z_Rinv_z)
      list(
        corrected = as.numeric(gls_vcov %*% crossprod(z_mat, R_inv_resid)),
        vcov = gls_vcov
      )
    }, error = function(e) {
      list(
        corrected = rep(NA_real_, ncol(z_mat)),
        vcov = matrix(NA_real_, nrow = ncol(z_mat), ncol = ncol(z_mat))
      )
    })

    out <- tibble::tibble(id = cluster_id)
    if (ncol(z_mat) == 1L) {
      out$corrected <- gls$corrected[[1]]
      out$corrected_intercept_full <- gls$corrected[[1]]
      out$gls_var11 <- gls$vcov[1, 1]
      out$ols_var11 <- gls$vcov[1, 1]
    } else {
      out$corrected_intercept_full <- gls$corrected[[1]]
      out$corrected_slope_full <- gls$corrected[[2]]
      out$gls_var11 <- gls$vcov[1, 1]
      out$gls_var12 <- gls$vcov[1, 2]
      out$gls_var22 <- gls$vcov[2, 2]
      out$ols_var11 <- gls$vcov[1, 1]
      out$ols_var12 <- gls$vcov[1, 2]
      out$ols_var22 <- gls$vcov[2, 2]
    }
    out
  })
}
