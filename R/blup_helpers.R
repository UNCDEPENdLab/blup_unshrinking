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

#' Compute likelihood-only cluster scores directly from residualized cluster OLS.
#'
#' @details
#' This helper computes corrected cluster scores without using posterior BLUP
#' covariance matrices. It first subtracts the fixed-effect prediction from each
#' observation within a cluster, then solves the cluster-level least-squares
#' problem for the random-effect design. For a random-intercept-only model, this
#' is the cluster mean residual. For a random-intercept/random-slope model, this
#' is the intercept and within-cluster slope from regressing residuals on
#' `(1, within_var)` inside each cluster.
#'
#' The function is useful for balanced comparisons where the likelihood-only
#' cluster score has a simple closed form and should not depend on the
#' prior-unweighting matrix algebra used by `get_corrected_scores()`.
#'
#' @param fit_obj Fitted `lme4` model whose fixed effects are used to residualize
#' the outcome.
#' @param data Original long-format data frame used for the mixed model.
#' @param cluster_var Character scalar naming the cluster/grouping column.
#' @param outcome_var Character scalar naming the outcome column.
#' @param within_var Optional character scalar naming the within-cluster
#' predictor. If `NULL`, an intercept-only cluster score is computed.
#'
#' @return
#' A tibble with one row per cluster and an `id` column. Intercept-only calls
#' return `corrected` and `corrected_intercept_full`; calls with `within_var`
#' return `corrected_intercept_full` and `corrected_slope_full`. Singular
#' within-cluster least-squares systems return `NA_real_` values for that
#' cluster's corrected scores.
get_closed_form_corrected_scores <- function(fit_obj, data, cluster_var, outcome_var, within_var = NULL) {
  cluster_ids <- unique(as.character(data[[cluster_var]]))
  beta_hat <- lme4::fixef(fit_obj)

  purrr::map_dfr(cluster_ids, function(cluster_id) {
    df_i <- data[data[[cluster_var]] == cluster_id, , drop = FALSE]

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
    # likelihood-only coefficients from the random-effect design.
    resid_i <- df_i[[outcome_var]] - as.numeric(x_mat %*% beta_vec)
    z_crossprod <- crossprod(z_mat)
    corrected <- tryCatch(
      as.numeric(solve(z_crossprod, crossprod(z_mat, resid_i))),
      error = function(e) rep(NA_real_, ncol(z_mat))
    )
    
    # Also calculate the OLS sampling variance: (Z'Z)^-1 * sigma^2
    sigma2_hat <- stats::sigma(fit_obj)^2
    ols_vcov <- tryCatch(
      solve(z_crossprod) * sigma2_hat,
      error = function(e) matrix(NA_real_, nrow = ncol(z_mat), ncol = ncol(z_mat))
    )

    out <- tibble::tibble(id = cluster_id)
    if (ncol(z_mat) == 1L) {
      out$corrected <- corrected[[1]]
      out$corrected_intercept_full <- corrected[[1]]
      out$ols_var11 <- ols_vcov[1, 1]
    } else {
      out$corrected_intercept_full <- corrected[[1]]
      out$corrected_slope_full <- corrected[[2]]
      out$ols_var11 <- ols_vcov[1, 1]
      out$ols_var12 <- ols_vcov[1, 2]
      out$ols_var22 <- ols_vcov[2, 2]
    }
    out
  })
}
