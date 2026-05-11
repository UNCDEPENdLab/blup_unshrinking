#' Legacy sandwich derivative helpers.
#'
#' These helpers were removed from `mlm_random_slope_blup_sandwich_coverage_sim.R`
#' because the active simulation path uses the derivative backend factory and
#' precomputed cluster objects instead. They are archived here for reference.

#' Calculate numerical gradient using central finite differences.
#' @param fn Function to differentiate.
#' @param x Parameter vector at which to evaluate the gradient.
#' @param eps Step size multiplier.
#' @return A numeric vector of the same length as x.
numeric_gradient <- function(fn, x, eps = 1e-6) {
  step <- eps * pmax(1, abs(x))
  out <- numeric(length(x))
  for (j in seq_along(x)) {
    x_plus <- x
    x_minus <- x
    x_plus[j] <- x_plus[j] + step[j]
    x_minus[j] <- x_minus[j] - step[j]
    out[j] <- (fn(x_plus) - fn(x_minus)) / (2 * step[j])
  }
  out
}

#' Calculate numerical Hessian using central finite differences.
#' @param fn Function to differentiate.
#' @param x Parameter vector at which to evaluate the Hessian.
#' @param eps Step size multiplier.
#' @return A symmetric numeric matrix of size length(x) by length(x).
numeric_hessian <- function(fn, x, eps = 1e-4) {
  step <- eps * pmax(1, abs(x))
  p <- length(x)
  out <- matrix(0, p, p)
  f0 <- fn(x)

  for (j in seq_len(p)) {
    x_plus <- x
    x_minus <- x
    x_plus[j] <- x_plus[j] + step[j]
    x_minus[j] <- x_minus[j] - step[j]
    out[j, j] <- (fn(x_plus) - 2 * f0 + fn(x_minus)) / (step[j]^2)
  }

  if (p > 1L) {
    for (j in seq_len(p - 1L)) {
      for (k in (j + 1L):p) {
        x_pp <- x
        x_pm <- x
        x_mp <- x
        x_mm <- x
        x_pp[c(j, k)] <- x_pp[c(j, k)] + c(step[j], step[k])
        x_pm[j] <- x_pm[j] + step[j]
        x_pm[k] <- x_pm[k] - step[k]
        x_mp[j] <- x_mp[j] - step[j]
        x_mp[k] <- x_mp[k] + step[k]
        x_mm[c(j, k)] <- x_mm[c(j, k)] - c(step[j], step[k])
        val <- (fn(x_pp) - fn(x_pm) - fn(x_mp) + fn(x_mm)) / (4 * step[j] * step[k])
        out[j, k] <- val
        out[k, j] <- val
      }
    }
  }

  out
}

#' Calculate the log-likelihood for a single cluster from a raw data frame.
#' @param cluster_df Data frame containing the observations for one cluster.
#' @param psi_vec The parameter vector on the unconstrained scale.
#' @return The scalar log-likelihood.
cluster_loglik <- function(cluster_df, psi_vec) {
  pars <- unpack_psi(psi_vec)
  x_mat <- model.matrix(~z, data = cluster_df)
  z_mat <- x_mat
  y_vec <- cluster_df$y
  resid_vec <- y_vec - drop(x_mat %*% pars$beta)
  vy <- z_mat %*% pars$G %*% t(z_mat) + (pars$sigma^2) * diag(nrow(cluster_df))
  chol_vy <- tryCatch(chol(vy), error = function(e) NULL)
  if (is.null(chol_vy)) {
    return(-Inf)
  }
  log_det <- 2 * sum(log(diag(chol_vy)))
  quad <- sum(backsolve(chol_vy, resid_vec, transpose = TRUE)^2)
  -0.5 * (nrow(cluster_df) * log(2 * pi) + log_det + quad)
}

#' Calculate the total log-likelihood across raw cluster data frames.
#' @param split_dat A list of cluster data frames.
#' @param psi_vec The parameter vector on the unconstrained scale.
#' @return The scalar total log-likelihood.
total_loglik <- function(split_dat, psi_vec) {
  sum(vapply(split_dat, cluster_loglik, numeric(1), psi_vec = psi_vec))
}

corrected_slope_from_psi <- function(cluster_df, psi_vec) {
  pars <- unpack_psi(psi_vec)
  x_mat <- model.matrix(~z, data = cluster_df)
  z_mat <- x_mat
  y_vec <- cluster_df$y
  resid_vec <- y_vec - drop(x_mat %*% pars$beta)
  # In this Gaussian model with R_i = sigma^2 I, sigma cancels from the likelihood-only score.
  bls <- solve(crossprod(z_mat), crossprod(z_mat, resid_vec))
  unname(bls[2])
}

corrected_slope_from_psi_var <- function(cluster_df, psi_vec) {
  pars <- unpack_psi_var(psi_vec)
  x_mat <- model.matrix(~z, data = cluster_df)
  z_mat <- x_mat
  y_vec <- cluster_df$y
  resid_vec <- y_vec - drop(x_mat %*% pars$beta)
  # For R_i = sigma^2 I, sigma^2 cancels from the GLS normal equations.
  bls <- solve(crossprod(z_mat), crossprod(z_mat, resid_vec))
  unname(bls[2])
}

corrected_component_gradient <- function(cluster_df, parameterization = c("sd", "var"), component = 2L) {
  parameterization <- match.arg(parameterization)
  x_mat <- model.matrix(~z, data = cluster_df)
  z_mat <- x_mat
  selector <- rep(0, ncol(z_mat))
  selector[component] <- 1

  # For the Gaussian likelihood-only score, d \tilde{b}_i / d beta' =
  # -e_component' (Z_i' Z_i)^{-1} Z_i' X_i. In this simulation X_i = Z_i.
  beta_grad <- -drop(selector %*% solve(crossprod(z_mat), crossprod(z_mat, x_mat)))
  c(beta_grad, rep(0, 4L))
}

corrected_component_gradient_precomputed <- function(cluster_obj, parameterization = c("sd", "var"), component = 2L) {
  parameterization <- match.arg(parameterization)
  selector <- rep(0, ncol(cluster_obj$Z))
  selector[component] <- 1
  beta_grad <- -drop(selector %*% solve(cluster_obj$ZtZ, cluster_obj$ZtX))
  c(beta_grad, rep(0, 4L))
}
