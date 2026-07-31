# Lai & Liu Study 1 data generation.

#' Simulate exactly the balanced Study 1 DGM from Lai and Liu's supplement.
#'
#' In particular, `x` is standardized with `sqrt(mean(x^2))`, so its sum of
#' squares is the cluster size.  Do not replace this with the VH reliability
#' time design, which uses a sample-SD convention and has sum of squares n - 1.
simulate_lai_study1_vh <- function(condition) {
  n_cluster <- as.integer(condition$num_clus[[1]])
  n_within <- as.integer(condition$clus_size[[1]])
  var_u1 <- as.numeric(condition$var_u1[[1]])
  icc <- as.numeric(condition$icc[[1]])
  rho <- as.numeric(condition$cor_u0_u1[[1]])
  sigma2 <- as.numeric(condition$sigma2[[1]])

  covu <- lai_study1_vh_dgm_covariance(icc, var_u1, rho)
  cid <- rep(seq_len(n_cluster), each = n_within)
  xj <- lai_study1_vh_within_design(n_within)[, "x"]
  x <- rep(xj, n_cluster)

  u <- MASS::mvrnorm(n_cluster, mu = c(0, 0), Sigma = covu)
  y <- lai_study1_vh_fixed_params$gamma0 +
    lai_study1_vh_fixed_params$gamma1 * x +
    rowSums(cbind(1, x) * u[cid, , drop = FALSE]) +
    stats::rnorm(n_cluster * n_within, sd = sqrt(sigma2))

  beta <- c(lai_study1_vh_fixed_params$beta_zu0, condition$beta_zu1[[1]])
  z_residual_variance <- 1 - drop(t(beta) %*% covu %*% beta)
  if (!is.finite(z_residual_variance) || z_residual_variance <= 0) {
    stop("The Lai Study 1 external-outcome residual variance must be positive.")
  }
  z <- lai_study1_vh_fixed_params$z_intercept + u %*% beta +
    stats::rnorm(n_cluster, sd = sqrt(z_residual_variance))

  list(
    lv1 = tibble::tibble(
      cid = factor(cid),
      cid_chr = as.character(cid),
      trial_index = ave(seq_along(cid), cid, FUN = seq_along),
      x = x,
      y = y
    ),
    lv2_true = tibble::tibble(
      id = as.character(seq_len(n_cluster)),
      z = as.numeric(z),
      true_u0 = u[, "u0"],
      true_u1 = u[, "u1"]
    ),
    covu = covu,
    z_residual_variance = z_residual_variance
  )
}
