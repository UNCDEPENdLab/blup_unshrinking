#' Lai Simulation Study 2: large unbalanced cluster sizes and heterogeneous X variance.

make_study2_cluster_sizes <- function(num_clus) {
  rep_len(
    c(seq(10, 75, length.out = 6), 80, 83, 86, 89, rep(91:100, each = 2)),
    length.out = num_clus
  )
}

simulate_study2 <- function(condition) {
  covu <- make_covu(condition)
  cluster_sizes <- make_study2_cluster_sizes(as.integer(condition$num_clus))
  cid <- rep(seq_len(as.integer(condition$num_clus)), cluster_sizes)
  x <- stats::rnorm(length(cid), sd = sqrt(cid / 20))

  u <- MASS::mvrnorm(as.integer(condition$num_clus), mu = c(0, 0), Sigma = covu)
  y <- fixed_params$gamma0 + fixed_params$gamma1 * x + rowSums(cbind(1, x) * u[cid, , drop = FALSE]) +
    draw_lai_level1_residuals(cluster_sizes, sigma = sqrt(condition$sigma2), condition = condition)
  ev_z <- 1 - drop(t(c(fixed_params$beta_zu0, condition$beta_zu1)) %*% covu %*% c(fixed_params$beta_zu0, condition$beta_zu1))
  z <- fixed_params$z_intercept + fixed_params$beta_zu0 * u[, 1] + condition$beta_zu1 * u[, 2] +
    stats::rnorm(as.integer(condition$num_clus), sd = sqrt(ev_z))

  list(
    lv1 = tibble::tibble(
      cid = factor(cid),
      cid_chr = as.character(cid),
      trial_index = ave(seq_along(cid), cid, FUN = seq_along),
      x = x,
      y = y
    ),
    lv2_true = tibble::tibble(
      id = as.character(seq_len(as.integer(condition$num_clus))),
      z = z,
      true_u0 = u[, 1],
      true_u1 = u[, 2]
    )
  )
}

run_study2_rep <- function(condition) {
  run_matched_outcome_rep(condition, simulate_study2(condition))
}
