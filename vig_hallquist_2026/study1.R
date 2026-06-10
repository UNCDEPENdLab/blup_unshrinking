#' Vig-Hallquist (2026) simulation study 1: BLUP as outcome
#' 
study1_methods <- function() {
  c(
    "oracle",
    "naive_blup",
    "closed_form",
    "fuller_closed_form",
    "fuller_alpha_stepdown_closed_form"
    "single_subject_ols",
    "lai_2spa",
    "direct_mlm"
  )
}

simulate_study1 <- function(condition) {
#   covu <- make_covu(condition)
#   cluster_sizes <- rep(as.integer(condition$clus_size), as.integer(condition$num_clus))
#   cid <- rep(seq_len(as.integer(condition$num_clus)), cluster_sizes)

#   xj <- seq(-1, 1, length.out = as.integer(condition$clus_size))
#   xj <- xj / sqrt(mean(xj^2))
#   x <- rep(xj, as.integer(condition$num_clus))

#   u <- MASS::mvrnorm(as.integer(condition$num_clus), mu = c(0, 0), Sigma = covu)
#   y <- fixed_params$gamma0 + fixed_params$gamma1 * x + rowSums(cbind(1, x) * u[cid, , drop = FALSE]) +
#     draw_lai_level1_residuals(cluster_sizes, sigma = sqrt(condition$sigma2), condition = condition)
#   ev_z <- 1 - drop(t(c(fixed_params$beta_zu0, condition$beta_zu1)) %*% covu %*% c(fixed_params$beta_zu0, condition$beta_zu1))
#   z <- fixed_params$z_intercept + fixed_params$beta_zu0 * u[, 1] + condition$beta_zu1 * u[, 2] +
#     stats::rnorm(as.integer(condition$num_clus), sd = sqrt(ev_z))

#   list(
#     lv1 = tibble::tibble(
#       cid = factor(cid),
#       cid_chr = as.character(cid),
#       trial_index = ave(seq_along(cid), cid, FUN = seq_along),
#       x = x,
#       y = y
#     ),
#     lv2_true = tibble::tibble(
#       id = as.character(seq_len(as.integer(condition$num_clus))),
#       z = z,
#       true_u0 = u[, 1],
#       true_u1 = u[, 2]
#     )
#   )
}

run_study1_rep <- function(condition) {
  # run_matched_outcome_rep(condition, simulate_study1(condition))
}
