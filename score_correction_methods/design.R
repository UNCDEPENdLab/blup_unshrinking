fixed_params <- list(
  gamma0 = 0,
  gamma1 = 0.5,
  beta_zu0 = 0.4,
  z_intercept = 1.5
)

make_score_correction_study_design <- function(max_conditions = Inf) {
  design <- tidyr::crossing(
    num_clus = c(30L, 100L, 500L),
    clus_size = c(3L, 10L, 25L),
    icc = c(0.05, 0.20, 0.50),
    vr_u1_u0 = c(0.5, 1.0, 2.0),
    cor_u0_u1 = c(-0.5, 0.0, 0.5),
    beta_zu1 = c(0.0, 0.4),
    r_structure = "iid"
  ) %>%
    dplyr::mutate(
      sigma2 = 1 - icc,
      var_u1 = vr_u1_u0 * icc,
      design_source = "supplement_script",
      condition_note = "From Lai study 1"
    )

  if (!is.na(max_conditions) && max_conditions > 0L) {
    design <- design %>% dplyr::slice_head(n = max_conditions)
  }
  
  design %>%
    dplyr::mutate(condition_id = seq_len(dplyr::n()))
}


simulate_score_correction <- function(condition) {
    covu <- make_covu(condition)
  cluster_sizes <- rep(as.integer(condition$clus_size), as.integer(condition$num_clus))
  cid <- rep(seq_len(as.integer(condition$num_clus)), cluster_sizes)

  xj <- seq(-1, 1, length.out = as.integer(condition$clus_size))
  xj <- xj / sqrt(mean(xj^2))
  x <- rep(xj, as.integer(condition$num_clus))

  u <- MASS::mvrnorm(as.integer(condition$num_clus), mu = c(0, 0), Sigma = covu)
  y <- fixed_params$gamma0 + fixed_params$gamma1 * x + rowSums(cbind(1, x) * u[cid, , drop = FALSE]) +
    draw_score_correction_level1_residuals(cluster_sizes, sigma = sqrt(condition$sigma2), condition = condition)
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
