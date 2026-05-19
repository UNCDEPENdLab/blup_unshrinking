#!/usr/bin/env Rscript

# Unit checks for shared Lai/OpenMx input construction helpers that are used by
# both the Lai replication and the sandwich coverage simulation.

suppressPackageStartupMessages({
  library(lme4)
})

source(file.path("R", "core_utils.R"), local = TRUE)
source(file.path("R", "lai_openmx_helpers.R"), local = TRUE)

cluster_df <- data.frame(z = c(-1, 0, 1))
stopifnot(
  isTRUE(all.equal(default_re_design(cluster_df), matrix(1, nrow = 3L, ncol = 1L))),
  isTRUE(all.equal(default_re_design(cluster_df, within_var = "z"), cbind(1, cluster_df$z)))
)

uni_row <- make_eb_output_row(
  id = "a",
  eb = 0.25,
  post_vcov = matrix(0.10, nrow = 1L),
  lambda = matrix(0.80, nrow = 1L),
  theta = matrix(0.20, nrow = 1L),
  prefix = "z_"
)
stopifnot(
  identical(names(uni_row), c("id", "z_u0_eb", "z_postvar11", "z_lambda11", "z_theta11")),
  isTRUE(all.equal(uni_row$z_u0_eb, 0.25))
)

bi_row <- make_eb_output_row(
  id = "b",
  eb = c(0.10, -0.20),
  post_vcov = matrix(c(0.40, 0.05, 0.05, 0.30), nrow = 2L),
  lambda = matrix(c(0.90, 0.10, 0.05, 0.85), nrow = 2L),
  theta = matrix(c(0.10, 0.02, 0.02, 0.15), nrow = 2L)
)
expected_bi_names <- c(
  "id", "u0_eb", "u1_eb", "postvar11", "postvar12", "postvar22",
  "lambda11", "lambda12", "lambda21", "lambda22", "theta11", "theta12", "theta22"
)
stopifnot(
  identical(names(bi_row), expected_bi_names),
  isTRUE(all.equal(bi_row$u1_eb, -0.20))
)

set.seed(2026)
n_id <- 10L
n_obs <- 5L
base_ids <- sprintf("cluster_%02d", seq_len(n_id))
stage2_ids <- rev(base_ids)
id_df <- data.frame(
  id = stage2_ids,
  x = seq(-1, 1, length.out = n_id)
)

sim_dat <- do.call(rbind, lapply(seq_along(base_ids), function(i) {
  z <- seq(-1, 1, length.out = n_obs)
  u0 <- rnorm(1, sd = 0.45)
  u1 <- rnorm(1, sd = 0.35)
  data.frame(
    id = base_ids[[i]],
    z = z,
    y = 0.2 + 0.6 * z + u0 + u1 * z + rnorm(n_obs, sd = 0.20)
  )
}))

fit <- lmer(y ~ 1 + z + (1 + z | id), data = sim_dat, REML = FALSE)
split_dat <- split(sim_dat, sim_dat$id)
out <- compute_lai_2spa_inputs(fit, split_dat, id_df)
diag_R_list <- stats::setNames(
  lapply(split_dat, function(df_i) stats::sigma(fit)^2 * diag(nrow(df_i))),
  names(split_dat)
)
out_diag_R <- compute_lai_2spa_inputs(fit, split_dat, id_df, R_list = diag_R_list)

required_cols <- c(
  "id", "x", "u0_eb", "u1_eb", "postvar11", "postvar12", "postvar22",
  "lambda11", "lambda12", "lambda21", "lambda22", "theta11", "theta12", "theta22"
)
missing_cols <- setdiff(required_cols, names(out))
if (length(missing_cols) > 0L) {
  stop("Missing expected Lai input columns: ", paste(missing_cols, collapse = ", "))
}

stopifnot(
  identical(out$id, id_df$id),
  isTRUE(all.equal(out$x, id_df$x)),
  all(is.finite(out$u0_eb)),
  all(is.finite(out$u1_eb)),
  all(is.finite(as.matrix(out[, setdiff(required_cols, c("id", "x")), drop = FALSE]))),
  isTRUE(all.equal(out$u0_eb, out_diag_R$u0_eb, tolerance = 1e-8)),
  isTRUE(all.equal(out$u1_eb, out_diag_R$u1_eb, tolerance = 1e-8)),
  isTRUE(all.equal(out$postvar11, out_diag_R$postvar11, tolerance = 1e-8)),
  isTRUE(all.equal(out$postvar12, out_diag_R$postvar12, tolerance = 1e-8)),
  isTRUE(all.equal(out$postvar22, out_diag_R$postvar22, tolerance = 1e-8)),
  isTRUE(all.equal(out$lambda11, out_diag_R$lambda11, tolerance = 1e-8)),
  isTRUE(all.equal(out$lambda12, out_diag_R$lambda12, tolerance = 1e-8)),
  isTRUE(all.equal(out$lambda21, out_diag_R$lambda21, tolerance = 1e-8)),
  isTRUE(all.equal(out$lambda22, out_diag_R$lambda22, tolerance = 1e-8)),
  isTRUE(all.equal(out$theta11, out_diag_R$theta11, tolerance = 1e-8)),
  isTRUE(all.equal(out$theta12, out_diag_R$theta12, tolerance = 1e-8)),
  isTRUE(all.equal(out$theta22, out_diag_R$theta22, tolerance = 1e-8))
)

rho_R <- 0.45
R_i <- stats::sigma(fit)^2 * outer(seq_len(n_obs), seq_len(n_obs), function(a, b) rho_R^abs(a - b))
ar_R_list <- stats::setNames(rep(list(R_i), n_id), base_ids)
out_ar_R <- compute_lai_2spa_inputs(fit, split_dat, id_df, R_list = ar_R_list)

first_id <- id_df$id[[1]]
df_1 <- split_dat[[first_id]]
Z_1 <- cbind(1, df_1$z)
G_hat <- as.matrix(lme4::VarCorr(fit)[["id"]])
beta_hat <- lme4::fixef(fit)
resid_1 <- df_1$y - drop(Z_1 %*% beta_hat[c("(Intercept)", "z")])
Sigma_1 <- Z_1 %*% G_hat %*% t(Z_1) + R_i
A_1 <- G_hat %*% t(Z_1) %*% solve(Sigma_1)
expected_eb_1 <- as.numeric(A_1 %*% resid_1)
expected_lambda_1 <- A_1 %*% Z_1
expected_theta_1 <- A_1 %*% R_i %*% t(A_1)
expected_post_1 <- solve(solve(G_hat) + crossprod(Z_1, solve(R_i, Z_1)))
row_1 <- out_ar_R[out_ar_R$id == first_id, , drop = FALSE]

stopifnot(
  isTRUE(all.equal(row_1$u0_eb, expected_eb_1[[1]], tolerance = 1e-8)),
  isTRUE(all.equal(row_1$u1_eb, expected_eb_1[[2]], tolerance = 1e-8)),
  isTRUE(all.equal(row_1$postvar11, expected_post_1[1, 1], tolerance = 1e-8)),
  isTRUE(all.equal(row_1$postvar12, expected_post_1[1, 2], tolerance = 1e-8)),
  isTRUE(all.equal(row_1$postvar22, expected_post_1[2, 2], tolerance = 1e-8)),
  isTRUE(all.equal(row_1$lambda11, expected_lambda_1[1, 1], tolerance = 1e-8)),
  isTRUE(all.equal(row_1$lambda12, expected_lambda_1[1, 2], tolerance = 1e-8)),
  isTRUE(all.equal(row_1$lambda21, expected_lambda_1[2, 1], tolerance = 1e-8)),
  isTRUE(all.equal(row_1$lambda22, expected_lambda_1[2, 2], tolerance = 1e-8)),
  isTRUE(all.equal(row_1$theta11, expected_theta_1[1, 1], tolerance = 1e-8)),
  isTRUE(all.equal(row_1$theta12, expected_theta_1[1, 2], tolerance = 1e-8)),
  isTRUE(all.equal(row_1$theta22, expected_theta_1[2, 2], tolerance = 1e-8))
)

cat("Lai/OpenMx input helper tests ok\n")
