#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(lme4)
})

source(file.path("R", "blup_helpers.R"), local = TRUE)

set.seed(20260519)
n_id <- 8L
n_obs <- 6L
dat <- do.call(rbind, lapply(seq_len(n_id), function(i) {
  z <- scale(seq(-1, 1, length.out = n_obs), center = TRUE, scale = TRUE)[, 1]
  u0 <- rnorm(1L, sd = 0.45)
  u1 <- rnorm(1L, sd = 0.35)
  data.frame(
    id = as.character(i),
    z = z,
    y = 0.3 + 0.8 * z + u0 + u1 * z + rnorm(n_obs, sd = 0.4)
  )
}))

fit <- lmer(y ~ 1 + z + (1 + z | id), data = dat, REML = FALSE)

iid_out <- get_closed_form_corrected_scores(
  fit_obj = fit,
  data = dat,
  cluster_var = "id",
  outcome_var = "y",
  within_var = "z"
)
diag_R_list <- stats::setNames(
  lapply(split(dat, dat$id), function(df_i) stats::sigma(fit)^2 * diag(nrow(df_i))),
  names(split(dat, dat$id))
)
diag_gls_out <- get_closed_form_corrected_scores(
  fit_obj = fit,
  data = dat,
  cluster_var = "id",
  outcome_var = "y",
  within_var = "z",
  R_list = diag_R_list
)

df_1 <- dat[dat$id == "1", , drop = FALSE]
Z_1 <- cbind(1, df_1$z)
resid_1 <- df_1$y - drop(Z_1 %*% lme4::fixef(fit)[c("(Intercept)", "z")])
expected_ols <- as.numeric(solve(crossprod(Z_1), crossprod(Z_1, resid_1)))
expected_ols_vcov <- solve(crossprod(Z_1)) * stats::sigma(fit)^2

stopifnot(
  isTRUE(all.equal(iid_out$corrected_intercept_full[[1]], expected_ols[[1]], tolerance = 1e-10)),
  isTRUE(all.equal(iid_out$corrected_slope_full[[1]], expected_ols[[2]], tolerance = 1e-10)),
  isTRUE(all.equal(iid_out$ols_var11[[1]], expected_ols_vcov[1, 1], tolerance = 1e-10)),
  isTRUE(all.equal(iid_out$ols_var12[[1]], expected_ols_vcov[1, 2], tolerance = 1e-10)),
  isTRUE(all.equal(iid_out$ols_var22[[1]], expected_ols_vcov[2, 2], tolerance = 1e-10)),
  isTRUE(all.equal(iid_out$gls_var22[[1]], iid_out$ols_var22[[1]], tolerance = 1e-12)),
  isTRUE(all.equal(iid_out$corrected_intercept_full, diag_gls_out$corrected_intercept_full, tolerance = 1e-10)),
  isTRUE(all.equal(iid_out$corrected_slope_full, diag_gls_out$corrected_slope_full, tolerance = 1e-10)),
  isTRUE(all.equal(iid_out$ols_var11, diag_gls_out$ols_var11, tolerance = 1e-10)),
  isTRUE(all.equal(iid_out$ols_var12, diag_gls_out$ols_var12, tolerance = 1e-10)),
  isTRUE(all.equal(iid_out$ols_var22, diag_gls_out$ols_var22, tolerance = 1e-10))
)

rho <- 0.55
R_i <- outer(seq_len(n_obs), seq_len(n_obs), function(a, b) rho^abs(a - b))
R_list <- stats::setNames(rep(list(R_i), n_id), as.character(seq_len(n_id)))
gls_out <- get_closed_form_corrected_scores(
  fit_obj = fit,
  data = dat,
  cluster_var = "id",
  outcome_var = "y",
  within_var = "z",
  R_list = R_list
)

Rinv_Z <- solve(R_i, Z_1)
Rinv_resid <- solve(R_i, resid_1)
expected_gls_vcov <- solve(crossprod(Z_1, Rinv_Z))
expected_gls <- as.numeric(expected_gls_vcov %*% crossprod(Z_1, Rinv_resid))

stopifnot(
  isTRUE(all.equal(gls_out$corrected_intercept_full[[1]], expected_gls[[1]], tolerance = 1e-10)),
  isTRUE(all.equal(gls_out$corrected_slope_full[[1]], expected_gls[[2]], tolerance = 1e-10)),
  isTRUE(all.equal(gls_out$gls_var11[[1]], expected_gls_vcov[1, 1], tolerance = 1e-10)),
  isTRUE(all.equal(gls_out$gls_var12[[1]], expected_gls_vcov[1, 2], tolerance = 1e-10)),
  isTRUE(all.equal(gls_out$gls_var22[[1]], expected_gls_vcov[2, 2], tolerance = 1e-10)),
  isTRUE(all.equal(gls_out$ols_var22[[1]], gls_out$gls_var22[[1]], tolerance = 1e-12))
)

cat("Closed-form GLS corrected-score tests ok\n")
