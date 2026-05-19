#!/usr/bin/env Rscript

# Unit checks for reusable stacked-sandwich helper functions. These exercise the
# parameter packing, cluster precomputation, likelihood scale equivalence, and
# small end-to-end sandwich assembly used by the random-slope simulation.

suppressPackageStartupMessages({
  library(lme4)
})

source(file.path("R", "derivative_backends.R"), local = TRUE)
source(file.path("R", "stacked_sandwich_helpers.R"), local = TRUE)

set.seed(1001)

n_id <- 18L
n_obs <- 6L
id_df <- data.frame(
  id = sprintf("id%02d", seq_len(n_id)),
  x = seq(-1.4, 1.4, length.out = n_id)
)

sim_dat <- do.call(rbind, lapply(seq_len(n_id), function(i) {
  z <- seq(-1.25, 1.25, length.out = n_obs)
  u0 <- 0.45 * id_df$x[[i]] + rnorm(1, sd = 0.35)
  u1 <- -0.25 * id_df$x[[i]] + rnorm(1, sd = 0.30)
  data.frame(
    id = id_df$id[[i]],
    x = id_df$x[[i]],
    z = z,
    y = 0.6 + 0.75 * z + u0 + u1 * z + rnorm(n_obs, sd = 0.25)
  )
}))

fit <- lmer(y ~ 1 + z + (1 + z | id), data = sim_dat, REML = FALSE)
split_dat <- split(sim_dat, sim_dat$id)

psi_sd <- pack_psi(fit)
psi_var <- pack_psi_var(fit)
pars_sd <- unpack_psi(psi_sd)
pars_var <- unpack_psi_var(psi_var)
numeric_matrix <- function(mat) matrix(as.numeric(mat), nrow = nrow(mat), ncol = ncol(mat))

stopifnot(
  length(psi_sd) == 6L,
  length(psi_var) == 6L,
  isTRUE(all.equal(pars_sd$beta, unname(fixef(fit)[c("(Intercept)", "z")]), tolerance = 1e-10)),
  isTRUE(all.equal(pars_var$beta, unname(fixef(fit)[c("(Intercept)", "z")]), tolerance = 1e-10)),
  isTRUE(all.equal(numeric_matrix(pars_sd$G), numeric_matrix(as.matrix(VarCorr(fit)$id)), tolerance = 1e-10)),
  isTRUE(all.equal(numeric_matrix(pars_var$G), numeric_matrix(as.matrix(VarCorr(fit)$id)), tolerance = 1e-10)),
  isTRUE(all.equal(pars_sd$sigma, sigma(fit), tolerance = 1e-10)),
  isTRUE(all.equal(pars_var$sigma, sigma(fit), tolerance = 1e-10))
)

cluster_objects <- prepare_cluster_objects(split_dat)
first_id <- names(split_dat)[[1]]
first_x <- model.matrix(~z, data = split_dat[[first_id]])
expected_ols <- drop(solve(crossprod(first_x), crossprod(first_x, split_dat[[first_id]]$y)))

stopifnot(
  identical(names(cluster_objects), names(split_dat)),
  isTRUE(all.equal(cluster_objects[[first_id]]$ols_coef, expected_ols, tolerance = 1e-12))
)

expected_corrected <- expected_ols[[2]] - fixef(fit)[["z"]]
stopifnot(
  isTRUE(all.equal(
    corrected_slope_from_precomputed(cluster_objects[[first_id]], psi_sd, parameterization = "sd"),
    expected_corrected,
    tolerance = 1e-12
  ))
)

ll_sd <- cluster_loglik_precomputed(cluster_objects[[first_id]], psi_sd, parameterization = "sd")
ll_var <- cluster_loglik_precomputed(cluster_objects[[first_id]], psi_var, parameterization = "var")
stopifnot(
  is.finite(ll_sd),
  is.finite(ll_var),
  isTRUE(all.equal(ll_sd, ll_var, tolerance = 1e-8))
)

sandwich_out <- stacked_sandwich_for_corrected_scores(
  split_dat = split_dat,
  id_df = id_df,
  fit_null = fit,
  psi_hat = psi_sd,
  derivative_backend = make_derivative_backend("handcoded")
)

stopifnot(
  length(sandwich_out$alpha_hat) == 2L,
  length(sandwich_out$corrected_scores) == n_id,
  all(names(sandwich_out$corrected_scores) == id_df$id)
)

for (vcov_name in paste0("vcov_hc", 0:3)) {
  vcov_mat <- sandwich_out[[vcov_name]]
  stopifnot(
    is.matrix(vcov_mat),
    identical(dim(vcov_mat), c(2L, 2L)),
    all(is.finite(vcov_mat))
  )
}

rho_R <- 0.4
R_i <- 0.25^2 * outer(seq_len(n_obs), seq_len(n_obs), function(a, b) rho_R^abs(a - b))
R_list <- stats::setNames(rep(list(R_i), n_id), id_df$id)
cluster_objects_R <- prepare_cluster_objects(split_dat, R_list = R_list)
first_R_inv_Z <- solve(R_i, first_x)
first_R_inv_y <- solve(R_i, split_dat[[first_id]]$y)
expected_gls_coef <- drop(solve(crossprod(first_x, first_R_inv_Z), crossprod(first_x, first_R_inv_y)))
expected_gls_score <- expected_gls_coef[[2]] - fixef(fit)[["z"]]

stopifnot(
  isTRUE(all.equal(cluster_objects_R[[first_id]]$ols_coef, expected_gls_coef, tolerance = 1e-10)),
  isTRUE(all.equal(
    corrected_slope_from_precomputed(
      cluster_objects_R[[first_id]],
      pack_psi_known_R(fit),
      parameterization = "known_R"
    ),
    expected_gls_score,
    tolerance = 1e-10
  ))
)

ll_known_R <- cluster_loglik_precomputed(
  cluster_objects_R[[first_id]],
  pack_psi_known_R(fit),
  parameterization = "known_R"
)
stopifnot(is.finite(ll_known_R))

sandwich_known_R <- stacked_sandwich_for_corrected_scores(
  split_dat = split_dat,
  id_df = id_df,
  fit_null = fit,
  psi_hat = psi_sd,
  derivative_backend = make_derivative_backend("handcoded"),
  R_list = R_list
)
known_R_pars <- unpack_psi_known_R(sandwich_known_R$psi_stage1)
expected_sandwich_gls_score <- expected_gls_coef[[2]] -
  drop(cluster_objects_R[[first_id]]$projection_X[2, ] %*% known_R_pars$beta)
known_R_score_sum <- colSums(t(vapply(cluster_objects_R, function(cluster_obj) {
  make_derivative_backend("handcoded")$gradient(
    fn = function(p) cluster_loglik_precomputed(cluster_obj, p, parameterization = "known_R"),
    x = sandwich_known_R$psi_stage1
  )
}, numeric(length(sandwich_known_R$psi_stage1)))))

stopifnot(
  length(sandwich_known_R$alpha_hat) == 2L,
  length(sandwich_known_R$corrected_scores) == n_id,
  all(names(sandwich_known_R$corrected_scores) == id_df$id),
  identical(sandwich_known_R$psi_parameterization, "known_R"),
  max(abs(known_R_score_sum)) < 1e-3,
  isTRUE(all.equal(
    sandwich_known_R$corrected_scores[[first_id]],
    expected_sandwich_gls_score,
    tolerance = 1e-8
  ))
)

for (vcov_name in paste0("vcov_hc", 0:3)) {
  vcov_mat <- sandwich_known_R[[vcov_name]]
  stopifnot(
    is.matrix(vcov_mat),
    identical(dim(vcov_mat), c(2L, 2L)),
    all(is.finite(vcov_mat))
  )
}

cat("stacked sandwich helper tests ok\n")
