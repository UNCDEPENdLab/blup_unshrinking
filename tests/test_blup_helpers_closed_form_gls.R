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
    trial_index = seq_len(n_obs),
    z = z,
    y = 0.3 + 0.8 * z + u0 + u1 * z + rnorm(n_obs, sd = 0.4)
  )
}))

fit <- lmer(y ~ 1 + z + (1 + z | id), data = dat, REML = FALSE)
manual_eb <- get_stage1_eb_components(
  fit_obj = fit,
  data = dat,
  cluster_var = "id",
  outcome_var = "y",
  within_var = "z"
)
lme4_re <- lme4::ranef(fit, condVar = TRUE)$id
lme4_post <- attr(lme4_re, "postVar")

stopifnot(
  isTRUE(all.equal(manual_eb$u0_eb, as.numeric(lme4_re[["(Intercept)"]]), tolerance = 1e-8)),
  isTRUE(all.equal(manual_eb$u1_eb, as.numeric(lme4_re[["z"]]), tolerance = 1e-8)),
  isTRUE(all.equal(manual_eb$postvar11, as.numeric(lme4_post[1, 1, ]), tolerance = 1e-8)),
  isTRUE(all.equal(manual_eb$postvar12, as.numeric(lme4_post[1, 2, ]), tolerance = 1e-8)),
  isTRUE(all.equal(manual_eb$postvar22, as.numeric(lme4_post[2, 2, ]), tolerance = 1e-8)),
  isTRUE(all(is.finite(manual_eb$lambda22))),
  isTRUE(all(is.finite(manual_eb$theta22))),
  isTRUE(all(is.finite(manual_eb$corrected_z))),
  isTRUE(all(is.finite(manual_eb$corrected_z_diag)))
)

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
mermod_R_error <- tryCatch({
  get_gls_corrected_scores(
    fit_obj = fit,
    data = dat,
    cluster_var = "id",
    outcome_var = "y",
    within_var = "z",
    R_list = R_list
  )
  FALSE
}, error = function(e) grepl("R_list.*merMod", conditionMessage(e)))

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
  isTRUE(all.equal(gls_out$ols_var22[[1]], gls_out$gls_var22[[1]], tolerance = 1e-12)),
  isTRUE(mermod_R_error)
)

singular_dat <- data.frame(
  id = rep(c("a", "b"), each = 4L),
  z = rep(c(-1.5, -0.5, 0.5, 1.5), 2L),
  y = c(0.1, 0.4, 0.7, 1.0, -0.2, 0.0, 0.2, 0.4)
)
singular_fit <- structure(list(), class = "singular_stage1_fit")
extract_stage1_components.singular_stage1_fit <- function(fit_obj, data, cluster_var, within_var = NULL,
                                                          R_list = NULL, group = NULL) {
  cluster_ids <- unique(as.character(data[[cluster_var]]))
  split_dat <- split(data, as.character(data[[cluster_var]]), drop = TRUE)[cluster_ids]
  list(
    beta_hat = c("(Intercept)" = 0.2, z = 0.3),
    G_hat = matrix(c(0.25, 0.10, 0.10, 0.04), nrow = 2L),
    R_list = stats::setNames(lapply(split_dat, function(df_i) 0.09 * diag(nrow(df_i))), cluster_ids),
    re_names_raw = c("(Intercept)", "z")
  )
}

singular_components <- get_stage1_eb_components(
  fit_obj = singular_fit,
  data = singular_dat,
  cluster_var = "id",
  outcome_var = "y",
  within_var = "z"
)

stopifnot(
  identical(singular_components$id, c("a", "b")),
  all(c("lambda11", "lambda12", "lambda21", "lambda22", "theta11", "theta12", "theta22") %in%
    names(singular_components)),
  all(is.finite(singular_components$u0_eb)),
  all(is.finite(singular_components$u1_eb)),
  all(is.finite(as.matrix(singular_components[, c(
    "postvar11", "postvar12", "postvar22",
    "lambda11", "lambda12", "lambda21", "lambda22",
    "theta11", "theta12", "theta22",
    "mle_z", "mle_z_var"
  ), drop = FALSE]))),
  all(is.na(singular_components$corrected_z)),
  all(is.na(singular_components$corrected_z_var))
)

if (requireNamespace("nlme", quietly = TRUE)) {
  fit_nlme <- tryCatch(
    nlme::lme(
      fixed = y ~ z,
      random = ~1 + z | id,
      correlation = nlme::corAR1(form = ~trial_index | id),
      data = dat,
      method = "ML",
      control = nlme::lmeControl(returnObject = TRUE, msMaxIter = 100L, opt = "optim")
    ),
    error = function(e) NULL
  )

  if (!is.null(fit_nlme)) {
    nlme_extractor_out <- get_gls_corrected_scores(
      fit_obj = fit_nlme,
      data = dat,
      cluster_var = "id",
      outcome_var = "y",
      within_var = "z",
      R_list = R_list
    )
    nlme_resid_1 <- df_1$y - drop(Z_1 %*% nlme::fixef(fit_nlme)[c("(Intercept)", "z")])
    expected_nlme_gls <- as.numeric(expected_gls_vcov %*% crossprod(Z_1, solve(R_i, nlme_resid_1)))

    nlme_fitted_R_out <- get_gls_corrected_scores(
      fit_obj = fit_nlme,
      data = dat,
      cluster_var = "id",
      outcome_var = "y",
      within_var = "z"
    )

    stopifnot(
      isTRUE(all.equal(nlme_extractor_out$mle_z[[1]], expected_nlme_gls[[2]], tolerance = 1e-10)),
      isTRUE(all.equal(nlme_extractor_out$mle_z_var[[1]], expected_gls_vcov[2, 2], tolerance = 1e-10)),
      isTRUE(all(is.finite(nlme_extractor_out$blup_intercept))),
      isTRUE(all(is.finite(nlme_extractor_out$blup_z))),
      isTRUE(all(is.finite(nlme_fitted_R_out$mle_z))),
      isTRUE(all(is.finite(nlme_fitted_R_out$corrected_z_var)))
    )
  }
}

cat("Closed-form GLS corrected-score tests ok\n")
