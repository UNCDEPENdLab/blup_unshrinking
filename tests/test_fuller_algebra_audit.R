#!/usr/bin/env Rscript

# Regression checks supporting the Fuller algebra audit and corrected
# production implementation.

source(file.path(
  "vig_hallquist_2026", "analysis", "fuller_algebra_audit_simulation.R"
))
source(file.path("R", "stage2_estimators.R"))

set.seed(20260817)
sim <- simulate_audit_replication(
  m = 200L,
  geometry = "regular",
  scenario = "outcome_homoskedastic"
)

dat <- data.frame(
  y = sim$y,
  x0 = sim$x[, 2L],
  x1 = sim$x[, 3L],
  meas11 = sim$omega_x[, 2L, 2L],
  meas12 = sim$omega_x[, 2L, 3L],
  meas22 = sim$omega_x[, 3L, 3L],
  measyy = sim$omega_y,
  measxy0 = sim$omega_xy[, 2L],
  measxy1 = sim$omega_xy[, 3L]
)

production <- fit_fuller_dual(
  dat,
  outcome = "y",
  predictor_u0 = "x0",
  predictor_u1 = "x1",
  meas11 = "meas11",
  meas12 = "meas12",
  meas22 = "meas22",
  outcome_meas_var = "measyy",
  predictor_outcome_meas_cov_u0 = "measxy0",
  predictor_outcome_meas_cov_u1 = "measxy1",
  skip_internal_scaling = TRUE
)
audit_current <- fit_fuller_audit(
  y = sim$y,
  x = sim$x,
  omega_x = sim$omega_x,
  omega_y = sim$omega_y,
  omega_xy = sim$omega_xy,
  preliminary = "current_recycled",
  se_bread = "current_final"
)
audit_cross_fixed <- fit_fuller_audit(
  y = sim$y,
  x = sim$x,
  omega_x = sim$omega_x,
  omega_y = sim$omega_y,
  omega_xy = sim$omega_xy,
  preliminary = "current_cross",
  se_bread = "current_final"
)

# The independently coded corrected branch must reproduce every focal result.
stopifnot(
  production$status_code[[1L]] == 0L,
  audit_cross_fixed$status[[1L]] == "ok",
  isTRUE(all.equal(
    production$estimate[[1L]], audit_cross_fixed$estimate[[1L]], tolerance = 1e-10
  )),
  isTRUE(all.equal(
    production$se[[1L]], audit_cross_fixed$se[[1L]], tolerance = 1e-9
  )),
  isTRUE(all.equal(
    production$fuller_lambda1[[1L]], audit_cross_fixed$lambda1[[1L]],
    tolerance = 1e-10
  )),
  isTRUE(all.equal(
    production$fuller_lambda2[[1L]], audit_cross_fixed$lambda2[[1L]],
    tolerance = 1e-10
  )),
  isTRUE(all.equal(
    production$fuller_sigma2[[1L]], audit_cross_fixed$sigma2[[1L]],
    tolerance = 1e-9
  )),
  identical(production$fuller_preliminary_moment[[1L]], "modified"),
  identical(production$fuller_variance_bread[[1L]], "modified"),
  identical(
    production$fuller_predictor_outcome_covariance_source[[1L]], "supplied"
  )
)

# With outcome error only and zero predictor-outcome error covariance, Fuller's
# preliminary cross moment is X'y. The corrected preliminary coefficient is
# therefore ordinary OLS; the production-recycled coefficient is not.
ols_gamma <- drop(solve(crossprod(sim$x), crossprod(sim$x, sim$y)))
stopifnot(
  isTRUE(all.equal(
    audit_cross_fixed$gamma0_target[[1L]], tail(ols_gamma, 1L),
    tolerance = 1e-10
  )),
  abs(audit_current$gamma0_target[[1L]] - tail(ols_gamma, 1L)) > 1e-3
)

# A nonzero predictor--outcome covariance must propagate through the complete
# covariance blocks, preliminary cross-moment, composite residual variance,
# weights, final point estimate, and sandwich meat.
set.seed(202608171)
sim_cross <- simulate_audit_replication(
  m = 180L,
  geometry = "regular",
  scenario = "both_homoskedastic"
)
sim_cross$omega_xy[, 2L] <- 0.08
sim_cross$omega_xy[, 3L] <- -0.04
cross_dat <- data.frame(
  y = sim_cross$y,
  x0 = sim_cross$x[, 2L],
  x1 = sim_cross$x[, 3L],
  meas11 = sim_cross$omega_x[, 2L, 2L],
  meas12 = sim_cross$omega_x[, 2L, 3L],
  meas22 = sim_cross$omega_x[, 3L, 3L],
  measyy = sim_cross$omega_y,
  measxy0 = sim_cross$omega_xy[, 2L],
  measxy1 = sim_cross$omega_xy[, 3L]
)
production_cross <- fit_fuller_dual(
  cross_dat,
  outcome = "y",
  predictor_u0 = "x0",
  predictor_u1 = "x1",
  meas11 = "meas11",
  meas12 = "meas12",
  meas22 = "meas22",
  outcome_meas_var = "measyy",
  predictor_outcome_meas_cov_u0 = "measxy0",
  predictor_outcome_meas_cov_u1 = "measxy1"
)
audit_cross <- fit_fuller_audit(
  y = sim_cross$y,
  x = sim_cross$x,
  omega_x = sim_cross$omega_x,
  omega_y = sim_cross$omega_y,
  omega_xy = sim_cross$omega_xy,
  preliminary = "current_cross",
  se_bread = "current_final"
)
stopifnot(
  production_cross$status_code[[1L]] == 0L,
  audit_cross$status[[1L]] == "ok",
  isTRUE(all.equal(
    production_cross$estimate[[1L]], audit_cross$estimate[[1L]], tolerance = 1e-10
  )),
  isTRUE(all.equal(
    production_cross$se[[1L]], audit_cross$se[[1L]], tolerance = 1e-9
  )),
  isTRUE(all.equal(
    production_cross$fuller_sigma2[[1L]], audit_cross$sigma2[[1L]],
    tolerance = 1e-9
  )),
  isTRUE(all.equal(
    production_cross$fuller_predictor_outcome_covariance_max_abs[[1L]],
    0.08,
    tolerance = 1e-12
  ))
)

# The same covariance algebra must hold for the Study-1-style single-predictor
# pathway; this also guards the slope coefficient's index in the composite
# measurement variance and tilde-covariance terms.
omega_x_single <- sim_cross$omega_x[, c(1L, 3L), c(1L, 3L), drop = FALSE]
dim(omega_x_single) <- c(nrow(cross_dat), 2L, 2L)
omega_xy_single <- sim_cross$omega_xy[, c(1L, 3L), drop = FALSE]
production_single_cross <- fit_fuller_dual(
  cross_dat,
  outcome = "y",
  predictor_u1 = "x1",
  meas22 = "meas22",
  outcome_meas_var = "measyy",
  predictor_outcome_meas_cov_u1 = "measxy1"
)
audit_single_cross <- fit_fuller_audit(
  y = sim_cross$y,
  x = cbind(1, sim_cross$x[, 3L]),
  omega_x = omega_x_single,
  omega_y = sim_cross$omega_y,
  omega_xy = omega_xy_single,
  preliminary = "current_cross",
  se_bread = "current_final"
)
stopifnot(
  production_single_cross$status_code[[1L]] == 0L,
  audit_single_cross$status[[1L]] == "ok",
  isTRUE(all.equal(
    production_single_cross$estimate[[1L]], audit_single_cross$estimate[[1L]],
    tolerance = 1e-10
  )),
  isTRUE(all.equal(
    production_single_cross$se[[1L]], audit_single_cross$se[[1L]],
    tolerance = 1e-9
  ))
)

# The production comparison helper is the corrected nested sequence used by
# VH Studies 1--4: stabilized, Fuller preliminary, then all Fuller equations.
production_variants <- fit_fuller_dual_variants(
  cross_dat,
  outcome = "y",
  predictor_u0 = "x0",
  predictor_u1 = "x1",
  meas11 = "meas11",
  meas12 = "meas12",
  meas22 = "meas22",
  outcome_meas_var = "measyy",
  predictor_outcome_meas_cov_u0 = "measxy0",
  predictor_outcome_meas_cov_u1 = "measxy1"
)
audit_book_preliminary <- fit_fuller_audit(
  sim_cross$y, sim_cross$x, sim_cross$omega_x, sim_cross$omega_y,
  sim_cross$omega_xy, preliminary = "documented", se_bread = "current_final"
)
audit_book <- fit_fuller_audit(
  sim_cross$y, sim_cross$x, sim_cross$omega_x, sim_cross$omega_y,
  sim_cross$omega_xy, preliminary = "documented", se_bread = "documented_full"
)
stopifnot(
  identical(
    production_variants$fuller_variant,
    c("stabilized", "fuller_preliminary", "fuller_equations")
  ),
  isTRUE(all.equal(
    production_variants$estimate[[2L]], audit_book_preliminary$estimate[[1L]],
    tolerance = 1e-10
  )),
  isTRUE(all.equal(
    production_variants$se[[2L]], audit_book_preliminary$se[[1L]],
    tolerance = 1e-9
  )),
  isTRUE(all.equal(
    production_variants$estimate[[3L]], audit_book$estimate[[1L]],
    tolerance = 1e-10
  )),
  isTRUE(all.equal(
    production_variants$se[[3L]], audit_book$se[[1L]], tolerance = 1e-9
  ))
)

# In the absence of measurement error all audited variants coincide.
set.seed(20260818)
no_error <- simulate_audit_replication(
  m = 100L,
  geometry = "weak_shape",
  scenario = "no_measurement_error"
)
no_error_fits <- lapply(seq_len(nrow(audit_variants)), function(i) {
  fit_fuller_audit(
    y = no_error$y,
    x = no_error$x,
    omega_x = no_error$omega_x,
    omega_y = no_error$omega_y,
    omega_xy = no_error$omega_xy,
    preliminary = audit_variants$preliminary[[i]],
    se_bread = audit_variants$se_bread[[i]]
  )
})
no_error_estimates <- vapply(no_error_fits, `[[`, numeric(1L), "estimate")
no_error_ses <- vapply(no_error_fits, `[[`, numeric(1L), "se")
stopifnot(
  max(abs(no_error_estimates - no_error_estimates[[1L]])) < 1e-10,
  max(abs(no_error_ses - no_error_ses[[1L]])) < 1e-10
)

# The book variance bread changes inference, not the final point estimator.
set.seed(20260819)
predictor_error <- simulate_audit_replication(
  m = 100L,
  geometry = "regular",
  scenario = "both_heteroskedastic"
)
documented_preliminary <- fit_fuller_audit(
  y = predictor_error$y,
  x = predictor_error$x,
  omega_x = predictor_error$omega_x,
  omega_y = predictor_error$omega_y,
  omega_xy = predictor_error$omega_xy,
  preliminary = "documented",
  se_bread = "current_final"
)
documented_full <- fit_fuller_audit(
  y = predictor_error$y,
  x = predictor_error$x,
  omega_x = predictor_error$omega_x,
  omega_y = predictor_error$omega_y,
  omega_xy = predictor_error$omega_xy,
  preliminary = "documented",
  se_bread = "documented_full"
)
stopifnot(
  isTRUE(all.equal(
    documented_preliminary$estimate[[1L]], documented_full$estimate[[1L]],
    tolerance = 1e-12
  )),
  abs(documented_preliminary$se[[1L]] - documented_full$se[[1L]]) > 1e-8
)

cat("Fuller algebra audit tests ok\n")
