source("R/stage2_estimators.R")

# simulate ill-conditioned EIV data
set.seed(42)
n <- 300
u0 <- rnorm(n)
u1 <- 0.5 * u0 + rnorm(n, sd = 0.5)
y <- 1.0 + 0.3 * u0 + 0.8 * u1 + rnorm(n)

# Extremely large measurement error (larger than true variance!)
meas_u0 <- rnorm(n, sd = 2.0)
meas_u1 <- rnorm(n, sd = 2.0)

x0_obs <- u0 + meas_u0
x1_obs <- u1 + meas_u1

dat <- data.frame(
  y_obs = y,
  x0_obs = x0_obs,
  x1_obs = x1_obs,
  meas11 = rep(4.0, n),
  meas12 = rep(0, n),
  meas22 = rep(4.0, n)
)

res <- fit_fuller_dual_stepdown(dat, "y_obs", "x0_obs", "x1_obs", "meas11", "meas12", "meas22")
print(res[, c("status_code", "fuller_measurement_weight_used", "fuller_auto_guard_reason", "fuller_auto_search_evaluations")])
