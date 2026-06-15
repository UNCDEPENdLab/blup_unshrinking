#!/usr/bin/env Rscript

# Small synthetic check for the Lai structural-slope OpenMx wrapper.

if (!requireNamespace("OpenMx", quietly = TRUE)) {
  cat("OpenMx not installed; skipping Lai OpenMx wrapper test\n")
  quit(status = 0L)
}

source(file.path("R", "core_utils.R"), local = TRUE)
source(file.path("R", "lai_openmx_helpers.R"), local = TRUE)

set.seed(3030)
n <- 30L
x <- seq(-1.5, 1.5, length.out = n)
u0_lat <- rnorm(n, sd = 0.55)
u1_lat <- 0.45 * x + rnorm(n, sd = 0.35)

stage2_df <- data.frame(
  x = x,
  u0_eb = 0.88 * u0_lat + 0.05 * u1_lat + rnorm(n, sd = 0.08),
  u1_eb = 0.04 * u0_lat + 0.86 * u1_lat + rnorm(n, sd = 0.08),
  lambda11 = rep(0.88, n),
  lambda12 = rep(0.05, n),
  lambda21 = rep(0.04, n),
  lambda22 = rep(0.86, n),
  theta11 = rep(0.03, n),
  theta12 = rep(0.005, n),
  theta22 = rep(0.03, n)
)

out <- fit_lai_2spa(stage2_df)
required_cols <- c(
  "estimate", "se", "ci_low", "ci_high", "status_code", "mx_status_msg",
  "mx_info_definite", "mx_condition_number", "mx_issue_class", "mx_issue_detail"
)
missing_cols <- setdiff(required_cols, names(out))
if (length(missing_cols) > 0L) {
  stop("Missing expected OpenMx output columns: ", paste(missing_cols, collapse = ", "))
}

stopifnot(
  nrow(out) == 1L,
  identical(out$mx_issue_class, "ok"),
  identical(as.integer(out$status_code), 0L),
  is.finite(out$estimate),
  is.finite(out$se),
  is.finite(out$ci_low),
  is.finite(out$ci_high)
)

set.seed(3031)
n_dual <- 180L
y0_lat <- rnorm(n_dual, sd = 0.8)
y1_lat <- 0.3 * y0_lat + rnorm(n_dual, sd = 0.6)
q0_lat <- rnorm(n_dual, sd = 0.7)
q1_lat <- 0.25 * y0_lat + 0.55 * y1_lat + rnorm(n_dual, sd = 0.5)
dual_df <- data.frame(
  u0_eb = y0_lat + rnorm(n_dual, sd = 0.1),
  u1_eb = y1_lat + rnorm(n_dual, sd = 0.1),
  q_u0_eb = q0_lat + rnorm(n_dual, sd = 0.1),
  q_u1_eb = q1_lat + rnorm(n_dual, sd = 0.1),
  lambda11 = 1,
  lambda12 = 0,
  lambda21 = 0,
  lambda22 = 1,
  theta11 = 0.01,
  theta12 = 0,
  theta22 = 0.01,
  q_lambda11 = 1,
  q_lambda12 = 0,
  q_lambda21 = 0,
  q_lambda22 = 1,
  q_theta11 = 0.01,
  q_theta12 = 0,
  q_theta22 = 0.01
)
dual_scale_1 <- fit_lai_2spa_dual_process(
  dual_df,
  theta0_start = 0.25,
  theta1_start = 0.55,
  reporting_scale = 1
)
dual_scale_2 <- fit_lai_2spa_dual_process(
  dual_df,
  theta0_start = 0.25,
  theta1_start = 0.55,
  reporting_scale = 2
)

stopifnot(
  identical(dual_scale_1$mx_issue_class, "ok"),
  identical(dual_scale_2$mx_issue_class, "ok"),
  isTRUE(all.equal(
    dual_scale_2$estimate,
    2 * dual_scale_1$estimate,
    tolerance = 1e-8
  )),
  isTRUE(all.equal(
    dual_scale_2$se,
    2 * dual_scale_1$se,
    tolerance = 1e-8
  ))
)

cat("Lai OpenMx wrapper test ok\n")
