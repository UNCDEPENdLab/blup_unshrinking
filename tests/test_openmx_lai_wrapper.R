#!/usr/bin/env Rscript

# Small synthetic check for the Lai structural-slope OpenMx wrapper.

if (!requireNamespace("OpenMx", quietly = TRUE)) {
  cat("OpenMx not installed; skipping Lai OpenMx wrapper test\n")
  quit(status = 0L)
}

source(file.path("R", "core_utils.R"), local = TRUE)
source(file.path("R", "lai_openmx_helpers.R"), local = TRUE)

# Failed OpenMx algebras sometimes evaluate to length zero. The diagnostics
# contract must still yield exactly one typed row so a 2S-PA method cannot
# disappear from a replication.
empty_diagnostics <- mx_diagnostics_tibble(list(
  mx_raw_focal_estimate = numeric(),
  mx_raw_focal_se = NULL,
  mx_latent_covariance_boundary = logical()
))
stopifnot(
  nrow(empty_diagnostics) == 1L,
  is.na(empty_diagnostics$mx_raw_focal_estimate),
  is.na(empty_diagnostics$mx_raw_focal_se),
  is.logical(empty_diagnostics$mx_latent_covariance_boundary),
  is.na(empty_diagnostics$mx_latent_covariance_boundary)
)

# `mxPath(connect = "unique.bivariate")` fills the RAM A matrix in
# column-major order. Use four distinct entries so this test fails if the
# Stage-1 row-major names (11, 12, 21, 22) are passed through unchanged and the
# off-diagonal loadings are silently transposed.
loading_probe <- c(
  lambda11 = 11,
  lambda12 = 12,
  lambda21 = 21,
  lambda22 = 22
)
openmx_loading_cols <- openmx_bivariate_loading_columns()
expected_loading_matrix <- matrix(
  c(11, 12, 21, 22),
  nrow = 2L,
  byrow = TRUE,
  dimnames = list(c("u0_eb", "u1_eb"), c("u0", "u1"))
)

loading_value_model <- OpenMx::mxModel(
  "loading_value_order_probe",
  type = "RAM",
  manifestVars = c("u0_eb", "u1_eb"),
  latentVars = c("u0", "u1"),
  OpenMx::mxPath(
    from = c("u0", "u1"),
    to = c("u0_eb", "u1_eb"),
    connect = "unique.bivariate",
    free = FALSE,
    values = unname(loading_probe[openmx_loading_cols])
  )
)
actual_loading_values <- loading_value_model$A$values[
  c("u0_eb", "u1_eb"),
  c("u0", "u1"),
  drop = FALSE
]

loading_label_model <- OpenMx::mxModel(
  "loading_label_order_probe",
  type = "RAM",
  manifestVars = c("u0_eb", "u1_eb"),
  latentVars = c("u0", "u1"),
  OpenMx::mxPath(
    from = c("u0", "u1"),
    to = c("u0_eb", "u1_eb"),
    connect = "unique.bivariate",
    free = FALSE,
    labels = paste0("data.", openmx_loading_cols)
  )
)
actual_loading_labels <- loading_label_model$A$labels[
  c("u0_eb", "u1_eb"),
  c("u0", "u1"),
  drop = FALSE
]
expected_loading_labels <- matrix(
  paste0("data.", c("lambda11", "lambda12", "lambda21", "lambda22")),
  nrow = 2L,
  byrow = TRUE,
  dimnames = dimnames(expected_loading_matrix)
)

stopifnot(
  identical(
    openmx_loading_cols,
    c("lambda11", "lambda21", "lambda12", "lambda22")
  ),
  isTRUE(all.equal(
    actual_loading_values,
    expected_loading_matrix,
    tolerance = 0
  )),
  identical(actual_loading_labels, expected_loading_labels)
)

# Inspect the A matrices actually assembled by the shared 2S-PA wrapper. Stub
# only model execution and result extraction; `fit_lai_2spa()` still constructs
# its complete OpenMx RAM model through the production code path.
wrapper_probe_df <- data.frame(
  x = c(-1.5, -0.5, 0.5, 1.5),
  u0_eb = c(-0.4, 0.2, -0.1, 0.3),
  u1_eb = c(0.1, -0.3, 0.4, -0.2),
  lambda11 = 11,
  lambda12 = 12,
  lambda21 = 21,
  lambda22 = 22,
  theta11 = 1,
  theta12 = 0,
  theta22 = 1
)
lai_helper_environment <- environment(fit_lai_2spa)
original_run_mx_safe <- get("run_mx_safe", envir = lai_helper_environment)
original_extract_mx_stats <- get(
  "extract_mx_stats",
  envir = lai_helper_environment
)
captured_wrapper_model <- NULL
capture_run_mx_safe <- function(mx_mod, ...) {
  captured_wrapper_model <<- mx_mod
  mx_mod
}
dummy_extract_mx_stats <- function(mx_fit, ...) {
  # Force the normally lazy first argument so the production call to
  # `run_mx_safe(mx_mod)` executes our model-capture stub.
  force(mx_fit)
  tibble::tibble(
    estimate = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
    status_code = NA_integer_
  )
}
assign("run_mx_safe", capture_run_mx_safe, envir = lai_helper_environment)
assign(
  "extract_mx_stats",
  dummy_extract_mx_stats,
  envir = lai_helper_environment
)

invisible(fit_lai_2spa(wrapper_probe_df, use_average = FALSE))
wrapper_definition_labels <- captured_wrapper_model$A$labels[
  c("u0_eb", "u1_eb"),
  c("u0", "u1"),
  drop = FALSE
]
invisible(fit_lai_2spa(wrapper_probe_df, use_average = TRUE))
wrapper_average_values <- captured_wrapper_model$A$values[
  c("u0_eb", "u1_eb"),
  c("u0", "u1"),
  drop = FALSE
]

assign("run_mx_safe", original_run_mx_safe, envir = lai_helper_environment)
assign(
  "extract_mx_stats",
  original_extract_mx_stats,
  envir = lai_helper_environment
)
rm(
  original_run_mx_safe,
  original_extract_mx_stats,
  captured_wrapper_model,
  capture_run_mx_safe,
  dummy_extract_mx_stats,
  lai_helper_environment
)

stopifnot(
  identical(wrapper_definition_labels, expected_loading_labels),
  isTRUE(all.equal(
    wrapper_average_values,
    expected_loading_matrix,
    tolerance = 0
  ))
)

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
  "mx_info_definite", "mx_condition_number", "mx_issue_class", "mx_issue_detail",
  "mx_raw_focal_estimate", "mx_raw_focal_se",
  "mx_latent_intercept_variance", "mx_latent_slope_variance",
  "mx_latent_intercept_slope_covariance",
  "mx_latent_intercept_slope_correlation",
  "mx_latent_covariance_min_eigenvalue", "mx_latent_covariance_boundary"
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
  is.finite(out$ci_high),
  is.finite(out$mx_raw_focal_estimate),
  is.finite(out$mx_raw_focal_se),
  is.finite(out$mx_latent_intercept_variance),
  is.finite(out$mx_latent_slope_variance),
  is.finite(out$mx_latent_covariance_min_eigenvalue),
  !out$mx_latent_covariance_boundary
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
  )),
  isTRUE(all.equal(
    dual_scale_2$mx_raw_focal_estimate,
    dual_scale_1$mx_raw_focal_estimate,
    tolerance = 1e-8
  )),
  isTRUE(all.equal(
    dual_scale_2$mx_raw_focal_se,
    dual_scale_1$mx_raw_focal_se,
    tolerance = 1e-8
  )),
  is.finite(dual_scale_1$mx_predictor_latent_slope_variance),
  is.finite(dual_scale_1$mx_outcome_latent_slope_variance)
)

cat("Lai OpenMx wrapper test ok\n")
