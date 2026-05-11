#!/usr/bin/env Rscript

# Unit checks for formatting stacked-sandwich outputs into simulation estimator
# rows.

source(file.path("R", "stage2_estimators.R"), local = TRUE)

sandwich_out <- list(
  alpha_hat = c(1.25, -0.40),
  vcov_hc0 = matrix(c(0.10, 0.01, 0.01, 0.04), nrow = 2L),
  vcov_hc1 = matrix(c(0.12, 0.01, 0.01, 0.09), nrow = 2L),
  vcov_hc2 = matrix(c(0.14, 0.01, 0.01, 0.16), nrow = 2L),
  vcov_hc3 = matrix(c(0.16, 0.01, 0.01, 0.25), nrow = 2L)
)

df <- 9L
out <- format_stacked_sandwich_rows(
  sandwich_out = sandwich_out,
  df = df,
  alpha_names = c("(Intercept)", "x"),
  term = "x",
  method_prefix = "unit_stacked"
)

expected_methods <- paste0("unit_stacked_hc", 0:3)
expected_se <- c(0.20, 0.30, 0.40, 0.50)
crit <- qt(0.975, df)

stopifnot(
  identical(out$method, expected_methods),
  isTRUE(all.equal(out$estimate, rep(-0.40, 4L), tolerance = 1e-12)),
  isTRUE(all.equal(out$se, expected_se, tolerance = 1e-12)),
  isTRUE(all.equal(out$ci_low, -0.40 - crit * expected_se, tolerance = 1e-12)),
  isTRUE(all.equal(out$ci_high, -0.40 + crit * expected_se, tolerance = 1e-12))
)

missing_term <- format_stacked_sandwich_rows(
  sandwich_out = sandwich_out,
  df = df,
  alpha_names = c("(Intercept)", "x"),
  term = "not_a_term",
  method_prefix = "unit_stacked"
)

stopifnot(
  identical(missing_term$method, expected_methods),
  all(is.na(missing_term$estimate)),
  all(is.na(missing_term$se)),
  all(is.na(missing_term$ci_low)),
  all(is.na(missing_term$ci_high))
)

cat("stage-2 estimator formatter tests ok\n")
