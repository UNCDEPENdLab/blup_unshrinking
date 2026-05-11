#!/usr/bin/env Rscript

# Smoke test that root helper files can be sourced in the order used by the
# simulation entry point and that key shared functions are available.

suppressPackageStartupMessages({
  library(lme4)
})

helper_files <- c(
  "R/core_utils.R",
  "R/derivative_backends.R",
  "R/tmb_stage1_helpers.R",
  "R/stacked_sandwich_helpers.R",
  "R/stats_helpers.R",
  "R/stage2_estimators.R",
  "R/blup_helpers.R",
  "R/lai_openmx_helpers.R"
)

for (helper_file in helper_files) {
  source(helper_file, local = TRUE)
}

expected_functions <- c(
  "safe_lmer",
  "make_derivative_backend",
  "make_tmb_stage1_data",
  "stacked_sandwich_for_corrected_scores",
  "extract_lm_stats",
  "format_stacked_sandwich_rows",
  "get_corrected_scores",
  "fit_lai_2spa"
)

missing_functions <- expected_functions[!vapply(expected_functions, exists, logical(1), mode = "function")]
if (length(missing_functions) > 0L) {
  stop("Missing expected helper functions: ", paste(missing_functions, collapse = ", "))
}

cat("helper source-order smoke test ok\n")
