#!/usr/bin/env Rscript

# Smoke test that root helper files can be sourced through the central source
# order used by simulation entry points and that key shared functions are
# available.

suppressPackageStartupMessages({
  library(lme4)
})

source(file.path("R", "source_helpers.R"), local = TRUE)

helper_files <- project_helper_source_order()
source_project_helpers(".", helper_files = helper_files)

expected_functions <- c(
  "project_helper_source_order",
  "source_project_helpers",
  "safe_lmer",
  "safe_lme",
  "slice_condition_chunk",
  "make_derivative_backend",
  "make_tmb_stage1_data",
  "stacked_sandwich_for_corrected_scores",
  "simulate_dataset",
  "extract_lm_stats",
  "format_stacked_sandwich_rows",
  "select_eiv_result_columns",
  "extract_stage1_components",
  "get_stage1_eb_components",
  "get_corrected_scores",
  "get_gls_corrected_scores",
  "get_stage1_diagnostics",
  "fit_lai_2spa"
)

missing_functions <- expected_functions[!vapply(expected_functions, exists, logical(1), mode = "function")]
if (length(missing_functions) > 0L) {
  stop("Missing expected helper functions: ", paste(missing_functions, collapse = ", "))
}

cat("helper source-order smoke test ok\n")
