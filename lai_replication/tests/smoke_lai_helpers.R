#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
})

find_repo_root <- function() {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  script_path <- if (length(script_arg) > 0L) {
    normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = FALSE)
  } else {
    normalizePath(file.path("lai_replication", "tests", "smoke_lai_helpers.R"), mustWork = FALSE)
  }

  candidates <- unique(normalizePath(c(
    getwd(),
    file.path(getwd(), ".."),
    file.path(dirname(script_path), "..", "..")
  ), mustWork = FALSE))

  for (candidate in candidates) {
    if (file.exists(file.path(candidate, "lai_replication", "mlm_random_slope_lai_apples_to_apples_sim.R"))) {
      return(candidate)
    }
  }

  getwd()
}

repo_root <- find_repo_root()
driver_path <- file.path(repo_root, "lai_replication", "mlm_random_slope_lai_apples_to_apples_sim.R")

if (!file.exists(driver_path)) {
  stop("Could not find Lai apples-to-apples driver from repo root: ", repo_root)
}

source(driver_path, local = TRUE)

out_dir <- file.path(tempdir(), paste0("lai_smoke_", Sys.getpid()))
res <- run_simulation(
  n_sim = 1L,
  study_arg = "1",
  out_dir = out_dir,
  n_cores = 1L,
  max_conditions = 1L,
  resume_existing = FALSE
)

required_cols <- c("method", "estimate", "se", "ci_low", "ci_high", "status_code", "truth")
missing_cols <- setdiff(required_cols, names(res$results))
if (length(missing_cols) > 0L) {
  stop("Missing expected result columns: ", paste(missing_cols, collapse = ", "))
}

required_methods <- c("oracle_dual", "naive_dual_eb", "corrected_dual", "lai_2spa")
missing_methods <- setdiff(required_methods, unique(res$results$method))
if (length(missing_methods) > 0L) {
  stop("Missing expected methods: ", paste(missing_methods, collapse = ", "))
}

if (nrow(res$summary) == 0L) {
  stop("Smoke run did not produce a summary.")
}

cat("lai helper smoke test ok\n")
