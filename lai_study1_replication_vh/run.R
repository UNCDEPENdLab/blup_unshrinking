#!/usr/bin/env Rscript

# Command-line entry point for the Lai Study 1 / VH-estimator replication.
#
# Usage:
# Rscript lai_study1_replication_vh/run.R 1000 outputs/lai_study1_vh 1 NA 1

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(MASS)
  library(lme4)
  library(OpenMx)
  library(sandwich)
  library(geigen)
  library(readr)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) normalizePath(sub("^--file=", "", script_arg[[1]])) else normalizePath("lai_study1_replication_vh/run.R")
repo_root <- normalizePath(file.path(dirname(script_path), ".."))
refresh_dir <- file.path(repo_root, "lai_study1_replication_vh")

source(file.path(repo_root, "R", "source_helpers.R"), local = TRUE)
source_project_helpers(repo_root)
source(file.path(refresh_dir, "designs.R"), local = TRUE)
source(file.path(refresh_dir, "study1.R"), local = TRUE)
source(file.path(refresh_dir, "estimators.R"), local = TRUE)
source(file.path(refresh_dir, "runner.R"), local = TRUE)

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  n_sim <- if (length(args) >= 1L) as.integer(args[[1]]) else 100L
  out_dir <- if (length(args) >= 2L) args[[2]] else file.path(refresh_dir, "outputs")
  n_cores <- if (length(args) >= 3L) as.integer(args[[3]]) else 1L
  max_conditions <- if (length(args) >= 4L && !args[[4]] %in% c("", "NA")) as.integer(args[[4]]) else NA_integer_
  resume_existing <- if (length(args) >= 5L) args[[5]] %in% c("1", "true", "TRUE") else TRUE
  methods <- if (length(args) >= 6L) strsplit(args[[6]], ",", fixed = TRUE)[[1]] else lai_study1_vh_methods()
  condition_ids <- if (length(args) >= 7L && !args[[7]] %in% c("", "NA")) {
    as.integer(strsplit(args[[7]], ",", fixed = TRUE)[[1]])
  } else {
    NULL
  }
  aggregate_only <- if (length(args) >= 8L) args[[8]] %in% c("1", "true", "TRUE") else FALSE
  write_aggregate <- if (length(args) >= 9L) {
    args[[9]] %in% c("1", "true", "TRUE")
  } else {
    is.null(condition_ids)
  }

  run_lai_study1_vh(
    n_sim = n_sim,
    out_dir = out_dir,
    n_cores = n_cores,
    max_conditions = max_conditions,
    methods = methods,
    resume_existing = resume_existing,
    condition_ids = condition_ids,
    aggregate_only = aggregate_only,
    write_aggregate = write_aggregate
  )
}
