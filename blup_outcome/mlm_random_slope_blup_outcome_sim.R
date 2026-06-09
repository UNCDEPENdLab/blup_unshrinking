#!/usr/bin/env Rscript

#' BLUP-as-outcome random-slope simulation entry point.

# Load every package used by the sourced modules up front. The simulation can
# run through `foreach` workers, so package availability needs to be consistent
# between the parent process and worker sessions.
suppressPackageStartupMessages({
  library(lme4)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(foreach)
  library(doParallel)
  library(OpenMx)
  library(sandwich)
})

# Keep BLAS/OpenMP and OpenMx from introducing nested parallelism. The runner
# controls parallelism explicitly through `n_cores`; allowing lower-level math
# libraries to spawn threads can oversubscribe SLURM jobs and make timings noisy.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
OpenMx::mxOption(NULL, "Number of Threads", 1L)

# Locate the repository root robustly for both direct CLI execution and the
# historical compatibility wrapper. Rscript exposes the executing file through
# a `--file=` argument, but sourced execution may not, so the current working
# directory is also considered.
locate_repo_root <- function() {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  script_path <- if (length(script_arg) > 0L) {
    normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = FALSE)
  } else {
    NA_character_
  }

  # Candidate order is intentional: prefer the caller's working directory for
  # normal repo-root invocations, then the script directory, then its parent.
  # The final `mlm_blups` fallback preserves a legacy local folder name used by
  # older scripts.
  candidates <- unique(normalizePath(c(
    getwd(),
    if (!is.na(script_path)) dirname(script_path) else character(),
    if (!is.na(script_path)) file.path(dirname(script_path), "..") else character(),
    "mlm_blups"
  ), mustWork = FALSE))

  # `R/sim_helpers.R` is a stable sentinel for this repository and avoids
  # accidentally treating a parent project or output directory as the root.
  roots <- candidates[file.exists(file.path(candidates, "R", "sim_helpers.R"))]
  if (length(roots) == 0L) {
    stop("Could not locate repository root containing shared R helpers.")
  }
  roots[[1]]
}

repo_root <- locate_repo_root()

source(file.path(repo_root, "R", "source_helpers.R"), local = TRUE)
source_project_helpers(repo_root)
source(file.path(repo_root, "blup_outcome", "designs.R"), local = TRUE)
source(file.path(repo_root, "blup_outcome", "study_common.R"), local = TRUE)
source(file.path(repo_root, "blup_outcome", "runner.R"), local = TRUE)

args <- commandArgs(trailingOnly = TRUE)

# In array/batch jobs, default to the scheduler-provided CPU allocation when
# the caller does not pass an explicit core count. Invalid or missing scheduler
# values fall back to a single-core run.
default_n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
if (is.na(default_n_cores) || default_n_cores < 1L) {
  default_n_cores <- 1L
}

# CLI argument contract:
#   1  n_sim              replications per selected condition
#   2  out_dir            output directory for condition and aggregate files
#   3  n_cores            foreach workers; defaults to SLURM_CPUS_PER_TASK
#   4  derivative_method  backend for stacked-sandwich derivatives
#   5  grid_mode          design grid selector
#   6  analysis_mode      "screen" for cheap checks, "full" for all estimators
#   7  chunk_index        one-based chunk id; defaults to SLURM_ARRAY_TASK_ID
#   8  chunk_size         conditions per chunk; can use BLUP_OUTCOME_CHUNK_SIZE
#   9  resume_existing    logical flag for skipping completed condition files
#   10 execution_mode     "run" condition files or "aggregate" completed files
#   11 max_conditions     optional pre-chunk cap for smoke/debug runs
n_sim <- if (length(args) >= 1L) as.integer(args[[1]]) else 100L
out_dir <- if (length(args) >= 2L) args[[2]] else file.path(repo_root, "outputs", "blup_outcome")
n_cores <- if (length(args) >= 3L) as.integer(args[[3]]) else default_n_cores
derivative_method <- if (length(args) >= 4L) args[[4]] else "handcoded"
grid_mode <- if (length(args) >= 5L) args[[5]] else "base"
analysis_mode <- if (length(args) >= 6L) args[[6]] else "full"

# Chunk defaults are scheduler-aware so SLURM arrays can call the same command
# for every task. Explicit CLI values always win over environment variables.
chunk_index <- if (length(args) >= 7L) parse_optional_integer_arg(args[[7]]) else parse_optional_integer_arg(Sys.getenv("SLURM_ARRAY_TASK_ID", NA_character_))
chunk_size <- if (length(args) >= 8L) parse_optional_integer_arg(args[[8]]) else parse_optional_integer_arg(Sys.getenv("BLUP_OUTCOME_CHUNK_SIZE", NA_character_))
resume_existing <- if (length(args) >= 9L) parse_logical_arg(args[[9]], default = TRUE) else TRUE

# `run` mode writes only per-condition outputs for chunk-safe distributed jobs.
# `aggregate` mode reconstructs top-level CSVs and plots after chunks finish.
execution_mode <- if (length(args) >= 10L) args[[10]] else Sys.getenv("BLUP_OUTCOME_EXECUTION_MODE", "run")
max_conditions <- if (length(args) >= 11L) parse_optional_integer_arg(args[[11]]) else NA_integer_

# Hand off to the documented runner. Keeping all simulation work inside the
# runner makes this script a thin, testable CLI adapter.
run_blup_outcome_simulation(
  n_sim = n_sim,
  out_dir = out_dir,
  n_cores = n_cores,
  derivative_method = derivative_method,
  grid_mode = grid_mode,
  analysis_mode = analysis_mode,
  chunk_index = chunk_index,
  chunk_size = chunk_size,
  resume_existing = resume_existing,
  execution_mode = execution_mode,
  max_conditions = max_conditions
)
