#!/usr/bin/env Rscript

#' ---
#' title: "Score Correction Method Comparison"
#' description: "Replicate the Lai & Liu random-slope predictor simulation designs
#'               using study-specific modules, then compare the performance of different 
#'               score correction methods."
#' ---

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
  library(MASS)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(foreach)
  library(doParallel)
  library(OpenMx)
  library(glmnet)
  library(sandwich)
  library(geigen)
})

current_script_path <- function() {
  frame_files <- vapply(sys.frames(), function(frame) {
    if (!is.null(frame$ofile)) frame$ofile else NA_character_
  }, character(1))
  frame_files <- stats::na.omit(frame_files)
  if (length(frame_files) > 0L) {
    return(normalizePath(frame_files[[length(frame_files)]], mustWork = FALSE))
  }

  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0L) {
    return(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE))
  }

  NA_character_
}

find_repo_root <- function() {
  script_path <- current_script_path()
  script_dir <- if (!is.na(script_path)) dirname(script_path) else getwd()
  candidates <- unique(normalizePath(c(
    getwd(),
    file.path(getwd(), ".."),
    script_dir,
    file.path(script_dir, ".."),
    file.path(script_dir, "..", "..")
  ), mustWork = FALSE))

  for (candidate in candidates) {
    if (file.exists(file.path(candidate, "R", "core_utils.R")) &&
      file.exists(file.path(candidate, "R", "blup_helpers.R"))) {
      return(candidate)
    }
  }

  stop("Could not locate repository root containing shared R/ helpers and score_correction_methods/ scripts.")
}

repo_root <- find_repo_root()
score_correction_dir <- file.path(repo_root, "score_correction_methods")

source(file.path(repo_root, "R", "core_utils.R"), local = TRUE)
source(file.path(repo_root, "R", "blup_helpers.R"), local = TRUE)
source(file.path(repo_root, "R", "sim_helpers.R"), local = TRUE)

source(file.path(score_correction_dir, "design.R"), local = TRUE)
source(file.path(score_correction_dir, "runner.R"), local = TRUE)


args <- commandArgs(trailingOnly = TRUE)

if (sys.nframe() == 0L) {
  n_sim <- if (length(args) >= 1L) as.integer(args[[1]]) else 100L
  out_dir <- if (length(args) >= 2L) args[[2]] else file.path(score_correction_dir, "outputs", "score_correction_comparison")
  n_cores <- if (length(args) >= 3L) as.integer(args[[3]]) else 1L
  max_conditions <- if (length(args) >= 4L) parse_optional_integer_arg(args[[4]]) else NA_integer_
  chunk_index <- if (length(args) >= 5L) parse_optional_integer_arg(args[[5]]) else NA_integer_
  chunk_size <- if (length(args) >= 6L) parse_optional_integer_arg(args[[6]]) else NA_integer_
  resume_existing <- if (length(args) >= 7L) parse_logical_arg(args[[7]], default = TRUE) else TRUE

  run_score_correction_simulation(
    n_sim = n_sim,
    out_dir = out_dir,
    n_cores = n_cores,
    max_conditions = max_conditions,
    chunk_index = chunk_index,
    chunk_size = chunk_size,
    resume_existing = resume_existing
  )
}