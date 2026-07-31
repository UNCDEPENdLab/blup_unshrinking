#!/usr/bin/env Rscript

# Usage:
# Rscript lai_study1_replication_vh/make_figure2_analogue.R \
#   outputs/lai_study1_vh outputs/lai_study1_vh/figure2_analogue

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[[1]]))
} else {
  normalizePath("lai_study1_replication_vh/make_figure2_analogue.R")
}
refresh_dir <- dirname(script_path)
source(file.path(refresh_dir, "analysis", "figure2_analogue.R"), local = TRUE)

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1L) {
    stop("Usage: Rscript make_figure2_analogue.R RESULTS_DIR [ANALYSIS_DIR]")
  }
  results_dir <- args[[1]]
  analysis_dir <- if (length(args) >= 2L) {
    args[[2]]
  } else {
    file.path(results_dir, "figure2_analogue")
  }
  run_lai_study1_vh_postestimation_figures(results_dir, analysis_dir)
}
