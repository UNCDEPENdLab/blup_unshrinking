#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
})

r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
for (f in r_files) {
  source(f, local = TRUE)
}
cat("All R/ files loaded successfully.\n")

lai_replication_files <- paste0("lai_replication/",
    c(
      "designs.R",
      "runner.R",
      "study_common.R",
      "study1.R",
      "study2.R",
      "study3.R",
      "study4.R"
    )
)
for (f in lai_replication_files) {
  source(f, local = TRUE)
}
cat("All lai_replication/ files loaded successfully.\n")

blup_outcome_files <- paste0("blup_outcome/",
    c(
      "designs.R",
      "runner.R",
      "study_common.R"
    )
)
for (f in blup_outcome_files) {
  source(f, local = TRUE)
}
cat("All blup_outcome/ files loaded successfully.\n")