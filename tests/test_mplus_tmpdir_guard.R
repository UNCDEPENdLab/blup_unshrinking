#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
})

source(file.path(
  "vig_hallquist_2026",
  "random_effects_structural_simulation.R"
))

condition <- select_design("5") %>%
  filter(
    calibration_arm == "shape_preserving_partial",
    mean_clus_size == 10L,
    marginal_rho == 0.5,
    standardized_beta_target == 0.2
  ) %>%
  slice(1L)

set.seed(20260814)
sim <- simulate_study5(condition)
old_tmpdir <- Sys.getenv("TMPDIR", unset = NA_character_)
scratch_before <- list.dirs(tempdir(), recursive = FALSE, full.names = TRUE)
fit <- fit_mplus_blup_predictor(
  level1_data = sim$lv1,
  level2_data = sim$lv2_true,
  outcome_variable = "y",
  within_component = "x",
  between_component = "z",
  cluster_id = "cid",
  reporting_scale = condition$tau1[[1]]
)
scratch_after <- list.dirs(tempdir(), recursive = FALSE, full.names = TRUE)

stopifnot(
  fit$status_code == 0L,
  fit$mx_issue_class == "ok",
  is.finite(fit$estimate),
  is.finite(fit$se),
  fit$se > 0,
  identical(Sys.getenv("TMPDIR", unset = NA_character_), old_tmpdir),
  setequal(scratch_before, scratch_after)
)

cat("Mplus writable-TMPDIR guard test ok\n")
