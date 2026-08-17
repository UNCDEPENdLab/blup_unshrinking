#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(dplyr))

source(file.path("R", "core_utils.R"), local = TRUE)
source(file.path("R", "stage2_estimators.R"), local = TRUE)
source(file.path("lai_replication", "study_common.R"), local = TRUE)
source(file.path("lai_replication", "runner.R"), local = TRUE)

rows <- tibble::tibble(
  method = c("oracle_dual", "naive_dual_eb", "fuller_alpha_stepdown", "lai_2spa", "ridge_dual_eb"),
  estimate = rep(0.4, 5L),
  se = c(0.1, 0.1, 0.1, 0.1, NA_real_),
  ci_low = c(0.2, 0.2, 0.2, 0.2, NA_real_),
  ci_high = c(0.6, 0.6, 0.6, 0.6, NA_real_),
  status_code = rep(0L, 5L),
  analysis_eligible = c(FALSE, FALSE, NA, NA, NA),
  analysis_exclusion_reason = c("stage2_near_collinear", "stage2_rank_deficient", NA, NA, NA),
  fuller_auto_guard_pass = c(NA, NA, FALSE, NA, NA),
  fuller_auto_guard_reason = c(NA, NA, "condition_cap", NA, NA),
  mx_issue_class = c(NA, NA, NA, "ok", NA),
  mx_info_definite = c(NA, NA, NA, FALSE, NA),
  mx_condition_number = rep(NA_real_, 5L)
)

classified <- add_matched_outcome_analysis_eligibility(rows)
stopifnot(
  identical(classified$analysis_eligible, c(FALSE, FALSE, FALSE, FALSE, TRUE)),
  identical(classified$interval_eligible, c(FALSE, FALSE, FALSE, FALSE, FALSE)),
  identical(
    classified$analysis_exclusion_reason,
    c(
      "stage2_near_collinear", "stage2_rank_deficient", "fuller_guard_condition_cap",
      "openmx_information_not_definite", NA_character_
    )
  ),
  identical(classified$interval_exclusion_reason[[5L]], "invalid_standard_error")
)

summary_input <- tibble::tibble(
  condition_id = 1L,
  study = "study4",
  method = "naive_dual_eb",
  num_clus = 50L,
  clus_size = 10L,
  icc = 0.1,
  vr_u1_u0 = 1,
  cor_u0_u1 = 0.2,
  beta_zu1 = 0.4,
  design_source = "unit",
  condition_note = "unit",
  estimate = c(1, 100, 2),
  truth = 1,
  se = c(0.1, 0.1, NA_real_),
  ci_low = c(0.8, 99.8, NA_real_),
  ci_high = c(1.2, 100.2, NA_real_),
  status_code = 0L,
  analysis_eligible = c(TRUE, FALSE, TRUE),
  interval_eligible = c(TRUE, FALSE, FALSE)
)

summary_out <- summarize_results_df(summary_input)
stopifnot(
  summary_out$n_success == 3L,
  summary_out$n_analysis_eligible == 2L,
  summary_out$n_interval_eligible == 1L,
  isTRUE(all.equal(summary_out$mean_estimate, 1.5)),
  isTRUE(all.equal(summary_out$rmse, sqrt(0.5))),
  isTRUE(all.equal(summary_out$coverage, 1))
)

stopifnot(all(c(
  "add_matched_outcome_analysis_eligibility",
  "assess_dual_ols_design"
) %in% lai_parallel_exports()))

cat("Lai analysis-eligibility tests ok\n")
