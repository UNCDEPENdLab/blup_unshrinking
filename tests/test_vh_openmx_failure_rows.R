#!/usr/bin/env Rscript

# Regression test for the exact Study-1v2 replication that exposed the
# length-zero OpenMx diagnostic bug. The 2S-PA fit is expected to fail in this
# difficult low-reliability sample, but it must remain one explicit failed
# method row and must not print the OpenMx `fitUnits`/`vcov()` error.

source(file.path(
  "vig_hallquist_2026",
  "random_effects_structural_simulation.R"
))

design <- select_design("v2")
condition <- design[design$condition_id == 2493L, , drop = FALSE]
stopifnot(nrow(condition) == 1L)

set.seed(vh_replication_seed(condition, 2L))
result <- run_study1_rep(condition)

stopifnot(
  nrow(result) == length(study1_methods()),
  all(table(result$method) == 1L),
  sum(result$method == "lai_2spa") == 1L
)
lai_row <- result[result$method == "lai_2spa", , drop = FALSE]
stopifnot(
  nrow(lai_row) == 1L,
  lai_row$status_code != 0L,
  !is.na(lai_row$mx_issue_class),
  lai_row$mx_issue_class != "ok",
  is.na(lai_row$estimate),
  is.na(lai_row$se)
)

cat("VH OpenMx failed-fit row regression test ok\n")
