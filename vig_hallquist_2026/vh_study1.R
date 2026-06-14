#' Vig-Hallquist (2026) simulation study 1: BLUP as outcome

study1_methods <- function() {
  c(
    "oracle",
    "naive_blup",
    "closed_form",
    "fuller_closed_form",
    "fuller_alpha_stepdown_closed_form",
    "single_subject_ols",
    "lai_2spa",
    "direct_mlm"
  )
}

run_study1_rep <- function(condition) {
  sim <- simulate_data_blup_as_outcome(condition)
  # run_matched_outcome_rep(condition, simulate_study1(condition))
}
