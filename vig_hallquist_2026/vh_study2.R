#' Vig-Hallquist (2026) simulation study 2: BLUP as predictor
#' 
study2_methods <- function() {
  c(
    "oracle",
    "naive_dual_blup",
    "closed_form",
    "naive_slope_blup",
    "fuller_closed_form",
    "fuller_alpha_stepdown_closed_form",
    "single_subject_ols",
    "lai_2spa",
    "msem"
  )
}

simulate_study2 <- function(condition) {
  
}

run_study2_rep <- function(condition) {
  sim <- simulate_data_blup_as_predictor(condition)
  # run_matched_outcome_rep(condition, simulate_study2(condition))
}