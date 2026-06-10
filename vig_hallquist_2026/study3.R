#' Vig-Hallquist (2026) simulation study 3: BLUP as predictor and outcome


study3_methods <- function() {
  c(
    "oracle",
    "naive_blup_on_blup",
    "closed_form_on_blup",
    "blup_on_closed_form",
    "closed_form_on_closed_form",
    "fuller_closed_form",
    "lai_2spa",
    "sem"
  )
}

simulate_study3 <- function(condition) {

}

run_study3_rep <- function(condition) {
  # run_matched_outcome_rep(condition, simulate_study3(condition))
}