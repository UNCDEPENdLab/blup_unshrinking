#' Condition grids for BLUP-as-stage-2-outcome simulations.

blup_outcome_balance_modes <- function() {
  c("balanced", "unbalanced", "informative_unbalanced")
}

balance_mode_to_sim_arg <- function(balance_mode) {
  switch(
    as.character(balance_mode[[1]]),
    balanced = "balanced",
    unbalanced = FALSE,
    informative_unbalanced = "highly_unbalanced",
    stop("Unsupported balance mode: ", balance_mode)
  )
}

make_blup_outcome_design <- function(grid_mode = "base", max_conditions = NA_integer_) {
  grid_mode <- as.character(grid_mode[[1]])

  design <- switch(
    grid_mode,
    smoke = tidyr::crossing(
      n_id = 40L,
      mean_n_trial = 8L,
      gamma_x_on_slope = 0.3,
      rho = 0.3,
      balance_mode = "unbalanced",
      min_n_trial = 4L,
      highly_unbalanced_min_n_trial = 3L,
      highly_unbalanced_power = 3,
      tau1 = 0.7,
      sigma = 1.0
    ),
    base = tidyr::crossing(
      n_id = c(50L, 200L),
      mean_n_trial = c(4L, 8L, 20L, 50L),
      gamma_x_on_slope = c(0.0, 0.1, 0.5),
      rho = c(-0.5, 0.0, 0.5),
      balance_mode = blup_outcome_balance_modes(),
      min_n_trial = 2L,
      highly_unbalanced_min_n_trial = 2L,
      highly_unbalanced_power = 3,
      tau1 = c(0.3, 1.2),
      sigma = c(0.5, 1.5)
    ),
    heteroscedastic = tidyr::crossing(
      n_id = c(50L, 100L, 200L),
      mean_n_trial = c(8L, 20L, 50L),
      gamma_x_on_slope = c(0.0, 0.5),
      rho = c(-0.5, 0.5),
      balance_mode = c("balanced", "informative_unbalanced"),
      min_n_trial = 2L,
      highly_unbalanced_min_n_trial = 2L,
      highly_unbalanced_power = c(2, 4),
      tau1 = c(0.7, 1.2),
      sigma = c(1.0, 1.5)
    ),
    heteroscedastic_sparse = tidyr::crossing(
      n_id = c(50L, 100L, 200L),
      mean_n_trial = c(3L, 4L, 6L, 8L, 12L, 20L),
      gamma_x_on_slope = c(0.0, 0.1, 0.5),
      rho = c(-0.5, 0.0, 0.5),
      balance_mode = c("balanced", "unbalanced", "informative_unbalanced"),
      min_n_trial = 2L,
      highly_unbalanced_min_n_trial = 2L,
      highly_unbalanced_power = 3,
      tau1 = c(0.3, 0.7, 1.2),
      sigma = c(0.5, 1.0, 1.5)
    ),
    expanded = tidyr::crossing(
      n_id = c(30L, 50L, 100L, 200L, 500L),
      mean_n_trial = c(3L, 4L, 6L, 8L, 12L, 20L, 50L, 100L),
      gamma_x_on_slope = c(0.0, 0.1, 0.3, 0.5),
      rho = c(-0.5, 0.0, 0.5, 0.8),
      balance_mode = blup_outcome_balance_modes(),
      min_n_trial = 2L,
      highly_unbalanced_min_n_trial = 2L,
      highly_unbalanced_power = c(2, 3, 4),
      tau1 = c(0.3, 0.7, 1.2),
      sigma = c(0.5, 1.0, 1.5)
    ),
    stop("`grid_mode` must be one of: smoke, base, heteroscedastic, heteroscedastic_sparse, expanded.")
  )

  design <- design %>%
    dplyr::mutate(
      condition_id = dplyr::row_number(),
      balanced = balance_mode,
      design_source = grid_mode
    ) %>%
    dplyr::relocate(condition_id)

  if (!is.na(max_conditions)) {
    design <- design %>% dplyr::slice(seq_len(min(max_conditions, nrow(design))))
  }

  design
}
