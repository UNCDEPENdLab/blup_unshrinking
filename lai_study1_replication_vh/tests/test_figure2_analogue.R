#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

find_repo_root <- function() {
  candidates <- unique(normalizePath(c(getwd(), file.path(getwd(), "..")), mustWork = FALSE))
  for (candidate in candidates) {
    if (file.exists(file.path(candidate, "lai_study1_replication_vh", "analysis", "figure2_analogue.R"))) {
      return(candidate)
    }
  }
  stop("Could not locate the repository root.")
}

repo_root <- find_repo_root()
source(file.path(repo_root, "lai_study1_replication_vh", "analysis", "figure2_analogue.R"), local = TRUE)

base_rows <- tibble::tibble(
  condition_id = 1L,
  method = rep(c("naive_dual_blup", "lai_2spa"), each = 2L),
  rep = rep(1:2, 2L),
  icc = 0.2,
  cor_u0_u1 = 0.5,
  beta_zu1 = 0.4,
  num_clus = 100L,
  clus_size = 10L,
  vr_u1_u0 = 1,
  var_u1 = 0.2,
  status_code = 0L,
  vh_analysis_eligible = c(TRUE, FALSE, TRUE, TRUE),
  lai_original_eligible = c(TRUE, TRUE, TRUE, FALSE)
)

raw_rows <- base_rows |>
  dplyr::mutate(reporting_scale = "raw", estimate = c(0.4, 0.4, 0.4, 0.4), truth = 0.4)
latent_rows <- base_rows |>
  dplyr::mutate(
    reporting_scale = "latent_sd", estimate = c(0.2, 0.4, 0.3, 0.5), truth = 0.4 * sqrt(0.2)
  )
historical_rows <- base_rows |>
  dplyr::filter(!(method == "lai_2spa" & rep == 2L)) |>
  dplyr::mutate(
    reporting_scale = "lai_original_standardized",
    estimate = c(0.1, 0.3, 0.5),
    truth = 0.4 * sqrt(0.2)
  )
results <- dplyr::bind_rows(raw_rows, latent_rows, historical_rows)

historical <- summarize_lai_study1_vh_figure2_analogue(results) |>
  dplyr::arrange(method)
primary <- summarize_lai_study1_vh_primary_figure(results) |>
  dplyr::arrange(method)

truth <- 0.4 * sqrt(0.2)
stopifnot(
  nrow(historical) == 2L,
  identical(historical$n_sim, c(2L, 2L)),
  identical(historical$n_lai_original_eligible, c(1L, 2L)),
  identical(historical$n_historical_scale_available, c(1L, 2L)),
  identical(historical$n_retained, c(1L, 2L)),
  isTRUE(all.equal(historical$robust_bias, c(0.5 - truth, 0.2 - truth), tolerance = 1e-12)),
  all(historical$reporting_scale == "lai_original_standardized"),
  all(historical$eligibility_rule == "lai_original_status_code_zero"),
  nrow(primary) == 2L,
  identical(primary$n_vh_analysis_eligible, c(2L, 1L)),
  identical(primary$n_retained, c(2L, 1L)),
  isTRUE(all.equal(primary$mean_bias, c(0.4 - truth, 0.2 - truth), tolerance = 1e-12)),
  all(primary$reporting_scale == "latent_sd")
)

cat("Lai Study 1 Figure 2 analogue tests ok\n")
