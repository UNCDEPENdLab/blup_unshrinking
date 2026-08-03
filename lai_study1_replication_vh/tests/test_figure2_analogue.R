#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
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

# The production postprocessor must cache the historical trimmed summaries,
# while deriving the VH companion from the compact aggregate summary.
cache_root <- tempfile("lai_study1_vh_figure_cache_")
cache_results_dir <- file.path(cache_root, "results")
cache_condition_dir <- file.path(cache_results_dir, "conditions")
cache_analysis_dir <- file.path(cache_root, "analysis")
dir.create(cache_condition_dir, recursive = TRUE)

cache_base <- tibble::tibble(
  condition_id = 1L,
  method = rep(c("naive_dual_blup", "lai_2spa", "msem"), each = 2L),
  rep = rep(1:2, 3L),
  icc = 0.2,
  cor_u0_u1 = 0.5,
  beta_zu1 = 0,
  num_clus = 100L,
  clus_size = 10L,
  vr_u1_u0 = 1,
  var_u1 = 0.2,
  status_code = 0L,
  vh_analysis_eligible = TRUE,
  lai_original_eligible = TRUE
)
cache_truth <- 0
cache_results <- dplyr::bind_rows(
  cache_base |> dplyr::mutate(reporting_scale = "raw", estimate = 0.4, truth = 0.4, ci_low = 0.2, ci_high = 0.6),
  cache_base |> dplyr::mutate(reporting_scale = "latent_sd", estimate = 0.45, truth = cache_truth, ci_low = 0.2, ci_high = 0.7),
  cache_base |> dplyr::mutate(reporting_scale = "lai_original_standardized", estimate = 0.45, truth = cache_truth, ci_low = 0.2, ci_high = 0.7)
)
cache_result_file <- file.path(cache_condition_dir, "condition_0001_replications.csv.gz")
readr::write_csv(cache_results, cache_result_file)

cache_primary <- summarize_lai_study1_vh_primary_figure(cache_results)
cache_aggregate <- cache_primary |>
  dplyr::transmute(
    condition_id, method, reporting_scale, icc, cor_u0_u1, beta_zu1,
    num_clus, clus_size, vr_u1_u0, var_u1, truth,
    n_rep = n_sim,
    n_vh_analysis_eligible,
    bias = mean_bias,
    vh_analysis_eligibility_rate
  )
readr::write_csv(cache_aggregate, file.path(cache_results_dir, "lai_study1_vh_summary.csv"))

cache_run <- run_lai_study1_vh_postestimation_figures(
  cache_results_dir,
  cache_analysis_dir
)
stopifnot(
  file.exists(file.path(cache_analysis_dir, "figure2_analogue_cache_manifest.csv")),
  file.exists(file.path(cache_analysis_dir, "historical_condition_cache", "condition_0001_historical_summary.csv")),
  file.exists(file.path(cache_analysis_dir, "historical_condition_cache", "condition_0001_historical_inference.csv")),
  file.exists(file.path(cache_analysis_dir, "figure3_coverage_analogue.png")),
  file.exists(file.path(cache_analysis_dir, "figure4_type1_analogue.png")),
  file.exists(file.path(cache_analysis_dir, "vh_primary_companion_cell_summary.csv")),
  nrow(cache_run$figure2_analogue) == 3L,
  nrow(cache_run$figure3_4_analogue) == 3L,
  nrow(cache_run$vh_primary) == 3L,
  identical(cache_run$vh_primary$n_sim, rep(2, 3L))
)

# A legacy combined summary can be promoted to the cache without rereading the
# compressed replication result, which is important after a pre-cache run.
bootstrap_analysis_dir <- file.path(cache_root, "bootstrap_analysis")
dir.create(bootstrap_analysis_dir, recursive = TRUE)
file.copy(
  lai_study1_vh_historical_summary_path(cache_analysis_dir),
  lai_study1_vh_historical_summary_path(bootstrap_analysis_dir)
)
bootstrapped <- lai_study1_vh_historical_summary_from_cache(cache_results_dir, bootstrap_analysis_dir)
stopifnot(
  nrow(bootstrapped) == 3L,
  file.exists(lai_study1_vh_historical_manifest_path(bootstrap_analysis_dir)),
  file.exists(file.path(bootstrap_analysis_dir, "historical_condition_cache", "condition_0001_historical_summary.csv"))
)

cat("Lai Study 1 Figure 2 analogue tests ok\n")
