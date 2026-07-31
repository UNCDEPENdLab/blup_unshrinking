# Condition and replication runner for the Lai Study 1 VH refresh.

lai_study1_vh_pipeline_version <- function() {
  "lai_study1_vh_raw_latent_sd_lai_original_scale_eligibility_array_v5"
}

#' Select explicit condition IDs for an array task or focused rerun.
select_lai_study1_vh_conditions <- function(design, condition_ids = NULL) {
  if (is.null(condition_ids)) {
    return(design)
  }
  condition_ids <- sort(unique(as.integer(condition_ids)))
  if (length(condition_ids) == 0L || any(!is.finite(condition_ids))) {
    stop("`condition_ids` must contain one or more finite condition IDs.")
  }
  missing_ids <- setdiff(condition_ids, design$condition_id)
  if (length(missing_ids) > 0L) {
    stop("Unknown Lai Study 1 VH condition ID(s): ", paste(missing_ids, collapse = ", "))
  }
  dplyr::filter(design, condition_id %in% condition_ids)
}

#' Determine whether a condition file is complete for the current pipeline.
#'
#' A pre-existing file is reusable only when it has the current pipeline
#' version, every requested method and replication ID, and exactly one raw and
#' common-latent-SD row for every method/replication pair.  Historical-Lai
#' rows are optional because failed 2S-PA/MSEM fits may not yield a fitted
#' method-specific multiplier; when present they must be nonduplicated and
#' limited to methods with a historical counterpart.
lai_study1_vh_condition_output_is_complete <- function(result_file, condition,
                                                        n_sim, methods) {
  if (!file.exists(result_file)) return(FALSE)
  results <- tryCatch(readr::read_csv(result_file, show_col_types = FALSE), error = function(e) NULL)
  required <- c("pipeline_version", "condition_id", "rep", "method", "reporting_scale")
  if (is.null(results) || !all(required %in% names(results))) return(FALSE)

  n_sim <- as.integer(n_sim)
  methods <- unique(as.character(methods))
  condition_id <- as.integer(condition$condition_id[[1]])
  if (!identical(unique(as.character(results$pipeline_version)), lai_study1_vh_pipeline_version()) ||
      !identical(sort(unique(as.integer(results$condition_id))), condition_id) ||
      !setequal(unique(as.character(results$method)), methods) ||
      !setequal(unique(as.integer(results$rep)), seq_len(n_sim))) {
    return(FALSE)
  }

  allowed_scales <- c("raw", "latent_sd", "lai_original_standardized")
  if (any(!as.character(results$reporting_scale) %in% allowed_scales)) return(FALSE)
  base <- results[as.character(results$reporting_scale) %in% c("raw", "latent_sd"), , drop = FALSE]
  expected_base_rows <- n_sim * length(methods) * 2L
  if (nrow(base) != expected_base_rows) return(FALSE)
  expected_keys <- as.vector(outer(
    as.vector(outer(seq_len(n_sim), methods, paste, sep = "\r")),
    c("raw", "latent_sd"),
    paste,
    sep = "\r"
  ))
  base_keys <- paste(base$rep, base$method, base$reporting_scale, sep = "\r")
  if (!setequal(base_keys, expected_keys) || anyDuplicated(base_keys)) return(FALSE)

  historical <- results[as.character(results$reporting_scale) == "lai_original_standardized", , drop = FALSE]
  historical_methods <- c("naive_dual_blup", "lai_2spa", "msem")
  historical_keys <- paste(historical$rep, historical$method, sep = "\r")
  !any(!as.character(historical$method) %in% historical_methods) && !anyDuplicated(historical_keys)
}

#' Rebuild aggregate summaries after all selected condition files are complete.
rebuild_lai_study1_vh_summary <- function(n_sim, out_dir,
                                          methods = lai_study1_vh_methods(),
                                          max_conditions = NA_integer_,
                                          condition_ids = NULL) {
  design <- make_lai_study1_vh_design(max_conditions = max_conditions) |>
    select_lai_study1_vh_conditions(condition_ids)
  condition_dir <- file.path(out_dir, "conditions")
  incomplete <- integer()
  summaries <- vector("list", nrow(design))
  for (i in seq_len(nrow(design))) {
    condition <- design[i, , drop = FALSE]
    result_file <- file.path(condition_dir, sprintf("condition_%04d_replications.csv.gz", condition$condition_id))
    if (!lai_study1_vh_condition_output_is_complete(result_file, condition, n_sim, methods)) {
      incomplete <- c(incomplete, condition$condition_id)
      next
    }
    results <- readr::read_csv(result_file, show_col_types = FALSE)
    summaries[[i]] <- summarize_lai_study1_vh_results(results)
    readr::write_csv(
      summaries[[i]],
      file.path(condition_dir, sprintf("condition_%04d_summary.csv", condition$condition_id))
    )
  }
  if (length(incomplete) > 0L) {
    stop(
      "Cannot rebuild Lai Study 1 VH aggregate summary; incomplete condition IDs: ",
      paste(incomplete, collapse = ", ")
    )
  }
  summary_df <- dplyr::bind_rows(summaries)
  readr::write_csv(summary_df, file.path(out_dir, "lai_study1_vh_summary.csv"))
  invisible(summary_df)
}

summarize_lai_study1_vh_results <- function(results) {
  results %>%
    dplyr::mutate(
      converged = is.finite(estimate) & (is.na(status_code) | status_code == 0L),
      analysis_ready = converged & analysis_eligible,
      interval_ready = analysis_ready & interval_eligible,
      bias = dplyr::if_else(analysis_ready, estimate - truth, NA_real_),
      sq_error = bias^2,
      covered = dplyr::if_else(interval_ready, ci_low <= truth & ci_high >= truth, NA)
    ) %>%
    dplyr::group_by(
      condition_id, study, reporting_scale, method, num_clus, clus_size, icc,
      vr_u1_u0, cor_u0_u1, beta_zu1, var_u1, sigma2,
      dgm_population_slope_sd, dgm_posterior_slope_reliability,
      dgm_source, dgm_commit
    ) %>%
    dplyr::summarise(
      truth = dplyr::first(truth),
      n_rep = dplyr::n(),
      n_converged = sum(converged),
      n_analysis_eligible = sum(analysis_ready),
      n_vh_analysis_eligible = sum(vh_analysis_eligible, na.rm = TRUE),
      n_lai_original_eligible = sum(lai_original_eligible, na.rm = TRUE),
      n_eligible_both = sum(vh_analysis_eligible & lai_original_eligible, na.rm = TRUE),
      n_lai_original_only = sum(lai_original_eligible & !vh_analysis_eligible, na.rm = TRUE),
      n_vh_only = sum(vh_analysis_eligible & !lai_original_eligible, na.rm = TRUE),
      convergence = mean(converged),
      vh_analysis_eligibility_rate = mean(vh_analysis_eligible, na.rm = TRUE),
      lai_original_eligibility_rate = mean(lai_original_eligible, na.rm = TRUE),
      mean_estimate = mean(estimate[analysis_ready], na.rm = TRUE),
      bias = mean(bias, na.rm = TRUE),
      mc_se_mean = if (sum(analysis_ready) > 1L) {
        stats::sd(estimate[analysis_ready]) / sqrt(sum(analysis_ready))
      } else NA_real_,
      rmse = sqrt(mean(sq_error, na.rm = TRUE)),
      coverage = mean(covered, na.rm = TRUE),
      mean_realized_true_slope_sd = mean(stage1_realized_true_slope_sd, na.rm = TRUE),
      mean_eb_slope_sd = mean(stage1_eb_slope_sd, na.rm = TRUE),
      mean_corrected_slope_sd = mean(stage1_corrected_slope_sd, na.rm = TRUE),
      mean_fitted_slope_sd = mean(stage1_fitted_slope_sd, na.rm = TRUE),
      mean_fitted_residual_sd = mean(stage1_fitted_residual_sd, na.rm = TRUE),
      mean_posterior_slope_variance = mean(stage1_mean_posterior_slope_variance, na.rm = TRUE),
      mean_lambda22 = mean(stage1_mean_lambda22, na.rm = TRUE),
      mean_fitted_posterior_slope_reliability = mean(stage1_fitted_posterior_slope_reliability, na.rm = TRUE),
      mean_eb_to_population_slope_sd = mean(stage1_eb_to_population_slope_sd, na.rm = TRUE),
      mean_corrected_to_population_slope_sd = mean(stage1_corrected_to_population_slope_sd, na.rm = TRUE),
      mean_fitted_to_population_slope_sd = mean(stage1_fitted_to_population_slope_sd, na.rm = TRUE),
      mean_lai_original_naive_eb_slope_sd = mean(lai_original_naive_eb_slope_sd, na.rm = TRUE),
      mean_lai_original_naive_eb_slope_sample_sd = mean(lai_original_naive_eb_slope_sample_sd, na.rm = TRUE),
      mean_lai_original_method_multiplier = mean(lai_original_method_multiplier, na.rm = TRUE),
      lai_original_reporting_available_rate = mean(lai_original_reporting_available, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      dplyr::across(
        where(is.numeric),
        ~ dplyr::if_else(is.nan(.x), NA_real_, .x)
      )
    )
}

run_lai_study1_vh_condition <- function(condition,
                                         n_sim,
                                         methods = lai_study1_vh_methods(),
                                         seed = 9308L,
                                         n_cores = 1L) {
  condition <- tibble::as_tibble(condition)
  if (nrow(condition) != 1L) stop("`condition` must have exactly one row.")
  n_sim <- as.integer(n_sim)
  if (!is.finite(n_sim) || n_sim < 1L) stop("`n_sim` must be a positive integer.")

  one_replication <- function(rep_id) {
    set.seed(as.integer(seed + condition$condition_id[[1]] * 100000L + rep_id))
    sim <- simulate_lai_study1_vh(condition)
    fit_lai_study1_vh_estimators(condition, sim, methods = methods) %>%
      add_lai_study1_vh_reporting_scales(condition) %>%
      dplyr::mutate(
        condition_id = condition$condition_id,
        study = condition$study,
        rep = rep_id,
        pipeline_version = lai_study1_vh_pipeline_version(),
        num_clus = condition$num_clus,
        clus_size = condition$clus_size,
        icc = condition$icc,
        vr_u1_u0 = condition$vr_u1_u0,
        cor_u0_u1 = condition$cor_u0_u1,
        beta_zu1 = condition$beta_zu1,
        var_u1 = condition$var_u1,
        sigma2 = condition$sigma2,
        dgm_population_slope_sd = condition$dgm_population_slope_sd,
        dgm_posterior_slope_reliability = condition$dgm_posterior_slope_reliability,
        dgm_source = condition$dgm_source,
        dgm_commit = condition$dgm_commit
      )
  }

  rep_ids <- seq_len(n_sim)
  results <- if (n_cores > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(rep_ids, one_replication, mc.cores = n_cores, mc.set.seed = FALSE)
  } else {
    lapply(rep_ids, one_replication)
  }
  dplyr::bind_rows(results)
}

run_lai_study1_vh <- function(n_sim = 1000L,
                               out_dir,
                               n_cores = 1L,
                               max_conditions = NA_integer_,
                               methods = lai_study1_vh_methods(),
                               resume_existing = TRUE,
                               condition_ids = NULL,
                               aggregate_only = FALSE,
                               write_aggregate = is.null(condition_ids)) {
  methods <- unique(as.character(methods))
  if (!setequal(methods, lai_study1_vh_methods()) &&
      length(setdiff(methods, lai_study1_vh_methods())) > 0L) {
    stop("Unknown Lai Study 1 VH method(s): ", paste(setdiff(methods, lai_study1_vh_methods()), collapse = ", "))
  }
  if (isTRUE(aggregate_only)) {
    return(invisible(rebuild_lai_study1_vh_summary(
      n_sim = n_sim,
      out_dir = out_dir,
      methods = methods,
      max_conditions = max_conditions,
      condition_ids = condition_ids
    )))
  }

  design <- make_lai_study1_vh_design(max_conditions = max_conditions) |>
    select_lai_study1_vh_conditions(condition_ids)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  condition_dir <- file.path(out_dir, "conditions")
  dir.create(condition_dir, recursive = TRUE, showWarnings = FALSE)

  for (i in seq_len(nrow(design))) {
    condition <- design[i, , drop = FALSE]
    result_file <- file.path(condition_dir, sprintf("condition_%04d_replications.csv.gz", condition$condition_id))
    summary_file <- file.path(condition_dir, sprintf("condition_%04d_summary.csv", condition$condition_id))
    if (isTRUE(resume_existing) &&
        lai_study1_vh_condition_output_is_complete(result_file, condition, n_sim, methods)) {
      results <- readr::read_csv(result_file, show_col_types = FALSE)
    } else {
      results <- run_lai_study1_vh_condition(
        condition, n_sim = n_sim, methods = methods, n_cores = n_cores
      )
      readr::write_csv(results, result_file)
    }
    condition_summary <- summarize_lai_study1_vh_results(results)
    readr::write_csv(condition_summary, summary_file)
  }

  if (isTRUE(write_aggregate)) {
    summary_df <- rebuild_lai_study1_vh_summary(
      n_sim = n_sim,
      out_dir = out_dir,
      methods = methods,
      max_conditions = max_conditions,
      condition_ids = condition_ids
    )
    return(invisible(list(design = design, summary = summary_df)))
  }
  invisible(list(design = design, summary = NULL))
}
