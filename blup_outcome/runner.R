#' Runner and aggregation for BLUP-as-outcome simulations.
#'
#' This module is the orchestration layer for the random-slope
#' BLUP/corrected-score-as-outcome simulation. Study-specific data generation
#' and estimator fitting live in `blup_outcome/study_common.R` and design-grid
#' construction lives in `blup_outcome/designs.R`; this file is responsible for
#' turning a selected design grid into per-condition replication files,
#' resumable progress logs, aggregate summaries, and diagnostic plots.

#' Build standard condition-level output paths.
#'
#' @details
#' BLUP-outcome runs write one replication file and one summary file per design
#' condition under `file.path(out_dir, "conditions")`. This wrapper fixes the
#' filename prefix and compression policy expected by this simulation while
#' delegating the path construction to the shared runner helper.
#'
#' @param out_dir Character scalar output directory.
#' @param condition_id Integer-like condition identifier.
#'
#' @return A named list of file paths for the condition replication, summary,
#'   issue-summary, and stage1-summary artifacts.
blup_outcome_condition_paths <- function(out_dir, condition_id) {
  condition_file_paths(out_dir, condition_id, prefix = "condition", compressed_replications = FALSE)
}

#' Return display labels for BLUP-outcome estimator method codes.
#'
#' @details
#' Simulation result rows use compact method codes that are stable for CSV
#' output and programmatic filtering. Reports and plots use these labels for
#' readability. The names of the returned vector are method codes and the values
#' are display labels.
#'
#' @return A named character vector mapping method codes to human-readable
#'   labels.
blup_outcome_method_labels <- function() {
  c(
    oracle = "Oracle latent slope",
    naive_blup = "Naive BLUP outcome",
    naive_blup_hc3 = "Naive BLUP outcome (HC3)",
    diag_corrected = "Diagonal-only corrected outcome",
    diag_corrected_hc3 = "Diagonal-only corrected outcome (HC3)",
    matrix_corrected = "Full-matrix corrected outcome",
    matrix_corrected_hc3 = "Full-matrix corrected outcome (HC3)",
    closed_form = "Closed-form score outcome",
    closed_form_hc3 = "Closed-form score outcome (HC3)",
    single_subject_ols = "Single-subject OLS slope",
    single_subject_ols_hc3 = "Single-subject OLS slope (HC3)",
    lai_2spa = "Lai 2S-PA",
    lai_2spaa = "Lai 2S-PAA",
    closed_form_stacked_hc0 = "Closed-form score + stacked sandwich",
    closed_form_stacked_hc1 = "Closed-form score + stacked sandwich (HC1)",
    closed_form_stacked_hc2 = "Closed-form score + stacked sandwich (HC2)",
    closed_form_stacked_hc3 = "Closed-form score + stacked sandwich (HC3)",
    direct_mlm = "Direct mixed model"
  )
}

#' Summarize BLUP-outcome replication-level results by condition and method.
#'
#' @details
#' This summary function computes convergence, bias, RMSE, coverage, empirical
#' standard deviation, Monte Carlo standard error of the mean estimate, OpenMx
#' status-10 failure rates, first-stage singularity/problem rates, and selected
#' score-diagnostic averages. It is intentionally tolerant of missing
#' diagnostic columns so that `screen` and `full` analysis modes can share the
#' same aggregation path.
#'
#' A result row is treated as converged when it has a finite estimate and is not
#' an OpenMx status-10 failure. Coverage is computed from the row's confidence
#' interval against the simulation truth.
#'
#' @param results Replication-level tibble/data frame produced by
#'   `run_blup_outcome_condition_replications()` or loaded from condition files.
#'
#' @return A condition-method summary tibble with method labels attached. Returns
#'   an empty tibble when `results` has zero rows.
summarize_blup_outcome_results <- function(results) {
  if (nrow(results) == 0L) {
    return(tibble::tibble())
  }

  # Only summarize diagnostics that are present in the current result set. This
  # lets lean screen-mode outputs and full outputs use the same function.
  diagnostic_cols <- intersect(
    c(
      "mean_realized_trials", "min_realized_trials", "prop_ids_leq_2_trials", "prop_ids_leq_3_trials",
      "mean_postvar_u1", "mean_theta_u1", "mean_lambda_u1",
      "blup_variance_ratio", "diag_corrected_variance_ratio", "matrix_corrected_variance_ratio",
      "closed_form_variance_ratio", "blup_true_cor", "matrix_corrected_true_cor",
      "closed_form_true_cor", "diag_corrected_failure_rate", "matrix_corrected_failure_rate",
      "closed_form_failure_rate", "cluster_size_x_cor", "stage1_re_corr", "stage1_eb_corr",
      "stage1_design_kappa"
    ),
    names(results)
  )

  summary_df <- results %>%
    dplyr::mutate(
      # Status 10 is the OpenMx broad failure code used by the Lai wrappers;
      # these rows should not be counted as usable estimates.
      status10_failure = !is.na(.data$status_code) & .data$status_code == 10L,
      converged = !.data$status10_failure & !is.na(.data$estimate),
      bias = .data$estimate - .data$truth,
      sq_error = (.data$estimate - .data$truth)^2,
      covered = .data$ci_low <= .data$truth & .data$ci_high >= .data$truth
    ) %>%
    dplyr::group_by(
      .data$condition_id, .data$method, .data$n_id, .data$mean_n_trial, .data$gamma_x_on_slope,
      .data$rho, .data$balance_mode, .data$tau1, .data$sigma, .data$min_n_trial,
      .data$highly_unbalanced_min_n_trial, .data$highly_unbalanced_power, .data$design_source
    ) %>%
    dplyr::summarise(
      truth = dplyr::first(.data$truth),
      n_rep = dplyr::n(),
      convergence = safe_mean(.data$converged),
      mean_estimate = safe_mean(.data$estimate),
      mean_se = safe_mean(.data$se),
      emp_sd = if (sum(!is.na(.data$estimate)) > 1L) stats::sd(.data$estimate, na.rm = TRUE) else NA_real_,
      mc_se_mean = if (sum(!is.na(.data$estimate)) > 1L) stats::sd(.data$estimate, na.rm = TRUE) / sqrt(sum(!is.na(.data$estimate))) else NA_real_,
      bias = safe_mean(.data$bias),
      rmse = if (all(is.na(.data$sq_error))) NA_real_ else sqrt(mean(.data$sq_error, na.rm = TRUE)),
      coverage = safe_mean(.data$covered),
      n_success = sum(.data$converged, na.rm = TRUE),
      n_status10_fail = sum(.data$status10_failure, na.rm = TRUE),
      prop_status10_fail = safe_mean(.data$status10_failure),
      stage1_singular_rate = safe_mean(.data$stage1_lmer_singular),
      stage1_problem_rate = safe_mean(.data$stage1_singular_problem),
      dplyr::across(dplyr::all_of(diagnostic_cols), safe_mean),
      .groups = "drop"
    )

  labels <- blup_outcome_method_labels()
  summary_df %>%
    dplyr::mutate(
      # Preserve unknown/new method codes instead of dropping them; this keeps
      # exploratory estimator rows visible even before labels are added here.
      method_label = labels[as.character(.data$method)],
      method_label = dplyr::if_else(is.na(.data$method_label), as.character(.data$method), .data$method_label),
      method_label = factor(.data$method_label, levels = labels)
    )
}

#' Summarize non-OK OpenMx issue classes for BLUP-outcome results.
#'
#' @details
#' The Lai 2S-PA/PAA estimators can fail for several qualitatively different
#' OpenMx reasons. This helper collapses replication-level issue details into a
#' condition/method/issue-class table, preserving a compact sample of distinct
#' diagnostic details for inspection.
#'
#' @param results Replication-level result tibble/data frame.
#'
#' @return A tibble with issue counts and sample details. Returns an empty tibble
#'   when no OpenMx issue columns are present or no rows are available.
summarize_blup_outcome_issues <- function(results) {
  if (nrow(results) == 0L || !("mx_issue_class" %in% names(results))) {
    return(tibble::tibble())
  }

  results %>%
    dplyr::filter(!is.na(.data$mx_issue_class), .data$mx_issue_class != "ok") %>%
    dplyr::group_by(
      .data$condition_id, .data$method, .data$mx_issue_class, .data$n_id, .data$mean_n_trial,
      .data$gamma_x_on_slope, .data$rho, .data$balance_mode, .data$tau1, .data$sigma
    ) %>%
    dplyr::summarise(
      n_rep = dplyr::n(),
      n_distinct_issue_details = dplyr::n_distinct(.data$mx_issue_detail[!is.na(.data$mx_issue_detail)]),
      sample_issue_detail = compact_message(.data$mx_issue_detail[!is.na(.data$mx_issue_detail)]),
      .groups = "drop"
    )
}

#' Write aggregate BLUP-outcome CSV outputs and diagnostic plots.
#'
#' @details
#' This function writes the aggregate replication file, method/condition summary,
#' and OpenMx issue summary. When summary rows are available, it also writes
#' coverage and bias plots for the main estimator subset used in reports.
#'
#' The plotting subset intentionally omits many diagnostic variants so the
#' default figures remain readable. The full method set remains available in the
#' CSV outputs.
#'
#' @param results Replication-level results from completed condition files.
#' @param out_dir Character scalar output directory.
#' @param prefix Character scalar prefix used for aggregate filenames.
#'
#' @return Invisibly returns a list with `summary` and `issue_summary` tibbles.
write_blup_outcome_aggregate_outputs <- function(results, out_dir, prefix = "blup_outcome") {
  summary_df <- summarize_blup_outcome_results(results)
  issue_df <- summarize_blup_outcome_issues(results)

  # Aggregate CSVs are the canonical artifacts; plots are derived convenience
  # outputs and are skipped when there is nothing to plot.
  write_csv_atomic(results, file.path(out_dir, sprintf("%s_replication_results.csv", prefix)))
  write_csv_atomic(summary_df, file.path(out_dir, sprintf("%s_summary.csv", prefix)))
  write_csv_atomic(issue_df, file.path(out_dir, sprintf("%s_issue_summary.csv", prefix)))

  if (nrow(summary_df) == 0L) {
    return(invisible(list(summary = summary_df, issue_summary = issue_df)))
  }

  plot_df <- summary_df %>%
    dplyr::filter(.data$method %in% c("oracle", "naive_blup", "matrix_corrected", "closed_form", "closed_form_stacked_hc3", "lai_2spa", "direct_mlm"))

  if (nrow(plot_df) > 0L) {
    # Lines require more than one x-value within a panel/group. Points are still
    # drawn for one-off cells, but dropping singletons from geom_line avoids
    # ggplot warnings and misleading single-point segments.
    line_df <- plot_df %>%
      dplyr::group_by(.data$method_label, .data$balance_mode, .data$tau1, .data$sigma, .data$n_id) %>%
      dplyr::filter(dplyr::n() > 1L) %>%
      dplyr::ungroup()

    p_cov <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$mean_n_trial, y = .data$coverage, color = .data$method_label, linetype = .data$balance_mode, shape = .data$balance_mode)) +
      ggplot2::geom_hline(yintercept = 0.95, linetype = 2, color = "grey40") +
      ggplot2::geom_line(data = line_df, linewidth = 0.8) +
      ggplot2::geom_point(size = 2) +
      ggplot2::facet_grid(tau1 + sigma ~ n_id, labeller = ggplot2::label_both) +
      ggplot2::scale_x_continuous(breaks = sort(unique(plot_df$mean_n_trial))) +
      ggplot2::coord_cartesian(ylim = c(0.5, 1)) +
      ggplot2::labs(
        title = "Coverage for x predicting recovered random-slope outcomes",
        x = "Target mean trials per subject",
        y = "Coverage",
        color = "Method",
        linetype = "Balance",
        shape = "Balance"
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(legend.position = "bottom", legend.box = "vertical")

    p_bias <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$mean_n_trial, y = .data$bias, color = .data$method_label, linetype = .data$balance_mode, shape = .data$balance_mode)) +
      ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
      ggplot2::geom_line(data = line_df, linewidth = 0.8) +
      ggplot2::geom_point(size = 2) +
      ggplot2::facet_grid(tau1 + sigma ~ n_id, labeller = ggplot2::label_both) +
      ggplot2::scale_x_continuous(breaks = sort(unique(plot_df$mean_n_trial))) +
      ggplot2::labs(
        title = "Bias for x predicting recovered random-slope outcomes",
        x = "Target mean trials per subject",
        y = "Bias",
        color = "Method",
        linetype = "Balance",
        shape = "Balance"
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(legend.position = "bottom", legend.box = "vertical")

    ggplot2::ggsave(file.path(out_dir, sprintf("%s_coverage.png", prefix)), p_cov, width = 10.8, height = 6.6, units = "in", dpi = 300)
    ggplot2::ggsave(file.path(out_dir, sprintf("%s_bias.png", prefix)), p_bias, width = 10.8, height = 6.6, units = "in", dpi = 300)
  }

  invisible(list(summary = summary_df, issue_summary = issue_df))
}

#' List objects exported to parallel BLUP-outcome workers.
#'
#' @details
#' `foreach` workers run in separate R sessions, so locally defined functions
#' sourced by the driver must be exported explicitly. This list covers the
#' replication wrapper, estimator helpers, OpenMx/Lai helpers, stacked-sandwich
#' helpers, derivative-backend helpers, and utility functions used transitively
#' by `run_blup_outcome_rep()`.
#'
#' @return Character vector of object names passed to `foreach(..., .export)`.
blup_outcome_parallel_exports <- function() {
  c(
    "run_blup_outcome_rep", "blup_outcome_methods", "empty_blup_outcome_result",
    "standardize_estimator_rows", "fit_score_outcome_ols", "empty_stacked_rows",
    "fit_direct_mlm_row", "make_blup_outcome_diagnostics", "balance_mode_to_sim_arg",
    "draw_random_effects", "simulate_dataset", "safe_lmer", "safe_mean",
    "compact_message", "get_stage1_diagnostics", "empty_stage1_diagnostics",
    "get_corrected_scores", "get_diagonal_corrected_scores", "get_closed_form_corrected_scores",
    "fit_observed_single", "finalize_ols_se_variants", "extract_lmer_stats",
    "default_re_design", "make_eb_output_row", "compute_eb_measurement_inputs",
    "compute_bivariate_eb_inputs", "compute_lai_2spa_inputs", "fit_lai_2spa",
    "run_mx_safe", "extract_mx_stats", "extract_mx_se_details", "classify_mx_issue",
    "pack_psi", "pack_psi_var", "unpack_psi", "unpack_psi_var",
    "prepare_cluster_objects", "corrected_slope_from_precomputed",
    "cluster_loglik_precomputed", "total_loglik_precomputed",
    "get_stage1_sandwich_inputs", "stacked_sandwich_for_corrected_scores",
    "format_stacked_sandwich_rows", "make_derivative_backend",
    "ensure_tmb_stage1_dll", "make_tmb_stage1_data", "make_tmb_stage1_object",
    "get_tmb_stage1_hessian", "project_to_pd"
  )
}

#' Run all replications for one BLUP-outcome design condition.
#'
#' @details
#' Each replication is seeded deterministically from a fixed base seed, the
#' condition id, and the replication id. This makes serial, parallel, chunked,
#' and resumed runs reproducible as long as the same condition ids and
#' replication counts are used.
#'
#' When `n_cores > 1`, the function assumes a parallel backend has already been
#' registered by `run_blup_outcome_simulation()` and uses `%dopar%`; otherwise it
#' runs serially with `purrr::map_dfr()`.
#'
#' @param condition One-row design condition tibble.
#' @param n_sim Positive integer number of replications.
#' @param n_cores Positive integer number of cores requested.
#' @param params List of fixed simulation parameters passed to
#'   `run_blup_outcome_rep()`.
#' @param derivative_backend Derivative backend object returned by
#'   `make_derivative_backend()`.
#' @param analysis_mode Character scalar, either `"screen"` or `"full"`.
#'
#' @return Replication-level result tibble with condition columns and a `rep`
#'   identifier appended.
run_blup_outcome_condition_replications <- function(condition, n_sim, n_cores, params, derivative_backend, analysis_mode) {
  rep_ids <- seq_len(n_sim)

  run_single_rep <- function(rep_id) {
    # Keep the seed formula independent of execution order so parallel and
    # serial runs produce identical draws for the same condition/rep pair.
    set.seed(20260423L + (as.integer(condition$condition_id[[1]]) * 100000L) + as.integer(rep_id))
    rep_out <- run_blup_outcome_rep(
      condition = condition,
      params = params,
      derivative_backend = derivative_backend,
      analysis_mode = analysis_mode
    )
    dplyr::bind_cols(
      rep_out,
      condition[rep(1L, nrow(rep_out)), , drop = FALSE],
      tibble::tibble(rep = rep_id)
    )
  }

  if (n_cores > 1L) {
    foreach::foreach(
      rep_id = rep_ids,
      .combine = dplyr::bind_rows,
      .inorder = FALSE,
      .packages = c("lme4", "MASS", "dplyr", "tidyr", "purrr", "tibble", "OpenMx", "sandwich", "glmnet"),
      .export = blup_outcome_parallel_exports()
    ) %dopar% {
      run_single_rep(rep_id)
    }
  } else {
    purrr::map_dfr(rep_ids, run_single_rep)
  }
}

#' Run or aggregate the BLUP-outcome simulation grid.
#'
#' @details
#' This is the main programmatic entry point for BLUP-as-outcome simulations. In
#' `execution_mode = "run"`, it selects the requested design grid, optionally
#' slices it to a chunk, runs each condition, and writes condition-level
#' replication/summary files plus a progress log. It deliberately does not
#' rewrite aggregate outputs in run mode, which makes concurrent chunk jobs safe
#' when sharing an output directory.
#'
#' In `execution_mode = "aggregate"`, it ignores chunk arguments, loads all
#' completed condition outputs for the selected grid, writes aggregate CSVs and
#' plots, prints the summary, and returns the aggregate objects.
#'
#' @param n_sim Positive integer number of replications per condition.
#' @param out_dir Character scalar output directory.
#' @param n_cores Positive integer number of worker cores. Values above 1 enable
#'   a `doParallel` backend.
#' @param derivative_method Character scalar derivative backend name accepted by
#'   `make_derivative_backend()`.
#' @param grid_mode Character scalar design-grid mode passed to
#'   `make_blup_outcome_design()`.
#' @param analysis_mode Character scalar. `"screen"` runs the lighter estimator
#'   set; `"full"` includes computationally heavier Lai and stacked-sandwich
#'   paths.
#' @param chunk_index Optional one-based chunk index for distributed runs.
#' @param chunk_size Optional number of conditions per chunk.
#' @param resume_existing Logical. If `TRUE`, skip condition files that already
#'   have both replication and summary artifacts.
#' @param execution_mode Character scalar, either `"run"` or `"aggregate"`.
#' @param max_conditions Optional cap on selected conditions before chunking.
#'
#' @return In run mode, invisibly returns the selected condition grid. In
#'   aggregate mode, invisibly returns the list from
#'   `write_blup_outcome_aggregate_outputs()`.
run_blup_outcome_simulation <- function(n_sim = 100L,
                                        out_dir,
                                        n_cores = 1L,
                                        derivative_method = "handcoded",
                                        grid_mode = "base",
                                        analysis_mode = "full",
                                        chunk_index = NA_integer_,
                                        chunk_size = NA_integer_,
                                        resume_existing = TRUE,
                                        execution_mode = "run",
                                        max_conditions = NA_integer_) {
  if (is.na(n_sim) || n_sim < 1L) {
    stop("`n_sim` must be a positive integer.")
  }
  if (is.na(n_cores) || n_cores < 1L) {
    stop("`n_cores` must be a positive integer.")
  }
  if (!(analysis_mode %in% c("full", "screen"))) {
    stop("`analysis_mode` must be one of: full, screen.")
  }
  if (!(execution_mode %in% c("run", "aggregate"))) {
    stop("`execution_mode` must be one of: run, aggregate.")
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  derivative_backend <- make_derivative_backend(method = derivative_method)
  full_condition_grid <- make_blup_outcome_design(grid_mode = grid_mode, max_conditions = max_conditions)

  # Run mode works on one chunk; aggregate mode always inspects the full
  # selected design so it can rebuild outputs after all chunks finish.
  condition_grid <- if (identical(execution_mode, "run")) {
    slice_condition_chunk(full_condition_grid, chunk_index = chunk_index, chunk_size = chunk_size)
  } else {
    slice_condition_chunk(full_condition_grid)
  }
  chunk_meta <- attr(condition_grid, "chunk_meta")
  chunk_label <- make_chunk_label(chunk_meta)
  file_prefix <- sprintf("blup_outcome_%s", chunk_label)

  write_csv_atomic(condition_grid, file.path(out_dir, sprintf("%s_manifest.csv", file_prefix)))
  progress_path <- file.path(out_dir, sprintf("%s_progress.csv", file_prefix))

  params <- list(
    beta_0 = 1.0,
    beta_z = 0.6,
    tau0 = 0.9
  )

  if (identical(execution_mode, "aggregate")) {
    # Aggregation is intentionally read-only with respect to condition files:
    # it reconstructs all top-level artifacts from completed per-condition CSVs.
    results <- load_completed_condition_results(
      full_condition_grid,
      out_dir,
      path_fun = blup_outcome_condition_paths
    )
    out <- write_blup_outcome_aggregate_outputs(results, out_dir, prefix = "mlm_random_slope_blup_outcome")
    message("Aggregated condition outputs from: ", normalizePath(out_dir))
    message("Conditions represented: ", if ("condition_id" %in% names(results)) dplyr::n_distinct(results$condition_id) else 0L)
    print(out$summary)
    return(invisible(out))
  }

  if (n_cores > 1L) {
    doParallel::registerDoParallel(cores = n_cores)
    on.exit(doParallel::stopImplicitCluster(), add = TRUE)
  }

  for (condition_idx in seq_len(nrow(condition_grid))) {
    condition <- condition_grid[condition_idx, , drop = FALSE]
    condition_id <- condition$condition_id[[1]]
    condition_paths <- blup_outcome_condition_paths(out_dir, condition_id)
    started_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

    # Resume requires both detailed replication rows and the condition summary;
    # if either is missing, rerun the condition and overwrite partial artifacts.
    if (isTRUE(resume_existing) && file.exists(condition_paths$replications) && file.exists(condition_paths$summary)) {
      write_progress_row(
        progress_path,
        tibble::tibble(
          condition_id = condition_id,
          chunk_label = chunk_label,
          status = "skipped_existing",
          started_at = started_at,
          finished_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
          elapsed_seconds = 0,
          analysis_mode = analysis_mode,
          replications_file = condition_paths$replications,
          summary_file = condition_paths$summary
        )
      )
      next
    }

    write_progress_row(
      progress_path,
      tibble::tibble(
        condition_id = condition_id,
        chunk_label = chunk_label,
        status = "started",
        started_at = started_at,
        finished_at = NA_character_,
        elapsed_seconds = NA_real_,
        analysis_mode = analysis_mode,
        replications_file = condition_paths$replications,
        summary_file = condition_paths$summary
      )
    )

    # The console message is deliberately concise but includes all design axes
    # needed to identify slow or fragile cells in long batch logs.
    message(
      sprintf(
        "Condition %d/%d: n_id=%s, mean_n_trial=%s, gamma=%s, rho=%s, balance=%s, tau1=%s, sigma=%s",
        condition_idx,
        nrow(condition_grid),
        condition$n_id[[1]],
        condition$mean_n_trial[[1]],
        condition$gamma_x_on_slope[[1]],
        condition$rho[[1]],
        condition$balance_mode[[1]],
        condition$tau1[[1]],
        condition$sigma[[1]]
      )
    )

    timing <- system.time({
      condition_results <- run_blup_outcome_condition_replications(
        condition = condition,
        n_sim = n_sim,
        n_cores = n_cores,
        params = params,
        derivative_backend = derivative_backend,
        analysis_mode = analysis_mode
      )
    })
    condition_summary <- summarize_blup_outcome_results(condition_results)

    # Condition-level outputs are written immediately so interrupted jobs can be
    # resumed without losing completed simulation cells.
    write_csv_atomic(condition_results, condition_paths$replications)
    write_csv_atomic(condition_summary, condition_paths$summary)

    write_progress_row(
      progress_path,
      tibble::tibble(
        condition_id = condition_id,
        chunk_label = chunk_label,
        status = "completed",
        started_at = started_at,
        finished_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        elapsed_seconds = unname(timing[["elapsed"]]),
        analysis_mode = analysis_mode,
        n_methods = nrow(condition_results),
        n_success = sum(!is.na(condition_results$estimate)),
        replications_file = condition_paths$replications,
        summary_file = condition_paths$summary
      )
    )
  }

  message("Saved condition outputs to: ", normalizePath(out_dir))
  message("Derivative backend: ", derivative_backend$name)
  message("Chunk label: ", chunk_label)
  message("Run aggregation separately with execution_mode='aggregate'.")
  invisible(condition_grid)
}
