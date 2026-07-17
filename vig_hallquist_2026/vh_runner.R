#' Generic runner and output aggregation for Vig-Hallquist (2026).
#'
#' Study-specific modules define the data generators and estimators; the runner
#' is responsible for selecting design conditions, executing replications,
#' writing per-condition artifacts, and rebuilding aggregate summaries from
#' completed condition files.

vh_pipeline_version <- function() {
  "fuller_study4_diagnostics_v2_20260717"
}

#' Parse an optional command-line integer argument.
#'
#' @details
#' Shell scripts and batch schedulers often pass optional values as empty
#' strings or literal missing-value tokens. This helper normalizes those
#' sentinel values to `NA_integer_` so downstream validation can treat omitted
#' and explicitly missing inputs the same way.
#'
#' @param x A scalar-like value, usually a character vector element from
#'   `commandArgs()`.
#'
#' @return An integer scalar, or `NA_integer_` when `x` is absent or one of the
#'   recognized missing-value tokens.
parse_optional_integer_arg <- function(x) {
  # Accept the common string encodings used by shell wrappers for omitted args.
  if (is.null(x) || length(x) == 0L || is.na(x) || x %in% c("", "NA", "NaN", "NULL", "null")) {
    return(NA_integer_)
  }
  as.integer(x)
}

#' Parse a command-line logical flag.
#'
#' @param x A scalar-like value from `commandArgs()`.
#' @param default Logical value returned for omitted or missing-value tokens.
#'
#' @return A logical scalar.
parse_logical_arg <- function(x, default = FALSE) {
  if (is.null(x) || length(x) == 0L || is.na(x) || x %in% c("", "NA", "NaN", "NULL", "null")) {
    return(default)
  }
  if (is.logical(x)) {
    return(isTRUE(x))
  }
  x_chr <- tolower(trimws(as.character(x[[1]])))
  if (x_chr %in% c("1", "true", "t", "yes", "y")) {
    return(TRUE)
  }
  if (x_chr %in% c("0", "false", "f", "no", "n")) {
    return(FALSE)
  }
  stop("Could not parse logical argument: ", x)
}

#' Restrict a simulation design to one chunk.
#'
#' @details
#' Chunking lets large condition grids be split across independent jobs without
#' changing the underlying design identifiers. The input design is sliced by
#' row position after `select_design()` has applied any study and
#' `max_conditions` filters. The returned data frame carries a `chunk_meta`
#' attribute used for file naming and progress reporting.
#'
#' If both chunk arguments are `NA`, the design is returned unchanged and the
#' metadata describes the full selected design. If one chunk argument is
#' supplied, the other must be supplied too.
#'
#' @param design A data frame or tibble returned by `select_design()`. Must
#'   contain `condition_id`.
#' @param chunk_index One-based chunk number to run, or `NA_integer_` for the
#'   full selected design.
#' @param chunk_size Number of design rows per chunk, or `NA_integer_` for the
#'   full selected design.
#'
#' @return A sliced design with a `chunk_meta` attribute containing
#'   `chunk_index`, `chunk_size`, `condition_start`, `condition_end`, and
#'   `n_conditions`.
slice_design_chunk <- function(design, chunk_index = NA_integer_, chunk_size = NA_integer_) {
  if (is.na(chunk_index) && is.na(chunk_size)) {
    attr(design, "chunk_meta") <- list(
      chunk_index = NA_integer_,
      chunk_size = NA_integer_,
      condition_start = min(design$condition_id),
      condition_end = max(design$condition_id),
      n_conditions = nrow(design)
    )
    return(design)
  }

  if (is.na(chunk_index) || is.na(chunk_size)) {
    stop("`chunk_index` and `chunk_size` must be supplied together.")
  }
  if (chunk_index < 1L || chunk_size < 1L) {
    stop("`chunk_index` and `chunk_size` must be positive integers.")
  }

  # Chunks are contiguous 1-based windows over the selected design, not over
  # the unfiltered full factorial grid.
  start_idx <- ((chunk_index - 1L) * chunk_size) + 1L
  end_idx <- min(nrow(design), chunk_index * chunk_size)
  if (start_idx > nrow(design)) {
    out <- design[0, , drop = FALSE]
    attr(out, "chunk_meta") <- list(
      chunk_index = chunk_index,
      chunk_size = chunk_size,
      condition_start = NA_integer_,
      condition_end = NA_integer_,
      n_conditions = 0L
    )
    return(out)
  }

  out <- design %>%
    dplyr::slice(start_idx:end_idx)

  attr(out, "chunk_meta") <- list(
    chunk_index = chunk_index,
    chunk_size = chunk_size,
    condition_start = min(out$condition_id),
    condition_end = max(out$condition_id),
    n_conditions = nrow(out)
  )

  out
}

#' Create the stable label used for chunk-specific aggregate files.
#'
#' @param chunk_meta The `chunk_meta` attribute produced by
#'   `slice_design_chunk()`.
#'
#' @return A character scalar. Full-design runs return `"full_selection"`;
#'   chunked runs return a label containing the zero-padded chunk number and
#'   first/last condition identifiers.
make_chunk_label <- function(chunk_meta) {
  if (is.null(chunk_meta) || is.na(chunk_meta$chunk_index)) {
    "full_selection"
  } else {
    sprintf(
      "chunk_%03d_conditions_%04d_%04d",
      chunk_meta$chunk_index,
      chunk_meta$condition_start,
      chunk_meta$condition_end
    )
  }
}

#' Build per-condition output paths.
#'
#' @details
#' Condition-level outputs live under `out_dir/conditions`. The replication
#' results are compressed because they can be large; summaries are small CSV
#' files intended for quick inspection.
#'
#' @param out_dir Root output directory for the simulation run.
#' @param condition_id Integer condition identifier.
#'
#' @return A named list with paths for `replications`, `summary`,
#'   `issue_summary`, and `stage1_summary`.
get_condition_file_paths <- function(out_dir, condition_id) {
  condition_dir <- file.path(out_dir, "conditions")
  dir.create(condition_dir, recursive = TRUE, showWarnings = FALSE)

  list(
    replications = file.path(condition_dir, sprintf("condition_%04d_replication_results.csv.gz", condition_id)),
    summary = file.path(condition_dir, sprintf("condition_%04d_summary.csv", condition_id)),
    issue_summary = file.path(condition_dir, sprintf("condition_%04d_issue_summary.csv", condition_id)),
    stage1_summary = file.path(condition_dir, sprintf("condition_%04d_stage1_problem_summary.csv", condition_id))
  )
}

#' Summarize replication-level estimator results.
#'
#' @details
#' The summary is computed at the condition-study-method level while preserving
#' the design parameters needed to compare rows across the study grids.
#' OpenMx status-code-10 failures are counted separately and are excluded from
#' the convergence flag even if an estimate column is present. Bias, RMSE, and
#' coverage are based on the standardized estimand stored in `truth`.
#'
#' @param results Replication-level result tibble produced by
#'   `run_condition_replications()` or loaded from per-condition files. Expected
#'   columns include `condition_id`, `study`, `method`, `estimate`, `truth`,
#'   confidence limits, status diagnostics, and the design descriptors.
#'
#' @return A tibble with one row per condition-study-method combination and
#'   Monte Carlo summary columns including convergence, mean estimate, MC SE of
#'   the mean, bias, coverage, RMSE, median and 95th-percentile absolute error,
#'   successful replications, and status-10 failure counts.
summarize_results_df <- function(results) {
  design_cols <- intersect(
    c(
      "condition_id", "study", "method", "method_role", "num_clus", "mean_clus_size",
      "target_reliability", "achieved_reliability", "marginal_rho",
      "standardized_beta_target", "structural_target", "structural_r2",
      "focal_unique_r2", "mean_clus_size_y", "mean_clus_size_q",
      "target_reliability_y", "target_reliability_q",
      "achieved_reliability_y", "achieved_reliability_q",
      "balance_mode", "r_structure", "r_rho", "sigma", "sigma_y", "sigma_q",
      "information_profile", "is_falsification_control", "information_matched",
      "profile_min_clus_size", "profile_max_clus_size",
      "reliability_sd", "reliability_iqr", "reliability_min",
      "reliability_max", "reliability_small", "reliability_large",
      "population_lambda22_mean", "population_lambda22_min",
      "population_lambda22_max",
      "population_lambda_matrix_frobenius_rms_dispersion",
      "population_theta22_mean", "population_theta22_min",
      "population_theta22_max",
      "population_theta_matrix_frobenius_rms_dispersion",
      "population_ols_var22_mean", "population_ols_var22_min",
      "population_ols_var22_max",
      "population_ols_cov_matrix_frobenius_rms_dispersion"
    ),
    names(results)
  )
  replication_diagnostic_cols <- intersect(
    c(
      "realized_reliability_mean", "realized_reliability_sd",
      "realized_reliability_iqr", "realized_reliability_min",
      "realized_reliability_max",
      "lambda22_mean", "lambda22_sd", "lambda22_min", "lambda22_max",
      "lambda_matrix_frobenius_rms_dispersion",
      "theta22_mean", "theta22_sd", "theta22_min", "theta22_max",
      "theta_matrix_frobenius_rms_dispersion",
      "ols_var22_mean", "ols_var22_sd", "ols_var22_min", "ols_var22_max",
      "ols_cov_matrix_frobenius_rms_dispersion",
      "blup_slope_bias_small", "blup_slope_rmse_small",
      "blup_slope_bias_large", "blup_slope_rmse_large",
      "corrected_slope_bias_small", "corrected_slope_rmse_small",
      "corrected_slope_bias_large", "corrected_slope_rmse_large",
      "average_measurement_slope_bias_small",
      "average_measurement_slope_rmse_small",
      "average_measurement_slope_bias_large",
      "average_measurement_slope_rmse_large",
      "fuller_measurement_weight_used",
      "fuller_alpha_step1_used", "fuller_alpha_step3_used",
      "fuller_alpha_scaling_used",
      "fuller_correction1", "fuller_correction_c",
      "fuller_correction_scaling",
      "fuller_sx1_star_condition", "fuller_sx1_star_min_eigen",
      "fuller_sx1_star_relative_min_eigen",
      "fuller_sx_star_condition", "fuller_sx_star_min_eigen",
      "fuller_sx_star_relative_min_eigen",
      "fuller_scaling_condition", "fuller_scaling_min_eigen",
      "fuller_scaling_relative_min_eigen",
      "fuller_auto_guard_pass", "fuller_auto_full_weight_guard_pass",
      "fuller_auto_search_evaluations"
    ),
    names(results)
  )
  results %>%
    dplyr::mutate(
      # Status 10 is an OpenMx optimizer failure class; keep it visible rather
      # than folding it into generic missing-estimate behavior.
      status10_failure = !is.na(status_code) & status_code == 10L,
      converged = !status10_failure & !is.na(estimate),
      bias = estimate - truth,
      abs_error = abs(estimate - truth),
      sq_error = (estimate - truth)^2,
      covered = ci_low <= truth & ci_high >= truth
    ) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(design_cols))) %>%
    dplyr::summarise(
      truth = dplyr::first(truth),
      n_rep = dplyr::n(),
      convergence = safe_mean(converged),
      mean_estimate = safe_mean(estimate),
      mc_se_mean = if (sum(!is.na(estimate)) > 1L) stats::sd(estimate, na.rm = TRUE) / sqrt(sum(!is.na(estimate))) else NA_real_,
      bias = safe_mean(bias),
      coverage = safe_mean(covered),
      rmse = if (all(is.na(sq_error))) NA_real_ else sqrt(mean(sq_error, na.rm = TRUE)),
      median_absolute_error = if (all(is.na(abs_error))) {
        NA_real_
      } else {
        stats::median(abs_error, na.rm = TRUE)
      },
      p95_absolute_error = if (all(is.na(abs_error))) {
        NA_real_
      } else {
        unname(stats::quantile(abs_error, probs = 0.95, na.rm = TRUE))
      },
      n_success = sum(converged, na.rm = TRUE),
      n_status10_fail = sum(status10_failure, na.rm = TRUE),
      prop_status10_fail = safe_mean(status10_failure),
      dplyr::across(
        dplyr::all_of(replication_diagnostic_cols),
        safe_mean,
        .names = "mean_{.col}"
      ),
      .groups = "drop"
    )
}

#' Summarize results by first-stage singularity/problem status.
#'
#' @details
#' Some simulation conditions can be dominated by first-stage mixed-model
#' singularity or design-matrix instability. When the replication results carry
#' those diagnostics, this helper stratifies the standard performance summaries
#' by `stage1_singular_problem` and adds aggregate first-stage diagnostic
#' measures. If the diagnostics are not present, it returns an empty tibble so
#' callers can always write a stable output file.
#'
#' @param results Replication-level result tibble.
#'
#' @return A tibble grouped by condition-study-method and
#'   `stage1_singular_problem`, or an empty tibble when first-stage problem
#'   diagnostics are unavailable.
summarize_stage1_problem_df <- function(results) {
  if (!("stage1_singular_problem" %in% names(results))) {
    return(tibble::tibble())
  }

  design_cols <- intersect(
    c(
      "condition_id", "study", "method", "method_role", "stage1_singular_problem",
      "num_clus", "mean_clus_size", "target_reliability",
      "achieved_reliability", "marginal_rho", "standardized_beta_target",
      "structural_target", "structural_r2", "focal_unique_r2",
      "mean_clus_size_y", "mean_clus_size_q",
      "target_reliability_y", "target_reliability_q",
      "achieved_reliability_y", "achieved_reliability_q",
      "balance_mode", "r_structure", "r_rho", "sigma", "sigma_y", "sigma_q",
      "information_profile", "is_falsification_control", "information_matched",
      "profile_min_clus_size", "profile_max_clus_size",
      "reliability_sd", "reliability_min", "reliability_max"
    ),
    names(results)
  )

  results %>%
    dplyr::mutate(
      status10_failure = !is.na(status_code) & status_code == 10L,
      converged = !status10_failure & !is.na(estimate),
      bias = estimate - truth,
      sq_error = (estimate - truth)^2,
      covered = ci_low <= truth & ci_high >= truth
    ) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(design_cols))) %>%
    dplyr::summarise(
      truth = dplyr::first(truth),
      n_rep = dplyr::n(),
      convergence = safe_mean(converged),
      mean_estimate = safe_mean(estimate),
      bias = safe_mean(bias),
      coverage = safe_mean(covered),
      rmse = if (all(is.na(sq_error))) NA_real_ else sqrt(mean(sq_error, na.rm = TRUE)),
      n_success = sum(converged, na.rm = TRUE),
      prop_status10_fail = safe_mean(status10_failure),
      prop_stage1_lmer_singular = safe_mean(stage1_lmer_singular),
      mean_stage1_re_corr = safe_mean(stage1_re_corr),
      mean_stage1_eb_corr = safe_mean(stage1_eb_corr),
      mean_stage1_design_kappa = safe_mean(stage1_design_kappa),
      sample_stage1_problem_detail = compact_message(stage1_problem_detail[!is.na(stage1_problem_detail)]),
      .groups = "drop"
    )
}

#' Summarize OpenMx and estimator issue classes.
#'
#' @details
#' Estimator wrappers classify OpenMx or standard-error problems into compact
#' `mx_issue_class` labels. This helper keeps only non-`ok` issue classes and
#' reports their frequency within each condition-study-method group, along with
#' a compact example detail string for inspection.
#'
#' @param results Replication-level result tibble.
#'
#' @return A tibble with issue-class counts and proportions, or an empty tibble
#'   when `mx_issue_class` is absent.
summarize_issue_df <- function(results) {
  if (!("mx_issue_class" %in% names(results))) {
    return(tibble::tibble())
  }

  design_cols <- intersect(
    c(
      "condition_id", "study", "method", "method_role", "mx_issue_class", "num_clus",
      "mean_clus_size", "target_reliability", "achieved_reliability",
      "marginal_rho", "standardized_beta_target", "structural_target",
      "structural_r2", "focal_unique_r2", "balance_mode", "r_structure",
      "r_rho", "sigma", "mean_clus_size_y", "mean_clus_size_q",
      "target_reliability_y", "target_reliability_q",
      "achieved_reliability_y", "achieved_reliability_q", "sigma_y", "sigma_q",
      "information_profile", "is_falsification_control", "information_matched",
      "profile_min_clus_size", "profile_max_clus_size",
      "reliability_sd", "reliability_min", "reliability_max"
    ),
    names(results)
  )

  results %>%
    dplyr::group_by(condition_id, study, method) %>%
    dplyr::mutate(n_rep_total = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(
      !is.na(mx_issue_class),
      nzchar(mx_issue_class),
      mx_issue_class != "ok"
    ) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(design_cols))) %>%
    dplyr::summarise(
      n_rep = dplyr::n(),
      n_rep_total = dplyr::first(n_rep_total),
      prop_rep = dplyr::n() / dplyr::first(n_rep_total),
      n_distinct_issue_details = dplyr::n_distinct(mx_issue_detail[!is.na(mx_issue_detail)]),
      sample_issue_detail = compact_message(mx_issue_detail[!is.na(mx_issue_detail)]),
      .groups = "drop"
    )
}

#' Append one row to a progress CSV.
#'
#' @details
#' The progress file is intentionally append-only so interrupted jobs retain a
#' record of completed and skipped conditions. The first write includes headers;
#' later writes append rows without column names.
#'
#' @param progress_path Path to the run-level progress CSV.
#' @param row_df One-row data frame or tibble with progress metadata.
#'
#' @return Invisibly returns the result of the underlying write call.
write_progress_row <- function(progress_path, row_df) {
  if (!file.exists(progress_path)) {
    utils::write.csv(row_df, file = progress_path, row.names = FALSE)
  } else {
    utils::write.table(
      row_df,
      file = progress_path,
      sep = ",",
      row.names = FALSE,
      col.names = FALSE,
      append = TRUE
    )
  }
}

#' Read and type-normalize a replication result file.
#'
#' @details
#' `data.table::fread()` is fast, but resumed aggregate builds can encounter
#' plain CSV files from older runs as well as current `.csv.gz` files. This
#' helper reads either form and coerces known numeric, integer, and logical
#' columns back to their expected types before row-binding across conditions.
#'
#' @param path Path to a per-condition replication result CSV or CSV.GZ file.
#'
#' @return A tibble containing the replication-level rows from `path`.
read_replication_results_file <- function(path) {
  out <- tibble::as_tibble(data.table::fread(path))
  fuller_logical_cols <- c(
    "fuller_auto_tempered", "fuller_auto_guard_pass",
    "fuller_auto_full_weight_guard_pass", "fuller_auto_search_nonmonotone"
  )
  fuller_integer_cols <- "fuller_auto_search_evaluations"
  fuller_character_cols <- c(
    "fuller_auto_guard_reason", "fuller_auto_full_weight_guard_reason"
  )
  fuller_numeric_cols <- setdiff(
    grep("^fuller_", names(out), value = TRUE),
    c(fuller_logical_cols, fuller_integer_cols, fuller_character_cols)
  )

  # Keep these casts centralized so aggregate rebuilds do not depend on the
  # exact type inference chosen by fread for any single condition file.
  numeric_cols <- unique(c(intersect(
    # TODO: update these columns
    c(
      "estimate", "se", "ci_low", "ci_high", "truth",
      "mx_condition_number",
      "stage1_re_corr", "stage1_eb_corr", "stage1_design_kappa",
      "stage1_y_re_corr", "stage1_y_eb_corr", "stage1_y_design_kappa",
      "stage1_q_re_corr", "stage1_q_eb_corr", "stage1_q_design_kappa",
      "target_reliability", "achieved_reliability", "marginal_rho",
      "standardized_beta_target", "structural_r2", "focal_unique_r2",
      "tau1", "beta1z", "beta2z", "outcome_residual_variance",
      "target_reliability_y", "target_reliability_q",
      "achieved_reliability_y", "achieved_reliability_q",
      "tau1_y", "tau1_q", "theta0", "theta1",
      "standardized_theta0", "slope_variance_marginal_y",
      "slope_variance_marginal_q", "slope_variance_residual_q",
      "tau1_residual_q", "rho_residual_q", "sigma_y", "sigma_q",
      "sigma2", "var_u1", "sigma_z", "fuller_lambda1", "fuller_lambda2",
      "fuller_sigma2", "fuller_weight_min", "fuller_weight_max",
      "fuller_correction_c", "profile_small_weight", "profile_large_weight",
      "reliability_sd", "reliability_iqr", "reliability_min", "reliability_max",
      "reliability_small", "reliability_large",
      "population_lambda22_mean", "population_lambda22_min",
      "population_lambda22_max",
      "population_lambda_matrix_frobenius_rms_dispersion",
      "population_theta22_mean", "population_theta22_min",
      "population_theta22_max",
      "population_theta_matrix_frobenius_rms_dispersion",
      "population_ols_var22_mean", "population_ols_var22_min",
      "population_ols_var22_max",
      "population_ols_cov_matrix_frobenius_rms_dispersion",
      "realized_reliability_mean", "realized_reliability_sd",
      "realized_reliability_iqr", "realized_reliability_min",
      "realized_reliability_max",
      "mean_realized_trials", "min_realized_trials", "max_realized_trials",
      "prop_ids_leq_3_trials",
      "lambda22_mean", "lambda22_sd", "lambda22_min", "lambda22_max",
      "lambda_matrix_frobenius_rms_dispersion",
      "theta22_mean", "theta22_sd", "theta22_min", "theta22_max",
      "theta_matrix_frobenius_rms_dispersion",
      "ols_var22_mean", "ols_var22_sd", "ols_var22_min", "ols_var22_max",
      "ols_cov_matrix_frobenius_rms_dispersion",
      "blup_slope_bias_small", "blup_slope_rmse_small",
      "blup_slope_bias_large", "blup_slope_rmse_large",
      "corrected_slope_bias_small", "corrected_slope_rmse_small",
      "corrected_slope_bias_large", "corrected_slope_rmse_large",
      "average_measurement_slope_bias_small",
      "average_measurement_slope_rmse_small",
      "average_measurement_slope_bias_large",
      "average_measurement_slope_rmse_large"
    ),
    names(out)
  ), fuller_numeric_cols))
  integer_cols <- intersect(
    c(
      "condition_id", "rep", "num_clus", "mean_clus_size",
      "mean_clus_size_y", "mean_clus_size_q", "status_code",
      "profile_min_clus_size", "profile_max_clus_size",
      "profile_small_clus_size", "profile_large_clus_size",
      "score_small_cluster_size", "score_large_cluster_size",
      fuller_integer_cols
    ),
    names(out)
  )
  logical_cols <- intersect(
    c(
      "stage1_singular_problem", "stage1_lmer_singular",
      "stage1_y_singular_problem", "stage1_y_lmer_singular",
      "stage1_q_singular_problem", "stage1_q_lmer_singular",
      "is_falsification_control", "information_matched",
      fuller_logical_cols
    ),
    names(out)
  )

  if (length(numeric_cols) > 0L) {
    out <- out %>% dplyr::mutate(dplyr::across(dplyr::all_of(numeric_cols), ~ suppressWarnings(as.numeric(.x))))
  }
  if (length(integer_cols) > 0L) {
    out <- out %>% dplyr::mutate(dplyr::across(dplyr::all_of(integer_cols), ~ suppressWarnings(as.integer(.x))))
  }
  if (length(logical_cols) > 0L) {
    out <- out %>% dplyr::mutate(dplyr::across(dplyr::all_of(logical_cols), as.logical))
  }

  out
}

#' Load all completed condition-level replication results for a design.
#'
#' @details
#' This is used after the condition loop to rebuild aggregate outputs from the
#' per-condition artifacts on disk. Rebuilding from disk makes `resume_existing`
#' runs deterministic: skipped conditions and newly completed conditions are
#' combined through the same code path.
#'
#' @param design Design rows that should be represented in the aggregate output.
#'   Must contain `condition_id`.
#' @param out_dir Root output directory containing the `conditions` subdirectory.
#'
#' @return A tibble formed by row-binding all discovered condition replication
#'   files, or an empty tibble if none are present.
load_completed_condition_results <- function(design, out_dir) {
  files <- purrr::map_chr(design$condition_id, function(condition_id) {
    gz_path <- get_condition_file_paths(out_dir, condition_id)$replications
    csv_path <- sub("\\.gz$", "", gz_path)
    if (file.exists(gz_path)) {
      gz_path
    } else if (file.exists(csv_path)) {
      csv_path
    } else {
      NA_character_
    }
  })
  existing <- stats::na.omit(files)
  if (length(existing) == 0L) {
    return(tibble::tibble())
  }

  purrr::map_dfr(existing, read_replication_results_file)
}

load_completed_condition_csv <- function(design, out_dir, artifact) {
  if (!(artifact %in% c("summary", "issue_summary", "stage1_summary"))) {
    stop("Unsupported condition artifact: ", artifact)
  }
  files <- purrr::map_chr(design$condition_id, function(condition_id) {
    path <- get_condition_file_paths(out_dir, condition_id)[[artifact]]
    if (file.exists(path)) path else NA_character_
  })
  existing <- stats::na.omit(files)
  if (length(existing) == 0L) {
    return(tibble::tibble())
  }

  purrr::map_dfr(existing, function(path) {
    tryCatch(
      tibble::as_tibble(utils::read.csv(path, check.names = FALSE)),
      error = function(e) tibble::tibble()
    )
  })
}

make_condition_replication_index <- function(design, out_dir) {
  purrr::map_dfr(design$condition_id, function(condition_id) {
    path <- get_condition_file_paths(out_dir, condition_id)$replications
    tibble::tibble(
      condition_id = condition_id,
      replication_file = path,
      file_exists = file.exists(path),
      file_size_bytes = if (file.exists(path)) file.info(path)$size else NA_real_
    )
  })
}

estimate_replication_result_rows <- function(design, n_sim) {
  method_counts <- vapply(
    seq_len(nrow(design)),
    function(i) length(study_methods_for_condition(design[i, , drop = FALSE])),
    integer(1)
  )
  as.double(n_sim) * sum(method_counts)
}

condition_output_is_complete <- function(
    path,
    expected_methods,
    expected_n_sim,
    expected_pipeline_version = vh_pipeline_version()) {
  if (!file.exists(path)) {
    return(FALSE)
  }
  out <- tryCatch(read_replication_results_file(path), error = function(e) NULL)
  required <- c("method", "rep", "pipeline_version")
  if (is.null(out) || !all(required %in% names(out))) {
    return(FALSE)
  }
  if (!all(out$pipeline_version == expected_pipeline_version)) {
    return(FALSE)
  }
  if (!setequal(unique(out$method), expected_methods) ||
      !setequal(unique(out$rep), seq_len(expected_n_sim))) {
    return(FALSE)
  }
  method_counts <- table(factor(out$method, levels = expected_methods))
  all(method_counts == expected_n_sim)
}

#' Dispatch one replication to the correct study module.
#'
#' @param condition One-row condition tibble with a `study` column equal to
#'  `"study0"`, `"study1"`, `"study2"`, `"study3"`, or `"study4"`.
#'
#' @return A replication-level result tibble from the corresponding
#'   `run_study*_rep()` function.
run_study_rep <- function(condition) {
  study_key <- as.character(condition$study[[1]])
  switch(
    study_key,
    study0 = run_study1_rep(condition),
    study1 = run_study1_rep(condition),
    study2 = run_study2_rep(condition),
    study3 = run_study3_rep(condition),
    study4 = run_study4_rep(condition),
    stop("Unsupported VH study key: ", study_key)
  )
}

#' List objects that must be exported to parallel workers.
#'
#' @details
#' `foreach` workers run in separate R processes and therefore need explicit
#' access to locally defined functions that are not package exports. This list
#' includes study simulators, estimator helpers, diagnostics, and OpenMx
#' utilities used anywhere below `run_study_rep()`.
#'
#' @return A character vector of object names passed to `foreach(..., .export)`.
vig_hallquist_parallel_exports <- function() {
  candidates <- c(
    "run_study_rep", "run_study1_rep", "run_study2_rep", "run_study3_rep", "run_study4_rep",
    "vh_pipeline_version",
    "add_study2_method_roles", "add_study4_method_roles",
    "rescale_fuller_to_population_sd",
    "prepare_fuller_average_measurement", "fit_fuller_average_measurement",
    "simulate_study1", "simulate_study2", "simulate_study3", "simulate_study4",
    "simulate_data_blup_as_outcome", "simulate_data_blup_as_predictor",
    "simulate_data_dual_blup", "simulate_data_study4",
    "make_failed_result", "add_study_result_context",
    "combine_dual_stage1_diagnostics",
    "fit_stage1", "condition_uses_non_iid_R", "condition_to_r_spec",
    "condition_to_nlme_correlation", "draw_level1_residuals",
    "balance_mode_to_sim_arg", "make_reliability_time_design",
    "study4_profile_spec", "study4_time_design", "make_study4_cluster_sizes",
    "study4_weighted_quantile", "study4_matrix_rms_dispersion",
    "study4_measurement_matrix_summary", "study4_measurement_diagnostics",
    "draw_random_effects",
    "matched_study_methods", "disparate_study_methods", "study_methods_for_condition", "tempered_eiv_methods",
    "condition_includes_tempered_eiv", "lai_condition_to_r_spec", "lai_condition_uses_non_iid_R",
    "lai_condition_to_nlme_correlation", "add_lai_trial_index", "draw_lai_level1_residuals",
    "fit_lai_stage1", "extract_centered_slope_eb", "fit_tempered_eiv_dual_set", "lai_truth", "make_covu",
    "make_study2_cluster_sizes", "make_study3_cluster_sizes", "fixed_params",
    "study1_methods", "study2_methods", "study3_methods", "study4_methods",
    "safe_lmer", "safe_lme", "empty_stage1_diagnostics", "get_stage1_diagnostics",
    "normalize_r_spec", "make_R_matrix", "draw_residuals_from_R",
    "make_random_effect_covariance", "posterior_random_effect_covariance",
    "get_closed_form_corrected_scores", "get_stage1_eb_components",
    "extract_stage1_components", "extract_stage1_components.merMod", "extract_stage1_components.lme",
    "extract_stage1_components.default", "normalize_R_list", "as_plain_vcov_matrix", "stage1_fixef",
    "format_stage1_eb_row", "default_re_design", "make_eb_output_row", "select_lai_measurement_columns",
    "compute_eb_measurement_inputs", "compute_bivariate_eb_inputs",
    "compute_univariate_eb_inputs", "compute_lai_2spa_inputs",
    "fit_observed_single", "fit_observed_dual",
    "finalize_ols_se_variants", "fit_eiv_dual", "fit_ridge_dual", "fit_fuller",
    "fit_fuller_dual",
    "fit_fuller_dual_stepdown", "fit_fuller_dual_alpha_stepdown",
    "fit_lai_2spa", "fit_lai_2spa_observed_outcome", "fit_lai_2spa_disparate",
    "fit_lai_2spa_dual_process",
    "run_mx_safe", "extract_mx_stats", "extract_mx_se_details",
    "classify_mx_issue", "compact_message", "project_to_pd", "fit_mplus_blup_predictor",
    "extract_mplus_stats"
  )
  candidates[vapply(candidates, exists, logical(1), mode = "function", inherits = TRUE) |
    candidates == "fixed_params"]
}

#' Run all replications for one design condition.
#'
#' @details
#' Each replication receives a deterministic seed derived from a fixed base
#' seed, the condition identifier, and the replication index. This keeps serial,
#' parallel, and resumed condition runs reproducible as long as the same
#' condition IDs and replication counts are used.
#'
#' Parallel execution uses `foreach` and expects a registered backend, which is
#' set up by `run_simulation()` when `n_cores > 1`.
#'
#' @param condition One-row condition tibble.
#' @param n_sim Number of replications to run for the condition.
#' @param n_cores Number of cores requested. Values greater than one select the
#'   `foreach` path.
#'
#' @return A tibble with one row per estimator/method result per replication,
#'   augmented with condition descriptors and the replication id.
run_condition_replications <- function(condition, n_sim, n_cores = 1L) {
  rep_ids <- seq_len(n_sim)

  run_single_rep <- function(rep_id) {
    # The large condition multiplier avoids seed collisions across conditions
    # even for high replication counts.
    set.seed(20260612 + (as.integer(condition$condition_id) * 100000L) + as.integer(rep_id))
    rep_out <- run_study_rep(condition)
    dplyr::bind_cols(
      rep_out,
      condition[rep(1L, nrow(rep_out)), , drop = FALSE] %>% dplyr::select(-study),
      tibble::tibble(
        rep = rep_id,
        pipeline_version = vh_pipeline_version()
      )
    )
  }

  if (n_cores > 1L) {
    foreach::foreach(
      rep_id = rep_ids,
      .combine = dplyr::bind_rows,
      # TODO: check these exports
      .packages = c("data.table", "lme4", "MASS", "dplyr", "tidyr", "purrr", "tibble", "OpenMx",
                    "glmnet", "sandwich", "geigen", "MplusAutomation", "glue"),
      .export = vig_hallquist_parallel_exports()
    ) %dopar% {
      run_single_rep(rep_id)
    }
  } else {
    purrr::map_dfr(rep_ids, run_single_rep)
  }
}

#' Run a Vig-Hallquist simulation batch.
#'
#' @details
#' This is the top-level entry point used by the replication script. It
#' selects the requested study design, optionally slices the design into a
#' chunk, runs each condition, writes condition-level artifacts immediately, and
#' finally rebuilds aggregate replication and summary files from the completed
#' condition outputs.
#'
#' The per-condition write pattern is deliberate: long simulation jobs can be
#' interrupted and later resumed with `resume_existing = TRUE`. A condition is
#' skipped only when both its replication file and main summary file already
#' exist. After the loop, aggregate outputs are rebuilt from all condition files
#' found for the selected design/chunk.
#'
#' Output files include:
#' - A manifest CSV with the selected design rows.
#' - An append-only progress CSV.
#' - Per-condition replication, summary, issue-summary, and first-stage
#'   problem-summary files under `conditions/`.
#' - Chunk-specific aggregate files. Full-selection runs also write legacy
#'   aggregate filenames without the chunk label.
#'
#' @param n_sim Positive integer number of replications per condition.
#' @param study_arg Study selector passed to `select_design()`, typically
#'   `"all"`, `"study0"`, `"study1"`, `"study2"`, `"study3"`, or `"study4"`.
#' @param out_dir Root output directory. Defaults to
#'   `file.path(vig_hallquist_dir, "outputs", "vig_hallquist")`.
#' @param n_cores Positive integer number of worker cores. Values greater than
#'   one register a `doParallel` backend and use the parallel replication path.
#' @param max_conditions Optional cap passed to `select_design()` after study
#'   selection.
#' @param chunk_index Optional one-based chunk index for distributed runs. Must
#'   be supplied with `chunk_size`.
#' @param chunk_size Optional number of selected design rows per chunk. Must be
#'   supplied with `chunk_index`.
#' @param resume_existing Logical. If `TRUE`, skip conditions whose
#'   per-condition replication and summary files already exist.
#' @param max_aggregate_replication_rows Maximum expected rows that may be
#'   loaded to materialize a combined replication table. Larger runs aggregate
#'   condition summaries and write a replication-file index instead.
#'
#' @return Invisibly returns a list with aggregate `results`, `summary`,
#'   `issue_summary`, and `stage1_summary` tibbles.
run_simulation <- function(n_sim = 100L,
                           study_arg = "all",
                           out_dir = file.path(vig_hallquist_dir, "outputs", "vig_hallquist"),
                           n_cores = 1L,
                           max_conditions = NA_integer_,
                           chunk_index = NA_integer_,
                           chunk_size = NA_integer_,
                           resume_existing = TRUE,
                           max_aggregate_replication_rows = 2e6) {
  if (is.na(n_sim) || n_sim < 1L) {
    stop("`n_sim` must be a positive integer.")
  }
  if (is.na(n_cores) || n_cores < 1L) {
    stop("`n_cores` must be a positive integer.")
  }
  if (is.na(max_aggregate_replication_rows) ||
      max_aggregate_replication_rows <= 0) {
    stop("`max_aggregate_replication_rows` must be positive or `Inf`.")
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  design <- select_design(study_arg = study_arg, max_conditions = max_conditions)
  design <- slice_design_chunk(design, chunk_index = chunk_index, chunk_size = chunk_size)
  if (nrow(design) == 0L) {
    message(
      "No conditions in requested chunk ", chunk_index,
      " for study selection `", study_arg, "`; exiting successfully."
    )
    return(invisible(list(
      results = tibble::tibble(),
      summary = tibble::tibble(),
      issue_summary = tibble::tibble(),
      stage1_summary = tibble::tibble()
    )))
  }
  chunk_meta <- attr(design, "chunk_meta")
  chunk_label <- make_chunk_label(chunk_meta)
  selection_label <- gsub(
    "[^a-z0-9]+",
    "_",
    tolower(as.character(study_arg[[1]]))
  )
  file_prefix <- sprintf(
    "vig_hallquist_%s_%s",
    selection_label,
    chunk_label
  )

  # Chunk labels keep concurrently run jobs from overwriting each other's
  # aggregate files while preserving condition-level paths shared by resume.
  manifest_path <- file.path(out_dir, sprintf("%s_manifest.csv", file_prefix))
  progress_path <- file.path(out_dir, sprintf("%s_progress.csv", file_prefix))
  aggregate_replications_path <- file.path(out_dir, sprintf("%s_replication_results.csv.gz", file_prefix))
  aggregate_summary_path <- file.path(out_dir, sprintf("%s_summary.csv", file_prefix))
  aggregate_issue_summary_path <- file.path(out_dir, sprintf("%s_issue_summary.csv", file_prefix))
  aggregate_stage1_summary_path <- file.path(out_dir, sprintf("%s_stage1_problem_summary.csv", file_prefix))
  replication_index_path <- file.path(
    out_dir,
    sprintf("%s_condition_replication_index.csv", file_prefix)
  )

  utils::write.csv(design, file = manifest_path, row.names = FALSE)

  if (n_cores > 1L) {
    doParallel::registerDoParallel(cores = n_cores)
    on.exit(doParallel::stopImplicitCluster(), add = TRUE)
  }

  for (i in seq_len(nrow(design))) {
    condition <- design[i, , drop = FALSE]
    condition_id <- condition$condition_id[[1]]
    paths <- get_condition_file_paths(out_dir, condition_id)
    started_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

    # A condition is resumable only when both the detailed rows and primary
    # summary exist; partial artifacts are overwritten by rerunning the condition.
    if (isTRUE(resume_existing) && file.exists(paths$replications) && file.exists(paths$summary)) {
      expected_methods <- study_methods_for_condition(condition)
      if (condition_output_is_complete(
          paths$replications,
          expected_methods = expected_methods,
          expected_n_sim = n_sim)) {
        write_progress_row(
          progress_path,
          tibble::tibble(
            condition_id = condition_id,
            chunk_label = chunk_label,
            status = "skipped_existing",
            started_at = started_at,
            finished_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            elapsed_seconds = 0,
            n_methods = NA_integer_,
            n_success = NA_integer_,
            replications_file = paths$replications,
            summary_file = paths$summary
          )
        )
        next
      }
    }

    timing <- system.time({
      condition_results <- run_condition_replications(condition, n_sim = n_sim, n_cores = n_cores)
    })
    condition_summary <- summarize_results_df(condition_results)
    condition_issue_summary <- summarize_issue_df(condition_results)
    condition_stage1_summary <- summarize_stage1_problem_df(condition_results)

    data.table::fwrite(condition_results, file = paths$replications)
    utils::write.csv(condition_summary, file = paths$summary, row.names = FALSE)
    utils::write.csv(condition_issue_summary, file = paths$issue_summary, row.names = FALSE)
    utils::write.csv(condition_stage1_summary, file = paths$stage1_summary, row.names = FALSE)

    write_progress_row(
      progress_path,
      tibble::tibble(
        condition_id = condition_id,
        chunk_label = chunk_label,
        status = "completed",
        started_at = started_at,
        finished_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        elapsed_seconds = unname(timing[["elapsed"]]),
        n_methods = nrow(condition_results),
        n_success = sum(!is.na(condition_results$estimate)),
        replications_file = paths$replications,
        summary_file = paths$summary
      )
    )
  }

  expected_result_rows <- estimate_replication_result_rows(design, n_sim)
  materialize_replications <- is.infinite(max_aggregate_replication_rows) ||
    expected_result_rows <= max_aggregate_replication_rows

  if (materialize_replications) {
    results <- load_completed_condition_results(design, out_dir)
    summary_df <- if (nrow(results) > 0L) summarize_results_df(results) else tibble::tibble()
    issue_summary_df <- if (nrow(results) > 0L) summarize_issue_df(results) else tibble::tibble()
    stage1_summary_df <- if (nrow(results) > 0L) summarize_stage1_problem_df(results) else tibble::tibble()
    data.table::fwrite(results, file = aggregate_replications_path)
  } else {
    results <- tibble::tibble()
    summary_df <- load_completed_condition_csv(design, out_dir, "summary")
    issue_summary_df <- load_completed_condition_csv(design, out_dir, "issue_summary")
    stage1_summary_df <- load_completed_condition_csv(design, out_dir, "stage1_summary")
    replication_index <- make_condition_replication_index(design, out_dir)
    utils::write.csv(replication_index, replication_index_path, row.names = FALSE)
    if (file.exists(aggregate_replications_path)) {
      unlink(aggregate_replications_path)
    }
    message(
      "Skipped in-memory replication aggregation for approximately ",
      format(expected_result_rows, big.mark = ",", scientific = FALSE),
      " rows. Condition-level replication files remain authoritative."
    )
  }

  utils::write.csv(summary_df, file = aggregate_summary_path, row.names = FALSE)
  utils::write.csv(issue_summary_df, file = aggregate_issue_summary_path, row.names = FALSE)
  utils::write.csv(stage1_summary_df, file = aggregate_stage1_summary_path, row.names = FALSE)
  # Preserve the original aggregate filenames for downstream scripts that do
  # not know about chunk-specific output naming.
  if (identical(chunk_label, "full_selection") &&
      identical(tolower(as.character(study_arg[[1]])), "all")) {
    if (materialize_replications) {
      data.table::fwrite(results, file = file.path(out_dir, "vig_hallquist_replication_results.csv.gz"))
    } else {
      legacy_replication_path <- file.path(
        out_dir,
        "vig_hallquist_replication_results.csv.gz"
      )
      if (file.exists(legacy_replication_path)) {
        unlink(legacy_replication_path)
      }
      utils::write.csv(
        replication_index,
        file = file.path(out_dir, "vig_hallquist_condition_replication_index.csv"),
        row.names = FALSE
      )
    }
    utils::write.csv(summary_df, file = file.path(out_dir, "vig_hallquist_summary.csv"), row.names = FALSE)
    utils::write.csv(issue_summary_df, file = file.path(out_dir, "vig_hallquist_issue_summary.csv"), row.names = FALSE)
    utils::write.csv(stage1_summary_df, file = file.path(out_dir, "vig_hallquist_stage1_problem_summary.csv"), row.names = FALSE)
  }

  message("Saved outputs to: ", normalizePath(out_dir))
  if (materialize_replications) {
    message("Aggregate replication results: ", aggregate_replications_path)
  } else {
    message("Condition replication index: ", replication_index_path)
  }
  message("Aggregate summary results: ", aggregate_summary_path)
  message("Aggregate issue summary results: ", aggregate_issue_summary_path)
  message("Aggregate stage1 problem summary results: ", aggregate_stage1_summary_path)

  invisible(list(results = results, summary = summary_df, issue_summary = issue_summary_df, stage1_summary = stage1_summary_df))
}
