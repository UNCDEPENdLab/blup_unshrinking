# Post-estimation Figure 2 analogue and VH-primary companion summaries.

lai_study1_vh_historical_figure2_methods <- function() {
  c("naive_dual_blup", "lai_2spa", "msem")
}

lai_study1_vh_historical_figure2_labels <- function() {
  c(
    naive_dual_blup = "Naive dual BLUP",
    lai_2spa = "2S-PA",
    msem = "MSEM"
  )
}

lai_study1_vh_primary_figure_labels <- function() {
  c(
    oracle_dual = "Oracle dual",
    naive_dual_blup = "Naive dual BLUP",
    closed_form_dual = "Closed-form dual",
    fuller_closed_form = "Fuller closed-form",
    fuller_alpha_stepdown_closed_form = "Fuller alpha-stepdown",
    lai_2spa = "2S-PA",
    msem = "MSEM"
  )
}

lai_study1_vh_safe_trimmed_mean <- function(x, trim = 0.1) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else mean(x, trim = trim)
}

lai_study1_vh_required_figure_columns <- function() {
  c(
    "condition_id", "method", "reporting_scale", "estimate", "truth",
    "status_code", "vh_analysis_eligible", "lai_original_eligible",
    "icc", "cor_u0_u1", "beta_zu1", "num_clus", "clus_size",
    "vr_u1_u0", "var_u1"
  )
}

lai_study1_vh_assert_figure_columns <- function(results) {
  missing <- setdiff(lai_study1_vh_required_figure_columns(), names(results))
  if (length(missing) > 0L) {
    stop(
      "Replication results are missing Figure 2 analysis columns: ",
      paste(missing, collapse = ", "),
      ". Rerun with the v4 Lai Study 1 VH pipeline."
    )
  }
  invisible(results)
}

#' Summarize the historical-style, restricted-method Figure 2 analogue.
#'
#' This intentionally mirrors the original Figure 2 aggregation: its outcome
#' is the 20%-trimmed Monte Carlo bias (10% trimmed from each tail), the
#' historical method-specific reporting scale, and Lai's `status_code == 0`
#' retention rule.  It is an analogue rather than an exact replication because
#' the VH bundle contains only three directly corresponding methods.
summarize_lai_study1_vh_figure2_analogue <- function(results) {
  lai_study1_vh_assert_figure_columns(results)
  methods <- lai_study1_vh_historical_figure2_methods()
  group_vars <- c(
    "condition_id", "method", "icc", "cor_u0_u1", "beta_zu1",
    "num_clus", "clus_size", "vr_u1_u0", "var_u1"
  )

  raw_population <- results |>
    dplyr::filter(reporting_scale == "raw", method %in% methods) |>
    dplyr::mutate(lai_original_eligible = dplyr::coalesce(lai_original_eligible, FALSE)) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarise(
      n_sim = dplyr::n(),
      n_lai_original_eligible = sum(lai_original_eligible),
      lai_original_eligibility_rate = mean(lai_original_eligible),
      historical_truth = dplyr::first(beta_zu1) * sqrt(dplyr::first(var_u1)),
      .groups = "drop"
    )

  historical_results <- results |>
    dplyr::filter(
      reporting_scale == "lai_original_standardized",
      method %in% methods
    ) |>
    dplyr::mutate(
      lai_original_eligible = dplyr::coalesce(lai_original_eligible, FALSE),
      retained = lai_original_eligible & is.finite(estimate)
    ) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarise(
      n_historical_scale_available = dplyr::n(),
      n_retained = sum(retained),
      historical_scale_availability_rate = mean(is.finite(estimate)),
      truth = dplyr::first(truth),
      mean_bias = mean(estimate[retained], na.rm = TRUE) - dplyr::first(truth),
      robust_bias = lai_study1_vh_safe_trimmed_mean(estimate[retained]) - dplyr::first(truth),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      dplyr::across(c(mean_bias, robust_bias), ~ dplyr::if_else(is.nan(.x), NA_real_, .x))
    )

  raw_population |>
    dplyr::left_join(historical_results, by = group_vars) |>
    dplyr::mutate(
      n_historical_scale_available = dplyr::coalesce(n_historical_scale_available, 0L),
      n_retained = dplyr::coalesce(n_retained, 0L),
      truth = dplyr::coalesce(truth, historical_truth),
      figure_type = "lai_figure2_analogue",
      reporting_scale = "lai_original_standardized",
      eligibility_rule = "lai_original_status_code_zero",
      bias_statistic = "20_percent_trimmed_mean_bias"
    )
}

#' Summarize the all-seven-method VH primary companion figure.
summarize_lai_study1_vh_primary_figure <- function(results) {
  lai_study1_vh_assert_figure_columns(results)
  group_vars <- c(
    "condition_id", "method", "icc", "cor_u0_u1", "beta_zu1",
    "num_clus", "clus_size", "vr_u1_u0", "var_u1"
  )

  results |>
    dplyr::filter(reporting_scale == "latent_sd") |>
    dplyr::mutate(
      vh_analysis_eligible = dplyr::coalesce(vh_analysis_eligible, FALSE),
      retained = vh_analysis_eligible & is.finite(estimate)
    ) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarise(
      truth = dplyr::first(truth),
      n_sim = dplyr::n(),
      n_vh_analysis_eligible = sum(vh_analysis_eligible),
      n_retained = sum(retained),
      vh_analysis_eligibility_rate = mean(vh_analysis_eligible),
      mean_bias = mean(estimate[retained], na.rm = TRUE) - dplyr::first(truth),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      mean_bias = dplyr::if_else(is.nan(mean_bias), NA_real_, mean_bias),
      figure_type = "vh_primary_companion",
      reporting_scale = "latent_sd",
      eligibility_rule = "vh_primary",
      bias_statistic = "mean_bias"
    )
}

lai_study1_vh_prepare_figure_data <- function(summary_df, labels) {
  known_methods <- intersect(names(labels), unique(summary_df$method))
  summary_df |>
    dplyr::filter(method %in% known_methods) |>
    dplyr::mutate(
      method_label = factor(method, levels = names(labels), labels = unname(labels)),
      icc_label = factor(icc, levels = sort(unique(icc))),
      cor_u0_u1_label = factor(
        paste0("Corr(u0, u1) = ", cor_u0_u1),
        levels = paste0("Corr(u0, u1) = ", sort(unique(cor_u0_u1)))
      ),
      beta_zu1_label = factor(
        paste0("beta[Z.u1] = ", beta_zu1),
        levels = paste0("beta[Z.u1] = ", sort(unique(beta_zu1)))
      ),
      num_clus_label = factor(
        paste0("J == ", num_clus),
        levels = paste0("J == ", sort(unique(num_clus)))
      )
    )
}

plot_lai_study1_vh_figure2_analogue <- function(summary_df) {
  plot_df <- lai_study1_vh_prepare_figure_data(
    summary_df,
    lai_study1_vh_historical_figure2_labels()
  )
  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = icc_label, y = robust_bias, colour = method_label)
  ) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey45") +
    ggplot2::geom_boxplot(na.rm = TRUE, outlier.alpha = 0.45) +
    ggplot2::facet_grid(cor_u0_u1_label ~ beta_zu1_label, labeller = ggplot2::label_parsed) +
    ggplot2::labs(
      title = "Lai Study 1 Figure 2 analogue",
      subtitle = "Historical Lai scale; status-code-zero retention; overlapping methods only",
      x = "ICC",
      y = "20%-trimmed bias",
      colour = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
}

plot_lai_study1_vh_primary_figure <- function(summary_df, central_limit = 0.2) {
  if (!is.finite(central_limit) || central_limit <= 0) {
    stop("`central_limit` must be positive and finite.")
  }
  plot_df <- lai_study1_vh_prepare_figure_data(
    summary_df,
    lai_study1_vh_primary_figure_labels()
  )
  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = icc_label, y = mean_bias, colour = method_label)
  ) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey45") +
    ggplot2::geom_boxplot(na.rm = TRUE, outlier.alpha = 0.35) +
    ggplot2::coord_cartesian(ylim = c(-central_limit, central_limit)) +
    ggplot2::facet_grid(cor_u0_u1_label ~ beta_zu1_label, labeller = ggplot2::label_parsed) +
    ggplot2::labs(
      title = "VH-primary companion: central bias range",
      subtitle = paste0(
        "All seven methods; common latent-SD scale; VH eligibility; displayed range ±",
        format(central_limit, trim = TRUE)
      ),
      x = "ICC",
      y = "Mean bias",
      colour = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
}

#' Show the full absolute-bias distribution omitted from the central-range plot.
#'
#' A log scale makes rare but consequential unstable conditions visible without
#' forcing the central comparison onto an uninformative multi-thousand-unit axis.
plot_lai_study1_vh_primary_extremes <- function(summary_df) {
  plot_df <- lai_study1_vh_prepare_figure_data(
    summary_df,
    lai_study1_vh_primary_figure_labels()
  ) |>
    dplyr::mutate(abs_mean_bias = pmax(abs(mean_bias), 1e-8))
  ggplot2::ggplot(plot_df, ggplot2::aes(x = method_label, y = abs_mean_bias)) +
    ggplot2::geom_hline(yintercept = c(0.1, 1), colour = "grey55", linetype = "dashed") +
    ggplot2::geom_boxplot(fill = "#d7eaf7", colour = "#1f4e79", outlier.alpha = 0.35) +
    ggplot2::scale_y_log10() +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "VH-primary companion: absolute-bias instability diagnostic",
      subtitle = "All conditions retained in the distribution; dashed references at |bias| = 0.1 and 1",
      x = NULL,
      y = "Absolute condition-level mean bias (log scale)"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")
}

lai_study1_vh_replication_files <- function(results_dir) {
  condition_dir <- file.path(results_dir, "conditions")
  if (!dir.exists(condition_dir)) {
    stop("Could not find condition replication outputs in: ", condition_dir)
  }
  files <- list.files(
    condition_dir,
    pattern = "^condition_[0-9]+_replications\\.csv\\.gz$",
    full.names = TRUE
  )
  if (length(files) == 0L) {
    stop("No condition replication CSV files were found in: ", condition_dir)
  }
  sort(files)
}

#' Cache format for the historical Figure 2 condition summaries.
#'
#' Increment when a change could alter the trimmed-bias calculation.  The
#' source-file signatures stored beside the cache prevent stale reuse after a
#' condition result is replaced.
lai_study1_vh_figure_cache_version <- function() {
  "lai_study1_vh_figure_cache_v1"
}

lai_study1_vh_historical_cache_dir <- function(analysis_dir) {
  file.path(analysis_dir, "historical_condition_cache")
}

lai_study1_vh_historical_manifest_path <- function(analysis_dir) {
  file.path(analysis_dir, "figure2_analogue_cache_manifest.csv")
}

lai_study1_vh_historical_summary_path <- function(analysis_dir) {
  file.path(analysis_dir, "figure2_analogue_cell_summary.csv")
}

lai_study1_vh_condition_id_from_result_path <- function(path) {
  match <- regexec("^condition_([0-9]+)_replications\\.csv\\.gz$", basename(path))
  pieces <- regmatches(basename(path), match)[[1]]
  if (length(pieces) != 2L) {
    stop("Could not extract a condition ID from: ", path)
  }
  as.integer(pieces[[2]])
}

lai_study1_vh_source_signature <- function(path) {
  info <- file.info(path)
  if (!isTRUE(file.exists(path)) || is.na(info$size) || is.na(info$mtime)) {
    stop("Could not obtain a source signature for: ", path)
  }
  tibble::tibble(
    condition_id = lai_study1_vh_condition_id_from_result_path(path),
    source_file = basename(path),
    source_size_bytes = as.numeric(info$size),
    source_mtime_unix = as.numeric(info$mtime)
  )
}

lai_study1_vh_source_signatures <- function(files) {
  dplyr::bind_rows(lapply(files, lai_study1_vh_source_signature)) |>
    dplyr::arrange(condition_id)
}

lai_study1_vh_historical_cache_columns <- function() {
  c(
    "condition_id", "method", "reporting_scale", "estimate", "truth",
    "status_code", "vh_analysis_eligible", "lai_original_eligible",
    "icc", "cor_u0_u1", "beta_zu1", "num_clus", "clus_size",
    "vr_u1_u0", "var_u1"
  )
}

lai_study1_vh_historical_methods_complete <- function(summary_df, condition_ids) {
  expected_methods <- lai_study1_vh_historical_figure2_methods()
  required <- c(
    "condition_id", "method", "robust_bias", "reporting_scale",
    "eligibility_rule", "bias_statistic"
  )
  if (!all(required %in% names(summary_df))) return(FALSE)
  keys <- summary_df |>
    dplyr::select(condition_id, method) |>
    dplyr::distinct()
  expected <- tidyr::crossing(
    condition_id = sort(unique(as.integer(condition_ids))),
    method = expected_methods
  )
  nrow(summary_df) == nrow(expected) &&
    nrow(keys) == nrow(expected) &&
    nrow(dplyr::anti_join(expected, keys, by = c("condition_id", "method"))) == 0L &&
    all(as.character(summary_df$reporting_scale) == "lai_original_standardized") &&
    all(as.character(summary_df$eligibility_rule) == "lai_original_status_code_zero") &&
    all(as.character(summary_df$bias_statistic) == "20_percent_trimmed_mean_bias")
}

lai_study1_vh_cache_metadata_matches <- function(cache_df, signature) {
  required <- c(
    "figure_cache_version", "source_file", "source_size_bytes",
    "source_mtime_unix"
  )
  if (!all(required %in% names(cache_df)) || nrow(cache_df) == 0L) return(FALSE)
  all(cache_df$figure_cache_version == lai_study1_vh_figure_cache_version()) &&
    all(cache_df$source_file == signature$source_file[[1]]) &&
    all(cache_df$source_size_bytes == signature$source_size_bytes[[1]]) &&
    isTRUE(all.equal(
      as.numeric(cache_df$source_mtime_unix),
      rep(signature$source_mtime_unix[[1]], nrow(cache_df)),
      tolerance = 0
    ))
}

lai_study1_vh_condition_cache_path <- function(analysis_dir, condition_id) {
  file.path(
    lai_study1_vh_historical_cache_dir(analysis_dir),
    sprintf("condition_%04d_historical_summary.csv", condition_id)
  )
}

lai_study1_vh_read_historical_condition_cache <- function(cache_path, signature) {
  if (!file.exists(cache_path)) return(NULL)
  cache_df <- tryCatch(
    readr::read_csv(cache_path, show_col_types = FALSE),
    error = function(e) NULL
  )
  if (is.null(cache_df) || !lai_study1_vh_cache_metadata_matches(cache_df, signature)) {
    return(NULL)
  }
  cache_df |>
    dplyr::select(-dplyr::any_of(c(
      "figure_cache_version", "source_file", "source_size_bytes", "source_mtime_unix"
    )))
}

lai_study1_vh_write_historical_condition_cache <- function(summary_df, signature, analysis_dir) {
  cache_dir <- lai_study1_vh_historical_cache_dir(analysis_dir)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cache_path <- lai_study1_vh_condition_cache_path(analysis_dir, signature$condition_id[[1]])
  cache_df <- summary_df |>
    dplyr::mutate(
      figure_cache_version = lai_study1_vh_figure_cache_version(),
      source_file = signature$source_file[[1]],
      source_size_bytes = signature$source_size_bytes[[1]],
      source_mtime_unix = signature$source_mtime_unix[[1]]
    )
  readr::write_csv(cache_df, cache_path)
  summary_df
}

lai_study1_vh_read_historical_condition <- function(path) {
  # gzip streams still have to be decompressed, but parsing only this small
  # column set materially reduces CPU and peak memory versus parsing every
  # diagnostic column in each replication result.
  results <- readr::read_csv(
    path,
    col_select = dplyr::all_of(lai_study1_vh_historical_cache_columns()),
    show_col_types = FALSE
  )
  summarize_lai_study1_vh_figure2_analogue(results)
}

lai_study1_vh_write_historical_manifest <- function(signatures, analysis_dir,
                                                     cache_origin = "computed") {
  manifest <- signatures |>
    dplyr::mutate(
      figure_cache_version = lai_study1_vh_figure_cache_version(),
      cache_origin = cache_origin,
      cached_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
    )
  readr::write_csv(manifest, lai_study1_vh_historical_manifest_path(analysis_dir))
  invisible(manifest)
}

lai_study1_vh_historical_manifest_matches <- function(signatures, analysis_dir) {
  manifest_path <- lai_study1_vh_historical_manifest_path(analysis_dir)
  if (!file.exists(manifest_path)) return(FALSE)
  manifest <- tryCatch(readr::read_csv(manifest_path, show_col_types = FALSE), error = function(e) NULL)
  required <- c(
    "condition_id", "source_file", "source_size_bytes", "source_mtime_unix",
    "figure_cache_version"
  )
  if (is.null(manifest) || !all(required %in% names(manifest))) return(FALSE)
  manifest <- manifest |>
    dplyr::select(dplyr::all_of(required)) |>
    dplyr::arrange(condition_id)
  expected <- signatures |>
    dplyr::mutate(figure_cache_version = lai_study1_vh_figure_cache_version()) |>
    dplyr::select(dplyr::all_of(required)) |>
    dplyr::arrange(condition_id)
  identical(manifest$condition_id, expected$condition_id) &&
    identical(manifest$source_file, expected$source_file) &&
    identical(manifest$source_size_bytes, expected$source_size_bytes) &&
    identical(manifest$figure_cache_version, expected$figure_cache_version) &&
    isTRUE(all.equal(manifest$source_mtime_unix, expected$source_mtime_unix, tolerance = 0))
}

lai_study1_vh_read_valid_historical_summary_cache <- function(signatures, analysis_dir) {
  summary_path <- lai_study1_vh_historical_summary_path(analysis_dir)
  if (!file.exists(summary_path) || !lai_study1_vh_historical_manifest_matches(signatures, analysis_dir)) {
    return(NULL)
  }
  summary_df <- tryCatch(readr::read_csv(summary_path, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(summary_df) || !lai_study1_vh_historical_methods_complete(summary_df, signatures$condition_id)) {
    return(NULL)
  }
  summary_df
}

lai_study1_vh_bootstrap_historical_cache <- function(signatures, analysis_dir) {
  # The pre-cache pipeline already writes this combined summary.  When its
  # structural keys are complete, promote it to the new cache rather than
  # forcing a second 486-file pass solely to establish a cache manifest.
  summary_path <- lai_study1_vh_historical_summary_path(analysis_dir)
  if (!file.exists(summary_path)) return(NULL)
  summary_df <- tryCatch(readr::read_csv(summary_path, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(summary_df) || !lai_study1_vh_historical_methods_complete(summary_df, signatures$condition_id)) {
    return(NULL)
  }
  for (i in seq_len(nrow(signatures))) {
    signature <- signatures[i, , drop = FALSE]
    cell_summary <- summary_df |>
      dplyr::filter(condition_id == signature$condition_id[[1]])
    lai_study1_vh_write_historical_condition_cache(cell_summary, signature, analysis_dir)
  }
  lai_study1_vh_write_historical_manifest(signatures, analysis_dir, cache_origin = "bootstrap_existing_summary")
  summary_df
}

lai_study1_vh_historical_summary_from_cache <- function(results_dir, analysis_dir,
                                                         force = FALSE) {
  files <- lai_study1_vh_replication_files(results_dir)
  signatures <- lai_study1_vh_source_signatures(files)
  if (!isTRUE(force)) {
    cached <- lai_study1_vh_read_valid_historical_summary_cache(signatures, analysis_dir)
    if (!is.null(cached)) return(cached)
    bootstrapped <- lai_study1_vh_bootstrap_historical_cache(signatures, analysis_dir)
    if (!is.null(bootstrapped)) return(bootstrapped)
  }

  cell_summaries <- vector("list", nrow(signatures))
  for (i in seq_len(nrow(signatures))) {
    signature <- signatures[i, , drop = FALSE]
    cache_path <- lai_study1_vh_condition_cache_path(analysis_dir, signature$condition_id[[1]])
    cell_summary <- if (!isTRUE(force)) {
      lai_study1_vh_read_historical_condition_cache(cache_path, signature)
    } else {
      NULL
    }
    if (is.null(cell_summary)) {
      path <- files[match(signature$condition_id[[1]], signatures$condition_id)]
      cell_summary <- lai_study1_vh_read_historical_condition(path)
      cell_summary <- lai_study1_vh_write_historical_condition_cache(cell_summary, signature, analysis_dir)
    }
    cell_summaries[[i]] <- cell_summary
  }
  historical_summary <- dplyr::bind_rows(cell_summaries)
  if (!lai_study1_vh_historical_methods_complete(historical_summary, signatures$condition_id)) {
    stop("Historical Figure 2 cache construction did not produce every expected condition-method row.")
  }
  readr::write_csv(historical_summary, lai_study1_vh_historical_summary_path(analysis_dir))
  lai_study1_vh_write_historical_manifest(signatures, analysis_dir)
  historical_summary
}

lai_study1_vh_historical_inference_cache_path <- function(analysis_dir, condition_id) {
  file.path(
    lai_study1_vh_historical_cache_dir(analysis_dir),
    sprintf("condition_%04d_historical_inference.csv", condition_id)
  )
}

lai_study1_vh_historical_inference_columns <- function() {
  c(
    "condition_id", "method", "reporting_scale", "ci_low", "ci_high", "truth",
    "lai_original_eligible", "icc", "cor_u0_u1", "beta_zu1", "num_clus",
    "clus_size", "vr_u1_u0", "var_u1"
  )
}

#' Summarize historical-style coverage and CI-based Type I error by condition.
#'
#' Lai's original code calculates coverage from confidence limits and Type I
#' error from p < .05.  This implementation has stored confidence intervals,
#' not estimator-specific p values, so the Type I analogue is the equivalent
#' two-sided CI-exclusion rule when those intervals are Wald 95% intervals.
summarize_lai_study1_vh_historical_inference <- function(results) {
  required <- lai_study1_vh_historical_inference_columns()
  missing <- setdiff(required, names(results))
  if (length(missing) > 0L) {
    stop("Replication results are missing historical inference columns: ", paste(missing, collapse = ", "))
  }
  group_vars <- c(
    "condition_id", "method", "icc", "cor_u0_u1", "beta_zu1",
    "num_clus", "clus_size", "vr_u1_u0", "var_u1"
  )
  results |>
    dplyr::filter(
      reporting_scale == "lai_original_standardized",
      method %in% lai_study1_vh_historical_figure2_methods()
    ) |>
    dplyr::mutate(
      lai_original_eligible = dplyr::coalesce(lai_original_eligible, FALSE),
      interval_available = is.finite(ci_low) & is.finite(ci_high) & is.finite(truth),
      retained = lai_original_eligible & interval_available,
      covered = ci_low < truth & ci_high > truth,
      rejected = ci_low >= truth | ci_high <= truth
    ) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarise(
      n_historical_scale_available = dplyr::n(),
      n_retained = sum(retained),
      coverage = mean(covered[retained], na.rm = TRUE),
      type1_error = if (dplyr::first(beta_zu1) == 0) mean(rejected[retained], na.rm = TRUE) else NA_real_,
      .groups = "drop"
    ) |>
    dplyr::mutate(
      dplyr::across(c(coverage, type1_error), ~ dplyr::if_else(is.nan(.x), NA_real_, .x)),
      figure_type = "lai_figure3_4_analogue",
      reporting_scale = "lai_original_standardized",
      eligibility_rule = "lai_original_status_code_zero",
      type1_definition = "two_sided_95_percent_ci_exclusion"
    )
}

lai_study1_vh_read_historical_inference_condition_cache <- function(cache_path, signature) {
  cache_df <- lai_study1_vh_read_historical_condition_cache(cache_path, signature)
  if (is.null(cache_df)) return(NULL)
  required <- c("coverage", "type1_error", "type1_definition")
  if (!all(required %in% names(cache_df))) return(NULL)
  cache_df
}

lai_study1_vh_write_historical_inference_condition_cache <- function(summary_df, signature, analysis_dir) {
  cache_dir <- lai_study1_vh_historical_cache_dir(analysis_dir)
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cache_path <- lai_study1_vh_historical_inference_cache_path(analysis_dir, signature$condition_id[[1]])
  cache_df <- summary_df |>
    dplyr::mutate(
      figure_cache_version = lai_study1_vh_figure_cache_version(),
      source_file = signature$source_file[[1]],
      source_size_bytes = signature$source_size_bytes[[1]],
      source_mtime_unix = signature$source_mtime_unix[[1]]
    )
  readr::write_csv(cache_df, cache_path)
  summary_df
}

lai_study1_vh_historical_inference_from_cache <- function(results_dir, analysis_dir,
                                                           force = FALSE) {
  files <- lai_study1_vh_replication_files(results_dir)
  signatures <- lai_study1_vh_source_signatures(files)
  cell_summaries <- vector("list", nrow(signatures))
  for (i in seq_len(nrow(signatures))) {
    signature <- signatures[i, , drop = FALSE]
    cache_path <- lai_study1_vh_historical_inference_cache_path(analysis_dir, signature$condition_id[[1]])
    cell_summary <- if (!isTRUE(force)) {
      lai_study1_vh_read_historical_inference_condition_cache(cache_path, signature)
    } else {
      NULL
    }
    if (is.null(cell_summary)) {
      path <- files[match(signature$condition_id[[1]], signatures$condition_id)]
      results <- readr::read_csv(
        path,
        col_select = dplyr::all_of(lai_study1_vh_historical_inference_columns()),
        show_col_types = FALSE
      )
      cell_summary <- summarize_lai_study1_vh_historical_inference(results)
      cell_summary <- lai_study1_vh_write_historical_inference_condition_cache(cell_summary, signature, analysis_dir)
    }
    cell_summaries[[i]] <- cell_summary
  }
  dplyr::bind_rows(cell_summaries)
}

plot_lai_study1_vh_coverage_analogue <- function(summary_df) {
  plot_df <- lai_study1_vh_prepare_figure_data(
    summary_df,
    lai_study1_vh_historical_figure2_labels()
  )
  ggplot2::ggplot(plot_df, ggplot2::aes(x = icc_label, y = coverage, colour = method_label)) +
    ggplot2::geom_hline(yintercept = 0.95, colour = "grey45") +
    ggplot2::geom_boxplot(na.rm = TRUE, outlier.alpha = 0.4) +
    ggplot2::facet_grid(cor_u0_u1_label ~ num_clus_label, labeller = ggplot2::label_parsed) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      title = "Lai Study 1 Figure 3 analogue: 95% coverage",
      subtitle = "Historical Lai scale; status-code-zero retention; overlapping methods only",
      x = "ICC", y = "Empirical 95% coverage", colour = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
}

plot_lai_study1_vh_type1_analogue <- function(summary_df) {
  plot_df <- lai_study1_vh_prepare_figure_data(
    dplyr::filter(summary_df, beta_zu1 == 0),
    lai_study1_vh_historical_figure2_labels()
  )
  ggplot2::ggplot(plot_df, ggplot2::aes(x = icc_label, y = type1_error, colour = method_label)) +
    ggplot2::geom_hline(yintercept = 0.05, colour = "grey45") +
    ggplot2::geom_boxplot(na.rm = TRUE, outlier.alpha = 0.4) +
    ggplot2::facet_grid(cor_u0_u1_label ~ num_clus_label, labeller = ggplot2::label_parsed) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      title = "Lai Study 1 Figure 4 analogue: Type I error",
      subtitle = "CI-exclusion test at beta[Z.u1] = 0; historical Lai scale and status-code-zero retention",
      x = "ICC", y = "Empirical Type I error", colour = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
}

#' Derive the all-method VH companion entirely from the aggregate summary.
#'
#' The old implementation recomputed these condition-level means from the raw
#' replications.  The aggregation already contains the same VH-primary mean
#' bias and eligibility counts, so rereading the 486 gzip files is unnecessary.
summarize_lai_study1_vh_primary_from_aggregate <- function(results_dir) {
  summary_path <- file.path(results_dir, "lai_study1_vh_summary.csv")
  required <- c(
    "condition_id", "method", "reporting_scale", "icc", "cor_u0_u1", "beta_zu1",
    "num_clus", "clus_size", "vr_u1_u0", "var_u1", "truth", "n_rep",
    "n_vh_analysis_eligible", "bias", "vh_analysis_eligibility_rate"
  )
  if (!file.exists(summary_path)) {
    stop("Could not find aggregate summary for VH companion: ", summary_path)
  }
  aggregate <- readr::read_csv(
    summary_path,
    col_select = dplyr::all_of(required),
    show_col_types = FALSE
  )
  aggregate |>
    dplyr::filter(
      reporting_scale == "latent_sd",
      method %in% names(lai_study1_vh_primary_figure_labels())
    ) |>
    dplyr::transmute(
      condition_id, method, icc, cor_u0_u1, beta_zu1, num_clus, clus_size,
      vr_u1_u0, var_u1, truth,
      n_sim = n_rep,
      n_vh_analysis_eligible,
      n_retained = n_vh_analysis_eligible,
      vh_analysis_eligibility_rate,
      mean_bias = bias,
      figure_type = "vh_primary_companion",
      reporting_scale = "latent_sd",
      eligibility_rule = "vh_primary",
      bias_statistic = "mean_bias"
    )
}

#' Build Figure 2 analogue and VH-primary outputs with a resumable cache.
#'
#' `force = FALSE` reuses historical condition summaries only when their
#' source-file signatures and cache version match.  The VH companion never
#' needs replication-level rereads after the production aggregate exists.
run_lai_study1_vh_postestimation_figures <- function(results_dir, analysis_dir,
                                                      force = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 package is required to create post-estimation figures.")
  }
  dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
  historical_summary <- lai_study1_vh_historical_summary_from_cache(
    results_dir = results_dir,
    analysis_dir = analysis_dir,
    force = force
  )
  historical_inference <- lai_study1_vh_historical_inference_from_cache(
    results_dir = results_dir,
    analysis_dir = analysis_dir,
    force = force
  )
  vh_primary_summary <- summarize_lai_study1_vh_primary_from_aggregate(results_dir)

  historical_csv <- lai_study1_vh_historical_summary_path(analysis_dir)
  historical_inference_csv <- file.path(analysis_dir, "figure3_4_analogue_cell_summary.csv")
  primary_csv <- file.path(analysis_dir, "vh_primary_companion_cell_summary.csv")
  readr::write_csv(historical_inference, historical_inference_csv)
  readr::write_csv(vh_primary_summary, primary_csv)

  historical_plot <- plot_lai_study1_vh_figure2_analogue(historical_summary)
  coverage_plot <- plot_lai_study1_vh_coverage_analogue(historical_inference)
  type1_plot <- plot_lai_study1_vh_type1_analogue(historical_inference)
  primary_plot <- plot_lai_study1_vh_primary_figure(vh_primary_summary)
  primary_extremes_plot <- plot_lai_study1_vh_primary_extremes(vh_primary_summary)
  plot_files <- c(
    figure2_analogue_png = file.path(analysis_dir, "figure2_analogue.png"),
    figure2_analogue_pdf = file.path(analysis_dir, "figure2_analogue.pdf"),
    figure3_coverage_analogue_png = file.path(analysis_dir, "figure3_coverage_analogue.png"),
    figure3_coverage_analogue_pdf = file.path(analysis_dir, "figure3_coverage_analogue.pdf"),
    figure4_type1_analogue_png = file.path(analysis_dir, "figure4_type1_analogue.png"),
    figure4_type1_analogue_pdf = file.path(analysis_dir, "figure4_type1_analogue.pdf"),
    vh_primary_png = file.path(analysis_dir, "vh_primary_companion.png"),
    vh_primary_pdf = file.path(analysis_dir, "vh_primary_companion.pdf"),
    vh_primary_extremes_png = file.path(analysis_dir, "vh_primary_extremes.png"),
    vh_primary_extremes_pdf = file.path(analysis_dir, "vh_primary_extremes.pdf")
  )
  ggplot2::ggsave(plot_files[["figure2_analogue_png"]], historical_plot, width = 9, height = 8, dpi = 300)
  ggplot2::ggsave(plot_files[["figure2_analogue_pdf"]], historical_plot, width = 9, height = 8)
  ggplot2::ggsave(plot_files[["figure3_coverage_analogue_png"]], coverage_plot, width = 10, height = 9, dpi = 300)
  ggplot2::ggsave(plot_files[["figure3_coverage_analogue_pdf"]], coverage_plot, width = 10, height = 9)
  ggplot2::ggsave(plot_files[["figure4_type1_analogue_png"]], type1_plot, width = 10, height = 9, dpi = 300)
  ggplot2::ggsave(plot_files[["figure4_type1_analogue_pdf"]], type1_plot, width = 10, height = 9)
  ggplot2::ggsave(plot_files[["vh_primary_png"]], primary_plot, width = 11, height = 8, dpi = 300)
  ggplot2::ggsave(plot_files[["vh_primary_pdf"]], primary_plot, width = 11, height = 8)
  ggplot2::ggsave(plot_files[["vh_primary_extremes_png"]], primary_extremes_plot, width = 9, height = 4.5, dpi = 300)
  ggplot2::ggsave(plot_files[["vh_primary_extremes_pdf"]], primary_extremes_plot, width = 9, height = 4.5)

  invisible(list(
    figure2_analogue = historical_summary,
    figure3_4_analogue = historical_inference,
    vh_primary = vh_primary_summary,
    files = c(
      figure2_analogue_summary = historical_csv,
      figure3_4_analogue_summary = historical_inference_csv,
      figure2_analogue_manifest = lai_study1_vh_historical_manifest_path(analysis_dir),
      vh_primary_summary = primary_csv,
      plot_files
    )
  ))
}
