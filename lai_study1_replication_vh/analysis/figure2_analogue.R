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
      cor_u0_u1_label = factor(cor_u0_u1, levels = sort(unique(cor_u0_u1))),
      beta_zu1_label = factor(beta_zu1, levels = sort(unique(beta_zu1)))
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
    ggplot2::facet_grid(cor_u0_u1_label ~ beta_zu1_label, labeller = ggplot2::label_both) +
    ggplot2::labs(
      title = "Lai Study 1 Figure 2 analogue",
      subtitle = "VH implementation; historical Lai scale, status-code-zero retention, and overlapping methods only",
      x = "ICC",
      y = "20%-trimmed bias",
      colour = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
}

plot_lai_study1_vh_primary_figure <- function(summary_df) {
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
    ggplot2::facet_grid(cor_u0_u1_label ~ beta_zu1_label, labeller = ggplot2::label_both) +
    ggplot2::labs(
      title = "Lai Study 1 VH primary companion",
      subtitle = "All seven VH methods; common latent-SD scale and VH primary eligibility",
      x = "ICC",
      y = "Mean bias",
      colour = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
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

#' Build Figure 2 analogue and VH-primary outputs without loading all cells at once.
run_lai_study1_vh_postestimation_figures <- function(results_dir, analysis_dir) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 package is required to create post-estimation figures.")
  }
  dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
  files <- lai_study1_vh_replication_files(results_dir)

  cell_summaries <- lapply(files, function(path) {
    results <- readr::read_csv(path, show_col_types = FALSE)
    list(
      historical = summarize_lai_study1_vh_figure2_analogue(results),
      vh_primary = summarize_lai_study1_vh_primary_figure(results)
    )
  })
  historical_summary <- dplyr::bind_rows(lapply(cell_summaries, `[[`, "historical"))
  vh_primary_summary <- dplyr::bind_rows(lapply(cell_summaries, `[[`, "vh_primary"))

  historical_csv <- file.path(analysis_dir, "figure2_analogue_cell_summary.csv")
  primary_csv <- file.path(analysis_dir, "vh_primary_companion_cell_summary.csv")
  readr::write_csv(historical_summary, historical_csv)
  readr::write_csv(vh_primary_summary, primary_csv)

  historical_plot <- plot_lai_study1_vh_figure2_analogue(historical_summary)
  primary_plot <- plot_lai_study1_vh_primary_figure(vh_primary_summary)
  plot_files <- c(
    figure2_analogue_png = file.path(analysis_dir, "figure2_analogue.png"),
    figure2_analogue_pdf = file.path(analysis_dir, "figure2_analogue.pdf"),
    vh_primary_png = file.path(analysis_dir, "vh_primary_companion.png"),
    vh_primary_pdf = file.path(analysis_dir, "vh_primary_companion.pdf")
  )
  ggplot2::ggsave(plot_files[["figure2_analogue_png"]], historical_plot, width = 8, height = 6, dpi = 300)
  ggplot2::ggsave(plot_files[["figure2_analogue_pdf"]], historical_plot, width = 8, height = 6)
  ggplot2::ggsave(plot_files[["vh_primary_png"]], primary_plot, width = 10, height = 6, dpi = 300)
  ggplot2::ggsave(plot_files[["vh_primary_pdf"]], primary_plot, width = 10, height = 6)

  invisible(list(
    figure2_analogue = historical_summary,
    vh_primary = vh_primary_summary,
    files = c(figure2_analogue_summary = historical_csv, vh_primary_summary = primary_csv, plot_files)
  ))
}
