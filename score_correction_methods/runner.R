make_covu <- function(condition) {
  matrix(
    c(
      condition$icc,
      condition$cor_u0_u1 * sqrt(condition$icc * condition$var_u1),
      condition$cor_u0_u1 * sqrt(condition$icc * condition$var_u1),
      condition$var_u1
    ),
    nrow = 2L,
    byrow = TRUE
  )
}

draw_score_correction_level1_residuals <- function(cluster_sizes, sigma, condition) {
  r_spec <- score_correction_condition_to_r_spec(condition)
  unlist(lapply(cluster_sizes, function(n_i) {
    draw_residuals_from_R(make_R_matrix(n_i, sigma = sigma, r_spec = r_spec))
  }), use.names = FALSE)
}

score_correction_condition_to_nlme_correlation <- function(condition, cluster_var, index_var = "trial_index") {
  r_spec <- normalize_r_spec(score_correction_condition_to_r_spec(condition))
  switch(
    r_spec$structure,
    iid = NULL,
    ar1 = {
      if (!requireNamespace("nlme", quietly = TRUE)) {
        stop("The `nlme` package is required for AR(1) Stage-1 residual covariance fits.")
      }
      nlme::corAR1(form = stats::as.formula(sprintf("~%s | %s", index_var, cluster_var)))
    },
    stop("No nlme residual-correlation adapter is implemented for `", r_spec$structure, "`.")
  )
}

score_correction_condition_to_r_spec <- function(condition) {
    r_structure <- if ("r_structure" %in% names(condition)) {
    as.character(condition$r_structure[[1]])
  } else if ("residual_structure" %in% names(condition)) {
    as.character(condition$residual_structure[[1]])
  } else {
    "iid"
  }

  switch(
    tolower(r_structure),
    iid = list(structure = "iid"),
    ar1 = {
      rho <- if ("r_rho" %in% names(condition)) condition$r_rho[[1]] else condition$residual_rho[[1]]
      list(structure = "ar1", rho = as.numeric(rho))
    },
    toeplitz = {
      correlations <- if ("r_correlations" %in% names(condition)) condition$r_correlations[[1]] else condition$residual_correlations[[1]]
      list(structure = "toeplitz", correlations = as.numeric(correlations))
    },
    stop("Unsupported residual structure for score correction replication: ", r_structure)
  )
}

score_correction_condition_uses_non_iid_R <- function(condition) {
  !identical(normalize_r_spec(score_correction_condition_to_r_spec(condition))$structure, "iid")
}

fit_score_correction_stage1 <- function(fixed, random, data, condition, cluster_var) {
  if (score_correction_condition_uses_non_iid_R(condition)) {
    if (!requireNamespace("nlme", quietly = TRUE)) {
      return(NULL)
    }
    return(safe_lme(
      fixed = fixed,
      random = random,
      data = data,
      correlation = score_correction_condition_to_nlme_correlation(condition, cluster_var = cluster_var),
      method = "REML",
      control = nlme::lmeControl(returnObject = TRUE, msMaxIter = 100L, opt = "optim")
    ))
  }

  fixed_txt <- paste(deparse(fixed), collapse = " ")
  random_txt <- sub("^~", "", paste(deparse(random), collapse = " "))
  safe_lmer(stats::as.formula(paste0(fixed_txt, " + (", random_txt, ")")), data = data)
}

add_score_correction_trial_index <- function(data, cluster_var, index_var = "trial_index") {
  cluster_ids <- as.character(data[[cluster_var]])
  data[[index_var]] <- ave(seq_len(nrow(data)), cluster_ids, FUN = seq_along)
  data
}

parse_optional_integer_arg <- function(x) {
  # Accept the common string encodings used by shell wrappers for omitted args.
  if (is.null(x) || length(x) == 0L || is.na(x) || x %in% c("", "NA", "NaN", "NULL", "null")) {
    return(NA_integer_)
  }
  as.integer(x)
}

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
    stop("Requested chunk starts after the end of the selected design.")
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

read_replication_results_file <- function(path) {
  out <- tibble::as_tibble(data.table::fread(path))

  # Keep these casts centralized so aggregate rebuilds do not depend on the
  # exact type inference chosen by fread for any single condition file.
  numeric_cols <- intersect(
    c(
      "score1", "score2", "icc", "vr_u1_u0", "cor_u0_u1", "beta_zu1",
      "sigma2", "var_u1"
    ),
    names(out)
  )
  integer_cols <- intersect(
    c("condition_id", "rep", "num_clus", "clus_size"),
    names(out)
  )
  logical_cols <- intersect(
    c("success"),
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

summarize_issue_df <- function(results) {
  results %>%
    dplyr::group_by(method) %>%
    dplyr::mutate(n_total = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!success) %>%
    dplyr::group_by(method, err_step) %>%
    summarize(
      n = dplyr::n(),
      n_total = dplyr::first(n_total),
      prop = n / dplyr::first(n_total),
      .groups = "drop"
    )
}

summarize_results_df <- function(results) {
  results %>%
    dplyr::group_by(
      method, condition_id, num_clus, clus_size, icc, vr_u1_u0, cor_u0_u1, beta_zu1,
      design_source, condition_note
    ) %>%
    dplyr::summarise(
      prop_success = safe_mean(success),
      mean_score1 = safe_mean(score1),
      mean_score2 = safe_mean(score2),
      .groups = "drop"
    )
}

condition_output_has_methods <- function(path, expected_methods) {
  if (!file.exists(path)) {
    return(FALSE)
  }
  out <- tryCatch(read_replication_results_file(path), error = function(e) NULL)
  if (is.null(out) || !("method" %in% names(out))) {
    return(FALSE)
  }
  all(expected_methods %in% out$method)
}

run_score_correction_rep <- function(condition, sim) {
  sim$lv1 <- add_score_correction_trial_index(sim$lv1, cluster_var = "cid")
  fit_y <- fit_score_correction_stage1(y ~ x, random = ~x | cid, data = sim$lv1, condition = condition, cluster_var = "cid")
  if (is.null(fit_y)) {
    return(tibble::tibble(
      method = score_correction_methods(),
      success = FALSE,
      err_step = "stage 1 fit",
      err_msg = "Stage-1 model fit failed or nlme is not available.",
      score1 = NA_real_,
      score2 = NA_real_
    ))
  }
  score_rows <- run_scoring_methods(fit_y, sim$lv1, cluster_var = "cid", outcome_var = "y", within_var = "x")
  return(score_rows)
}

score_correction_methods <- function() {
  c("closed_form", "unweight_full", "unweight_diag", "g_r_conversion", "g_sigma_conversion")
}


score_correction_columns <- function() {
  c("method", "success", "err_step", "err_msg", "score1", "score2")
}


score_correction_methods_exports <- function() {
  c(
    "simulate_score_correction", "fit_score_correction_stage1", "add_score_correction_trial_index",
    "fixed_params", "make_score_correction_study_design",
    "score_correction_methods", "score_correction_columns",
    "run_condition_replications", "run_one_score_correction_rep", "run_score_correction_rep",
    "summarize_results_df", "summarize_issue_df", "condition_output_has_methods",
    "slice_design_chunk", "make_chunk_label", "write_progress_row", "read_replication_results_file", 
    "load_completed_condition_results", "get_condition_file_paths", "simulate_study1",
    "draw_score_correction_level1_residuals", "score_correction_condition_to_nlme_correlation", "score_correction_condition_to_r_spec",
    "fit_score_correction_stage1", "score_correction_condition_uses_non_iid_R", "safe_mean",
    "make_covu", "safe_lmer", "safe_lme", "normalize_r_spec", "make_R_matrix", "draw_residuals_from_R",
    "extract_stage1_components", "extract_stage1_components.merMod", "extract_stage1_components.lme",
    "extract_stage1_components.default", "normalize_R_list", "as_plain_vcov_matrix",
    "try_scoring_methods", "try_scoring_closed_form", "try_scoring_unweight_full", "try_scoring_unweight_diag",
    "try_scoring_g_r_conversion", "try_scoring_g_sigma_conversion", "run_scoring_methods", "sanitize_re_names"
  )
}

run_one_score_correction_rep <- function(condition) {
  sim <- simulate_score_correction(condition)
  run_score_correction_rep(condition, sim)
}


run_condition_replications <- function(condition, n_sim, n_cores = 1L) {
  rep_ids <- seq_len(n_sim)

  run_single_rep <- function(rep_id) {
    # The large condition multiplier avoids seed collisions across conditions
    # even for high replication counts.
    set.seed(20260419 + (as.integer(condition$condition_id) * 100000L) + as.integer(rep_id))
    rep_out <- run_one_score_correction_rep(condition)
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
      .packages = c("data.table", "lme4", "MASS", "dplyr", "tidyr", "purrr", "tibble", "OpenMx", "glmnet", "sandwich", "geigen"),
      .export = score_correction_methods_exports()
    ) %dopar% {
      run_single_rep(rep_id)
    }
  } else {
    purrr::map_dfr(rep_ids, run_single_rep)
  }
}


run_score_correction_simulation <- function(n_sim = 100L,
                                            out_dir = file.path(score_correction_dir, "outputs", "score_correction_comparison"),
                                            n_cores = 1L,
                                            max_conditions = NA_integer_,
                                            chunk_index = NA_integer_,
                                            chunk_size = NA_integer_,
                                            resume_existing = TRUE) {
  if (is.na(n_sim) || n_sim < 1L) {
    stop("`n_sim` must be a positive integer.")
  }
  if (is.na(n_cores) || n_cores < 1L) {
    stop("`n_cores` must be a positive integer.")
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  design <- make_score_correction_study_design(max_conditions = max_conditions)
  design <- slice_design_chunk(design, chunk_index = chunk_index, chunk_size = chunk_size)
  chunk_meta <- attr(design, "chunk_meta")
  chunk_label <- make_chunk_label(chunk_meta)
  file_prefix <- sprintf("score_correction_comparison_%s", chunk_label)

  # Chunk labels keep concurrently run jobs from overwriting each other's
  # aggregate files while preserving condition-level paths shared by resume.
  manifest_path <- file.path(out_dir, sprintf("%s_manifest.csv", file_prefix))
  progress_path <- file.path(out_dir, sprintf("%s_progress.csv", file_prefix))
  aggregate_replication_path <- file.path(out_dir, sprintf("%s_replications.csv.gz", file_prefix))
  aggregate_summary_path <- file.path(out_dir, sprintf("%s_summary.csv", file_prefix))
  aggregate_issue_summary_path <- file.path(out_dir, sprintf("%s_issue_summary.csv", file_prefix))

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
      expected_methods <- score_correction_methods()
      if (condition_output_has_methods(paths$replications, expected_methods)) {
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

    write_progress_row(
      progress_path,
      tibble::tibble(
        condition_id = condition_id,
        chunk_label = chunk_label,
        status = "started",
        started_at = started_at,
        finished_at = NA_character_,
        n_methods = NA_integer_,
        n_success = NA_integer_,
        replications_file = paths$replications,
        summary_file = paths$summary
      )
    )

    timing <- system.time({
      condition_results <- run_condition_replications(condition, n_sim = n_sim, n_cores = n_cores)
    })
    condition_summary <- summarize_results_df(condition_results)
    condition_issue_summary <- summarize_issue_df(condition_results)

    data.table::fwrite(condition_results, file = paths$replications)
    utils::write.csv(condition_summary, file = paths$summary, row.names = FALSE)
    utils::write.csv(condition_issue_summary, file = paths$issue_summary, row.names = FALSE)

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
        n_success = sum(condition_results$success, na.rm = TRUE),
        replications_file = paths$replications,
        summary_file = paths$summary
      )
    )
  }

  # Rebuild aggregates from disk so skipped and newly completed conditions are
  # treated identically, and so a resumed run repairs missing aggregate files.
  results <- load_completed_condition_results(design, out_dir)
  summary_df <- if (nrow(results) > 0L) summarize_results_df(results) else tibble::tibble()
  issue_summary_df <- if (nrow(results) > 0L) summarize_issue_df(results) else tibble::tibble()

  data.table::fwrite(results, file = aggregate_replication_path)
  utils::write.csv(summary_df, file = aggregate_summary_path, row.names = FALSE)
  utils::write.csv(issue_summary_df, file = aggregate_issue_summary_path, row.names = FALSE)
  # Preserve the original aggregate filenames for downstream scripts that do
  # not know about chunk-specific output naming.
  if (identical(chunk_label, "full_selection")) {
    data.table::fwrite(results, file = file.path(out_dir, "score_correction_comparison_results.csv.gz"))
    utils::write.csv(summary_df, file = file.path(out_dir, "score_correction_comparison_summary.csv"), row.names = FALSE)
    utils::write.csv(issue_summary_df, file = file.path(out_dir, "score_correction_comparison_issue_summary.csv"), row.names = FALSE)
  }

  message("Saved outputs to: ", normalizePath(out_dir))
  message("Aggregate replication results: ", aggregate_replication_path)
  message("Aggregate summary results: ", aggregate_summary_path)
  message("Aggregate issue summary results: ", aggregate_issue_summary_path)

  invisible(list(results = results, summary = summary_df, issue_summary = issue_summary_df))

}