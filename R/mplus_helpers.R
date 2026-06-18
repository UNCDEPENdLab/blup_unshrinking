extract_mplus_stats <- function(output_file, yvar, xvar, ci_multiplier = stats::qnorm(0.975),
                                reporting_scale = NULL, type = "reg") {
  if (type != "reg") {
    stop("`extract_mplus_stats` is only set up for regressions currently")
  }
  
  results <- tryCatch(
    MplusAutomation::readModels(output_file),
    error = function(e) return(NULL)
  )
  if (is.null(results)) {
    return(tibble::tibble(
      estimate = NA_real_,
      se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      status_code = 2L,
      mx_issue_class = "mplus_null_output_file",
      mx_issue_detail = NA_character_
    ))
  }
  
  errors <- results$errors
  if (length(errors)) {
    parsed_errors <- unlist(lapply(errors, paste, collapse = " "))
    return(tibble::tibble(
      estimate = NA_real_,
      se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      status_code = 1L,
      mx_issue_class = "mplus_error_thrown",
      mx_issue_detail = paste(parsed_errors, collapse = "\n")
    ))
  } else {
    
    pars <- results$parameters$unstandardized
    yvar <- toupper(yvar)
    xvar <- toupper(xvar)
    
    scale_u1 <- if (is.null(reporting_scale)) {
      stop("no default reporting scale for mplus models currently")
      # sqrt(with(pars, est[paramHeader == "Variances" & param == xvar]))
    } else {
      reporting_scale <- as.numeric(reporting_scale[[1]])
      if (!is.finite(reporting_scale) || reporting_scale <= 0) {
        stop("`reporting_scale` must be finite and positive.")
      }
      reporting_scale
    }
    
    operator <- switch(
      type,
      reg = ".ON",
      var = ".WITH"
    )
    rows <- which(pars$paramHeader == paste0(yvar, operator) & pars$param == xvar)
    est <- pars[rows, "est"] * scale_u1
    se <- pars[rows, "se"] * scale_u1
    
    return(tibble::tibble(
      estimate = est,
      se = se,
      ci_low = est - ci_multiplier * se,
      ci_high = est + ci_multiplier * se,
      status_code = 0L,
      mx_issue_class = "ok",
      mx_issue_detail = "ok"
    ))
  }
}

fit_mplus_blup_predictor <- function(level1_data, level2_data, join_by = dplyr::join_by("cid_chr" == "id"),
                                     title = "blup_predictor_msem", outcome_variable = "y", slope_only = FALSE,
                                     within_component = "x", between_component = "z", cluster_id = "cid",
                                     random_slope_label = "s", model_file = tempfile(title, fileext = ".inp"),
                                     output_file = gsub(".inp", ".out", model_file), reporting_scale = NULL) {
  
  if (!grepl(".inp$", model_file)) {
    stop("`model_file` must end in .inp")
  }
  if (model_file != gsub(".out", ".inp", model_file)) {
    stop("`model_file` and `output_file` must have the same name")
  }
  if (!(join_by$x %in% names(level1_data) && join_by$y %in% names(level2_data))) {
    stop("invalid join_by statement -- data missing columns")
  }
  
  data <- dplyr::left_join(level1_data, level2_data, by = join_by)
  
  model <- MplusAutomation::mplusObject(
    TITLE = title,
    ANALYSIS = "TYPE = TWOLEVEL RANDOM",
    VARIABLE = glue::glue(
      "WITHIN = {within_component};
       BETWEEN = {between_component};
       CLUSTER IS {cluster_id};"
    ),
    MODEL = glue::glue(
      "%WITHIN%
       {random_slope_label} | {outcome_variable} ON {within_component};
       %BETWEEN%
       {between_component} ON {if (slope_only) '' else outcome_variable} {random_slope_label};
       {outcome_variable} WITH {random_slope_label};"
    ),
    rdata = data
  )
  
  fit <- suppressWarnings({
    tryCatch(
      MplusAutomation::mplusModeler(model, run = 1L, modelout = model_file),
      error = function(e) return(NULL)
    )
  })
  
  if (is.null(fit)) {
    out <- tibble::tibble(
      estimate = NA_real_,
      se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      status_code = 2L,
      mx_issue_class = "mplus_null_fit",
      mx_issue_detail = NA_character_
    )
  } else {
    out <- extract_mplus_stats(
      output_file,
      xvar = random_slope_label,
      yvar = between_component,
      reporting_scale = reporting_scale
    )
  }
  
  return(out)
  
}

fit_mplus_dual_process <- function(proc1_data, proc2_data, title = "dual_process_sem", outcome1_var = "y", outcome2_var = "q",
                                   cluster_id = "cid", time_prefix = "t", time_index_var = "trial_index", round_time = 5,
                                   time_value_var = "x", slope_only = FALSE, model_file = tempfile(title, fileext = ".inp"),
                                   output_file = gsub(".inp", ".out", model_file), reporting_scale = NULL) {
  
  if (is.null(time_value_var)) time_value_var <- time_index_var
  
  prepare_process_data <- function(data, outcome, cluster, idx, time_val, time_prefix) {
    
    cols <- c(outcome, cluster, idx, time_val)
    
    if (!all(cols%in% colnames(data))) {
      stop("dual process data missing specified columns: ", 
           setdiff(cols, colnames(data)))
    }
    
    wide <- data %>%
      dplyr::select(!!cluster, !!idx, !!outcome) %>%
      dplyr::arrange(!!cluster, !!idx) %>%
      tidyr::pivot_wider(
        names_from = !!idx, 
        values_from = !!outcome, 
        names_prefix = paste0(time_prefix, outcome)
      )
    
    time_map <- data %>%
      dplyr::select(!!idx, !!time_val) %>%
      dplyr::distinct() %>%
      dplyr::arrange(!!idx)
    
    col_names <- paste0(time_prefix, outcome, time_map[[idx]])
    at_values  <- round(time_map[[time_val]], round_time)
    
    model_string <- paste0(col_names, "@", at_values)
    time_var_range <- paste(col_names[1], col_names[length(col_names)], sep = "-")
    
    return(list(wide = wide, model_string = model_string, time_var_range = time_var_range))
    
  }
  
  y_out <- prepare_process_data(
    data = proc1_data, 
    outcome= outcome1_var, 
    cluster = cluster_id, 
    idx = time_index_var,
    time_val = time_value_var,
    time_prefix = time_prefix
  )
  
  q_out <- prepare_process_data(
    data = proc2_data, 
    outcome= outcome2_var, 
    cluster = cluster_id, 
    idx = time_index_var,
    time_val = time_value_var,
    time_prefix = time_prefix
  )
  
  data <- dplyr::left_join(y_out$wide, q_out$wide, by = cluster_id)
  
  model <- MplusAutomation::mplusObject(
    TITLE = title,
    MODEL = glue::glue(
      "i1 s1 | {paste(y_out$model_string, collapse = '\n')};
       i2 s2 | {paste(q_out$model_string, collapse = '\n')};
       s2 ON {if (slope_only) 's1' else paste0(c('i1', 's1'), collapse = ' ')};
       {paste(y_out$time_var_range)} (1);
       {paste(q_out$time_var_range)} (2);
       i1 WITH s1;
       i2 WITH s2;"
    ),
    rdata = data
  )
  
  fit <- suppressWarnings({
    tryCatch(
      MplusAutomation::mplusModeler(model, run = 1L, modelout = model_file),
      error = function(e) return(NULL)
    )
  })
  
  if (is.null(fit)) {
    out <- tibble::tibble(
      estimate = NA_real_,
      se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      status_code = 2L,
      mx_issue_class = "mplus_null_fit",
      mx_issue_detail = NA_character_
    )
  } else {
    out <- extract_mplus_stats(
      output_file,
      xvar = "s1",
      yvar = "s2",
      reporting_scale = reporting_scale
    )
  }
  
  return(out)
  
}
