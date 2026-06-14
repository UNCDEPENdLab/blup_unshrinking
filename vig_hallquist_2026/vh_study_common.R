#' Shared study machinery that is not part of the general score helpers.

# vig_hallquist_truth <- function(condition) {
#   condition$gamma * sqrt(condition$tau1)
# }

simulate_data_blup_as_outcome <- function(condition) {
  sim <- simulate_dataset(
    n_id = condition$num_clus,
    mean_n_trial = condition$mean_clus_size,
    params = list(
      tau0 = fixed_params$tau0,
      rho = condition$rho,
      gamma_x_on_slope = condition$beta1w,
      beta_0 = fixed_params$gamma0_outcome,
      beta_z = fixed_params$gamma1_outcome
    ),
    tau1 = condition$tau1,
    sigma = condition$sigma,
    has_random_slope = TRUE,
    balanced = balance_mode_to_sim_arg(condition$balance_mode),
    min_n_trial = condition$min_clus_size,
    highly_unbalanced_min_n_trial = condition$highly_unbalanced_min_clus_size,
    highly_unbalanced_power = condition$highly_unbalanced_power,
    r_spec = condition_to_r_spec(condition)
  )
  
  list(
    lv1 = tibble::tibble(
      cid = sim$dat$id,
      cid_chr = as.character(sim$dat$id),
      trial_index = sim$dat$trial_index,
      x = sim$dat$z, # level 1 predictor is now "x" across the board
      y = sim$dat$y
    ),
    lv2_true = tibble::tibble(
      id = sim$id_df$id,
      w = sim$id_df$x, # level 2 predictor of slope is now "w" across the board
      z = NA_real_,
      true_u0 = sim$id_df$true_intercept_dev,
      true_u1 = sim$id_df$true_slope_dev
    ), 
    R_list = sim$R_list,
    r_spec = sim$r_spec,
    mean_realized_trials = sim$mean_realized_trials,
    min_realized_trials = sim$min_realized_trials,
    prop_ids_leq_2_trials = sim$prop_ids_leq_2_trials,
    prop_ids_leq_3_trials = sim$prop_ids_leq_3_trials,
    study_structure = "w",
    balance_mode = condition$balance_mode
  )
}

simulate_data_blup_as_predictor <- function(condition) {
  covu <- make_random_effect_covariance(
    intercept_variance = fixed_params$tau0,
    slope_variance = condition$tau1,
    intercept_slope_correlation = condition$marginal_rho
  )
  
  r_spec <- normalize_r_spec(r_spec = condition_to_r_spec(condition))
  
  if (!is.finite(condition$min_clus_size) || condition$min_clus_size < 1L) {
    stop("`min_clus_size` must be an integer >= 1.")
  }
  
  balanced <- balance_mode_to_sim_arg(condition$balance_mode)
  if (isFALSE(balanced) || (is.character(balanced) && balanced != "balanced")) {
    stop("`balance_mode != 'balanced' not yet implemented for study structure 'z'")
  }
  
  cluster_sizes <- rep(as.integer(condition$mean_clus_size), as.integer(condition$num_clus))
  cid <- rep(seq_len(as.integer(condition$num_clus)), cluster_sizes)
  
  xj <- seq(-1, 1, length.out = as.integer(condition$mean_clus_size))
  xj <- xj / sqrt(mean(xj^2))
  x <- rep(xj, as.integer(condition$num_clus))
  
  u <- draw_random_effects(n_id = condition$num_clus, tau0 = fixed_params$tau0, tau1 = condition$tau1, rho = fixed_params$marginal_rho)
  y <- fixed_params$gamma0_predictor + fixed_params$gamma1_predictor * x + rowSums(cbind(1, x) * u[cid, , drop = FALSE]) +
    draw_level1_residuals(cluster_sizes, sigma = condition$sigma, condition = condition)
  sigma_e <- fixed_params$z_variance - drop(t(c(fixed_params$beta1z, condition$beta2z)) %*% covu %*% c(fixed_params$beta1z, condition$beta2z))
  z <- fixed_params$beta0z + fixed_params$beta1z * u[, 1] + condition$beta2z * u[, 2] +
    stats::rnorm(as.integer(condition$num_clus), sd = sqrt(sigma_e))

  R_list <- stats::setNames(
    lapply(cluster_sizes, make_R_matrix, sigma = condition$sigma, r_spec = r_spec),
    unique(cid)
  )
  
  list(
    lv1 = tibble::tibble(
      cid = cid,
      cid_chr = as.character(cid),
      trial_index = ave(seq_along(cid), cid, FUN = seq_along),
      x = x,
      y = y
    ),
    lv2_true = tibble::tibble(
      id = as.character(seq_len(as.integer(condition$num_clus))),
      w = NA_real_,
      z = z,
      true_u0 = u[, 1],
      true_u1 = u[, 2]
    ), 
    R_list = R_list,
    r_spec = r_spec,
    mean_realized_trials = NA,
    min_realized_trials = NA,
    prop_ids_leq_2_trials = NA,
    prop_ids_leq_3_trials = NA,
    study_structure = "z",
    balance_mode = condition$balance_mode
  )
}

simulate_data_dual_blup <- function(condition) {
  stop("Dual BLUP simulation not yet implemented")
}

condition_to_r_spec <- function(condition) {
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
    stop("Unsupported residual structure: ", r_structure)
  )
}


condition_to_nlme_correlation <- function(condition) {
  r_structure <- if ("r_structure" %in% names(condition)) {
    as.character(condition$r_structure[[1]])
  } else {
    "iid"
  }
  
  switch(
    r_structure,
    iid = NULL,
    ar1 = {
      if (!requireNamespace("nlme", quietly = TRUE)) {
        stop("The `nlme` package is required for AR(1) Stage-1 residual covariance fits.")
      }
      nlme::corAR1(form = ~trial_index | id)
    },
    stop("Unsupported residual structure for nlme Stage-1 fit: ", r_structure)
  )
}

condition_uses_non_iid_R <- function(condition) {
  !identical(normalize_r_spec(condition_to_r_spec(condition))$structure, "iid")
}


add_trial_index <- function(data, cluster_var, index_var = "trial_index") {
  cluster_ids <- as.character(data[[cluster_var]])
  data[[index_var]] <- ave(seq_len(nrow(data)), cluster_ids, FUN = seq_along)
  data
}

draw_level1_residuals <- function(cluster_sizes, sigma, condition) {
  r_spec <- list(structure = "iid")
  unlist(lapply(cluster_sizes, function(n_i) {
    draw_residuals_from_R(make_R_matrix(n_i, sigma = sigma, r_spec = r_spec))
  }), use.names = FALSE)
}

# TODO: this is only for BLUPs as predictors right now (see README for models)
fit_stage1 <- function(fixed, random, data, condition, cluster_var) {
  if (condition_uses_non_iid_R(condition)) {
    if (!requireNamespace("nlme", quietly = TRUE)) {
      return(NULL)
    }
    return(safe_lme(
      fixed = fixed,
      random = random,
      data = data,
      correlation = condition_to_nlme_correlation(condition, cluster_var = cluster_var),
      method = "REML",
      control = nlme::lmeControl(returnObject = TRUE, msMaxIter = 100L, opt = "optim")
    ))
  }

  fixed_txt <- paste(deparse(fixed), collapse = " ")
  random_txt <- sub("^~", "", paste(deparse(random), collapse = " "))
  safe_lmer(stats::as.formula(paste0(fixed_txt, " + (", random_txt, ")")), data = data)
}

extract_centered_slope_eb <- function(fit_obj, ordered_ids) {
  if (is.null(fit_obj)) {
    return(tibble::tibble(id = ordered_ids, centered_u1_eb = NA_real_))
  }

  re_df <- tryCatch({
    if (inherits(fit_obj, "lme")) {
      nlme::ranef(fit_obj)
    } else {
      lme4::ranef(fit_obj)[[1]]
    }
  }, error = function(e) NULL)

  if (is.null(re_df)) {
    return(tibble::tibble(id = ordered_ids, centered_u1_eb = NA_real_))
  }

  tibble::tibble(id = rownames(re_df), centered_u1_eb = as.numeric(re_df[[ncol(re_df)]]))
}

# study method vectors in study specific files
study_methods_for_condition <- function(condition) {
  study_key <- as.character(condition$study[[1]])
  switch(
    study_key,
    study1 = study1_methods(),
    study2 = study2_methods(),
    study3 = study3_methods(),
    study4 = study4_methods(),
    stop("Unsupported study key: ", study_key)
  )
}

make_failed_result <- function(condition, methods, truth) {
  dplyr::bind_cols(
    tibble::tibble(
      study = condition$study,
      method = methods,
      estimate = NA_real_,
      se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      status_code = NA_integer_,
      truth = truth
    ),
    empty_stage1_diagnostics()[rep(1L, length(methods)), , drop = FALSE]
  )
}

add_study_result_context <- function(results, condition, truth) {
  dplyr::mutate(results, study = condition$study, truth = truth)
}

# TODO: is this matched across studies?
run_matched_outcome_rep <- function(condition, sim) {
  # truth <- vig_hallquist_truth(condition)
  sim$lv1 <- add_trial_index(sim$lv1, cluster_var = "cid")
  fit_y <- fit_stage1(y ~ x, random = ~x | cid, data = sim$lv1, condition = condition, cluster_var = "cid")
  if (is.null(fit_y)) {
    return(make_failed_result(condition, matched_study_methods(), truth))
  }

  ordered_ids <- sim$lv2_true$id
  # TODO: how are scoring failures handled in stage 2
  stage1_y <- get_stage1_eb_components(
    fit_obj = fit_y,
    data = sim$lv1,
    cluster_var = "cid",
    outcome_var = "y",
    within_var = "x"
  )
  corrected_y <- get_closed_form_corrected_scores(
    fit_obj = fit_y,
    data = sim$lv1,
    cluster_var = "cid",
    outcome_var = "y",
    within_var = "x"
  )

  yc <- sim$lv1$y - ave(sim$lv1$y, sim$lv1$cid)
  xc <- sim$lv1$x - ave(sim$lv1$x, sim$lv1$cid)
  centered_dat <- dplyr::mutate(sim$lv1, yc = yc, xc = xc)
  fit_centered <- fit_stage1(yc ~ 0 + xc, random = ~0 + xc | cid, data = centered_dat, condition = condition, cluster_var = "cid")
  centered_u1 <- extract_centered_slope_eb(fit_centered, ordered_ids = ordered_ids)

  stage2_df <- sim$lv2_true %>%
    dplyr::left_join(stage1_y, by = "id") %>%
    dplyr::left_join(corrected_y %>% dplyr::select(id, corrected_intercept_full, corrected_slope_full, ols_var11, ols_var12, ols_var22), by = "id") %>%
    dplyr::left_join(centered_u1, by = "id")
  stage1_diag <- get_stage1_diagnostics(fit_y, stage2_df)

  # TODO: update these methods
  results <- dplyr::bind_rows(
    fit_observed_dual(stage2_df, outcome = "z", predictor_u0 = "true_u0", predictor_u1 = "true_u1") %>%
      dplyr::filter(se_type == "naive") %>%
      dplyr::transmute(method = "oracle_dual", estimate, se, ci_low, ci_high, status_code),
    finalize_ols_se_variants(fit_observed_single(stage2_df, outcome = "z", predictor = "u1_eb"), "naive_slope_only"),
    finalize_ols_se_variants(fit_observed_single(stage2_df, outcome = "z", predictor = "centered_u1_eb"), "centered_slope_only"),
    finalize_ols_se_variants(fit_observed_dual(stage2_df, outcome = "z", predictor_u0 = "u0_eb", predictor_u1 = "u1_eb"), "naive_dual_eb"),
    finalize_ols_se_variants(fit_observed_single(stage2_df, outcome = "z", predictor = "corrected_slope_full"), "corrected_slope_only"),
    finalize_ols_se_variants(fit_observed_dual(stage2_df, outcome = "z", predictor_u0 = "corrected_intercept_full", predictor_u1 = "corrected_slope_full"), "corrected_dual"),
    fit_lai_2spa_observed_outcome(stage2_df, use_average = FALSE) %>%
      dplyr::mutate(method = "lai_2spa") %>%
      dplyr::select(method, dplyr::everything()),
    fit_lai_2spa_observed_outcome(stage2_df, use_average = TRUE) %>%
      dplyr::mutate(method = "lai_2spaa") %>%
      dplyr::select(method, dplyr::everything()),
    fit_fuller_dual(
      stage2_df,
      outcome = "z",
      predictor_u0 = "corrected_intercept_full",
      predictor_u1 = "corrected_slope_full",
      meas11 = "ols_var11",
      meas12 = "ols_var12",
      meas22 = "ols_var22"
    ) %>%
      dplyr::mutate(method = "fuller") %>%
      dplyr::select(method, dplyr::everything()),
    fit_fuller_dual_stepdown(
      stage2_df,
      outcome = "z",
      predictor_u0 = "corrected_intercept_full",
      predictor_u1 = "corrected_slope_full",
      meas11 = "ols_var11",
      meas12 = "ols_var12",
      meas22 = "ols_var22"
    ) %>%
      dplyr::mutate(method = "fuller_stepdown") %>%
      dplyr::select(method, dplyr::everything()),
    fit_fuller_dual_alpha_stepdown(
      stage2_df,
      outcome = "z",
      predictor_u0 = "corrected_intercept_full",
      predictor_u1 = "corrected_slope_full",
      meas11 = "ols_var11",
      meas12 = "ols_var12",
      meas22 = "ols_var22"
    ) %>%
      dplyr::mutate(method = "fuller_alpha_stepdown") %>%
      dplyr::select(method, dplyr::everything())
  )

  dplyr::bind_cols(results, stage1_diag[rep(1L, nrow(results)), , drop = FALSE]) %>%
    add_study_result_context(condition, truth)
}
