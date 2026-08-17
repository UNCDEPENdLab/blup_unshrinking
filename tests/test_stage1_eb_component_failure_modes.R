#!/usr/bin/env Rscript

# Diagnostic check for `get_stage1_eb_components()` failure modes.
#
# This script uses matched Lai Study 1 and Study 4 conditions from the shared
# design grid. Study 1 provides an iid-residual control; Study 4 provides the
# AR(1) non-diagonal residual structure that has been more likely to trigger
# Stage-1 failures in practice. The local classifier below separates the two
# matrix inversions inside the Stage-1 loop so we can tell whether the problem
# is happening in the `R_i` solve, the `Sigma_i` solve, or both.

suppressPackageStartupMessages({
  library(dplyr)
})

find_repo_root <- function() {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  script_path <- if (length(script_arg) > 0L) {
    normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = FALSE)
  } else {
    normalizePath(file.path("tests", "test_stage1_eb_component_failure_modes.R"), mustWork = FALSE)
  }

  candidates <- unique(normalizePath(c(
    getwd(),
    file.path(getwd(), ".."),
    file.path(dirname(script_path), "..")
  ), mustWork = FALSE))

  for (candidate in candidates) {
    if (file.exists(file.path(candidate, "R", "blup_helpers.R"))) {
      return(candidate)
    }
  }

  getwd()
}

repo_root <- find_repo_root()
source(file.path(repo_root, "R", "source_helpers.R"), local = TRUE)
source_project_helpers(repo_root)
source(file.path(repo_root, "lai_replication", "designs.R"), local = TRUE)
source(file.path(repo_root, "lai_replication", "study_common.R"), local = TRUE)
source(file.path(repo_root, "lai_replication", "study1.R"), local = TRUE)
source(file.path(repo_root, "lai_replication", "study4.R"), local = TRUE)

pick_condition <- function(design, spec) {
  match_idx <- rep(TRUE, nrow(design))
  for (nm in names(spec)) {
    if (!nm %in% names(design)) {
      next
    }
    value <- spec[[nm]]
    if (is.na(value)) {
      match_idx <- match_idx & is.na(design[[nm]])
    } else {
      match_idx <- match_idx & design[[nm]] == value
    }
  }
  match_rows <- which(match_idx)
  if (length(match_rows) == 0L) {
    stop("Could not find the requested Lai design condition for ", spec$study, ".")
  }
  design[match_rows[[1L]], , drop = FALSE]
}

diagnose_stage1_inversions <- function(fit_obj, data, cluster_var, outcome_var, within_var, R_list = NULL, group = NULL) {
  cluster_ids <- unique(as.character(data[[cluster_var]]))
  split_dat <- split(data, as.character(data[[cluster_var]]), drop = TRUE)[cluster_ids]

  stage1 <- extract_stage1_components(
    fit_obj = fit_obj,
    data = data,
    cluster_var = cluster_var,
    within_var = within_var,
    R_list = R_list,
    group = group
  )
  beta_hat <- stage1$beta_hat
  g_hat <- stage1$G_hat
  R_list <- normalize_R_list(stage1$R_list, cluster_ids)

  out <- lapply(cluster_ids, function(cluster_id) {
    df_i <- split_dat[[cluster_id]]
    if (is.null(within_var)) {
      z_mat <- matrix(1, nrow = nrow(df_i), ncol = 1L)
      x_mat <- z_mat
      beta_vec <- beta_hat[[1L]]
    } else {
      z_vec <- df_i[[within_var]]
      z_mat <- cbind(1, z_vec)
      x_mat <- z_mat
      beta_vec <- c(beta_hat[[1L]], beta_hat[[within_var]])
    }

    resid_i <- df_i[[outcome_var]] - as.numeric(x_mat %*% beta_vec)
    R_i <- as.matrix(R_list[[cluster_id]])

    r_inverse_error <- tryCatch({
      solve(R_i, z_mat)
      solve(R_i, resid_i)
      NULL
    }, error = function(e) conditionMessage(e))

    sigma_y_i <- z_mat %*% g_hat %*% t(z_mat) + R_i
    sigma_y_error <- tryCatch({
      solve(sigma_y_i)
      NULL
    }, error = function(e) conditionMessage(e))

    data.frame(
      id = cluster_id,
      r_inverse_failed = !is.null(r_inverse_error),
      sigma_y_failed = !is.null(sigma_y_error),
      failure_mode = if (!is.null(r_inverse_error) && !is.null(sigma_y_error)) {
        "both"
      } else if (!is.null(r_inverse_error)) {
        "R_inverse"
      } else if (!is.null(sigma_y_error)) {
        "sigma_y_inverse"
      } else {
        "none"
      },
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}

summarize_stage1_diagnostics <- function(label, condition, sim_fun, seeds) {
  rows <- lapply(seeds, function(seed) {
    set.seed(seed)
    sim <- sim_fun(condition)
    sim$lv1 <- add_lai_trial_index(sim$lv1, cluster_var = "cid")

    fit_obj <- fit_lai_stage1(
      y ~ x,
      random = ~x | cid,
      data = sim$lv1,
      condition = condition,
      cluster_var = "cid"
    )
    if (is.null(fit_obj)) {
      stop("Stage-1 fit failed for diagnostic condition: ", label)
    }

    stage1 <- get_stage1_eb_components(
      fit_obj = fit_obj,
      data = sim$lv1,
      cluster_var = "cid",
      outcome_var = "y",
      within_var = "x"
    )
    diag_rows <- diagnose_stage1_inversions(
      fit_obj = fit_obj,
      data = sim$lv1,
      cluster_var = "cid",
      outcome_var = "y",
      within_var = "x"
    )

    stage1_failed <- vapply(seq_len(nrow(stage1)), function(i) {
      all(is.na(stage1[i, setdiff(names(stage1), "id"), drop = FALSE]))
    }, logical(1))

    if (!identical(stage1_failed, diag_rows$r_inverse_failed | diag_rows$sigma_y_failed)) {
      stop("Stage-1 helper failure rows did not match the low-level inversion classifier for ", label)
    }

    transform(
      diag_rows,
      condition = label,
      r_structure = as.character(condition$r_structure[[1]]),
      r_rho = as.numeric(condition$r_rho[[1]]),
      r_rho_key = if (is.na(condition$r_rho[[1]])) -1 else as.numeric(condition$r_rho[[1]]),
      seed = seed,
      helper_failed = stage1_failed
    )
  })

  do.call(rbind, rows)
}

study1_design <- select_design("1")
study4_design <- select_design("4")

diagnostic_specs <- list(
  study1_iid_control = list(
    condition = pick_condition(study1_design, list(
      study = "study1",
      num_clus = 30L,
      clus_size = 3L,
      icc = 0.05,
      vr_u1_u0 = 2.0,
      cor_u0_u1 = 0.5,
      beta_zu1 = 0.4
    )),
    sim_fun = simulate_study1
  ),
  study4_ar1_small_clusters_low_icc = list(
    condition = pick_condition(study4_design, list(
      study = "study4",
      num_clus = 30L,
      clus_size = 3L,
      icc = 0.05,
      vr_u1_u0 = 2.0,
      cor_u0_u1 = 0.5,
      beta_zu1 = 0,
      r_rho = 0.6
    )),
    sim_fun = simulate_study4
    )
)

diagnostics <- do.call(rbind, lapply(names(diagnostic_specs), function(label) {
  spec <- diagnostic_specs[[label]]
  summarize_stage1_diagnostics(label, spec$condition, spec$sim_fun, c(20260527L))
}))

summary_rows <- do.call(rbind, lapply(split(diagnostics, diagnostics$condition), function(df_i) {
  data.frame(
    condition = df_i$condition[[1L]],
    r_structure = df_i$r_structure[[1L]],
    r_rho_key = df_i$r_rho_key[[1L]],
    r_inverse_failed = sum(df_i$r_inverse_failed),
    sigma_y_failed = sum(df_i$sigma_y_failed),
    helper_failed = sum(df_i$helper_failed),
    stringsAsFactors = FALSE
  )
}))

summary_rows$mode <- ifelse(
  summary_rows$r_inverse_failed > 0 & summary_rows$sigma_y_failed > 0,
  "both",
  ifelse(
    summary_rows$r_inverse_failed > 0,
    "R_inverse",
    ifelse(summary_rows$sigma_y_failed > 0, "sigma_y_inverse", "none")
  )
)

print(summary_rows, row.names = FALSE)

stopifnot(
  nrow(diagnostics) > 0L,
  all(diagnostics$failure_mode %in% c("none", "R_inverse", "sigma_y_inverse", "both")),
  all(diagnostics$helper_failed == (diagnostics$r_inverse_failed | diagnostics$sigma_y_failed)),
  any(summary_rows$condition == "study1_iid_control"),
  any(summary_rows$condition == "study4_ar1_small_clusters_low_icc"),
  all(summary_rows$r_structure[summary_rows$condition == "study1_iid_control"] == "iid"),
  all(summary_rows$r_structure[summary_rows$condition != "study1_iid_control"] == "ar1"),
  all(summary_rows$r_rho_key[summary_rows$condition == "study1_iid_control"] == -1)
)

cat("Stage-1 inversion failure-mode diagnostic ok\n")