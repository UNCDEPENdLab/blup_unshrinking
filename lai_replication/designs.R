#' Study-specific condition grids for Lai apples-to-apples replications.

fixed_params <- list(
  gamma0 = 0,
  gamma1 = 0.5,
  beta_zu0 = 0.4,
  z_intercept = 1.5
)

make_study1_design <- function() {
  tidyr::crossing(
    study = "study1",
    num_clus = c(30L, 100L, 500L),
    clus_size = c(3L, 10L, 25L),
    icc = c(0.05, 0.20, 0.50),
    vr_u1_u0 = c(0.5, 1.0, 2.0),
    cor_u0_u1 = c(-0.5, 0.0, 0.5),
    beta_zu1 = c(0.0, 0.4)
  ) %>%
    dplyr::mutate(
      sigma2 = 1 - icc,
      var_u1 = vr_u1_u0 * icc,
      design_source = "supplement_script",
      condition_note = "Supplement sim1.R uses clus_size = 3, 10, 25; this differs from the paper text."
    )
}

make_study2_design <- function() {
  tidyr::crossing(
    study = "study2",
    num_clus = c(30L, 150L),
    clus_size = 100L,
    icc = c(0.05, 0.50),
    vr_u1_u0 = c(0.5, 2.0),
    cor_u0_u1 = c(-0.5, 0.0, 0.5),
    beta_zu1 = c(0.0, -0.3)
  ) %>%
    dplyr::mutate(
      sigma2 = 1 - icc,
      var_u1 = vr_u1_u0 * icc,
      design_source = "supplement_script",
      condition_note = "Study 2 follows the supplement scripts and notebooks."
    )
}

make_study3_design <- function() {
  tidyr::crossing(
    study = "study3",
    num_clus = c(30L, 150L),
    clus_size = 20L,
    sigma_z = 2.0,
    icc = c(0.05, 0.50),
    vr_u1_u0 = c(0.5, 2.0),
    cor_u0_u1 = c(-0.5, 0.0, 0.5),
    beta_zu1 = c(0.0, -0.3)
  ) %>%
    dplyr::mutate(
      sigma2 = 1 - icc,
      var_u1 = vr_u1_u0 * icc,
      design_source = "supplement_script",
      condition_note = paste(
        "Study 3 uses simulation_scripts/sim3_revise.R as the operational definition.",
        "The paper, notebook text, and script disagree on the exact cluster-size profile."
      )
    )
}

all_designs <- dplyr::bind_rows(
  make_study1_design(),
  make_study2_design(),
  make_study3_design()
)

select_design <- function(study_arg = "all", max_conditions = NA_integer_) {
  study_arg <- tolower(study_arg)
  design <- if (identical(study_arg, "all")) {
    all_designs
  } else {
    all_designs %>%
      dplyr::filter(study %in% paste0("study", unlist(strsplit(study_arg, ","))))
  }

  if (!is.na(max_conditions) && max_conditions > 0L) {
    design <- design %>% dplyr::slice_head(n = max_conditions)
  }

  if (nrow(design) == 0L) {
    stop("No study conditions selected. Use `all`, `1`, `2`, `3`, or a comma-separated combination like `1,2`.")
  }

  design %>%
    dplyr::mutate(condition_id = seq_len(dplyr::n()))
}
