#' Study-specific condition grids for Vig-Hallquist (2026)

# fixed_params <- list(
#   gamma0 = 0,
#   gamma1 = 0.5,
#   beta_zu0 = 0.4,
#   z_intercept = 1.5
# )

make_study1_design <- function() {
  tidyr::crossing(
    study = "study1",
    num_clus = c(30L, 50L, 100L, 150L, 300L),
    clus_size = c(3L, 5L, 10L, 25L),
    exp_post_rel = c(0.2, 0.5, 0.8)
    cor_u0_u1 = c(-0.5, 0.0, 0.5),
    effect_size = c(0.0, 0.2, 0.4, 0.6),
    r_structure = "iid"
    ) %>%
    dplyr::mutate(
    # sigma2 =
    # var_u1 =
    study_label = "BLUP as Outcome"
    )
}

make_study2_design <- function() {
  tidyr::crossing(
    study = "study1",
    num_clus = c(30L, 50L, 100L, 150L, 300L),
    clus_size = c(3L, 5L, 10L, 25L),
    exp_post_rel = c(0.2, 0.5, 0.8)
    cor_u0_u1 = c(-0.5, 0.0, 0.5),
    effect_size = c(0.0, 0.2, 0.4, 0.6),
    struc_tar = c("slope_only", "intercept_slope"),
    r_structure = "iid"
    ) %>%
    dplyr::mutate(
    # sigma2 =
    # var_u1 =
    study_label = "BLUP as Outcome"
    )
}

make_study3_design <- function() {
  # ...
}

make_study4_design <- function() {
  # ...
}

all_designs <- dplyr::bind_rows(
  make_study1_design(),
  make_study2_design(),
  make_study3_design(),
  make_study4_design()
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
    stop("No study conditions selected. Use `all`, `1`, `2`, `3`, `4`, or a comma-separated combination like `1,2`.")
  }

  design %>%
    dplyr::mutate(condition_id = seq_len(dplyr::n()))
}
