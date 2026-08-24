#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

source(file.path(
  "vig_hallquist_2026",
  "random_effects_structural_simulation.R"
))

# Legacy selectors and identifiers remain stable.
legacy1 <- select_design("1")
legacy5 <- select_design("5")
stopifnot(
  nrow(legacy1) == 720L,
  identical(range(legacy1$condition_id), c(1L, 720L)),
  all(legacy1$study == "study1"),
  all(legacy1$study_version == "legacy"),
  all(legacy1$calibration_version == "legacy_solve_tau1_sq"),
  identical(range(legacy5$condition_id), c(2417L, 2452L)),
  identical(
    nrow(make_study1_legacy_design()),
    nrow(make_study1_design())
  )
)

v2 <- select_design("allv2")
study1 <- dplyr::filter(v2, study == "study1v2")
study2 <- dplyr::filter(v2, study == "study2v2")
study3 <- dplyr::filter(v2, study == "study3v2")
study4 <- dplyr::filter(v2, study == "study4v2")

stopifnot(
  nrow(v2) == 2415L,
  nrow(study1) == 720L,
  nrow(study2) == 720L,
  nrow(study3) == 864L,
  nrow(study4) == 111L,
  identical(range(v2$condition_id), c(2453L, 4867L)),
  all(v2$study_version == "v2"),
  all(v2$calibration_version == "shape_preserving_v2"),
  all(v2$covariance_shape_fixed),
  all(v2$calibration_metric == "marginal_slope"),
  all(abs(v2$slope_intercept_variance_ratio - 1) < 1e-12),
  all(abs(v2$calibration_tau0_sq - 0.81) < 1e-12),
  all(abs(v2$calibration_tau1_sq - 0.81) < 1e-12),
  all(v2$G_condition_number <= 3 + 1e-10)
)

stopifnot(
  all(abs(study1$achieved_reliability - study1$target_reliability) < 1e-8),
  all(abs(
    study1$slope_variance_residual + study1$beta1w^2 -
      study1$slope_variance_marginal
  ) < 1e-12),
  all(abs(
    study1$rho_residual * fixed_params$tau0 * study1$tau1 -
      study1$marginal_rho * fixed_params$tau0 *
        sqrt(study1$slope_variance_marginal)
  ) < 1e-12),
  all(study1$structural_residual_G_min_eigen > 0),
  all(study2$outcome_residual_variance > 0),
  all(abs(study2$achieved_reliability - study2$target_reliability) < 1e-8),
  all(abs(study2$tau1^2 - 0.81) < 1e-12),
  all(abs(study3$achieved_reliability_y - study3$target_reliability_y) < 1e-8),
  all(abs(study3$achieved_reliability_q - study3$target_reliability_q) < 1e-8),
  all(study3$structural_residual_G_min_eigen > 0),
  all(study3$structural_joint_G_min_eigen > 0),
  all(abs(study3$tau1_y^2 - 0.81) < 1e-12),
  all(abs(study3$tau1_q^2 - 0.81) < 1e-12),
  all(abs(study4$achieved_reliability - study4$target_reliability) < 1e-8),
  all(study4$outcome_residual_variance > 0),
  all(abs(study4$tau1^2 - 0.81) < 1e-12)
)

rho_zero_1 <- dplyr::filter(study1, marginal_rho == 0)
rho_zero_2 <- dplyr::filter(study2, marginal_rho == 0)
stopifnot(
  all(abs(
    rho_zero_1$achieved_reliability -
      rho_zero_1$achieved_partial_reliability
  ) < 1e-10),
  all(abs(
    rho_zero_2$achieved_reliability -
      rho_zero_2$achieved_partial_reliability
  ) < 1e-10)
)

# The legacy one-coordinate solve is still callable and still visibly differs
# from the fixed-shape v2 DGM in the difficult low-R correlated cell.
legacy_extreme <- legacy1 %>%
  filter(
    mean_clus_size == 25L,
    target_reliability == 0.25,
    marginal_rho == 0.5
  ) %>%
  slice(1L)
stopifnot(
  legacy_extreme$slope_variance_marginal / fixed_params$tau0^2 < 0.01,
  legacy_extreme$calibration_version == "legacy_solve_tau1_sq"
)

# The ICC bridge is paired exactly at m=10, then deliberately holds different
# quantities constant away from that reference design.
bridge <- select_design("iccbridge")
bridge_first_stage <- bridge %>%
  distinct(
    calibration_arm, mean_clus_size, marginal_rho,
    posterior_reliability_anchor, icc_anchor, sigma,
    achieved_reliability, intercept_icc
  )
paired_m10 <- bridge_first_stage %>%
  filter(mean_clus_size == 10L) %>%
  group_by(marginal_rho, posterior_reliability_anchor) %>%
  summarise(
    sigma_range = diff(range(sigma)),
    reliability_range = diff(range(achieved_reliability)),
    icc_range = diff(range(intercept_icc)),
    .groups = "drop"
  )
posterior_arm <- bridge_first_stage %>%
  filter(calibration_arm == "posterior_reliability_targeted") %>%
  group_by(marginal_rho, posterior_reliability_anchor) %>%
  summarise(reliability_range = diff(range(achieved_reliability)), .groups = "drop")
icc_arm <- bridge_first_stage %>%
  filter(calibration_arm == "icc_anchored_at_m10") %>%
  group_by(marginal_rho, posterior_reliability_anchor) %>%
  summarise(icc_range = diff(range(intercept_icc)), .groups = "drop")
bridge_pairs <- bridge %>%
  group_by(bridge_pair_id, bridge_pair_label, simulation_seed_group) %>%
  summarise(
    n_arms = n_distinct(calibration_arm),
    n_rows = n(),
    .groups = "drop"
  )

stopifnot(
  nrow(bridge) == 108L,
  identical(range(bridge$condition_id), c(4868L, 4975L)),
  all(bridge$study == "iccbridge"),
  all(bridge$covariance_shape_fixed),
  nrow(bridge_pairs) == 54L,
  all(bridge_pairs$n_arms == 2L),
  all(bridge_pairs$n_rows == 2L),
  all(bridge$simulation_seed_group == 10000L + bridge$bridge_pair_id),
  all(abs(bridge$slope_intercept_variance_ratio - 1) < 1e-12),
  all(paired_m10$sigma_range < 1e-10),
  all(paired_m10$reliability_range < 1e-10),
  all(paired_m10$icc_range < 1e-10),
  all(posterior_arm$reliability_range < 1e-10),
  all(icc_arm$icc_range < 1e-10),
  any(abs(
    bridge_first_stage$achieved_reliability -
      bridge_first_stage$posterior_reliability_anchor
  ) > 0.05)
)

# The paired seed is an executable design contract, not just a manifest label.
# At the m=10 anchor the arm-specific condition IDs differ but the derived
# seeds and generated data are exactly identical.
anchor_pair <- bridge %>%
  filter(
    mean_clus_size == 10L,
    marginal_rho == -0.5,
    posterior_reliability_anchor == 0.25,
    standardized_beta_target == 0.4
  ) %>%
  arrange(calibration_arm)
stopifnot(
  nrow(anchor_pair) == 2L,
  anchor_pair$condition_id[[1L]] != anchor_pair$condition_id[[2L]],
  vh_replication_seed(anchor_pair[1L, ], 7L) ==
    vh_replication_seed(anchor_pair[2L, ], 7L)
)
set.seed(vh_replication_seed(anchor_pair[1L, ], 7L))
anchor_sim_icc <- simulate_study2(anchor_pair[1L, ])
set.seed(vh_replication_seed(anchor_pair[2L, ], 7L))
anchor_sim_posterior <- simulate_study2(anchor_pair[2L, ])
stopifnot(
  identical(anchor_sim_icc$lv1, anchor_sim_posterior$lv1),
  identical(anchor_sim_icc$lv2_true, anchor_sim_posterior$lv2_true)
)

crosswalk <- make_icc_posterior_reliability_crosswalk()
stopifnot(
  nrow(crosswalk) == 189L,
  all(crosswalk$G_condition_number < 5),
  all(abs(
    crosswalk$intercept_icc[
      crosswalk$calibration_metric == "intercept_icc"
    ] -
      crosswalk$calibration_target_value[
        crosswalk$calibration_metric == "intercept_icc"
      ]
  ) < 1e-12),
  all(abs(
    crosswalk$achieved_marginal_reliability[
      crosswalk$calibration_metric == "marginal_slope"
    ] -
      crosswalk$calibration_target_value[
        crosswalk$calibration_metric == "marginal_slope"
      ]
  ) < 1e-8)
)

# Generator smoke tests ensure each existing pathway accepts its v2 wrapper.
set.seed(20260817)
s1_sim <- simulate_study1(study1 %>% slice(1L))
s2_sim <- simulate_study2(study2 %>% slice(1L))
s3_sim <- simulate_study3(study3 %>% slice(1L))
s4_sim <- simulate_study4(study4 %>% slice(1L))
stopifnot(
  nrow(s1_sim$lv2_true) == study1$num_clus[[1]],
  nrow(s2_sim$lv2_true) == study2$num_clus[[1]],
  nrow(s3_sim$lv2_true) == study3$num_clus[[1]],
  nrow(s4_sim$lv2_true) == study4$num_clus[[1]],
  setequal(study_methods_for_condition(study1 %>% slice(1L)), study1_methods()),
  setequal(study_methods_for_condition(study2 %>% slice(1L)), study2_methods()),
  setequal(study_methods_for_condition(study3 %>% slice(1L)), study3_methods()),
  setequal(study_methods_for_condition(study4 %>% slice(1L)), study4_methods()),
  setequal(study_methods_for_condition(bridge %>% slice(1L)), study2_methods())
)

cat("VH amended v2 design and ICC bridge tests ok\n")
