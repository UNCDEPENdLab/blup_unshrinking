# Lai & Liu Study 1 condition manifest for the Vig--Hallquist estimator refresh.

lai_study1_vh_fixed_params <- list(
  gamma0 = 0,
  gamma1 = 0.5,
  beta_zu0 = 0.4,
  z_intercept = 1.5
)

#' Return Lai Study 1's balanced within-cluster random-effect design.
#'
#' The supplement scales the equally spaced values from -1 to 1 by their RMS,
#' rather than their sample SD.  Consequently, `crossprod(Z)[2, 2]` equals
#' the cluster size.  Keeping this in one helper makes the DGM and the
#' condition-level reliability diagnostic use exactly the same design.
lai_study1_vh_within_design <- function(clus_size) {
  clus_size <- as.integer(clus_size[[1]])
  if (!is.finite(clus_size) || clus_size < 2L) {
    stop("`clus_size` must be an integer of at least two.")
  }
  x <- seq(-1, 1, length.out = clus_size)
  x <- x / sqrt(mean(x^2))
  cbind(`(Intercept)` = 1, x = x)
}

#' Construct the true random-effect covariance for one Lai Study 1 condition.
lai_study1_vh_dgm_covariance <- function(icc, var_u1, cor_u0_u1) {
  icc <- as.numeric(icc[[1]])
  var_u1 <- as.numeric(var_u1[[1]])
  cor_u0_u1 <- as.numeric(cor_u0_u1[[1]])
  if (!is.finite(icc) || !is.finite(var_u1) || !is.finite(cor_u0_u1) ||
      icc <= 0 || var_u1 <= 0 || abs(cor_u0_u1) >= 1) {
    stop("Lai Study 1 random-effect covariance parameters are invalid.")
  }
  matrix(
    c(
      icc, cor_u0_u1 * sqrt(icc * var_u1),
      cor_u0_u1 * sqrt(icc * var_u1), var_u1
    ),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(c("u0", "u1"), c("u0", "u1"))
  )
}

#' Compute the true posterior slope reliability implied by a Lai Study 1 cell.
#'
#' This is a DGM property, not a calibrated input.  It uses the full
#' intercept--slope covariance matrix, so it remains valid when the random
#' effects are correlated.  In the special independent case it reduces to
#' `var_u1 * clus_size / (sigma2 + var_u1 * clus_size)` under Lai's RMS x
#' scaling.
lai_study1_vh_dgm_posterior_slope_reliability <- function(icc, var_u1,
                                                            cor_u0_u1, sigma2,
                                                            clus_size) {
  sigma2 <- as.numeric(sigma2[[1]])
  if (!is.finite(sigma2) || sigma2 <= 0) {
    stop("`sigma2` must be finite and positive.")
  }
  G <- lai_study1_vh_dgm_covariance(icc, var_u1, cor_u0_u1)
  Z <- lai_study1_vh_within_design(clus_size)
  V_post <- posterior_random_effect_covariance(
    G = G,
    Z = Z,
    R = sigma2 * diag(nrow(Z))
  )
  1 - V_post[2L, 2L] / G[2L, 2L]
}

#' Construct the original Lai Study 1 factorial DGM.
#'
#' The values and the `x` scaling are defined by the original supplement's
#' `simulation_scripts/sim1.R`.  Reliability is deliberately not calibrated in
#' this design: it remains an emergent property of ICC, slope variance, and
#' cluster size.
make_lai_study1_vh_design <- function(max_conditions = NA_integer_) {
  design <- tidyr::crossing(
    study = "lai_study1_vh",
    num_clus = c(30L, 100L, 500L),
    clus_size = c(3L, 10L, 25L),
    icc = c(0.05, 0.20, 0.50),
    vr_u1_u0 = c(0.5, 1.0, 2.0),
    cor_u0_u1 = c(-0.5, 0, 0.5),
    beta_zu1 = c(0, 0.4)
  ) %>%
    dplyr::mutate(
      sigma2 = 1 - icc,
      var_u1 = vr_u1_u0 * icc
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      dgm_population_slope_sd = sqrt(var_u1),
      dgm_posterior_slope_reliability = lai_study1_vh_dgm_posterior_slope_reliability(
        icc, var_u1, cor_u0_u1, sigma2, clus_size
      ),
      dgm_source = "Lai Liu supplement sim1.R",
      dgm_commit = "9ffe53168f6bb04e13ef977dc19a8d953d0bf29d",
      reporting_scales = "raw,latent_sd"
    ) %>%
    dplyr::ungroup()

  if (!is.na(max_conditions) && max_conditions > 0L) {
    design <- dplyr::slice_head(design, n = max_conditions)
  }

  dplyr::mutate(design, condition_id = dplyr::row_number())
}
