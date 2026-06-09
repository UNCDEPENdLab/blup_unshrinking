# Posterior-reliability calibration for random-intercept/random-slope designs.
#
# This module converts design-level targets into the population parameters used
# by `simulate_dataset()`. The calibration has two distinct stages:
#
# 1. Solve for the marginal random-slope variance that produces a requested
#    expected posterior reliability under the planned Z and R matrices.
# 2. Decompose that marginal slope into `gamma * x + u1` for a requested
#    structural R-squared, preserving the marginal intercept-slope covariance.
#
# The target first-stage covariance and the covariance used to draw residual
# random effects are therefore different objects:
#
#   G_marginal = Var(b0, b1), where b1 = gamma * x + u1
#   G_residual = Var(b0, u1)
#
# In particular, the marginal correlation belongs in the reliability
# calculation, while the residual correlation belongs in `simulate_dataset()`.
# Keeping both explicit prevents the structural effect from silently changing
# the intended first-stage covariance structure.

#' Construct the within-cluster random-effect design used by simulate_dataset().
#'
#' @details
#' The simulation generator creates equally spaced time values on `[-1, 1]`,
#' centers them, and by default divides by their sample standard deviation.
#' This helper reproduces that transformation exactly so calibration uses the
#' same random-effect design matrix as data generation.
#'
#' With `standardize = TRUE`, the slope column satisfies
#' `sum(z^2) = n - 1`. Centering also makes the intercept and slope columns
#' orthogonal in balanced iid conditions. The returned matrix is the
#' random-effect design `Z_i`, not the fixed-effect design `X_i`.
#'
#' @param n Integer cluster size, at least 2.
#' @param standardize Logical; standardize the centered time scores to unit
#'   sample standard deviation. This should normally remain `TRUE` to match
#'   `simulate_dataset()`.
#'
#' @return A numeric `n x 2` matrix with columns named `intercept` and `slope`.
make_reliability_time_design <- function(n, standardize = TRUE) {
  n <- as.integer(n[[1]])
  if (!is.finite(n) || n < 2L) {
    stop("`n` must be one integer >= 2.")
  }

  # Match the time-score construction in `simulate_dataset()` rather than using
  # an abstract information approximation based only on cluster size.
  z <- seq(-1, 1, length.out = n)
  z <- z - mean(z)
  if (isTRUE(standardize)) {
    # `scale(..., scale = TRUE)` uses the sample SD, so S_zz equals n - 1.
    z <- as.numeric(scale(z, center = FALSE, scale = TRUE))
  }
  cbind(intercept = 1, slope = z)
}

#' Construct an iid, AR(1), or Toeplitz residual covariance matrix.
#'
#' This standalone constructor mirrors `make_R_matrix()` without requiring the
#' simulation helper module to have been sourced first. It supports the same
#' residual structures used by the simulation machinery:
#'
#' - `"iid"`: identity correlation matrix;
#' - `"ar1"`: correlation `rho^abs(j - k)`;
#' - `"toeplitz"`: user-supplied lag correlations, padded with zero at
#'   unspecified higher lags.
#'
#' The returned matrix is `R_i = sigma^2 C_i`. A Cholesky factorization is used
#' as a final positive-definiteness check because the posterior calculation
#' requires an invertible residual covariance matrix.
#'
#' @param n Integer cluster size.
#' @param sigma Positive residual standard deviation.
#' @param r_spec Residual covariance specification. May be `NULL`, a character
#'   scalar, or a list containing `structure`. AR(1) specifications require
#'   `rho`; Toeplitz specifications require `correlations` or use `rho` as a
#'   backward-compatible alias.
#'
#' @return A positive-definite numeric `n x n` residual covariance matrix.
make_reliability_residual_covariance <- function(n, sigma = 1, r_spec = NULL) {
  n <- as.integer(n[[1]])
  sigma <- as.numeric(sigma[[1]])
  if (!is.finite(n) || n < 1L || !is.finite(sigma) || sigma <= 0) {
    stop("`n` must be positive and `sigma` must be finite and positive.")
  }

  if (is.null(r_spec)) {
    r_spec <- list(structure = "iid")
  } else if (is.character(r_spec) && length(r_spec) == 1L) {
    r_spec <- list(structure = r_spec)
  }
  if (!is.list(r_spec)) {
    stop("`r_spec` must be NULL, a character scalar, or a list.")
  }

  structure_name <- r_spec$structure
  if (is.null(structure_name)) {
    structure_name <- r_spec$type
  }
  if (is.null(structure_name) || length(structure_name) != 1L ||
      !is.character(structure_name)) {
    stop("`r_spec` must include one character `structure` field.")
  }
  structure_name <- tolower(as.character(structure_name[[1]]))

  # Build a correlation matrix first, then apply the residual variance scale.
  # This separation makes the role of sigma in reliability calibration clear.
  cor_mat <- switch(
    structure_name,
    iid = diag(n),
    ar1 = {
      rho <- as.numeric(r_spec$rho[[1]])
      if (!is.finite(rho) || abs(rho) >= 1) {
        stop("AR(1) `r_spec` requires `rho` with absolute value < 1.")
      }
      idx <- seq_len(n)
      outer(idx, idx, function(a, b) rho^abs(a - b))
    },
    toeplitz = {
      correlations <- r_spec$correlations
      if (is.null(correlations)) {
        correlations <- r_spec$rho
      }
      correlations <- as.numeric(correlations)
      if (length(correlations) < 1L || any(!is.finite(correlations)) ||
          any(abs(correlations) >= 1)) {
        stop("Toeplitz `r_spec` requires finite lag correlations with absolute values < 1.")
      }
      lag_cor <- numeric(n)
      lag_cor[[1]] <- 1
      n_lag <- min(n - 1L, length(correlations))
      if (n_lag > 0L) {
        # Unspecified lags remain zero, matching `make_R_matrix()`.
        lag_cor[seq_len(n_lag) + 1L] <- correlations[seq_len(n_lag)]
      }
      stats::toeplitz(lag_cor)
    },
    stop("Unsupported residual covariance structure: ", structure_name)
  )

  R <- sigma^2 * cor_mat
  # Fail during condition construction rather than deep inside a simulation
  # replication if the requested residual structure is inadmissible.
  tryCatch(
    chol(R),
    error = function(e) stop("The requested residual covariance matrix is not positive definite.")
  )
  R
}

#' Construct a bivariate marginal random-effects covariance matrix.
#'
#' @details
#' Given marginal intercept variance `tau0^2`, marginal slope variance
#' `tau1^2`, and correlation `rho`, this function constructs
#'
#' `G = [[tau0^2, rho*tau0*tau1], [rho*tau0*tau1, tau1^2]]`.
#'
#' In posterior-reliability calibration this is the covariance of the total
#' random intercept and total random slope seen by the first-stage model. It is
#' not necessarily the covariance used to draw residual random effects after a
#' structural predictor is added to the slope.
#'
#' @param intercept_variance Positive marginal random-intercept variance.
#' @param slope_variance Positive marginal random-slope variance.
#' @param intercept_slope_correlation Correlation between the marginal random
#'   intercept and marginal random slope; must have absolute value below one.
#'
#' @return A symmetric `2 x 2` covariance matrix with intercept/slope dimnames.
make_random_effect_covariance <- function(
    intercept_variance,
    slope_variance,
    intercept_slope_correlation) {
  intercept_variance <- as.numeric(intercept_variance[[1]])
  slope_variance <- as.numeric(slope_variance[[1]])
  rho <- as.numeric(intercept_slope_correlation[[1]])

  if (!is.finite(intercept_variance) || intercept_variance <= 0 ||
      !is.finite(slope_variance) || slope_variance <= 0) {
    stop("Random-effect variances must be finite and positive.")
  }
  if (!is.finite(rho) || abs(rho) >= 1) {
    stop("`intercept_slope_correlation` must have absolute value < 1.")
  }

  # Convert the scale-free correlation into the covariance needed by G.
  covariance <- rho * sqrt(intercept_variance * slope_variance)
  matrix(
    c(intercept_variance, covariance, covariance, slope_variance),
    nrow = 2L,
    dimnames = list(c("intercept", "slope"), c("intercept", "slope"))
  )
}

#' Compute the conditional posterior covariance for one cluster.
#'
#' @details
#' For the Gaussian mixed model with prior covariance `G`, random-effect design
#' `Z`, and residual covariance `R`, the conditional covariance is
#'
#' `V_post = (G^{-1} + Z' R^{-1} Z)^{-1}`.
#'
#' This covariance depends only on the design and variance components, not on
#' the observed response values. That property allows reliability to be
#' calibrated before any Monte Carlo outcomes are generated.
#'
#' @param G Positive-definite `2 x 2` marginal random-effects covariance matrix.
#' @param Z Numeric `n x 2` random-effect design matrix.
#' @param R Positive-definite `n x n` residual covariance matrix.
#'
#' @return The symmetric `2 x 2` posterior covariance matrix for one cluster.
posterior_random_effect_covariance <- function(G, Z, R) {
  if (!all(dim(G) == c(2L, 2L)) || ncol(Z) != 2L ||
      nrow(R) != nrow(Z) || ncol(R) != nrow(Z)) {
    stop("Dimensions of `G`, `Z`, and `R` are incompatible.")
  }
  # `solve(R, Z)` avoids explicitly constructing R^{-1}; `crossprod` then
  # forms the likelihood information Z'R^{-1}Z.
  likelihood_information <- crossprod(Z, solve(R, Z))
  solve(solve(G) + likelihood_information)
}

#' Compute expected posterior slope reliability across reference clusters.
#'
#' @details
#' Slope reliability for cluster type `i` is defined as
#'
#' `1 - V_post_i[2, 2] / G[2, 2]`.
#'
#' This function averages the posterior-variance ratio before subtracting it
#' from one. When `weights` are supplied, `Z_list` and `R_list` may contain only
#' unique cluster types and the weights record how often each type occurs in a
#' deterministic reference population.
#'
#' The target is an expected design-level reliability. Individual clusters can
#' have different reliabilities when cluster sizes or residual structures vary.
#'
#' @param G Positive-definite `2 x 2` marginal random-effects covariance matrix.
#' @param Z_list Nonempty list of `n_i x 2` random-effect design matrices.
#' @param R_list Nonempty list of matching `n_i x n_i` residual covariance
#'   matrices.
#' @param weights Optional nonnegative frequency weights, one per list element.
#'   Defaults to equal weights.
#'
#' @return A numeric scalar in principle between zero and one.
expected_slope_reliability <- function(G, Z_list, R_list, weights = NULL) {
  if (length(Z_list) == 0L || length(Z_list) != length(R_list)) {
    stop("`Z_list` and `R_list` must be nonempty lists of equal length.")
  }
  if (is.null(weights)) {
    weights <- rep(1, length(Z_list))
  }
  weights <- as.numeric(weights)
  if (length(weights) != length(Z_list) || any(!is.finite(weights)) ||
      any(weights < 0) || sum(weights) <= 0) {
    stop("`weights` must be nonnegative, finite, and match `Z_list`.")
  }

  # Only the slope diagonal is needed for the scalar reliability target, but
  # each posterior is computed from the full 2 x 2 covariance structure.
  posterior_slope_variances <- vapply(
    seq_along(Z_list),
    function(i) posterior_random_effect_covariance(
      G,
      Z_list[[i]],
      R_list[[i]]
    )[2L, 2L],
    numeric(1L)
  )
  # Average the remaining uncertainty fraction, then convert to reliability.
  1 - stats::weighted.mean(
    posterior_slope_variances / G[2L, 2L],
    w = weights
  )
}

#' Solve for marginal slope variance from a posterior-reliability target.
#'
#' @details
#' Reliability alone does not identify the full random-effects covariance
#' matrix. This function treats the intercept variance, marginal
#' intercept-slope correlation, `Z_i`, and `R_i` as fixed, and solves only for
#' the marginal slope variance.
#'
#' The root problem is
#'
#' `expected_slope_reliability(G(tau1^2)) - target_reliability = 0`.
#'
#' The search is performed on `log(tau1^2)`. This guarantees a positive
#' candidate variance and provides a stable search over many orders of
#' magnitude. Reliability at the two search bounds is also used to detect
#' infeasible targets before calling `uniroot()`.
#'
#' @param target_reliability Desired expected posterior slope reliability,
#'   strictly between zero and one.
#' @param Z_list Nonempty list of random-effect design matrices.
#' @param R_list Nonempty list of matching residual covariance matrices.
#' @param weights Optional nonnegative frequency weights for unique reference
#'   cluster types.
#' @param intercept_variance Fixed positive marginal intercept variance.
#' @param intercept_slope_correlation Fixed marginal intercept-slope
#'   correlation.
#' @param log_variance_bounds Length-two numeric interval for
#'   `log(slope_variance)`. The defaults cover an intentionally wide range.
#' @param tolerance Numerical tolerance passed to `uniroot()` and used for the
#'   attainable-range check.
#'
#' @return A list containing the target and achieved reliability, calibrated
#'   marginal slope variance, marginal `G` matrix, posterior covariance for each
#'   supplied cluster type, and reliability values at the search bounds.
calibrate_slope_variance <- function(
    target_reliability,
    Z_list,
    R_list,
    weights = NULL,
    intercept_variance = 1,
    intercept_slope_correlation = 0,
    log_variance_bounds = c(-30, 30),
    tolerance = 1e-10) {
  target_reliability <- as.numeric(target_reliability[[1]])
  if (!is.finite(target_reliability) ||
      target_reliability <= 0 || target_reliability >= 1) {
    stop("`target_reliability` must be strictly between 0 and 1.")
  }

  # Work on the log-variance scale so every candidate is strictly positive.
  reliability_at_log_variance <- function(log_slope_variance) {
    G <- make_random_effect_covariance(
      intercept_variance = intercept_variance,
      slope_variance = exp(log_slope_variance),
      intercept_slope_correlation = intercept_slope_correlation
    )
    expected_slope_reliability(G, Z_list, R_list, weights = weights)
  }

  # Correlated intercepts and slopes can impose a nonzero lower bound on slope
  # reliability. Evaluate the requested search interval before root finding so
  # infeasible design cells fail with an interpretable message.
  attainable <- vapply(
    log_variance_bounds,
    reliability_at_log_variance,
    numeric(1L)
  )
  if (target_reliability < attainable[[1L]] - tolerance ||
      target_reliability > attainable[[2L]] + tolerance) {
    stop(sprintf(
      paste0(
        "Target reliability %.4f is outside the attainable interval ",
        "[%.4f, %.4f] for this design and fixed covariance structure."
      ),
      target_reliability, attainable[[1L]], attainable[[2L]]
    ))
  }

  # Reliability is monotone in the marginal slope variance for the supported
  # Gaussian designs, making a one-dimensional bracketed root solve sufficient.
  root <- stats::uniroot(
    function(log_variance) {
      reliability_at_log_variance(log_variance) - target_reliability
    },
    interval = log_variance_bounds,
    tol = tolerance
  )

  # Convert the solved log variance back to the scale used in G and simulation.
  slope_variance <- exp(root$root)
  G <- make_random_effect_covariance(
    intercept_variance = intercept_variance,
    slope_variance = slope_variance,
    intercept_slope_correlation = intercept_slope_correlation
  )

  list(
    target_reliability = target_reliability,
    achieved_reliability = expected_slope_reliability(
      G,
      Z_list,
      R_list,
      weights = weights
    ),
    intercept_variance = as.numeric(intercept_variance[[1]]),
    slope_variance_marginal = slope_variance,
    intercept_slope_correlation_marginal = as.numeric(intercept_slope_correlation[[1]]),
    G_marginal = G,
    # Retain these matrices for diagnostics and independent verification of the
    # calibrated condition.
    posterior_covariances = Map(
      function(Z, R) posterior_random_effect_covariance(G, Z, R),
      Z_list,
      R_list
    ),
    attainable_reliability = stats::setNames(
      attainable,
      c("lower_bound", "upper_bound")
    )
  )
}

#' Decompose a calibrated marginal slope using a structural R-squared.
#'
#' @details
#' The first-stage model sees a total slope
#'
#' `b1_i = gamma * x_i + u1_i`,
#'
#' where `Var(x) = 1` and `x` is independent of `(b0, u1)`. If structural
#' `R^2` is the fraction of marginal slope variance explained by `x`, then
#'
#' `gamma^2 = R^2 * Var(b1)`
#'
#' and
#'
#' `Var(u1) = (1 - R^2) * Var(b1)`.
#'
#' The marginal covariance `Cov(b0, b1)` is preserved as `Cov(b0, u1)`.
#' Because the residual slope variance is smaller than the marginal slope
#' variance, the residual correlation passed to `simulate_dataset()` generally
#' differs from the marginal correlation used during reliability calibration.
#'
#' @param calibration Output from `calibrate_slope_variance()`.
#' @param structural_r_squared Fraction of marginal slope variance explained by
#'   the standardized structural predictor; must be in `[0, 1)`.
#' @param effect_sign Sign of the structural coefficient. Positive values yield
#'   a positive coefficient and negative values yield a negative coefficient.
#'
#' @return A list containing the raw structural coefficient, standardized beta,
#'   residual slope variance and standard deviation, preserved covariance,
#'   residual intercept-slope correlation, and residual `G` matrix used to draw
#'   `(b0, u1)`.
decompose_structural_slope <- function(
    calibration,
    structural_r_squared = 0,
    effect_sign = 1) {
  structural_r_squared <- as.numeric(structural_r_squared[[1]])
  effect_sign <- sign(as.numeric(effect_sign[[1]]))
  if (!is.finite(structural_r_squared) ||
      structural_r_squared < 0 || structural_r_squared >= 1) {
    stop("`structural_r_squared` must be in [0, 1).")
  }
  if (effect_sign == 0) {
    stop("`effect_sign` must be positive or negative.")
  }

  G_marginal <- calibration$G_marginal
  total_slope_variance <- G_marginal[2L, 2L]
  covariance_01 <- G_marginal[1L, 2L]
  # With Var(x) = 1, gamma^2 is exactly the explained slope variance.
  gamma <- effect_sign * sqrt(structural_r_squared * total_slope_variance)
  residual_slope_variance <- total_slope_variance - gamma^2

  # Preserve the marginal intercept-slope covariance. Since x is independent of
  # b0, Cov(b0, b1) = Cov(b0, u1).
  G_residual <- matrix(
    c(
      G_marginal[1L, 1L], covariance_01,
      covariance_01, residual_slope_variance
    ),
    nrow = 2L,
    dimnames = dimnames(G_marginal)
  )
  # A large structural R-squared can leave too little residual slope variance
  # to support the requested covariance. Detect that before simulation.
  eigenvalues <- eigen(G_residual, symmetric = TRUE, only.values = TRUE)$values
  if (min(eigenvalues) <= sqrt(.Machine$double.eps)) {
    stop(sprintf(
      paste0(
        "The requested structural R-squared %.3f is incompatible with the ",
        "calibrated intercept-slope covariance. Residual G is not positive definite."
      ),
      structural_r_squared
    ))
  }

  # This is the correlation expected by draw_random_effects() and
  # simulate_dataset(); it is not the target marginal correlation.
  residual_correlation <- covariance_01 / sqrt(
    G_residual[1L, 1L] * G_residual[2L, 2L]
  )

  list(
    gamma_x_on_slope = gamma,
    standardized_beta = effect_sign * sqrt(structural_r_squared),
    structural_r_squared = structural_r_squared,
    slope_variance_marginal = total_slope_variance,
    slope_variance_residual = residual_slope_variance,
    slope_sd_residual = sqrt(residual_slope_variance),
    intercept_slope_covariance_residual = covariance_01,
    intercept_slope_correlation_residual = residual_correlation,
    G_residual = G_residual
  )
}

#' Build a deterministic cluster-size profile for condition-level calibration.
#'
#' The profile approximates the population distribution of cluster sizes without
#' using replication-specific random draws. Calibrating against it once when the
#' condition grid is built keeps population parameters fixed across replications.
#'
#' @details
#' The returned counts reproduce the cluster-size mechanisms in
#' `simulate_dataset()` using deterministic midpoint quantiles instead of random
#' draws:
#'
#' - `"balanced"` returns one cluster size equal to `mean_n_trial`;
#' - `"unbalanced"` maps uniform midpoint quantiles through the simulation's
#'   `[0.6, 1.4] * mean_n_trial` rule;
#' - `"informative_unbalanced"` maps deterministic rank percentiles through the
#'   same decreasing power function used to relate cluster size to `x`.
#'
#' For unbalanced conditions, increasing `n_reference` gives a finer
#' approximation to the population distribution. These counts are calibration
#' devices only; each Monte Carlo replication still generates its own realized
#' cluster sizes using `simulate_dataset()`.
#'
#' @param mean_n_trial Target mean cluster size used by the simulation.
#' @param balance_mode One of `"balanced"`, `"unbalanced"`, or
#'   `"informative_unbalanced"`.
#' @param min_n_trial Lower bound for random unbalanced counts.
#' @param highly_unbalanced_min_n_trial Lower bound for informative unbalanced
#'   counts.
#' @param highly_unbalanced_power Positive power controlling the informative
#'   imbalance profile.
#' @param n_reference Number of deterministic midpoint quantiles for unbalanced
#'   profiles. Ignored in effect for balanced profiles, which need one type.
#'
#' @return An integer vector of deterministic reference cluster sizes.
reference_trial_counts <- function(
    mean_n_trial,
    balance_mode = "balanced",
    min_n_trial = 2L,
    highly_unbalanced_min_n_trial = 2L,
    highly_unbalanced_power = 3,
    n_reference = 1001L) {
  mean_n_trial <- as.integer(mean_n_trial[[1]])
  min_n_trial <- as.integer(min_n_trial[[1]])
  highly_unbalanced_min_n_trial <- as.integer(highly_unbalanced_min_n_trial[[1]])
  highly_unbalanced_power <- as.numeric(highly_unbalanced_power[[1]])
  n_reference <- as.integer(n_reference[[1]])
  balance_mode <- as.character(balance_mode[[1]])

  if (mean_n_trial < 2L || n_reference < 2L) {
    stop("`mean_n_trial` and `n_reference` must be at least 2.")
  }

  if (identical(balance_mode, "balanced")) {
    return(mean_n_trial)
  }

  # Midpoint quantiles avoid random-number generation and exclude exact 0/1
  # endpoints, making the profile deterministic and stable across runs.
  probabilities <- (seq_len(n_reference) - 0.5) / n_reference
  if (identical(balance_mode, "unbalanced")) {
    # Deterministic analogue of runif(n, 0.6, 1.4) in simulate_dataset().
    return(pmax(
      min_n_trial,
      as.integer(round((0.6 + 0.8 * probabilities) * mean_n_trial))
    ))
  }

  if (identical(balance_mode, "informative_unbalanced")) {
    if (mean_n_trial < highly_unbalanced_min_n_trial) {
      stop("`mean_n_trial` must be >= `highly_unbalanced_min_n_trial`.")
    }
    # The generator gives high-x subjects fewer observations. Reliability
    # calibration needs the resulting marginal count distribution; the
    # structural association with x is retained during actual simulation.
    reverse_percentile <- 1 - probabilities
    amplitude <- (mean_n_trial - highly_unbalanced_min_n_trial) *
      (highly_unbalanced_power + 1)
    return(pmax(
      highly_unbalanced_min_n_trial,
      as.integer(round(
        highly_unbalanced_min_n_trial +
          amplitude * reverse_percentile^highly_unbalanced_power
      ))
    ))
  }

  stop("Unsupported `balance_mode`: ", balance_mode)
}

#' Build deterministic Z and R reference lists for one simulation condition.
#'
#' @details
#' This function converts the deterministic count profile into the design
#' objects needed by `calibrate_slope_variance()`. Repeated cluster sizes are
#' compressed into unique `Z_i`/`R_i` pairs plus integer frequency weights.
#' Compression substantially reduces repeated matrix inversions while producing
#' the same weighted expected reliability as the full reference profile.
#'
#' @param mean_n_trial Target mean cluster size.
#' @param sigma Positive residual standard deviation.
#' @param balance_mode Cluster-size mechanism; see `reference_trial_counts()`.
#' @param min_n_trial Lower bound for random unbalanced counts.
#' @param highly_unbalanced_min_n_trial Lower bound for informative unbalanced
#'   counts.
#' @param highly_unbalanced_power Power controlling informative imbalance.
#' @param r_spec Residual covariance specification passed to
#'   `make_reliability_residual_covariance()`.
#' @param n_reference Number of deterministic reference quantiles.
#'
#' @return A list containing the full reference count profile, unique counts,
#'   frequency weights, matching `Z_list` and `R_list`, and count summaries.
make_reliability_reference_design <- function(
    mean_n_trial,
    sigma,
    balance_mode = "balanced",
    min_n_trial = 2L,
    highly_unbalanced_min_n_trial = 2L,
    highly_unbalanced_power = 3,
    r_spec = NULL,
    n_reference = 1001L) {
  counts <- reference_trial_counts(
    mean_n_trial = mean_n_trial,
    balance_mode = balance_mode,
    min_n_trial = min_n_trial,
    highly_unbalanced_min_n_trial = highly_unbalanced_min_n_trial,
    highly_unbalanced_power = highly_unbalanced_power,
    n_reference = n_reference
  )

  # Store one matrix pair per unique cluster size and carry multiplicities as
  # weights. This is exact for the deterministic reference profile.
  unique_counts <- sort(unique(counts))
  count_weights <- tabulate(match(counts, unique_counts), nbins = length(unique_counts))
  list(
    trial_counts = counts,
    unique_trial_counts = unique_counts,
    count_weights = count_weights,
    Z_list = lapply(unique_counts, make_reliability_time_design),
    R_list = lapply(
      unique_counts,
      make_reliability_residual_covariance,
      sigma = sigma,
      r_spec = r_spec
    ),
    mean_trial_count = mean(counts),
    min_trial_count = min(counts),
    max_trial_count = max(counts)
  )
}

#' Calibrate all simulation inputs for one reliability/R-squared condition.
#'
#' @details
#' This is the high-level condition-calibration entry point used by simulation
#' design builders. It:
#'
#' 1. constructs a deterministic reference distribution of cluster designs;
#' 2. solves for the marginal slope variance that attains the reliability target;
#' 3. decomposes the marginal slope according to structural `R^2`;
#' 4. returns the residual slope SD and residual correlation expected by
#'    `simulate_dataset()`.
#'
#' Call this function once when a condition manifest is constructed, not once
#' per Monte Carlo replication. Storing its scalar outputs in the manifest keeps
#' the population parameters fixed while realized subjects, cluster sizes, and
#' outcomes vary across replications.
#'
#' @param target_reliability Desired expected posterior slope reliability.
#' @param structural_r_squared Fraction of marginal slope variance explained by
#'   the standardized cluster-level predictor.
#' @param marginal_rho Target correlation between the total random intercept and
#'   total random slope in the first-stage marginal `G`.
#' @param tau0 Positive marginal random-intercept standard deviation.
#' @param mean_n_trial Target mean cluster size.
#' @param sigma Positive residual standard deviation.
#' @param balance_mode Cluster-size mechanism.
#' @param min_n_trial Lower bound for random unbalanced counts.
#' @param highly_unbalanced_min_n_trial Lower bound for informative unbalanced
#'   counts.
#' @param highly_unbalanced_power Power controlling informative imbalance.
#' @param r_spec Residual covariance specification.
#' @param n_reference Number of deterministic reference quantiles used for
#'   unbalanced calibration.
#' @param effect_sign Sign of the structural coefficient.
#'
#' @return A list containing:
#'
#' - target and achieved posterior reliability;
#' - structural `R^2` and standardized beta;
#' - marginal and residual intercept-slope correlations;
#' - cell-specific `gamma_x_on_slope`;
#' - marginal and residual slope variances;
#' - residual slope SD (`tau1_residual`);
#' - marginal and residual `G` matrices;
#' - summaries of the deterministic reference count profile.
calibrate_random_slope_condition <- function(
    target_reliability,
    structural_r_squared,
    marginal_rho,
    tau0,
    mean_n_trial,
    sigma,
    balance_mode = "balanced",
    min_n_trial = 2L,
    highly_unbalanced_min_n_trial = 2L,
    highly_unbalanced_power = 3,
    r_spec = NULL,
    n_reference = 1001L,
    effect_sign = 1) {
  # The reference profile is deterministic, so repeated calls for the same
  # condition return identical population parameters.
  reference <- make_reliability_reference_design(
    mean_n_trial = mean_n_trial,
    sigma = sigma,
    balance_mode = balance_mode,
    min_n_trial = min_n_trial,
    highly_unbalanced_min_n_trial = highly_unbalanced_min_n_trial,
    highly_unbalanced_power = highly_unbalanced_power,
    r_spec = r_spec,
    n_reference = n_reference
  )
  # Reliability is defined using the marginal G seen by a Stage-1 model that
  # omits the structural predictor x.
  calibration <- calibrate_slope_variance(
    target_reliability = target_reliability,
    Z_list = reference$Z_list,
    R_list = reference$R_list,
    weights = reference$count_weights,
    intercept_variance = tau0^2,
    intercept_slope_correlation = marginal_rho
  )
  # Convert the marginal slope into the residual random-effect parameters that
  # the existing data generator already accepts.
  structural <- decompose_structural_slope(
    calibration,
    structural_r_squared = structural_r_squared,
    effect_sign = effect_sign
  )

  list(
    target_reliability = calibration$target_reliability,
    achieved_reliability = calibration$achieved_reliability,
    structural_r_squared = structural$structural_r_squared,
    standardized_beta = structural$standardized_beta,
    marginal_rho = as.numeric(marginal_rho[[1]]),
    rho_residual = structural$intercept_slope_correlation_residual,
    gamma_x_on_slope = structural$gamma_x_on_slope,
    slope_variance_marginal = structural$slope_variance_marginal,
    slope_variance_residual = structural$slope_variance_residual,
    tau1_residual = structural$slope_sd_residual,
    G_marginal = calibration$G_marginal,
    G_residual = structural$G_residual,
    reference_mean_n_trial = reference$mean_trial_count,
    reference_min_n_trial = reference$min_trial_count,
    reference_max_n_trial = reference$max_trial_count,
    calibration_reference_n = length(reference$trial_counts)
  )
}
