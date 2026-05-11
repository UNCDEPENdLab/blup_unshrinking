#!/usr/bin/env Rscript

# Unit check for the stage-2 HC2/HC3 leverage adjustment used by
# `stacked_sandwich_for_corrected_scores()`. This replaces the old
# `masters_sim/test_sim_edit3.R` scratch script with assertions.

make_stage2_hc_scores <- function(x_stage2, residuals, leverage_cap = 0.999) {
  h_mat <- x_stage2 %*% solve(crossprod(x_stage2)) %*% t(x_stage2)
  h_ii <- pmin(diag(h_mat), leverage_cap)

  list(
    h_ii = h_ii,
    hc2 = (residuals / sqrt(1 - h_ii)) * x_stage2,
    hc3 = (residuals / (1 - h_ii)) * x_stage2
  )
}

x <- c(-2, -1, 0, 1, 2, 8)
x_stage2 <- cbind(1, x)
residuals <- c(-0.4, 0.1, 0.3, -0.2, 0.5, -0.1)

out <- make_stage2_hc_scores(x_stage2, residuals)
expected_h <- diag(x_stage2 %*% solve(crossprod(x_stage2)) %*% t(x_stage2))
expected_h <- pmin(expected_h, 0.999)

stopifnot(
  isTRUE(all.equal(out$h_ii, expected_h, tolerance = 1e-12)),
  all(out$h_ii >= 0),
  all(out$h_ii < 1)
)

expected_hc2 <- (residuals / sqrt(1 - expected_h)) * x_stage2
expected_hc3 <- (residuals / (1 - expected_h)) * x_stage2

stopifnot(
  isTRUE(all.equal(out$hc2, expected_hc2, tolerance = 1e-12)),
  isTRUE(all.equal(out$hc3, expected_hc3, tolerance = 1e-12))
)

hc0_meat <- crossprod(residuals * x_stage2)
hc2_meat <- crossprod(out$hc2)
hc3_meat <- crossprod(out$hc3)

stopifnot(
  all(diag(hc2_meat) >= diag(hc0_meat)),
  all(diag(hc3_meat) >= diag(hc2_meat))
)

cap_out <- make_stage2_hc_scores(
  x_stage2 = diag(3),
  residuals = c(1, -2, 3),
  leverage_cap = 0.999
)

stopifnot(
  all(cap_out$h_ii == 0.999),
  all(is.finite(cap_out$hc2)),
  all(is.finite(cap_out$hc3))
)

cat("HC2/HC3 leverage adjustment test ok\n")
