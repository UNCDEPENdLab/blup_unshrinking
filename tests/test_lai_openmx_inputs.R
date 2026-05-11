#!/usr/bin/env Rscript

# Unit checks for shared Lai/OpenMx input construction helpers that are used by
# both the Lai replication and the sandwich coverage simulation.

suppressPackageStartupMessages({
  library(lme4)
})

source(file.path("R", "core_utils.R"), local = TRUE)
source(file.path("R", "lai_openmx_helpers.R"), local = TRUE)

cluster_df <- data.frame(z = c(-1, 0, 1))
stopifnot(
  isTRUE(all.equal(default_re_design(cluster_df), matrix(1, nrow = 3L, ncol = 1L))),
  isTRUE(all.equal(default_re_design(cluster_df, within_var = "z"), cbind(1, cluster_df$z)))
)

uni_row <- make_eb_output_row(
  id = "a",
  eb = 0.25,
  post_vcov = matrix(0.10, nrow = 1L),
  lambda = matrix(0.80, nrow = 1L),
  theta = matrix(0.20, nrow = 1L),
  prefix = "z_"
)
stopifnot(
  identical(names(uni_row), c("id", "z_u0_eb", "z_postvar11", "z_lambda11", "z_theta11")),
  isTRUE(all.equal(uni_row$z_u0_eb, 0.25))
)

bi_row <- make_eb_output_row(
  id = "b",
  eb = c(0.10, -0.20),
  post_vcov = matrix(c(0.40, 0.05, 0.05, 0.30), nrow = 2L),
  lambda = matrix(c(0.90, 0.10, 0.05, 0.85), nrow = 2L),
  theta = matrix(c(0.10, 0.02, 0.02, 0.15), nrow = 2L)
)
expected_bi_names <- c(
  "id", "u0_eb", "u1_eb", "postvar11", "postvar12", "postvar22",
  "lambda11", "lambda12", "lambda21", "lambda22", "theta11", "theta12", "theta22"
)
stopifnot(
  identical(names(bi_row), expected_bi_names),
  isTRUE(all.equal(bi_row$u1_eb, -0.20))
)

set.seed(2026)
n_id <- 10L
n_obs <- 5L
base_ids <- sprintf("cluster_%02d", seq_len(n_id))
stage2_ids <- rev(base_ids)
id_df <- data.frame(
  id = stage2_ids,
  x = seq(-1, 1, length.out = n_id)
)

sim_dat <- do.call(rbind, lapply(seq_along(base_ids), function(i) {
  z <- seq(-1, 1, length.out = n_obs)
  u0 <- rnorm(1, sd = 0.45)
  u1 <- rnorm(1, sd = 0.35)
  data.frame(
    id = base_ids[[i]],
    z = z,
    y = 0.2 + 0.6 * z + u0 + u1 * z + rnorm(n_obs, sd = 0.20)
  )
}))

fit <- lmer(y ~ 1 + z + (1 + z | id), data = sim_dat, REML = FALSE)
split_dat <- split(sim_dat, sim_dat$id)
out <- compute_lai_2spa_inputs(fit, split_dat, id_df)

required_cols <- c(
  "id", "x", "u0_eb", "u1_eb", "postvar11", "postvar12", "postvar22",
  "lambda11", "lambda12", "lambda21", "lambda22", "theta11", "theta12", "theta22"
)
missing_cols <- setdiff(required_cols, names(out))
if (length(missing_cols) > 0L) {
  stop("Missing expected Lai input columns: ", paste(missing_cols, collapse = ", "))
}

stopifnot(
  identical(out$id, id_df$id),
  isTRUE(all.equal(out$x, id_df$x)),
  all(is.finite(out$u0_eb)),
  all(is.finite(out$u1_eb)),
  all(is.finite(as.matrix(out[, setdiff(required_cols, c("id", "x")), drop = FALSE])))
)

cat("Lai/OpenMx input helper tests ok\n")
