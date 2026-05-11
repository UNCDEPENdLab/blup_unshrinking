# MLM BLUP Simulations

This repository contains reusable R helpers, simulation drivers, and technical notes for work on empirical-Bayes/BLUP score correction in multilevel models. The code is script-oriented rather than packaged as an installed R package: simulation drivers source helper files from the repository-level `R/` directory.

## Repository Layout

- `R/`: shared helper modules used by multiple simulations.
- `lai_replication/`: Lai and Liu apples-to-apples replication modules and runner.
- `documentation/`: technical notes, manuscripts, references, and rendered PDFs.
- `archive/`: removed legacy or experimental helpers retained for reference.
- `src/`: TMB C++ source used by the TMB derivative backend.
- `tests/`: lightweight checks for shared calculations.

Top-level simulation and analysis scripts include:

- `mlm_blup_correction_sim.R`: simple random-intercept BLUP correction simulation.
- `mlm_random_slope_blup_correction_sim.R`: random-slope correction simulation.
- `mlm_random_slope_blup_predictor_comparison_sim.R`: predictor comparison simulation.
- `mlm_random_slope_blup_sandwich_coverage_sim.R`: random-slope corrected-score sandwich coverage simulation.
- `plot_blup_correction_summaries.R`, `plot_zoomed.R`, and `anova_analysis.R`: plotting and summary analyses for generated outputs.

## Shared Helpers

- `R/core_utils.R`: small shared utilities such as NA-stable means, compact diagnostics, quiet `lmer()` fitting, and matrix regularization.
- `R/blup_helpers.R`: EB/BLUP extraction, Vig-style prior unweighting, and closed-form corrected cluster scores.
- `R/sim_helpers.R`: common simulation data-generation helpers.
- `R/sim_diagnostics.R`: first-stage singularity, EB collinearity, and design-condition diagnostics.
- `R/stats_helpers.R`: extraction of estimates, standard errors, and intervals from `lm` and `lmer` fits.
- `R/stage2_estimators.R`: observed-score, HC3, ridge, EIV, and stacked-sandwich result-row formatting helpers.
- `R/derivative_backends.R`: finite-difference, `numDeriv`, `merDeriv`, TMB, and analytical derivative backend selection.
- `R/tmb_stage1_helpers.R`: TMB model compilation/object/hessian helpers for the stage-1 Gaussian model.
- `R/stacked_sandwich_helpers.R`: parameter packing, cluster likelihoods, corrected-score extraction, and stacked M-estimation sandwich assembly.
- `R/lai_openmx_helpers.R`: shared Lai-style EB measurement-model inputs, OpenMx retry/status extraction, and 2S-PA/2S-PAA wrappers.

## Lai Replication

`lai_replication/` contains the Lai apples-to-apples simulation driver, study-specific modules, condition runner, and smoke test. Shared EB measurement, OpenMx, score calculation, and stage-2 estimator code lives under root `R/` so non-Lai simulation scripts can reuse it. See `lai_replication/README.md` for the file-by-file structure.

The command-line entry point is `lai_replication/mlm_random_slope_lai_apples_to_apples_sim.R`.

Small smoke run:

```sh
Rscript lai_replication/mlm_random_slope_lai_apples_to_apples_sim.R 1 1 /private/tmp/lai_smoke 1 1
```

Arguments are:

1. number of replications,
2. study selector (`all`, `1`, `2`, `3`, or comma-separated),
3. output directory,
4. number of cores,
5. maximum number of conditions,
6. chunk index,
7. chunk size,
8. resume existing outputs (`1`/`0`; default is resume).

## Sandwich Coverage Simulation

The command-line entry point is `mlm_random_slope_blup_sandwich_coverage_sim.R`. It uses the shared derivative, TMB, Lai/OpenMx, and stacked-sandwich helper modules in `R/`.

```sh
Rscript mlm_random_slope_blup_sandwich_coverage_sim.R 1 /private/tmp/blup_sandwich_smoke 1 handcoded base screen
```

Arguments are:

1. number of replications,
2. output directory,
3. number of cores,
4. derivative backend (`handcoded`, `numDeriv`, `merDeriv`, `tmb`, or `analytical`),
5. grid mode (`base`, `heteroscedastic`, `heteroscedastic_sparse`, or `expanded`),
6. analysis mode (`screen` or `full`).

`screen` mode skips the full stacked sandwich/OpenMx path and is useful for quick checks. `full` mode includes Lai 2S-PA and stacked sandwich HC0-HC3 estimators.

## Dependencies

The scripts assume an R environment with the packages used by the drivers, including `lme4`, `dplyr`, `tidyr`, `purrr`, `tibble`, `ggplot2`, `foreach`, `doParallel`, `OpenMx`, `sandwich`, `MASS`, `glmnet`, and `data.table`. Optional derivative backends require their corresponding packages (`numDeriv`, `merDeriv`, and `TMB`).

## Tests

The tests are standalone `Rscript` files rather than a `testthat` suite. Run
them from the repository root.

```sh
Rscript tests/test_hc2_hc3_leverage.R
Rscript tests/test_helper_source_order.R
Rscript tests/test_lai_openmx_inputs.R
Rscript tests/test_openmx_lai_wrapper.R
Rscript tests/test_stacked_sandwich_helpers.R
Rscript tests/test_stage2_estimators.R
Rscript lai_replication/tests/smoke_lai_helpers.R
```

Coverage is intentionally focused on shared pipeline pieces:

- `test_hc2_hc3_leverage.R`: leverage scaling used by HC2/HC3 stacked-sandwich variants.
- `test_helper_source_order.R`: root helper source order and availability of key shared functions.
- `test_lai_openmx_inputs.R`: Lai/OpenMx EB measurement input construction, output naming, and subject ordering.
- `test_openmx_lai_wrapper.R`: small synthetic check for the Lai structural-slope OpenMx wrapper.
- `test_stacked_sandwich_helpers.R`: parameter packing, cluster precomputation, likelihood equivalence, corrected scores, and stacked-sandwich covariance output shape.
- `test_stage2_estimators.R`: formatting of stacked-sandwich HC0-HC3 result rows.
- `lai_replication/tests/smoke_lai_helpers.R`: one-replication Study 1 smoke test through the Lai runner.
