# MLM BLUP Simulations

This repository contains reusable R helpers, simulation drivers, and technical notes for work on empirical-Bayes/BLUP score correction in multilevel models. The code is script-oriented rather than packaged as an installed R package: simulation drivers source helper files from the repository-level `R/` directory.

## Repository Layout

- `R/`: shared helper modules used by multiple simulations.
- `blup_outcome/`: focused BLUP-as-stage-2-outcome simulation modules.
- `lai_replication/`: Lai and Liu apples-to-apples replication modules and runner.
- `documentation/`: technical notes, manuscripts, references, and rendered PDFs.
- `archive/`: removed legacy or experimental helpers retained for reference.
- `src/`: TMB C++ source used by the TMB derivative backend.
- `tests/`: lightweight checks for shared calculations.

Simulation and analysis scripts include:

- `archive/mlm_blup_correction_sim.R`: archived scalar random-intercept BLUP correction demo.
- `archive/mlm_random_slope_blup_correction_sim.R`: archived random-slope full-matrix correction demo.
- `mlm_random_slope_blup_predictor_comparison_sim.R`: predictor comparison simulation.
- `blup_outcome/mlm_random_slope_blup_outcome_sim.R`: random-slope BLUP/corrected-score-as-outcome simulation.
- `plot_blup_correction_summaries.R`, `plot_zoomed.R`, and `anova_analysis.R`: plotting and summary analyses for generated outputs.

## Shared Helpers

- `R/core_utils.R`: small shared utilities such as NA-stable means, compact diagnostics, quiet `lmer()` fitting, and matrix regularization.
- `R/source_helpers.R`: canonical source order for repository-level helpers and a helper loader for script entry points.
- `R/simulation_runner_helpers.R`: shared chunking, progress-file, and atomic CSV helpers for condition-grid simulations.
- `R/blup_helpers.R`: EB/BLUP extraction, Vig-style prior unweighting, and closed-form corrected cluster scores.
- `R/sim_helpers.R`: common simulation data-generation helpers.
- `R/reliability_calibration.R`: posterior-reliability calibration,
  standardized-beta conversion, and deterministic reference cluster-size
  profiles.
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
8. resume existing outputs (`1`/`0`; default is resume),
9. include tempered EIV sensitivity rows (`1`/`0`; default is core methods only).

The core Lai simulation excludes the tempered EIV regularization path. Passing
argument 9 as `1` appends `tempered_eiv_dual_corrected_l25/l50/l75` plus their
`_hc0` and `_hc3` EIV standard-error variants. These use lambda-weighted
covariance subtraction as sensitivity checks rather than as the classical full
EIV correction.

For SLURM clusters, `lai_replication/slurm/lai_condition_array.sbatch` wraps
the same entry point as a job array over condition chunks. It defaults to
`LAI_CHUNK_SIZE=5` and `#SBATCH --array=1-117`, which covers the current
582-condition all-studies Lai grid. See `lai_replication/README.md` for the
override variables and study-specific array ranges.

## BLUP-Outcome Simulation

The command-line entry point is `blup_outcome/mlm_random_slope_blup_outcome_sim.R`.

```sh
Rscript blup_outcome/mlm_random_slope_blup_outcome_sim.R 1 /private/tmp/blup_outcome_smoke 1 handcoded smoke screen
```

Arguments are:

1. number of replications,
2. output directory,
3. number of cores,
4. derivative backend (`handcoded`, `numDeriv`, `merDeriv`, `tmb`, or `analytical`),
5. grid mode (`posterior_reliability_smoke`, `posterior_reliability`, `smoke`,
   `base`, `heteroscedastic`, `heteroscedastic_sparse`, or `expanded`),
6. analysis mode (`screen` or `full`),
7. chunk index,
8. chunk size,
9. resume existing outputs (`1`/`0`; default is resume),
10. execution mode (`run` or `aggregate`; default is `run`),
11. maximum number of conditions to select before chunking.

`screen` mode skips the full stacked sandwich/OpenMx path and is useful for quick checks. `full` mode includes Lai 2S-PA/2S-PAA and stacked sandwich HC0-HC3 estimators.

The BLUP-outcome grids now explicitly test:

- raw EB/BLUPs as stage-2 outcomes alongside diagonal, full-matrix, and closed-form corrected outcomes;
- null, small, and moderate `x -> random slope` effects;
- negative, zero, and positive random intercept/slope correlations;
- balanced, randomly unbalanced, and informative cluster-size imbalance;
- sparse through high-trial first-stage regimes;
- fully crossed random-slope variance (`tau1`) and level-1 noise (`sigma`).

The opt-in `posterior_reliability` grid uses the amended design directly:

- posterior slope reliability: `.25`, `.50`, `.80`;
- standardized beta: `0`, `.20`, `.40`, `.60`;
- marginal intercept-slope correlation: `-.50`, `0`, `.50`;
- cluster size: `3`, `5`, `10`, `25`;
- number of clusters: `30`, `50`, `100`, `150`, `300`.

These conditions are calibrated once when the manifest is built. Each row
stores both `marginal_rho`, which defines the target first-stage G matrix, and
`rho_residual`, which is passed to `simulate_dataset()` when drawing residual
random effects. The calibrated `gamma_x_on_slope`, `tau1_residual`, and
`rho_residual` remain fixed across replications. Structural R-squared is stored
only as the derived square of the standardized-beta target.

Chunked runs mirror the Lai runner model. In `run` mode, the script writes only
condition-level files under `conditions/` plus a chunk-specific manifest and
progress file. It does not rewrite aggregate outputs, which makes it safe for
SLURM array jobs that share an output directory. If `chunk_index` is omitted,
the script uses `SLURM_ARRAY_TASK_ID` when available; if `n_cores` is omitted,
it uses `SLURM_CPUS_PER_TASK` when available. `chunk_size` can also be supplied
through `BLUP_OUTCOME_CHUNK_SIZE`.

Example chunk run:

```sh
Rscript blup_outcome/mlm_random_slope_blup_outcome_sim.R 100 /path/to/out 1 handcoded base full 3 5 1 run
```

After all chunks finish, run aggregation separately:

```sh
Rscript blup_outcome/mlm_random_slope_blup_outcome_sim.R 100 /path/to/out 1 handcoded base full NA NA 1 aggregate
```

## Dependencies

The scripts assume an R environment with the packages used by the drivers, including `lme4`, `dplyr`, `tidyr`, `purrr`, `tibble`, `ggplot2`, `foreach`, `doParallel`, `OpenMx`, `sandwich`, `MASS`, `glmnet`, and `data.table`. Optional derivative backends require their corresponding packages (`numDeriv`, `merDeriv`, and `TMB`).

## Tests

The tests are standalone `Rscript` files rather than a `testthat` suite. Run
them from the repository root.

```sh
Rscript tests/test_hc2_hc3_leverage.R
Rscript tests/test_helper_source_order.R
Rscript tests/test_blup_outcome_designs.R
Rscript tests/test_lai_openmx_inputs.R
Rscript tests/test_openmx_lai_wrapper.R
Rscript tests/test_stacked_sandwich_helpers.R
Rscript tests/test_stage2_estimators.R
Rscript lai_replication/tests/smoke_lai_helpers.R
```

Coverage is intentionally focused on shared pipeline pieces:

- `test_hc2_hc3_leverage.R`: leverage scaling used by HC2/HC3 stacked-sandwich variants.
- `test_helper_source_order.R`: root helper source order and availability of key shared functions.
- `test_blup_outcome_designs.R`: BLUP-outcome condition-grid coverage for null effects, sparse trials, balance modes, signed correlations, and fully crossed `tau1`/`sigma`.
- `test_lai_openmx_inputs.R`: Lai/OpenMx EB measurement input construction, output naming, and subject ordering.
- `test_openmx_lai_wrapper.R`: small synthetic check for the Lai structural-slope OpenMx wrapper.
- `test_stacked_sandwich_helpers.R`: parameter packing, cluster precomputation, likelihood equivalence, corrected scores, and stacked-sandwich covariance output shape.
- `test_stage2_estimators.R`: formatting of stacked-sandwich HC0-HC3 result rows.
- `lai_replication/tests/smoke_lai_helpers.R`: one-replication Study 1 smoke test through the Lai runner.
