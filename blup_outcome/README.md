# BLUP-Outcome Simulation

This folder contains the random-slope BLUP/corrected-score-as-outcome simulation. Study-specific design, data generation, and replication orchestration live here. Shared EB measurement-model construction, corrected-score recovery, diagnostics, and second-stage estimators live in the repository-level `R/` folder.

## Entry Point

- `mlm_random_slope_blup_outcome_sim.R`: command-line entry point. It loads packages, locates the repository root, sources shared helpers from `../R/`, sources the BLUP-outcome modules in this folder, and calls `run_blup_outcome_simulation()`.

Small smoke run from the repository root:

```sh
Rscript blup_outcome/mlm_random_slope_blup_outcome_sim.R 1 /private/tmp/blup_outcome_smoke 1 handcoded smoke screen
```

Arguments are:

1. number of replications,
2. output directory,
3. number of cores,
4. derivative backend (`handcoded`, `numDeriv`, `merDeriv`, `tmb`, or `analytical`),
5. grid mode (`smoke`, `residual_ar1`, `base`, `heteroscedastic`, `heteroscedastic_sparse`, or `expanded`),
6. analysis mode (`screen` or `full`),
7. chunk index,
8. chunk size,
9. resume existing outputs (`1`/`0`; default is resume),
10. execution mode (`run` or `aggregate`; default is `run`),
11. maximum number of conditions to select before chunking.

`screen` mode skips the full Lai/stacked-sandwich path and is useful for quick checks. `full` mode includes Lai 2S-PA/2S-PAA, Fuller variants (including alpha stepdown), and stacked-sandwich HC0-HC3 estimators.

Outputs are written under the requested output directory. For each selected condition, the runner writes condition-level replication, summary, issue-summary, and stage-1-diagnostic files under `conditions/`. Aggregate replication and summary files are written at the top-level output directory after chunk completion (or when `execution_mode` is `aggregate`).

## Modules

- `designs.R`: fixed population parameters, design grids, and `make_blup_outcome_design()`.
- `study_common.R`: BLUP-outcome utilities and the per-replication estimator pipeline.
- `runner.R`: condition chunking, replication loops, result aggregation, progress writing, summaries, and plotting.

## SLURM Array Runs

`slurm/blup_outcome_condition_array.sbatch` runs the BLUP-outcome entry point as a SLURM array job. The script defaults to the base grid with `BLUP_OUTCOME_CHUNK_SIZE=5` and a 173-task array (`base` has 864 conditions), so a typical launch is:

```sh
sbatch blup_outcome/slurm/blup_outcome_condition_array.sbatch
```

Common overrides can be supplied as environment variables:

```sh
sbatch \
  --export=ALL,BLUP_OUTCOME_N_SIM=1000,BLUP_OUTCOME_GRID_MODE=base,BLUP_OUTCOME_OUT_DIR=/path/to/blup_outcome,BLUP_OUTCOME_CHUNK_SIZE=5 \
  blup_outcome/slurm/blup_outcome_condition_array.sbatch
```

The wrapper accepts these variables:

- `BLUP_OUTCOME_REPO_ROOT`: repository checkout path; defaults to the wrapper's parent repository.
- `BLUP_OUTCOME_RSCRIPT`: Rscript executable; defaults to `Rscript`.
- `BLUP_OUTCOME_N_SIM`: replications per condition; defaults to `100`.
- `BLUP_OUTCOME_OUT_DIR`: shared output directory; defaults to `outputs/blup_outcome_slurm` under the repository.
- `BLUP_OUTCOME_N_CORES`: cores passed to the R runner; defaults to `SLURM_CPUS_PER_TASK`.
- `BLUP_OUTCOME_DERIVATIVE_METHOD`: `handcoded`, `numDeriv`, `merDeriv`, `tmb`, or `analytical`; defaults to `handcoded`.
- `BLUP_OUTCOME_GRID_MODE`: `smoke`, `residual_ar1`, `base`, `heteroscedastic`, `heteroscedastic_sparse`, or `expanded`; defaults to `base`.
- `BLUP_OUTCOME_ANALYSIS_MODE`: `screen` or `full`; defaults to `full`.
- `BLUP_OUTCOME_MAX_CONDITIONS`: optional pre-chunk condition cap; defaults to `NA`.
- `BLUP_OUTCOME_CHUNK_INDEX`: explicit chunk index; defaults to `SLURM_ARRAY_TASK_ID`.
- `BLUP_OUTCOME_CHUNK_SIZE`: conditions per array task; defaults to `5`.
- `BLUP_OUTCOME_RESUME_EXISTING`: skip completed condition outputs; defaults to `1`.
- `BLUP_OUTCOME_EXECUTION_MODE`: `run` or `aggregate`; defaults to `run`.

Each array task writes `blup_outcome_condition_array_<job>_<chunk>.out`,
`blup_outcome_condition_array_<job>_<chunk>.err`, and a combined
`blup_outcome_condition_array_<job>_<chunk>.log` under `${BLUP_OUTCOME_OUT_DIR}/logs`.

After all array jobs finish, rebuild the full aggregate files by running the
standard entry point once with no chunk arguments and resume enabled:

```sh
Rscript blup_outcome/mlm_random_slope_blup_outcome_sim.R \
  100 /path/to/blup_outcome 1 handcoded base full NA NA 1 aggregate NA
```

## Tests

Run these from the repository root:

```sh
Rscript tests/test_blup_outcome_designs.R
Rscript tests/test_helper_source_order.R
Rscript tests/test_stacked_sandwich_helpers.R
```
