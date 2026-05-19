# Lai Replication

This folder contains the Lai and Liu apples-to-apples simulation replication. Study-specific design, data generation, and replication orchestration live here. Shared EB measurement-model construction, OpenMx 2S-PA/2S-PAA wrappers, score recovery, diagnostics, and second-stage estimators live in the repository-level `R/` folder.

## Entry Point

- `mlm_random_slope_lai_apples_to_apples_sim.R`: command-line entry point. It loads packages, locates the repository root, sources shared helpers from `../R/`, sources the Lai modules in this folder, and calls `run_simulation()`.

Small smoke run from the repository root:

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

By default, the runner excludes the tempered EIV regularization path from the
core method set. Set argument 9 to `1` to append
`tempered_eiv_dual_corrected_l25/l50/l75` plus their `_hc0` and `_hc3` EIV
standard-error variants for sensitivity analyses.

Outputs are written under the requested output directory. For each selected condition, the runner writes condition-level replication, summary, issue-summary, and stage-1-diagnostic files under `conditions/`. It also writes aggregate files with the prefix `lai_apples_to_apples_<selection>`; for the full unchunked selection it additionally writes compatibility filenames without the selection suffix.

## Study Modules

- `designs.R`: fixed population parameters, Study 1/2/3 condition grids, and `select_design()`.
- `study1.R`: Study 1 data generation and replication wrapper for the balanced matched-clustering design.
- `study2.R`: Study 2 data generation and replication wrapper for the larger unbalanced matched-clustering design.
- `study3.R`: Study 3 data generation and replication wrapper for the disparate-clustering design where `Y` and `Z` have separate first-stage models.
- `study_common.R`: Lai-only utilities shared by the study modules, including covariance construction, failure rows, truth calculation, and the matched-outcome analysis pipeline used by Studies 1 and 2.

## Runner

- `runner.R`: generic condition chunking, replication loops, result aggregation, progress writing, summaries, and study dispatch.

## SLURM Array Runs

`slurm/lai_condition_array.sbatch` runs the existing Lai command-line entry
point as a SLURM array job. The script defaults to the full current Lai grid
(`all` studies, 582 conditions) with `LAI_CHUNK_SIZE=5`, so the default array
range is 117 tasks:

```sh
sbatch lai_replication/slurm/lai_condition_array.sbatch
```

Common overrides can be supplied as environment variables:

```sh
sbatch \
  --export=ALL,LAI_N_SIM=1000,LAI_STUDY=all,LAI_OUT_DIR=/path/to/lai_outputs,LAI_CHUNK_SIZE=5 \
  lai_replication/slurm/lai_condition_array.sbatch
```

For study-specific runs with `LAI_CHUNK_SIZE=5`, override the array range:

```sh
sbatch --array=1-98  --export=ALL,LAI_STUDY=1 lai_replication/slurm/lai_condition_array.sbatch
sbatch --array=1-10  --export=ALL,LAI_STUDY=2 lai_replication/slurm/lai_condition_array.sbatch
sbatch --array=1-10  --export=ALL,LAI_STUDY=3 lai_replication/slurm/lai_condition_array.sbatch
```

The wrapper accepts these variables:

- `LAI_REPO_ROOT`: repository checkout path; defaults to the wrapper's parent
  repository.
- `LAI_RSCRIPT`: Rscript executable; defaults to `Rscript`.
- `LAI_N_SIM`: replications per condition; defaults to `1000`.
- `LAI_STUDY`: `all`, `1`, `2`, `3`, or comma-separated selectors; defaults
  to `all`.
- `LAI_OUT_DIR`: shared output directory; defaults to
  `outputs/lai_apples_to_apples_slurm` under the repository.
- `LAI_N_CORES`: cores passed to the R runner; defaults to
  `SLURM_CPUS_PER_TASK`.
- `LAI_MAX_CONDITIONS`: optional pre-chunk condition cap; defaults to `NA`.
- `LAI_CHUNK_INDEX`: explicit chunk index; defaults to `SLURM_ARRAY_TASK_ID`.
- `LAI_CHUNK_SIZE`: conditions per array task; defaults to `5`.
- `LAI_RESUME_EXISTING`: skip completed condition outputs; defaults to `1`.
- `LAI_INCLUDE_TEMPERED_EIV`: append tempered EIV sensitivity rows; defaults
  to `0`.

Each array task writes `lai_condition_array_<job>_<chunk>.out`,
`lai_condition_array_<job>_<chunk>.err`, and a combined
`lai_condition_array_<job>_<chunk>.log` under `${LAI_OUT_DIR}/logs`.

After all array jobs finish, rebuild the full aggregate files by running the
standard entry point once with no chunk arguments and resume enabled:

```sh
Rscript lai_replication/mlm_random_slope_lai_apples_to_apples_sim.R \
  1000 all /path/to/lai_outputs 1 NA NA NA 1 0
```

## Shared Code Outside This Folder

These repository-level helpers are intentionally shared with non-Lai simulation scripts:

- `../R/source_helpers.R`: canonical source order for repository-level helpers and a helper loader for script entry points.
- `../R/core_utils.R`: small shared utilities such as `safe_lmer()`, `safe_mean()`, compact diagnostics, and positive-definite matrix projection.
- `../R/blup_helpers.R`: EB/BLUP extraction, Vig-style prior unweighting, and closed-form corrected cluster score calculation.
- `../R/stage2_estimators.R`: observed-score, HC3, ridge, EIV, and stacked-sandwich result-row formatting helpers.
- `../R/sim_diagnostics.R`: first-stage singularity and EB collinearity diagnostics.
- `../R/lai_openmx_helpers.R`: shared EB measurement-model input builders, OpenMx retry/status extraction, and 2S-PA/2S-PAA model wrappers.

## Tests

- `tests/smoke_lai_helpers.R`: one-replication smoke test for the Study 1 path through the full Lai runner.
- `../tests/test_lai_openmx_inputs.R`: unit checks for shared Lai/OpenMx EB measurement input construction.
- `../tests/test_openmx_lai_wrapper.R`: small synthetic check for the shared Lai structural-slope OpenMx wrapper.
- `../tests/test_helper_source_order.R`: smoke check that the repository-level helpers used by the Lai runner source cleanly in script order.

Run them from the repository root:

```sh
Rscript tests/test_lai_openmx_inputs.R
Rscript tests/test_openmx_lai_wrapper.R
Rscript tests/test_helper_source_order.R
Rscript lai_replication/tests/smoke_lai_helpers.R
```
