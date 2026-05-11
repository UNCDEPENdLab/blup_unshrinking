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
8. resume existing outputs (`1`/`0`; default is resume).

Outputs are written under the requested output directory. For each selected condition, the runner writes condition-level replication, summary, issue-summary, and stage-1-diagnostic files under `conditions/`. It also writes aggregate files with the prefix `lai_apples_to_apples_<selection>`; for the full unchunked selection it additionally writes compatibility filenames without the selection suffix.

## Study Modules

- `designs.R`: fixed population parameters, Study 1/2/3 condition grids, and `select_design()`.
- `study1.R`: Study 1 data generation and replication wrapper for the balanced matched-clustering design.
- `study2.R`: Study 2 data generation and replication wrapper for the larger unbalanced matched-clustering design.
- `study3.R`: Study 3 data generation and replication wrapper for the disparate-clustering design where `Y` and `Z` have separate first-stage models.
- `study_common.R`: Lai-only utilities shared by the study modules, including covariance construction, failure rows, truth calculation, and the matched-outcome analysis pipeline used by Studies 1 and 2.

## Runner

- `runner.R`: generic condition chunking, replication loops, result aggregation, progress writing, summaries, and study dispatch.

## Shared Code Outside This Folder

These repository-level helpers are intentionally shared with non-Lai simulation scripts:

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
