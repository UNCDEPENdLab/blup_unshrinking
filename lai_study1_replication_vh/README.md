# Lai Study 1 replication with Vig--Hallquist estimators

This directory will hold a fresh replication of Lai and Liu Study 1.  The
replication should preserve the original data-generating mechanism (DGM) while
running the seven-method Vig--Hallquist Study 2 estimator bundle.

## Scale-alignment issue in the original Study 1 code

The original supplementary code does **not** compare every method on one common
standardized scale.  Its `sim1.R` code constructs a method-specific
standardized estimate named `stdbeta`:

- The naive EB regressions multiply their raw path coefficient by the *observed
  EB slope standard deviation*, `sqrt(S[EB slope, EB slope])`.
- 2S-PA and MSEM multiply their raw path coefficient by the estimated *latent
  random-slope* standard deviation.
- The population value used for all bias, RMSE, and coverage calculations is
  `beta_zu1 * sqrt(vr_u1_u0 * icc)`: the raw population path coefficient
  multiplied by the **true** random-slope standard deviation.

Consequently, an EB method's reported standardized estimate can be on a
different scale from its benchmark.  In low-information conditions, EB slopes
are more strongly shrunken and their observed standard deviation is smaller
than the true random-slope standard deviation.  This alone can produce a
negative component in the reported bias, even apart from bias in the raw
second-stage path coefficient.  It is therefore unsafe to interpret the
original `raw_bias` field as raw-coefficient bias: it is the untrimmed Monte
Carlo bias for this method-specifically standardized estimand.  The original
notebook primarily plots the 20% trimmed version (`robust_bias2`).

This convention is part of the original result and should be reproduced in a
strict historical replication.  It should not, however, be the only analysis
in the refreshed study.

## Required estimand views

Each replication should retain the raw estimator output and produce three
explicit reporting views:

1. `original_method_standardized`: reproduce the original supplement's
   method-specific `stdbeta` convention.  This is the view for matching the
   published Study 1 figure.
2. `raw_path`: use scale 1 for every method and compare all estimates with
   `beta_zu1`.  This isolates raw path-coefficient bias.
3. `common_latent_sd_standardized`: multiply every method and its truth by the
   known population `sqrt(var_u1)`.  This makes standardized effects comparable
   without allowing EB shrinkage to change the reporting scale.

The third view is the direct bridge to the current Vig--Hallquist Study 2
workflow, whose estimators are explicitly passed a common latent-slope
reporting scale.  Fuller estimates are raw by default and must only be
rescaled for the third view; they should not be combined with the original
method-specific standardized rows without an explicit conversion.

## Interpretation implications

- A negative bias pattern for naive EB OLS in low-reliability cells may partly
  reflect the original standardization convention, rather than only a change
  in the raw OLS coefficient.
- Bias that remains in `raw_path` is evidence about the estimator itself.
- A difference between `raw_path` and `common_latent_sd_standardized` is only
  a fixed, condition-level unit conversion.  A difference involving
  `original_method_standardized` can additionally reflect EB shrinkage.
- Replication figures should report the exact retention rule used.  The
  original notebook computes parameter summaries among converged replications;
  modern Vig--Hallquist eligibility and conditioning screens should be supplied
  as a separate diagnostic analysis, not silently substituted into the
  historical figure.

## Source references

The original repository was inspected at commit
`9ffe53168f6bb04e13ef977dc19a8d953d0bf29d` (2026-01-16).

- DGM and Study 1 condition grid: [`simulation_scripts/sim1.R`, lines 29--72](https://github.com/marklhc/2spa-random-slopes-supp/blob/9ffe53168f6bb04e13ef977dc19a8d953d0bf29d/simulation_scripts/sim1.R#L29-L72).
- Naive EB `stdbeta` definitions: [`sim1.R`, lines 135--183](https://github.com/marklhc/2spa-random-slopes-supp/blob/9ffe53168f6bb04e13ef977dc19a8d953d0bf29d/simulation_scripts/sim1.R#L135-L183).
- 2S-PA and MSEM `stdbeta` definitions: [`sim1.R`, lines 185--280](https://github.com/marklhc/2spa-random-slopes-supp/blob/9ffe53168f6bb04e13ef977dc19a8d953d0bf29d/simulation_scripts/sim1.R#L185-L280).
- Population target and the `raw_bias` / `robust_bias` calculations: [`sim1.R`, lines 289--325](https://github.com/marklhc/2spa-random-slopes-supp/blob/9ffe53168f6bb04e13ef977dc19a8d953d0bf29d/simulation_scripts/sim1.R#L289-L325).
- Historical Figure 2 bias, Figure 3 coverage, and Figure 4 Type I-error
  calculations: [`notebooks/_sim1_results_supp.qmd`](https://github.com/marklhc/2spa-random-slopes-supp/blob/9ffe53168f6bb04e13ef977dc19a8d953d0bf29d/notebooks/_sim1_results_supp.qmd).

Relevant local code:

- [`lai_replication/study1.R`](/proj/mnhallqlab/projects/blup_unshrinking/lai_replication/study1.R) preserves the Lai Study 1 DGM.
- [`lai_replication/study_common.R`](/proj/mnhallqlab/projects/blup_unshrinking/lai_replication/study_common.R) is the current local Stage 1/Stage 2 path.
- [`R/stage2_estimators.R`](/proj/mnhallqlab/projects/blup_unshrinking/R/stage2_estimators.R) contains the current OLS and Fuller estimators.

## Current implementation

The refresh code in this directory implements the original Study 1 DGM and
fits the seven-method VH Study 2 bundle on the **raw** coefficient scale:

1. oracle dual regression;
2. naive dual-BLUP regression;
3. closed-form-score dual regression;
4. Fuller closed-form correction;
5. Fuller alpha-stepdown closed-form correction;
6. Lai 2S-PA; and
7. MSEM.

It then creates `raw` and `latent_sd` reporting rows from the same fitted
result:

```sh
Rscript lai_study1_replication_vh/run.R \
  1000 outputs/lai_study1_vh 1 NA 1
```

Arguments are, in order: replications per condition, output directory, cores,
an optional condition cap, resume flag, and an optional comma-separated method
list.  For a no-Mplus smoke run, for example:

```sh
Rscript lai_study1_replication_vh/run.R \
  2 /tmp/lai_study1_vh_smoke 1 1 0 \
  oracle_dual,naive_dual_blup,closed_form_dual,fuller_closed_form,lai_2spa
```

The MSEM helper is included in the default bundle but can be excluded from
smoke runs because it requires Mplus.

In addition, `lai_original_standardized` is a historical-comparison view for
the three methods that correspond directly to the original Study 1 bundle:
`naive_dual_blup`, `lai_2spa`, and `msem`.  It uses the original method-specific
multiplier while retaining the common population benchmark
`beta_zu1 * sqrt(var_u1)`.  The output also records
`lai_original_naive_eb_slope_sd` (the raw-data-ML EB-slope SD used by the
original naïve model), its `J - 1` sample-SD counterpart, each method's
historical multiplier, and its source.  Oracle and Fuller methods do not have
an original-Lai counterpart, so they deliberately have no historical-view row.

## Eligibility conventions

Each result row retains two eligibility decisions and their reasons:

- `vh_analysis_eligible` / `vh_analysis_exclusion_reason` is the primary
  analysis rule and matches VH Study 2: valid estimate and SE, OLS design
  checks, 2S-PA OpenMx diagnostics, and MSEM warning/target-parameter checks.
  The legacy public fields `analysis_eligible` and
  `analysis_exclusion_reason` are aliases for this primary rule.
- `lai_original_eligible` / `lai_original_exclusion_reason` implements the
  original Study 1 summary rule exactly: retain a method replication when its
  recorded optimizer status `code == 0`.  No extra numerical or inferential
  screening is imposed.

Condition summaries report counts and rates for each rule, their overlap, and
each rule's exclusive set.  This permits a direct post-estimation comparison;
the primary bias and coverage summaries continue to use the VH rule.
`eligibility_comparison` labels individual rows as `both`,
`lai_original_only`, `vh_only`, or `neither`.

## Post-estimation Figure 2--4 analogue outputs

After the model runs complete, generate the historical-style analogue and its
VH-primary companion without rerunning any models:

```sh
Rscript lai_study1_replication_vh/make_figure2_analogue.R \
  outputs/lai_study1_vh outputs/lai_study1_vh/figure2_analogue
```

The command reads each condition replication file one at a time and writes:

- `figure2_analogue_cell_summary.csv`, plus PNG/PDF plots: only naïve dual
  BLUP, 2S-PA, and MSEM; `lai_original_standardized`; `status_code == 0`;
  20%-trimmed bias. This is labeled a Figure 2 **analogue**, not a full
  reproduction, because the original `u1`-only EB and 2S-PAB methods are not
  in the VH seven-method bundle.
- `vh_primary_companion_cell_summary.csv`, plus PNG/PDF plots: all seven VH
  methods; common `latent_sd` reporting scale; VH primary eligibility; mean
  bias.
- `figure3_4_analogue_cell_summary.csv`, plus `figure3_coverage_analogue` and
  `figure4_type1_analogue` PNG/PDF plots: the same three overlapping methods,
  `lai_original_standardized` scale, and `status_code == 0` retention.  Figure
  3 is empirical coverage of the stored 95% confidence interval.  Figure 4 is
  limited to `beta_zu1 == 0` and defines rejection as exclusion of the truth
  from that interval.  Lai's original Figure 4 uses `p < .05`; the CI rule is
  equivalent for ordinary two-sided Wald intervals but is explicitly an
  analogue when an estimator's interval and p-value constructions differ.

### Cached figure rebuilds

The historical Figure 2 analogue needs replication-level estimates because its
20%-trimmed bias cannot be recovered from the aggregate summary.  Figures 3--4
also need replication-level intervals.  On their first run, the postprocessor
reduces each condition to small historical summaries and writes
signature-validated caches in
`figure2_analogue/historical_condition_cache/`, with a
`figure2_analogue_cache_manifest.csv` recording the source result-file size and
modification time.  Figure 2 summaries and Figure 3--4 inference summaries
are cached separately, so later runs reuse both caches and only redraw plots.

The all-seven-method VH companion is now derived directly from
`lai_study1_vh_summary.csv`; it does not reread the replication files.
Existing pre-cache `figure2_analogue_cell_summary.csv` output is promoted to
the cache when its expected condition-method rows are complete, avoiding an
unnecessary second full pass after an earlier successful figure job.

Force a full historical cache rebuild only after intentionally changing the
figure calculation itself:

```sh
Rscript lai_study1_replication_vh/make_figure2_analogue.R \
  outputs/lai_study1_vh outputs/lai_study1_vh/figure2_analogue true
```

Changed or replaced condition result files invalidate only their corresponding
historical cache entries; the combined historical summaries and plots are then
rebuilt from the cached and refreshed cells.

For the first Figures 3--4 build after a production run, submit the dedicated
post-estimation job (it makes one interval-only pass through the condition
files and then caches the results):

```sh
sbatch --export=ALL,LAI_VH_OUT_DIR=outputs/lai_study1_vh_production \
  lai_study1_replication_vh/slurm/lai_study1_vh_figures.sbatch
```

## Production Slurm execution

The production launcher runs one DGM condition per array task, keeps each task
single-core for Mplus safety, and validates existing condition files before
skipping them.  After a one-condition MSEM-inclusive preflight, submit the full
Lai grid (486 cells, 5,000 replications each) with:

```sh
sbatch --array=1-486%10 \
  --export=ALL,LAI_VH_N_SIM=5000,LAI_VH_OUT_DIR=outputs/lai_study1_vh_production \
  lai_study1_replication_vh/slurm/lai_study1_vh_condition_array.sbatch
```

The `%10` throttle controls simultaneous Mplus jobs; reduce it if the license
or cluster policy requires fewer concurrent processes.  Once every array task
finishes, rebuild summaries without rerunning any complete condition:

Submit from the repository root (as above), or set `LAI_VH_REPO_ROOT` to that
absolute path when submitting from elsewhere.

```sh
sbatch --export=ALL,LAI_VH_N_SIM=5000,LAI_VH_OUT_DIR=outputs/lai_study1_vh_production \
  lai_study1_replication_vh/slurm/rebuild_lai_study1_vh_summary.sbatch
```

The rebuild job fails and identifies any incomplete condition rather than
silently producing a partial aggregate summary.  To test one array-style cell
before the full submission, override the array and replication count:

```sh
sbatch --array=1-1%1 \
  --export=ALL,LAI_VH_N_SIM=10,LAI_VH_OUT_DIR=outputs/lai_study1_vh_preflight \
  lai_study1_replication_vh/slurm/lai_study1_vh_condition_array.sbatch
```

## Slope-SD and reliability diagnostics

Each condition records `dgm_population_slope_sd` and
`dgm_posterior_slope_reliability`.  The latter is the true posterior
reliability, computed using the full two-dimensional random-effect covariance,
the exact Lai RMS-scaled within-cluster design, and `sigma2 I` residual
covariance.  It is an emergent DGM property, not a reliability-calibration
input.

Each replication also records the realized true, EB, corrected-score, and
fitted slope SDs; the fitted residual SD; mean posterior slope variance;
`lambda22`; and fitted posterior slope reliability.  The SD ratios are all
relative to the population slope SD.  `lambda22` is deliberately reported as
a loading diagnostic only: it is not treated as slope reliability in cells
with correlated random effects.  Summary files average these replication-level
diagnostics within each condition, method, and reporting scale.
