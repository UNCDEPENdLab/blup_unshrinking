# Standardized-Beta Pipeline Changes for VH Studies 1-3

Date: June 15, 2026

## Purpose

This document records the mathematical and code changes made to use a
standardized regression coefficient as the primary structural effect-size
metric throughout the Vig-Hallquist Studies 1-3 pipeline.

The earlier refactor beginning around commit
`66bdd6242e0eb44ba149e3e34c590001b57b2ffc` made structural \(R^2\) the
manipulated effect-size quantity in parts of the reliability-calibration
pathway. That parameterization is straightforward when there is one structural
predictor, but it does not isolate a focal random-slope coefficient when a
random intercept is also included with a nonzero coefficient. The revised
pipeline therefore uses standardized beta as the manipulated quantity and
retains \(R^2\) only as a derived diagnostic.

## Common Estimand Convention

The design manifests now use `standardized_beta_target` as the manipulated
effect-size field.

The general definition is

$$
\beta_{\mathrm{std}}
=
\beta_{\mathrm{raw}}
\frac{\operatorname{SD}(X)}
     {\operatorname{SD}(Y)}.
$$

All estimators within a study are reported using the same calibrated
population standard deviations. They are not standardized using the observed
SD of each proxy or a method-specific estimate of a latent variance. This keeps
the estimand fixed across oracle, BLUP, corrected-score, Fuller, Lai, and direct
model estimators.

### Why calibrated reporting scales are passed through the pipeline

The simulation target is a fixed population standardized coefficient, not a
coefficient standardized separately by each estimator. Once a condition has
been calibrated, the population standard deviations that define its estimand
must therefore be passed to every Stage 2 estimator and any other reporting
path that transforms a raw structural coefficient.

If each method instead used the observed SD of its own predictor or a fitted
latent SD, method \(m\) would report

$$
\widehat\beta_{m,\mathrm{std}}
=
\widehat\beta_{m,\mathrm{raw}}
\frac{\widehat{\operatorname{SD}}_m(X)}
     {\widehat{\operatorname{SD}}_m(Y)}.
$$

This scale is random and method-specific. BLUPs generally have smaller
variances because of shrinkage, corrected scores generally have larger
variances, and Fuller and Lai models may estimate their own latent variances.
The resulting estimates would consequently target different finite-sample
quantities. Comparisons among methods would combine accuracy of the structural
coefficient with differences in proxy variance, latent-variance estimation,
and replication-specific sampling variation.

This issue is particularly important for naive BLUP regression. Typically,

$$
\operatorname{SD}(\widehat b_1)
<
\operatorname{SD}(b_1).
$$

Multiplying the naive coefficient by the shrunken BLUP SD does not place it on
the calibrated latent population scale. It can partially conceal or otherwise
alter the apparent attenuation that the simulation is intended to measure.

The fixed reporting multipliers are:

1. Study 1:

   $$
   c_1=\frac{1}{\sqrt{G_{22}}}.
   $$

2. Study 2:

   $$
   c_2=\frac{\sqrt{G_{22}}}{\sigma_z}.
   $$

3. Study 3:

   $$
   c_3=
   \sqrt{\frac{G_{y,22}}{G_{q,22}}}.
   $$

Each method first estimates its raw structural coefficient and then applies
the appropriate fixed multiplier. The same linear transformation must be
applied to the estimate, standard error, and confidence limits:

$$
\widehat\beta_{\mathrm{std}}=c\widehat\beta_{\mathrm{raw}},
\qquad
\operatorname{SE}(\widehat\beta_{\mathrm{std}})
=
c\operatorname{SE}(\widehat\beta_{\mathrm{raw}}).
$$

Using an estimated latent SD can be appropriate when the scientific estimand is
an estimator-specific or sample-standardized coefficient. It may also converge
to the same population quantity when the variance estimator is consistent.
That is not the estimand used in these simulations. Here, passing the calibrated
scale through the pipeline isolates the intended comparison: how accurately
each method estimates the same population structural effect. It also makes
bias, RMSE, standard-error calibration, and confidence-interval coverage
directly comparable across methods.

Structural \(R^2\) remains in manifests and summaries where useful, but it is
computed after beta calibration. It is not used to choose the focal effect.

## Study 1: Random Slope as Outcome

### Model

The random slope is

$$
b_{1i} = \gamma x_i + u_{1i},
\qquad
\operatorname{Var}(x_i)=1.
$$

Let the reliability-calibrated marginal random-effect covariance be

$$
\mathbf G =
\begin{bmatrix}
G_{11} & G_{12}\\
G_{12} & G_{22}
\end{bmatrix}.
$$

The target effect is

$$
\beta_{\mathrm{std}}
=
\frac{\gamma}{\sqrt{G_{22}}},
$$

so the raw coefficient is

$$
\gamma
=
\beta_{\mathrm{std}}\sqrt{G_{22}}.
$$

The derived structural variance fraction is

$$
R^2_{\mathrm{struct}}
=
\frac{\gamma^2}{G_{22}}
=
\beta_{\mathrm{std}}^2.
$$

### Residual random-effect covariance

Because \(G_{22}\) is the marginal variance of the total slope, the residual
slope variance must be

$$
\operatorname{Var}(u_{1i})
=
G_{22}-\gamma^2
=
G_{22}(1-\beta_{\mathrm{std}}^2).
$$

The marginal intercept-slope covariance is preserved:

$$
\operatorname{Cov}(b_{0i},u_{1i})=G_{12}.
$$

Therefore the covariance used to draw the residual random effects is

$$
\mathbf G_{\mathrm{residual}}
=
\begin{bmatrix}
G_{11} & G_{12}\\
G_{12} & G_{22}-\gamma^2
\end{bmatrix}.
$$

The residual correlation passed to the simulator is

$$
\rho_{\mathrm{residual}}
=
\frac{G_{12}}
{\sqrt{G_{11}(G_{22}-\gamma^2)}}.
$$

It is generally not equal to the requested marginal correlation. Positive
definiteness requires

$$
\beta_{\mathrm{std}}^2 < 1-\rho_{\mathrm{marginal}}^2.
$$

### Reporting scale

Raw estimates of \(\gamma\), including their SEs and confidence limits, are
multiplied by

$$
\frac{1}{\sqrt{G_{22}}}
$$

before comparison with `standardized_beta_target`.

### Implemented estimator set

The executable Study 1 roster contains eight method rows:

- Benchmark: `oracle`.
- Score regressions: `naive_blup`, `closed_form`, and
  `single_subject_ols`.
- Measurement-error corrections: `fuller_closed_form` and
  `fuller_alpha_stepdown_closed_form`.
- Latent-score and joint-model estimators: `lai_2spa` and `direct_mlm`.

## Study 2: Random Intercept and Slope as Predictors

### Model

The external outcome is

$$
z_i
=
\beta_0+\beta_1 b_{0i}+\beta_2 b_{1i}+e_i,
\qquad
\operatorname{Var}(z_i)=\sigma_z^2.
$$

The focal standardized slope coefficient is

$$
\beta_{2,\mathrm{std}}
=
\beta_2\frac{\sqrt{G_{22}}}{\sigma_z}.
$$

Thus

$$
\beta_2
=
\beta_{2,\mathrm{std}}
\frac{\sigma_z}{\sqrt{G_{22}}}.
$$

For `slope_only`, \(\beta_1=0\). For `intercept_slope`, \(\beta_1\) is the
fixed nuisance intercept coefficient.

The structural predictor variance is

$$
V_{\mathrm{pred}}
=
\beta_1^2G_{11}
+2\beta_1\beta_2G_{12}
+\beta_2^2G_{22},
$$

and the outcome residual variance is

$$
\sigma_e^2=\sigma_z^2-V_{\mathrm{pred}}.
$$

The total structural \(R^2\) is derived as

$$
R^2_{\mathrm{total}}
=
\frac{V_{\mathrm{pred}}}{\sigma_z^2}.
$$

For the intercept-plus-slope model, the focal slope's unique variance
contribution is also recorded:

$$
R^2_{\mathrm{focal,unique}}
=
\frac{
  \beta_2^2
  \operatorname{Var}(b_{1i}\mid b_{0i})
}{
  \sigma_z^2
},
$$

where

$$
\operatorname{Var}(b_{1i}\mid b_{0i})
=
G_{22}-\frac{G_{12}^2}{G_{11}}.
$$

### Why total \(R^2\) was removed as the target

When \(\beta_1\ne0\), total \(R^2\) combines:

1. the nuisance random-intercept effect;
2. the focal random-slope effect; and
3. the covariance contribution between the two predictors.

Consequently, a fixed nonzero \(\beta_1\) creates a positive lower bound on
attainable total \(R^2\). Null and small focal slope effects can then be
infeasible even though \(\beta_2=0\) is a valid and scientifically important
condition.

The previous quadratic `calibrate_beta2` proposal is mathematically valid if
the scientific target is total model \(R^2\). It is not appropriate for this
study's focal-effect grid because changing \(G_{12}\), reliability, or the
nuisance intercept effect changes the resulting focal \(\beta_2\).

### Reporting scale

Raw Study 2 slope coefficients are multiplied by

$$
\frac{\sqrt{G_{22}}}{\sigma_z}.
$$

The current design fixes \(\sigma_z^2=1\), so this reduces to the calibrated
latent slope SD.

Primary Study 2 models always include both the random-intercept and
random-slope predictors. Slope-only fitted models remain as diagnostics for
omitted-intercept contamination.

### Implemented estimator set

The executable Study 2 roster contains 15 method rows:

- Benchmark: `oracle_dual`.
- Primary dual-predictor estimators: `naive_dual_blup`,
  `naive_dual_blup_hc3`, `closed_form_dual`,
  `closed_form_dual_hc3`, `fuller_closed_form`,
  `fuller_alpha_stepdown_closed_form`, `lai_2spa`, and `lai_2spaa`.
- Slope-only diagnostics: `naive_slope_blup`,
  `naive_slope_blup_hc3`, `centered_slope_blup`,
  `centered_slope_blup_hc3`, `closed_form_slope`, and
  `closed_form_slope_hc3`.

The `_hc3` rows use the same point estimate as their corresponding OLS row and
change only the estimated sampling variance and confidence interval. The
`method_role` field distinguishes `benchmark`, `primary_dual`, and
`slope_only_diagnostic` rows so slope-only omitted-variable diagnostics are not
pooled with the primary estimators.

The former generic names `oracle` and `closed_form` were renamed
`oracle_dual` and `closed_form_dual` to make the fitted predictor set explicit.
`single_subject_ols` is not a separate Study 2 method because the full
closed-form score is the cluster-specific OLS score under iid residuals and
its GLS analogue under non-iid residuals. `msem` was removed from the roster
because no MSEM estimator was implemented; retaining the label would have
created a nominal method with no executable definition.

## Study 3: Random Effects as Predictors and Outcome

### Model

Study 3 contains two random-intercept/random-slope processes:

$$
\mathbf b_{yi}
=
\begin{bmatrix}
b_{y0i}\\
b_{y1i}
\end{bmatrix}
\sim N(\mathbf 0,\mathbf G_y)
$$

and

$$
\mathbf b_{qi}
=
\begin{bmatrix}
b_{q0i}\\
b_{q1i}
\end{bmatrix}
\sim N(\mathbf 0,\mathbf G_q).
$$

The Q slope is generated from the Y random effects:

$$
b_{q1i}
=
\theta_0b_{y0i}
+\theta_1b_{y1i}
+v_i.
$$

The focal target is the standardized \(\theta_1\):

$$
\theta_{1,\mathrm{std}}
=
\theta_1
\frac{\operatorname{SD}(b_{y1i})}
     {\operatorname{SD}(b_{q1i})}.
$$

Therefore

$$
\theta_1
=
\theta_{1,\mathrm{std}}
\frac{\operatorname{SD}(b_{q1i})}
     {\operatorname{SD}(b_{y1i})}.
$$

For the intercept-plus-slope condition, the nuisance coefficient is separately
defined on the same standardized scale:

$$
\theta_0
=
\theta_{0,\mathrm{std}}
\frac{\operatorname{SD}(b_{q1i})}
     {\operatorname{SD}(b_{y0i})}.
$$

For the slope-only condition, \(\theta_0=0\). A null focal condition always
sets \(\theta_1=0\), including when \(\theta_0\ne0\).

### Reliability and covariance calibration

\(\mathbf G_y\) and \(\mathbf G_q\) are separately calibrated from the Y and Q
posterior slope-reliability targets before either structural coefficient is
computed.

Let

$$
\mathbf c=
\begin{bmatrix}
\theta_0\\
\theta_1
\end{bmatrix}.
$$

The variance explained in the Q slope is

$$
V_{\mathrm{struct}}
=
\mathbf c^\top\mathbf G_y\mathbf c.
$$

The residual Q-slope variance is

$$
\operatorname{Var}(v_i)
=
G_{q,22}-V_{\mathrm{struct}}.
$$

Assuming \(b_{q0i}\) is independent of the Y random effects, preservation of
the target marginal Q intercept-slope covariance requires

$$
\operatorname{Cov}(b_{q0i},v_i)=G_{q,12}.
$$

Thus the residual Q block is

$$
\mathbf G_{q,\mathrm{residual}}
=
\begin{bmatrix}
G_{q,11} & G_{q,12}\\
G_{q,12} & G_{q,22}-V_{\mathrm{struct}}
\end{bmatrix}.
$$

The implied marginal covariance for
\((b_{y0i},b_{y1i},b_{q0i},b_{q1i})^\top\) is

$$
\mathbf G_{\mathrm{joint}}
=
\begin{bmatrix}
\mathbf G_y & \mathbf 0 & \mathbf G_y\mathbf c\\
\mathbf 0^\top & G_{q,11} & G_{q,12}\\
\mathbf c^\top\mathbf G_y & G_{q,12} & G_{q,22}
\end{bmatrix}.
$$

Both \(\mathbf G_{q,\mathrm{residual}}\) and
\(\mathbf G_{\mathrm{joint}}\) are checked for positive definiteness.

The derived total structural \(R^2\) is

$$
R^2_{\mathrm{total}}
=
\frac{\mathbf c^\top\mathbf G_y\mathbf c}{G_{q,22}}.
$$

When both Y predictors are included, the focal unique contribution is

$$
R^2_{\mathrm{focal,unique}}
=
\frac{
  \theta_1^2
  \operatorname{Var}(b_{y1i}\mid b_{y0i})
}{
  G_{q,22}
}.
$$

This derived unique \(R^2\) is not the manipulated effect. The manipulated
effect remains \(\theta_{1,\mathrm{std}}\).

### Reporting scale

Raw \(\theta_1\) estimates are multiplied by

$$
\frac{\operatorname{SD}(b_{y1i})}
     {\operatorname{SD}(b_{q1i})}
$$

for all Study 3 estimators.

### Implemented estimator set

The executable Study 3 roster contains 13 method rows:

- Benchmark: `oracle_dual`.
- OLS combinations, each with naive and HC3 inference:
  `naive_blup_on_blup`, `closed_form_on_blup`,
  `blup_on_closed_form`, and `closed_form_on_closed_form`.
- Measurement-error corrections: `fuller_closed_form` and
  `fuller_alpha_stepdown_closed_form`.
- Latent-score models: `lai_2spa` and `lai_2spaa`.

For names of the form `outcome_on_predictor`, the first component identifies
the Q-slope outcome score and the second identifies the Y-process predictor
scores. Thus `closed_form_on_blup` uses a closed-form Q outcome with Y BLUP
predictors, whereas `blup_on_closed_form` uses a Q BLUP outcome with
closed-form Y predictors.

## Major Code Changes

### Shared reliability calibration

File: `R/reliability_calibration.R`

- Replaced `structural_r_squared` and `effect_sign` inputs in
  `decompose_structural_slope()` with a signed `standardized_beta`.
- Simplified `calibrate_random_slope_condition()` to calibrate Study 1-style
  random-slope outcomes only. The old study-specific `"w"`/`"z"` branching was
  removed.
- Added `calibrate_blup_predictor_effect()` and
  `calibrate_blup_predictor_condition()` for Study 2.
- Added `calibrate_dual_process_effect()` for Study 3.
- Added explicit residual-variance and positive-definiteness checks.
- Kept total and focal unique \(R^2\) values as derived outputs.

### Design manifests

Files:

- `blup_outcome/designs.R`
- `vig_hallquist_2026/vh_designs.R`

Changes:

- Replaced manipulated grids of `structural_r2 = c(0, .04, .16, .36)` with
  `standardized_beta_target = c(0, .2, .4, .6)` for Studies 1 and 2.
- Added Study 3 targets `c(0, .2, .5)`.
- Calibrated first-stage reliability before structural coefficients in Studies
  2 and 3.
- Stored raw coefficients, population scales, residual variances, residual
  correlations, and derived \(R^2\) values in each condition.
- Implemented canonical designs containing 720 Study 1 conditions, 1,440
  Study 2 conditions, and 1,728 Study 3 conditions.
- Assigned stable global condition IDs independent of the selected study:
  Study 1 uses 1-720, Study 2 uses 721-2160, and Study 3 uses 2161-3888.
- Retained `.25`, `.50`, and `.80` as the Study 1/2 reliability levels.
- Used `.25` and `.80` for Study 3. A `.20` target is not attainable in some
  sparse, correlated conditions under the current posterior-reliability
  definition.

### Study runners and simulation

Files:

- `blup_outcome/study_common.R`
- `vig_hallquist_2026/vh_study_common.R`
- `vig_hallquist_2026/vh_study1.R`
- `vig_hallquist_2026/vh_study2.R`
- `vig_hallquist_2026/vh_study3.R`
- `vig_hallquist_2026/vh_runner.R`

Changes:

- Converted truth, estimates, SEs, and confidence limits to the common
  population standardized-beta scale.
- Completed the Study 1 estimator runner.
- Completed the Study 2 dual-predictor runner and classified primary
  dual-predictor methods separately from slope-only diagnostics.
- Implemented the complete dual-process Study 3 simulator and estimator runner.
- Added both combined and process-specific Y/Q first-stage diagnostics to
  Study 3 result rows.
- Added Study 3 condition fields to result summaries and result-file type
  normalization.
- Added all Study 3 functions to parallel worker exports.
- Added empirical standardized-beta diagnostics to the standalone BLUP-outcome
  pathway.

### Stage-2 estimator scaling

File: `R/stage2_estimators.R`

- Added an optional fixed `reporting_scale` to observed single- and
  dual-predictor regressions.
- Preserved the old observed-predictor-SD behavior when no fixed scale is
  supplied.
- Used fixed calibrated population scales in VH Studies 1-3 so proxy-specific
  variances do not redefine the estimand.

### Lai/OpenMx models

File: `R/lai_openmx_helpers.R`

- Added prefix support for all bivariate measurement-model columns, allowing Y
  and Q measurement inputs to coexist in one Study 3 data frame.
- Allowed the Study 2 Lai model to receive its nuisance intercept-path start
  and fixed reporting scale explicitly.
- Added a dual-process 2S-PA/2S-PAA model for Study 3.
- Defined the Study 3 reported algebra using the fixed calibrated population
  scale as

  $$
  \widehat\theta_1
  \frac{\operatorname{SD}_{\mathrm{cal}}(b_{y1})}
       {\operatorname{SD}_{\mathrm{cal}}(b_{q1})}.
  $$

- Included both structural paths and their covariance contribution when
  reconstructing the marginal Q-slope variance in the optional fitted-scale
  OpenMx algebra.

### Production execution and output integrity

Files:

- `vig_hallquist_2026/vh_runner.R`
- `vig_hallquist_2026/random_effects_structural_simulation.R`
- `vig_hallquist_2026/slurm/vig_hallquist_condition_array.sbatch`
- `vig_hallquist_2026/README.Rmd`

Changes:

- Added `pipeline_version = "standardized_beta_v1_20260615"` to every
  replication row.
- Resume mode now skips a condition only when its replication file has the
  current pipeline version, the exact expected method set, the requested
  replication-ID set, and the expected row count for every method. Files from
  an older pipeline or a different `n_sim` are rerun.
- Aggregate filenames include the study selection and chunk label so
  concurrently running study-specific or chunked jobs do not overwrite one
  another.
- Out-of-range array chunks return successfully with no conditions rather than
  failing the job.
- The canonical combined design has 3,888 conditions. With five conditions per
  task, the bundled SLURM script uses `#SBATCH --array=1-778%50`.
  Study-specific task counts are 144, 288, and 346 for Studies 1, 2, and 3.
- At 1,000 replications, Studies 1-3 produce 49,824,000 method rows. When the
  expected aggregate exceeds the default two-million-row threshold, the runner
  combines condition summaries without materializing a monolithic replication
  table and writes `*_condition_replication_index.csv` instead. The compressed
  per-condition files remain authoritative.
- The aggregation threshold can be changed with
  `VIG_HALLQUIST_MAX_AGGREGATE_ROWS`; setting it to `Inf` explicitly requests a
  full in-memory replication aggregate.

## Corrected Mathematical and Code Errors

### 1. Total \(R^2\) did not isolate the focal Study 2 coefficient

The earlier total-\(R^2\) calibration mixed the random-intercept effect,
random-slope effect, and their covariance. It also imposed a nonzero minimum
attainable \(R^2\) when the intercept coefficient was fixed. This made valid
null and small focal slope conditions impossible.

Correction: calibrate the signed standardized slope coefficient directly and
derive total \(R^2\) afterward.

### 2. The original `calibrate_beta2` quadratic solved a different problem

The quadratic equation itself was valid for selecting \(\beta_2\) to attain a
target total model \(R^2\). The error was treating that solution as if it held
the focal slope effect constant across conditions. It does not.

Correction: retain the quadratic only as an optional total-\(R^2\) method in
the technical note, not as the Study 2 design calibration.

### 3. Study 1 used the marginal correlation when drawing residual effects

After removing the explained slope variance, the residual slope SD changes.
Passing the original marginal correlation to the simulator no longer preserves
the requested marginal covariance.

Correction: compute and pass `rho_residual`, while retaining `marginal_rho` for
reliability calibration and reporting.

### 4. Study 2 passed standard deviations where covariance construction expected variances

The Study 2 simulator supplied `tau0` and `tau1` directly to
`make_random_effect_covariance()`, whose arguments are variances.

Correction: pass `tau0^2` and `tau1^2`.

### 5. Study 2 used the wrong intercept-slope correlation source

The simulator referenced `fixed_params$marginal_rho`, which was not the
condition-specific calibrated correlation.

Correction: use `condition$marginal_rho`.

### 6. Study 2 residual-outcome variance used the wrong intercept coefficient

The residual variance calculation and outcome generation used the global
`fixed_params$beta1z` even in `slope_only` conditions, where the condition's
intercept coefficient must be zero.

Correction: use `condition$beta1z` consistently in both the variance equation
and data generation, and verify the result against the stored calibrated
residual variance.

### 7. Stage-2 methods were not reporting the same estimand

Several estimators defaulted to scaling by the observed SD of their own BLUP or
corrected-score predictor. Because shrinkage and correction change those SDs,
methods could report different standardized quantities.

Correction: pass a common calibrated population reporting scale to every
method and rescale Fuller estimates, SEs, and confidence intervals explicitly.

### 8. Lai Study 2 used a hard-coded nuisance path and method-specific scale

The OpenMx wrapper expected a `fixed_params$beta_zu0` value used by another
simulation design and standardized the focal path using the fitted latent
slope variance.

Correction: pass the condition-specific nuisance start and the fixed
population reporting scale.

### 9. Bivariate Lai measurement columns could not be prefixed

The helper supported prefixes for univariate measurement inputs but ignored
them for bivariate inputs. This caused naming collisions when Y and Q each
contributed intercept/slope measurement models.

Correction: prefix all non-ID bivariate columns.

### 10. Study 3 lacked a correct focal-\(\theta_1\) calibration

Study 3 had no implemented dual-process calibration or simulator. A direct
reuse of total \(R^2\) would have confounded \(\theta_0\), \(\theta_1\), and
the Y intercept-slope covariance.

Correction: standardize \(\theta_0\) and \(\theta_1\) separately, preserve the
calibrated marginal Q covariance through the residual block, and reject
non-positive-definite conditions.

### 11. Study 3 needed the marginal Q-slope variance, not only residual variance, for standardization

The denominator of the focal standardized coefficient is
\(\operatorname{SD}(b_{q1})\), including structural and residual variation.
Using \(\operatorname{SD}(v)\) would change the target as the effect changed.

Correction: convert and report \(\theta_1\) using the calibrated marginal
\(\sqrt{G_{q,22}}\).

### 12. Residual covariance generation ignored non-iid condition settings

`draw_level1_residuals()` hard-coded iid residuals even when a condition
specified AR(1) or another supported residual covariance.

Correction: obtain the residual specification from `condition_to_r_spec()`.

### 13. The nlme AR(1) helper hard-coded the cluster variable

The AR(1) formula always grouped by `id`, while the VH data use `cid` in
several pathways.

Correction: accept `cluster_var` and construct the correlation formula from
the actual grouping column.

### 14. Study 1, Study 2, and Study 3 runners were incomplete

Study 1 and Study 2 were scaffolds, and the dual-process simulator explicitly
stopped as unimplemented.

Correction: implement complete simulation, fitting, correction, reporting, and
failure-row paths for all three studies.

### 15. Documentation still labeled beta-squared levels as the manipulated effect

Some documentation called `0`, `.04`, `.16`, and `.36` the structural-effect
grid even after the code had moved toward standardized beta.

Correction: document `0`, `.20`, `.40`, and `.60` as standardized-beta targets
and beta squared as a derived Study 1 diagnostic.

### 16. The documented low reliability level disagreed with the feasible design

Some study notes still listed `.20`, while the implemented posterior
reliability grids used `.25`. The `.20` target is infeasible for some planned
sparse and correlated conditions.

Correction: use and document `.25` as the low reliability level.

### 17. The calibration demonstration created duplicate beta columns

The example grid and calibration output both used `standardized_beta`, causing
duplicate names after column binding.

Correction: name the manipulated grid field `standardized_beta_target` and
retain `standardized_beta` for the calibrated value.

### 18. Study-specific design selection could change condition identities

Condition IDs based only on the selected rows would assign different seeds and
output paths to the same condition depending on whether a study was run alone
or as part of the combined design.

Correction: assign each study a canonical global ID range before applying
study selection or chunking.

### 19. Resume mode could silently reuse stale or incomplete results

Checking only for the existence of output files could reuse results generated
with older estimator definitions, a different replication count, or an
incomplete method set.

Correction: validate pipeline version, method membership, per-method row
counts, and replication IDs before skipping a condition.

### 20. Full aggregation was not feasible for the production design

The 1,000-replication combined run produces 49,824,000 method rows. Loading all
condition files into one table creates unnecessary memory and I/O risk.

Correction: aggregate condition-level summaries directly above a configurable
row threshold and write an index of authoritative replication files.

### 21. Study 3 first-stage diagnostics did not identify the affected process

A single unqualified diagnostic block could not distinguish failures in the Y
predictor process from failures in the Q outcome process.

Correction: retain combined diagnostics and add `stage1_y_*` and `stage1_q_*`
fields for singularity, fitted random-effect correlation, empirical-Bayes
score correlation, and Stage 2 design conditioning.

## Verification Added

The following checks now cover the beta-first pathway:

- Exact identity between the target beta and raw-coefficient conversion.
- Study 1 identity
  \(R^2_{\mathrm{struct}}=\beta_{\mathrm{std}}^2\).
- Exact attainment of posterior reliability targets.
- Positive-definiteness of residual and joint covariance matrices.
- Positive Study 2 external-outcome residual variance.
- Exact null focal coefficients when the target beta is zero.
- Separate nuisance \(\theta_0\) and focal \(\theta_1\) behavior in Study 3.
- Empirical recovery of the Study 3 focal standardized beta.
- Empirical recovery of the calibrated marginal Q covariance.
- End-to-end oracle recovery for Studies 1-3.
- Complete expected method sets for each study.
- Study 2 method-role classification and HC3 method expansion.
- Canonical condition counts, stable cross-selector condition IDs, and
  non-overlapping study ID ranges.
- Resume rejection for stale pipeline versions and mismatched replication
  counts.
- Process-specific Y/Q first-stage diagnostics in Study 3.
- Serial and parallel Study 3 execution.
- Small-run aggregate materialization and large-run condition-index
  aggregation.
- Backward compatibility for legacy non-calibrated BLUP-outcome grids.

Relevant tests include:

- `tests/test_stage2_estimators.R`
- `tests/test_blup_helpers_closed_form_gls.R`
- `tests/test_reliability_calibration.R`
- `tests/test_blup_outcome_designs.R`
- `tests/test_vh_blup_predictor_simulation.R`
- `tests/test_vh_study2_pathway.R`
- `tests/test_vh_standardized_beta_pathways.R`
- `tests/test_lai_openmx_inputs.R`
- `tests/test_openmx_lai_wrapper.R`

## Final Design Interpretation

The effect-size factor now has one interpretation across all three studies:

- Study 1: standardized effect of an observed cluster predictor on a latent
  random slope.
- Study 2: standardized effect of a latent random slope on an observed
  cluster outcome, adjusted for the latent random intercept in the primary
  model.
- Study 3: standardized effect of one process's latent random slope on another
  process's latent random slope, adjusted for the first process's latent random
  intercept.

Reliability determines the marginal random-effect variances. Standardized beta
then determines the focal raw coefficient. Any total or unique \(R^2\) is a
consequence of those choices rather than the quantity used to define them.
