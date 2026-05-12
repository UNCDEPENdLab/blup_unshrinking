# Technical Note: `fit_eiv_dual`

## Purpose

`fit_eiv_dual()` implements a dual-predictor errors-in-variables (EIV)
correction for Stage-2 regressions that use extracted random-effect scores as
predictors. In the current simulation workflow, the target regression is a
corrected outcome score regressed on two corrected predictor scores:
$$
y_i = \beta_0 + \beta_1 u_{0i} + \beta_2 u_{1i} + e_i,
$$

where \(u_{0i}\) is the intercept-like random-effect score and \(u_{1i}\) is
the slope-like random-effect score. The reported estimand is \(\beta_2\) scaled
by the estimated latent standard deviation of \(u_1\), so it is comparable to
the other one-SD Stage-2 estimator rows.

## Errors-in-Variables Usage

The EIV correction is applied to unshrunk or corrected-score predictors, not to
BLUPs. This distinction matters because classical EIV assumes an additive
measurement model:

$$
\mathbf{w}_i = \mathbf{u}_i + \boldsymbol{\eta}_i,\qquad
\operatorname{E}(\boldsymbol{\eta}_i \mid \mathbf{u}_i) = \mathbf{0},
\qquad
\operatorname{Var}(\boldsymbol{\eta}_i) = \mathbf{S}_i.
$$

Under that model, the observed predictor covariance is inflated:

$$
\operatorname{Var}(\mathbf{w}) =
\operatorname{Var}(\mathbf{u}) +
\operatorname{E}(\mathbf{S}_i).
$$

Therefore the latent predictor covariance is estimated by subtracting the known
or estimated sampling-error covariance. Applying this subtraction directly to
BLUPs would be conceptually wrong because BLUPs are already shrunken; their
variance is not an additively inflated version of the latent random-effect
variance.

## Corrected Normal Equations

For each complete case, define

$$
\mathbf{x}_i =
\begin{bmatrix}
1 \\
w_{0i} \\
w_{1i}
\end{bmatrix},
\qquad
\mathbf{M}_i =
\begin{bmatrix}
0 & 0 & 0 \\
0 & s_{11,i} & s_{12,i} \\
0 & s_{12,i} & s_{22,i}
\end{bmatrix}.
$$

The intercept is treated as measured without error. The supplied columns
`meas11`, `meas12`, and `meas22` populate the measurement-error block for the
two predictors. The estimator solves the corrected moment equation

$$
\left[
\sum_i \left(\mathbf{x}_i\mathbf{x}_i^\top -
\lambda\mathbf{M}_i\right)
\right]\hat{\boldsymbol{\beta}}
=
\sum_i \mathbf{x}_i y_i,
$$

where \(\lambda\) is `measurement_weight`. The default \(\lambda = 1\) is the
full classical EIV correction. Values \(0 < \lambda < 1\) deliberately subtract
only part of the measurement-error covariance and should be interpreted as
tempered or regularized EIV sensitivity analyses, not as the classical
correction.

The latent covariance diagnostic used for scaling and admissibility is

$$
\widehat{\boldsymbol{\Sigma}}_u(\lambda)
=
\frac{
\sum_i
(\mathbf{w}_i - \bar{\mathbf{w}})
(\mathbf{w}_i - \bar{\mathbf{w}})^\top
-
\lambda\sum_i \mathbf{S}_i
}{n - 1}.
$$

The reported estimate is

$$
\widehat{\theta}
=
\hat{\beta}_2
\sqrt{\widehat{\Sigma}_{u,22}(\lambda)}.
$$

## Standard Error Variants

`fit_eiv_dual()` now returns three conditional EIV standard-error variants:

- `naive`: a conventional homoskedastic model-based variance for the corrected
  normal equations;
- `hc0`: the empirical sandwich variance based on the corrected EIV estimating
  contributions;
- `hc3`: an HC3-style leverage-adjusted version of the same EIV sandwich.

The base method name, such as `eiv_dual_corrected`, corresponds to the `naive`
row. The robust variants append `_hc0` and `_hc3`. Earlier one-row versions of
`fit_eiv_dual()` were closest to the new `_hc0` row.

All three variants start from the same corrected estimating equation. Without
additional stabilization, the per-row estimating contribution is

$$
\boldsymbol{\psi}_i(\boldsymbol{\beta})
=
\mathbf{x}_i y_i
-
\mathbf{x}_i\mathbf{x}_i^\top\boldsymbol{\beta}
+
\lambda \mathbf{M}_i\boldsymbol{\beta}.
$$

This form is just the row-level version of the corrected normal equations:

$$
\sum_i \boldsymbol{\psi}_i(\hat{\boldsymbol{\beta}}) = \mathbf{0}.
$$

Equivalently, the usual OLS contribution
\(\mathbf{x}_i(y_i - \mathbf{x}_i^\top\boldsymbol{\beta})\) is augmented by
\(\lambda \mathbf{M}_i\boldsymbol{\beta}\), which restores the measurement-error
covariance subtracted from the predictor cross-products. Because
\(\mathbf{M}_i\) has zeros in the intercept row and column, this correction only
affects the two measured predictor equations.

Let

$$
\mathbf{A}
=
\sum_i \left(\mathbf{x}_i\mathbf{x}_i^\top - \lambda\mathbf{M}_i\right),
\qquad
\mathbf{B}
=
\sum_i
\boldsymbol{\psi}_i(\hat{\boldsymbol{\beta}})
\boldsymbol{\psi}_i(\hat{\boldsymbol{\beta}})^\top.
$$

The HC0 covariance estimator is

$$
\widehat{\operatorname{Var}}_{HC0}(\hat{\boldsymbol{\beta}})
=
\mathbf{A}^{-1}\mathbf{B}\mathbf{A}^{-1}.
$$

This is the unnormalized sandwich form corresponding to the summed estimating
equation. It is heteroskedasticity-robust at the Stage-2 row level because the
meat matrix \(\mathbf{B}\) is the empirical cross-product of the corrected
estimating contributions rather than a homoskedastic residual variance
multiplied by \(\mathbf{A}^{-1}\).

In this application, each Stage-2 row corresponds to a Stage-1 cluster. HC0
therefore lets clusters contribute to the variance according to the size and
direction of their corrected EIV estimating contribution. A cluster with a large
residual, extreme corrected scores, or a large measurement-error correction
block \(\mathbf{M}_i\) contributes more to the meat matrix than a cluster whose
corrected estimating contribution is small. This is useful because Stage-1
clusters are not equally informative after score extraction: some clusters have
more unstable corrected scores or stronger influence on the corrected normal
equations.

The naive/model-based row instead uses

$$
\widehat{\operatorname{Var}}_{naive}(\hat{\boldsymbol{\beta}})
=
\hat{\sigma}^2
\mathbf{A}^{-1}
\left(\sum_i \mathbf{x}_i\mathbf{x}_i^\top\right)
\mathbf{A}^{-1},
$$

with

$$
\hat{\sigma}^2
=
\frac{
\sum_i (y_i - \mathbf{x}_i^\top\hat{\boldsymbol{\beta}})^2
}{n - p}.
$$

This reduces to the usual OLS model-based variance when no EIV subtraction is
applied. With EIV subtraction, it treats the corrected normal-equation matrix as
fixed and uses a homoskedastic model for $$\operatorname{Var}(\sum_i
\mathbf{x}_i y_i)$$.

The HC3 row uses the same corrected estimating contributions as HC0 but inflates
each row by a leverage factor:
$$
\boldsymbol{\psi}_{i,HC3}
=
\frac{\boldsymbol{\psi}_i}{1 - h_i},
\qquad
h_i =
\mathbf{x}_i^\top \mathbf{A}^{-1}\mathbf{x}_i.
$$

The implementation caps \(h_i\) to \([0, 0.999]\) for numerical stability. The
HC3 meat is then

$$
\mathbf{B}_{HC3}
=
\sum_i
\boldsymbol{\psi}_{i,HC3}
\boldsymbol{\psi}_{i,HC3}^\top,
$$

and the covariance is

$$
\widehat{\operatorname{Var}}_{HC3}(\hat{\boldsymbol{\beta}})
=
\mathbf{A}^{-1}\mathbf{B}_{HC3}\mathbf{A}^{-1}.
$$

The HC3 adjustment is especially relevant for clusters with high leverage in
the corrected EIV design. Here leverage is not just ordinary OLS leverage from
the observed score values; it is computed against the EIV bread
\(\mathbf{A}^{-1}\), so it reflects the geometry after subtracting the
measurement-error covariance. A cluster can be high leverage because its
corrected intercept/slope score pair is unusual, because the corrected
predictor covariance is close to singular, or because the EIV correction makes
one predictor direction weakly identified. Dividing by \(1 - h_i\) increases
that cluster's contribution to the meat matrix:

$$
\boldsymbol{\psi}_{i,HC3}\boldsymbol{\psi}_{i,HC3}^\top
=
\frac{
\boldsymbol{\psi}_i\boldsymbol{\psi}_i^\top
}{(1 - h_i)^2}.
$$

Thus HC3 does not reweight the point estimate or give high-leverage clusters
more pull in solving the corrected normal equations. Instead, it gives those
clusters more weight in the uncertainty calculation, guarding against the usual
small-sample problem where high-leverage observations can make residual-based
variance estimates too optimistic.

If `stabilize_a_mat = TRUE` or `ridge_predictor_block = TRUE`, the function
solves a modified matrix \(\mathbf{A}^\ast\) rather than the raw corrected
matrix \(\mathbf{A}\). In that case, each estimating contribution is adjusted by

$$
\frac{\mathbf{A}^\ast - \mathbf{A}}{n}\hat{\boldsymbol{\beta}},
$$

so that the empirical contributions sum to zero for the same equation that was
actually solved:

$$
\boldsymbol{\psi}_i^\ast
=
\boldsymbol{\psi}_i
-
\frac{\mathbf{A}^\ast - \mathbf{A}}{n}\hat{\boldsymbol{\beta}},
\qquad
\sum_i \boldsymbol{\psi}_i^\ast = \mathbf{0}.
$$

The sandwich then uses \(\mathbf{A}^{\ast -1}\) as the bread and
\(\sum_i \boldsymbol{\psi}_i^\ast\boldsymbol{\psi}_i^{\ast\top}\) as the meat.
If direct inversion fails, the implementation falls back to a generalized
inverse for the bread. The same stabilized bread is used for all three
calculations.

Finally, the raw standard error for the target coefficient is extracted from
the \((3,3)\) element of the coefficient covariance matrix:

$$
\widehat{\operatorname{SE}}(\hat{\beta}_2)
=
\sqrt{
\widehat{\operatorname{Var}}(\hat{\boldsymbol{\beta}})_{33}
}.
$$

It is then put on the same one-SD scale as the reported estimate:

$$
\widehat{\operatorname{SE}}(\widehat{\theta})
=
\widehat{\operatorname{SE}}(\hat{\beta}_2)
\sqrt{\widehat{\Sigma}_{u,22}(\lambda)}.
$$

This standard error is conditional on the supplied corrected scores and
measurement-error covariance columns. It reflects empirical variability in the
corrected Stage-2 estimating equation and uses the supplied
measurement-error covariance in the EIV moment correction. It should not be
confused with the separate HC0-HC3 stacked-sandwich implementation in
`stacked_sandwich_for_corrected_scores()`.

The stacked-sandwich estimator explicitly stacks the Stage-1 mixed-model score
equations with the Stage-2 regression equations. In block form, it uses

$$
\mathbf{g}_i(\boldsymbol{\psi}, \boldsymbol{\alpha})
=
\begin{bmatrix}
\mathbf{s}_{1i}(\boldsymbol{\psi}) \\
\mathbf{s}_{2i}(\boldsymbol{\psi}, \boldsymbol{\alpha})
\end{bmatrix},
\qquad
\mathbf{A}
=
\begin{bmatrix}
\mathbf{A}_{11} & \mathbf{0} \\
\mathbf{A}_{21} & \mathbf{A}_{22}
\end{bmatrix},
$$

where \(\mathbf{s}_{1i}\) is the cluster-level Stage-1 likelihood score,
\(\mathbf{A}_{11}\) is the Stage-1 bread, \(\mathbf{A}_{22}\) is the Stage-2
bread, and \(\mathbf{A}_{21}\) captures how the corrected Stage-2 scores change
as the Stage-1 mixed-model parameters change. That \(\mathbf{A}_{21}\) block is
the piece that propagates Stage-1 parameter uncertainty into the Stage-2
coefficient variance.

By contrast, `fit_eiv_dual()` does not currently build this stacked block
system. Its sandwich uses only the corrected EIV Stage-2 estimating
contributions after the corrected scores and their OLS sampling covariance have
already been supplied. Thus, the EIV standard error accounts for row-level
variability in the corrected EIV estimating equation, but it is not the same as
the full Stage-1-plus-Stage-2 HC3 stacked sandwich.

## Non-Positive-Definite Checks

Before solving the corrected normal equations, `fit_eiv_dual()` checks whether
\(\widehat{\boldsymbol{\Sigma}}_u(\lambda)\) is positive definite. The main diagnostics are:

- `eiv_latent_cov_min_eigen`: the minimum eigenvalue of the corrected latent
  predictor covariance;
- `eiv_latent_cov_condition_number`: the ratio of maximum to minimum eigenvalue
  when the minimum eigenvalue is positive;
- `mx_issue_class`: `"ok"` on success or
  `"corrected_predictor_cov_not_pd"` when the corrected covariance is not
  admissible.

If `regularize = FALSE`, a non-positive-definite corrected predictor covariance
returns an unavailable estimate with `status_code = 3L`. This is intentional:
the full EIV subtraction has implied a latent predictor covariance that cannot
represent a valid covariance matrix.

## Regularization Options

`fit_eiv_dual()` has three stabilization paths, each with a different
interpretation.

1. `regularize = TRUE` adaptively reduces `measurement_weight` by binary search
   until \(\widehat{\boldsymbol{\Sigma}}_u(\lambda)\) is positive definite. This
   preserves the EIV form but changes the estimand to a tempered correction.
   The result records both `eiv_measurement_weight_requested` and
   `eiv_measurement_weight_used`.

2. `stabilize_a_mat = TRUE` projects the full corrected normal-equation matrix
   to positive definite by flooring eigenvalues at `min_eigen`. This is a
   numerical near-PD stabilization of the equation being solved, not a pure EIV
   correction.

3. `ridge_predictor_block = TRUE` adds ridge only to the corrected two-predictor
   block, targeting `ridge_min_eigen`, while leaving the intercept coupling in
   place. This is useful as a diagnostic for near-singular predictor geometry.

These options are best treated as sensitivity analyses. If the full correction
frequently needs them, the substantive conclusion is that the available
Stage-1 information is not strong enough to support an unregularized
two-predictor EIV deconvolution.

## Practical Reading

A successful `fit_eiv_dual()` row means the corrected latent predictor
covariance was admissible, the corrected normal equations were solvable, and
the scaled estimate and sandwich standard error were finite. A failure with
`corrected_predictor_cov_not_pd` is not just a numerical nuisance; it is a
substantive warning that the requested measurement-error subtraction is too
large relative to the observed predictor covariance in that replication or
dataset.
