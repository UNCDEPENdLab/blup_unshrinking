# From Posterior Reliability to a Random-Effects Structure

Consider the random-intercept/random-slope model

$$
\mathbf{y}_i = \mathbf{X}_i\boldsymbol{\beta}
  + \mathbf{Z}_i\mathbf{b}_i + \mathbf{e}_i,
\qquad
\mathbf{b}_i \sim N(\mathbf{0},\mathbf{G}),
\qquad
\mathbf{e}_i \sim N(\mathbf{0},\mathbf{R}_i).
$$

We want to choose the simulation parameters so that the random slope has a
target expected posterior reliability.

The revised simulation design uses:

$$
r_{\text{slope}}\in\{.25,.50,.80\}
$$

for low, medium, and high reliability, and

$$
R_{\text{struct}}^2\in\{0,.04,.16,.36\}
$$

for null, small, moderate, and large structural effects.

## 1. Define the target

For a given $\mathbf{G}$, design matrix $\mathbf{Z}_i$, and residual covariance
$\mathbf{R}_i$, the posterior covariance of the random effects is

$$
\mathbf{V}_i
=
\operatorname{Var}(\mathbf{b}_i\mid\mathbf{y}_i)
=
\left(
\mathbf{G}^{-1}
+\mathbf{Z}_i^\top\mathbf{R}_i^{-1}\mathbf{Z}_i
\right)^{-1}.
$$

If the slope is the second random effect, its expected posterior reliability is

$$
r_{\text{slope}}
=
1-\frac{1}{N}\sum_{i=1}^N
\frac{V_{i,22}}{G_{22}}.
$$

Thus, reliability is the expected proportion of marginal slope variance
remaining after posterior uncertainty is removed.

## 2. Fix the quantities reliability cannot identify

A reliability target does not uniquely determine every element of
$\mathbf{G}$. We therefore fix:

- the residual scale and structure, $\mathbf{R}_i$;
- the marginal random-intercept variance, $\tau_0^2$;
- the marginal intercept-slope correlation, $\rho$;
- each planned cluster's actual time design, $\mathbf{Z}_i$.

We then solve only for the marginal slope variance $\tau_1^2$.

## 3. Construct and calibrate the marginal G matrix

For any candidate $\tau_1^2$, define

$$
\mathbf{G}(\tau_1^2)
=
\begin{bmatrix}
\tau_0^2 & \rho\tau_0\tau_1 \\
\rho\tau_0\tau_1 & \tau_1^2
\end{bmatrix}.
$$

Compute

$$
f(\tau_1^2)
=
1-\frac{1}{N}\sum_{i=1}^N
\frac{
\left[
\left\{
\mathbf{G}(\tau_1^2)^{-1}
+\mathbf{Z}_i^\top\mathbf{R}_i^{-1}\mathbf{Z}_i
\right\}^{-1}
\right]_{22}
}{\tau_1^2}
-r_{\text{target}}.
$$

Numerically solve

$$
f(\tau_1^2)=0.
$$

The resulting matrix $\mathbf{G}_{\text{marginal}}$ is the covariance structure
that the first-stage model should see.

For the special case of iid residuals, centered time, and $\rho=0$,

$$
r_{\text{slope}}
=
\frac{\tau_1^2 S_{zz}}
{\sigma^2+\tau_1^2 S_{zz}},
\qquad
S_{zz}=\sum_j z_{ij}^2,
$$

so the solution is available directly:

$$
\frac{\tau_1^2}{\sigma^2}
=
\frac{r_{\text{target}}}
{(1-r_{\text{target}})S_{zz}}.
$$

## 4. Set the minimum reliability to .25

When intercepts and slopes are correlated, information about the intercept also
provides information about the slope. Consequently, slope reliability cannot
always approach zero as $\tau_1^2$ approaches zero.

For balanced clusters, iid residuals, centered time, and
$x=m\tau_0^2/\sigma^2$, the lower reliability bound is

$$
r_{\min}
=
\frac{\rho^2x}{1+x}
=
\rho^2
\frac{m\tau_0^2}{\sigma^2+m\tau_0^2}.
$$

Under the planned values $\tau_0=.9$, $\sigma=1$, and $|\rho|=.5$, this gives:

| Cluster size $m$ | Minimum attainable reliability |
|---:|---:|
| 3 | .177 |
| 5 | .201 |
| 10 | .223 |
| 25 | .238 |

The original low target of $.20$ is therefore infeasible for
$m\in\{5,10,25\}$ when $|\rho|=.5$. Raising the low target to $.25$ places it
above the lower bound for every planned cluster size and correlation condition.
The amended reliability grid is therefore

$$
\boxed{r_{\text{slope}}\in\{.25,.50,.80\}}.
$$

This amendment addresses feasibility of the posterior-reliability target. It
does not by itself solve the scaling problem created by holding a raw
structural coefficient constant across cells.

## 5. Account for a structural predictor of the slope

Suppose the data-generating model is

$$
b_{1i}=\gamma x_i+u_{1i},
\qquad
\operatorname{Var}(x_i)=1,
$$

where $b_{1i}$ is the total slope seen by a first-stage model that omits $x_i$.
Then

$$
\operatorname{Var}(b_{1i})
=
\gamma^2+\operatorname{Var}(u_{1i}).
$$

The standardized structural coefficient is

$$
\beta_{\text{std}}
=
\frac{\gamma}
{\sqrt{\operatorname{Var}(b_{1i})}}.
$$

Posterior-reliability calibration changes
$\operatorname{Var}(b_{1i})=G_{\text{marginal},22}$ across cluster sizes,
reliability levels, and correlations. Therefore, holding a raw $\gamma$
constant does not hold the standardized effect constant:

$$
\beta_{\text{std},c}
=
\frac{\gamma}
{\sqrt{G_{\text{marginal},22,c}}},
$$

where $c$ indexes simulation cells. A raw value such as $\gamma=.4$ can
represent a moderate effect in one cell, an extremely large effect in another,
or an impossible effect if $\gamma^2>G_{\text{marginal},22,c}$.

Instead, define the structural effect by the proportion of total slope variance
explained by $x_i$:

$$
R_{\text{struct}}^2
=
\frac{\gamma^2}
{G_{\text{marginal},22}}.
$$

With one standardized predictor,

$$
\beta_{\text{std}}
=
\operatorname{sign}(\gamma)\sqrt{R_{\text{struct}}^2}.
$$

The original standardized-effect targets
$\{0,.2,.4,.6\}$ therefore correspond exactly to

| Standardized effect $\beta_{\text{std}}$ | Structural $R^2$ |
|---:|---:|
| 0 | 0 |
| .20 | .04 |
| .40 | .16 |
| .60 | .36 |

The amended effect-size grid is

$$
\boxed{R_{\text{struct}}^2\in\{0,.04,.16,.36\}}.
$$

This preserves the original intended standardized effects while allowing the
raw coefficient to adapt to the slope variance in each reliability cell.

Calibrate the **marginal** slope variance first. For each cell, set

$$
\gamma_c
=
\operatorname{sign}(\beta_{\text{std}})
\sqrt{R_{\text{struct}}^2\,G_{\text{marginal},22}},
$$

and

$$
\operatorname{Var}(u_{1i})
=
(1-R_{\text{struct}}^2)G_{\text{marginal},22}.
$$

Thus, every $R_{\text{struct}}^2=.16$ cell has a standardized effect of $.40$,
even though its raw $\gamma_c$ differs. The fraction of slope variance explained
by $x_i$ is held constant rather than the arbitrary raw scale.

If $x_i$ is independent of the intercept and residual random effects, preserve
the calibrated covariance:

$$
\operatorname{Cov}(b_{0i},u_{1i})
=
G_{\text{marginal},12}.
$$

Therefore, the covariance matrix used to draw the residual random effects is

$$
\mathbf{G}_{\text{residual}}
=
\begin{bmatrix}
G_{\text{marginal},11}
&
G_{\text{marginal},12}
\\
G_{\text{marginal},12}
&
(1-R_{\text{struct}}^2)G_{\text{marginal},22}
\end{bmatrix}.
$$

The simulation draws $(b_{0i},u_{1i})^\top$ from this matrix and then constructs
$b_{1i}=\gamma_c x_i+u_{1i}$.

## 6. Check feasibility

Both covariance matrices must be positive definite:

$$
\lambda_{\min}(\mathbf{G}_{\text{marginal}})>0,
\qquad
\lambda_{\min}(\mathbf{G}_{\text{residual}})>0.
$$

Some low-reliability targets are impossible when the intercept is measured
precisely and strongly correlated with the slope. Such cells should be detected
and reported rather than silently modified.

When the marginal intercept-slope correlation is $\rho$, positive definiteness
of $\mathbf{G}_{\text{residual}}$ additionally requires

$$
R_{\text{struct}}^2 < 1-\rho^2.
$$

The largest planned structural effect is $.36$. With $|\rho|=.5$,
$1-\rho^2=.75$, so all four amended structural-effect levels satisfy this
condition.

## 7. Revised calibration workflow

For each cluster-size, time-design, residual-structure, reliability, and
correlation cell:

1. Set $r_{\text{target}}\in\{.25,.50,.80\}$.
2. Fix $\mathbf{Z}_i$, $\mathbf{R}_i$, $\tau_0^2$, and the marginal $\rho$.
3. Solve for $G_{\text{marginal},22}$.
4. Select $R_{\text{struct}}^2\in\{0,.04,.16,.36\}$.
5. Compute the cell-specific $\gamma_c$ and
   $\mathbf{G}_{\text{residual}}$.
6. Draw $(b_{0i},u_{1i})^\top$ from
   $\mathbf{G}_{\text{residual}}$ and construct
   $b_{1i}=\gamma_cx_i+u_{1i}$.

## 8. Why calibration occurs once per condition

Each simulation condition is intended to represent one fixed population
data-generating process. Its population parameters include

$$
\mathbf{G}_{\text{marginal}},\quad
\mathbf{G}_{\text{residual}},\quad
\gamma_c,\quad
\sigma,\quad
\rho_{\text{marginal}},\quad
\rho_{\text{residual}}.
$$

These quantities should be calibrated once when the condition is constructed
and then held fixed across all Monte Carlo replications of that condition.
Replications should differ because new subjects, covariates, random effects,
residuals, and possibly cluster sizes are sampled, not because the underlying
population changes.

For an unbalanced condition, let $\mathcal{D}_c$ denote the planned population
distribution of cluster designs in condition $c$. The calibrated reliability is

$$
r_c
=
1-
\mathbb{E}_{(\mathbf{Z}_i,\mathbf{R}_i)\sim\mathcal{D}_c}
\left[
\frac{V_{i,22}}{G_{c,22}}
\right].
$$

The deterministic reference profile approximates this expectation and produces
one value of $G_{c,22}$ for the condition. A particular replication will have a
finite realized set of cluster sizes, so its design-implied reliability,

$$
\widehat r_{cr}
=
1-\frac{1}{N_c}\sum_{i=1}^{N_c}
\frac{V_{cri,22}}{G_{c,22}},
$$

can fluctuate around $r_c$. That fluctuation is part of ordinary Monte Carlo
sampling variation and should remain visible.

Recalibrating within every replication would instead solve for a different
$G_{cr,22}$ using that replication's realized cluster sizes:

$$
\widehat r_{cr}(G_{cr,22})=r_c.
$$

This would create several problems:

1. The slope variance, raw $\gamma$, residual slope variance, and residual
   correlation would change from replication to replication.
2. Replications labeled as the same condition would no longer come from the
   same population.
3. Sampling variability in cluster sizes would be artificially removed from
   the reliability dimension of the simulation.
4. Bias, RMSE, and coverage would average over changing data-generating
   parameters rather than evaluate an estimator under one defined parameter
   point.
5. Reproducibility and interpretation would become harder because the true
   coefficient would itself depend on each realized design.

Condition-level calibration therefore separates two sources of variation:

- **Between conditions:** deliberate changes in target reliability,
  structural $R^2$, cluster-size mechanism, residual structure, and
  correlation.
- **Within a condition:** random sampling of subjects, realized cluster sizes,
  random effects, and residual outcomes under one fixed set of population
  parameters.

The condition manifest stores the calibrated population parameters before any
replications run. Every replication reads the same stored values, including
`gamma_x_on_slope`, `tau1_residual`, and `rho_residual`.

## 9. Verify each design cell

For every calibrated cell:

1. Recompute $r_{\text{slope}}$ from
   $\mathbf{G}_{\text{marginal}}$, $\mathbf{Z}_i$, and $\mathbf{R}_i$.
2. Confirm that it equals $r_{\text{target}}$ within numerical tolerance.
3. Simulate the residual random effects and reconstruct the total slopes.
4. Confirm that their empirical covariance approaches
   $\mathbf{G}_{\text{marginal}}$.
5. Confirm that the empirical structural $R^2$ approaches its target and that
   $\mathbf{G}_{\text{residual}}$ is positive definite.

The accompanying implementation is
[`posterior_reliability_calibration_demo.R`](../posterior_reliability_calibration_demo.R).

The reusable functions live in
[`R/reliability_calibration.R`](../R/reliability_calibration.R). The
BLUP-outcome simulation exposes two opt-in grid modes:

- `posterior_reliability_smoke`: one calibrated integration-test condition;
- `posterior_reliability`: the full amended 720-condition balanced grid.

Calibration occurs once when `make_blup_outcome_design()` constructs the
condition manifest. Each calibrated condition stores:

- `marginal_rho`: correlation in the target
  $\mathbf{G}_{\text{marginal}}$;
- `rho_residual`: correlation passed to `simulate_dataset()` when drawing
  $(b_{0i},u_{1i})$;
- `gamma_x_on_slope`: cell-specific raw coefficient implied by structural
  $R^2$;
- `tau1_residual`: residual slope standard deviation passed to
  `simulate_dataset()`.

Legacy grids do not contain these calibration fields and continue to pass their
existing `rho`, `tau1`, and `gamma_x_on_slope` values unchanged.
