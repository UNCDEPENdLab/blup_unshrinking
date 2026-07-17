# Study 4 Design: Heterogeneous Cluster Information and Reliability

**Date:** 17Jul2026

Michael Hallquist

## Purpose

Studies 1–3 test the performance of MLEs/unshrunk BLUPs with Fuller in three extracted random effects contexts:

1. a BLUP or corrected score used as an outcome;
2. a BLUP or corrected score used as a predictor;
3. scores from two first-stage models related to one another.

Those studies manipulate average cluster size and average posterior slope reliability between simulation conditions. Within each implemented condition, however, cluster size and the time design are identical, so the first-stage information is homogeneous across clusters. This is an important remaining gap because the main distinction between Lai's 2S-PA and 2S-PAA is whether the measurement loading and error matrices are allowed to vary by cluster. Thus, we need to address whether Zach's approach handles heterogeneous reliabilities as well as 2S-PA.

The primary question of Study 4 is: When clusters have the same average slope reliability but substantially different individual reliabilities, do row-specific unshrinking plus Fuller correction and Lai 2S-PA retain good point estimation, inference, and numerical stability, while their average-measurement counterparts (Fuller average measurement [see below] and 2S-PAA) deteriorate?

This is intended as a focused test of heterogeneous measurement information, not as a general robustness study. Informative missingness, residual autocorrelation, nonnormal random effects, and first-stage misspecification are outside the primary design because they introduce different  problems and would make it difficult to attribute method differences to reliability heterogeneity.

## Study 4 should build on study 2

Study 4 will use the Study 2 structural model: an observed cluster-level outcome is predicted by a latent random intercept and random slope from one first-stage mixed model. This extension is preferable to alternatives I considered because:

1. It isolates the predictor-measurement problem using only one first-stage model (relative to Study 3).
2. It is the setting in which unshrinking without an errors-in-variables correction is insufficient: unshrunk predictor scores are conditionally unbiased but noisy, so ordinary regression remains attenuated.
3. It permits a direct comparison among row-specific Fuller, row-specific Lai 2S-PA, and average-loading Lai 2S-PAA without simultaneously varying outcome-side score reliability.

## Important Predictor-Side Identity

Let the full-vector BLUP be

$$
\mathbf{m}_i = E(\mathbf{b}_i \mid D_i),
$$

where $D_i$ denotes the first-stage data. The posterior-mean projection
identity implies

$$
E(\mathbf{m}_i\mathbf{b}_i^\top)
=
E(\mathbf{m}_i\mathbf{m}_i^\top).
$$

Consequently, if

$$
z_i = \beta_0 + \boldsymbol{\beta}^\top\mathbf{b}_i + r_i
$$

and the full relevant random-effect vector is included, population regression of $z_i$ on $\mathbf{m}_i$ can recover $\boldsymbol{\beta}$ under the ideal Gaussian model. Thus, naive full-vector BLUP regression should not be treated as a guaranteed point-estimate failure. Extending Skrodal and Laake, the naive OLS standard errors are likely to still be wrong because score uncertainty is ignored.

This distinction makes the Study 4 comparison particularly informative:

- naive BLUP regression tests the posterior-mean cancellation result;
- unshrunk-score OLS demonstrates predictor measurement-error attenuation;
- unshrinking plus Fuller tests the complete proposed correction;
- average-measurement Fuller tests the same common-loading approximation as
  2S-PAA within the Fuller estimation framework;
- 2S-PA tests the corresponding row-specific latent measurement model;
- 2S-PAA tests the effect of replacing row-specific measurement inputs by
  averages.

## Data-Generating Model (extension of Study 2)

### First-stage repeated-measures model

For cluster $i$ and observation $j$,

$$
y_{ij}
=
\gamma_0 + \gamma_1 x_{ij}
+ u_{0i} + u_{1i}x_{ij} + e_{ij},
$$

with

$$
\begin{bmatrix}u_{0i}\\u_{1i}\end{bmatrix}
\sim
N\left(
\begin{bmatrix}0\\0\end{bmatrix},
\mathbf{G}
\right),
\qquad
e_{ij}\sim N(0,\sigma^2).
$$

As in Study 2, we use iid Gaussian residuals and the correctly specified random-intercept/random-slope model. Likeewise, the within-cluster time scores are centered and standardized. For a cluster of size $m_i$, this gives

$$
\sum_j x_{ij}=0,
\qquad
\sum_j x_{ij}^2=m_i-1.
$$

Cluster size therefore changes the slope information while preserving a common interpretation of the within-cluster predictor.

### Second-stage structural model

Generate the observed cluster-level outcome as

$$
z_i
=
\beta_{0z}
+ \beta_{1z}u_{0i}
+ \beta_{2z}u_{1i}
+ r_i,
\qquad
r_i\sim N(0,\sigma_z^2).
$$

The focal estimand is the standardized slope effect

$$
\beta_{2,\mathrm{std}}
=
\beta_{2z}
\frac{\operatorname{SD}(u_1)}{\operatorname{SD}(z)}.
$$

The nuisance intercept effect should remain fixed at the Study 2 value. The residual variance of $z_i$ is adjusted in each condition so that $\operatorname{Var}(z_i)=1$ and the requested standardized focal effect is preserved.

## Cluster-Information Profiles

Related to Lai's approach, we will vary cluster-level information (heterogeneity) by manipulating cluster sizes without introducing systematic biases between cluster conditions. All primary profiles have mean cluster size 10. Cluster-size assignments are generated independently of $u_{0i}$, $u_{1i}$, $z_i$, and all residuals. The
population profile is fixed within a condition, while cluster labels are randomly permuted in each replication.

| Profile | Cluster-size distribution | Purpose |
|---|---:|---|
| Homogeneous | all $m_i=10$ | Baseline with essentially common measurement inputs |
| Moderate heterogeneity | 50% $m_i=5$, 50% $m_i=15$ | Realistic unequal information at fixed mean cluster size |
| Severe heterogeneity | 50% $m_i=3$, 50% $m_i=17$ | Strong test of row-specific correction and numerical stability |

For odd group counts, allocate the extra cluster randomly to one size group so that the mechanism remains independent of the latent variables. Calibration
uses the exact population weights of one half for each group.

### Expected reliability separation

We try to retain the same approach to slope reliability targets as before (though we can omit the negative $\rho$ condition since it doesn't add anything new here. Using the current standardized time design and calibrating mean slope reliability to .50 gives approximately:

| Profile | $\rho(u_0,u_1)$ | Small-cluster reliability | Large-cluster reliability |
|---|---:|---:|---:|
| Moderate, $m=5/15$ | 0 | .35 | .65 |
| Severe, $m=3/17$ | 0 | .26 | .74 |
| Severe, $m=3/17$ | .50 | .31 | .69 |

These values provide substantial within-condition heterogeneity without introducing a misspecified residual model.

## Reliability Calibration

For any candidate marginal slope variance, define the cluster-type posterior covariance
$$
\mathbf{V}_g
=
\left(
\mathbf{G}^{-1}
+ \mathbf{Z}_g^\top\mathbf{R}_g^{-1}\mathbf{Z}_g
\right)^{-1},
$$

where $g$ indexes the cluster-size group. Calibrate $G_{22}$ so that the weighted mean posterior slope reliability equals the target:

$$
\bar\rho_{\mathrm{slope}}
=
\sum_g p_g
\left(
1-\frac{V_{g,22}}{G_{22}}
\right)
=
\rho_{\mathrm{target}}.
$$

In Study 2, every cluster within a condition has the same size and time design, so calibration sets the reliability of that single common design equal to the target. Study 4 instead solves for one common $G_{22}$ using the weighted mixture of cluster-size groups, so only the profile-average reliability equals the target while group-specific reliabilities deliberately differ. This is necessary to compare information profiles at the same mean reliability without calibrating away the reliability heterogeneity that Study 4 is designed to test.

Calibration must be performed separately for every combination of information profile, target mean reliability, and intercept-slope correlation. The standardized structural coefficient is then converted to raw units using the calibrated marginal slope variance, as in Studies 1–3.

The manifest should store both the mean target and distributional diagnostics:

- weighted achieved reliability;
- reliability in each cluster-size group;
- SD, range, and IQR of cluster-specific reliability;
- mean and range of `lambda22`;
- dispersion of the full $\boldsymbol{\Lambda}_i$ matrices;
- mean and range of `ols_var22`;
- dispersion of the full score-error covariance matrices.

## Primary Design Grid

| Factor | Levels |
|---|---|
| Number of clusters, $N$ | 50, 150, 300 |
| Mean posterior slope reliability | .25, .50, .80 |
| Information profile | homogeneous, moderate, severe |
| Intercept-slope correlation | 0, .50 |
| Standardized focal effect | 0, .40 |
| Structural target | intercept plus slope |
| Residual structure | Gaussian iid, correctly specified |

The primary grid contains

$$
3\times3\times3\times2\times2=108
$$

conditions. The information-matched control adds three targeted conditions (see below for details).

Use 1,000 Monte Carlo replications per final condition.

## Estimators

### Primary methods

1. **Oracle dual regression:** regress $z_i$ on the true $u_{0i}$ and $u_{1i}$.
2. **Naive full-vector BLUP regression:** regress $z_i$ on both BLUPs.
3. **Unshrunk full-vector OLS:** regress $z_i$ on the full corrected-score
   vector without predictor measurement-error correction.
4. **Unshrinking plus Fuller:** the main proposed method, using the
   cluster-specific full score-error covariance matrix.
5. **Average-measurement Fuller:** a Fuller-side analogue of 2S-PAA, using a
   common average loading matrix and average BLUP-error covariance matrix.
6. **Lai 2S-PA:** use cluster-specific $\boldsymbol{\Lambda}_i$ and $\boldsymbol{\Theta}_i$.
7. **Lai 2S-PAA:** replace the cluster-specific loading and error matrices by their sample averages.
8. **Direct MSEM:** fit a joint latent/multilevel model as a one-stage benchmark.
9. **Unshrinking plus Fuller alpha-stepdown:** retain as a stability estimator, especially for severe heterogeneity, low reliability, and small $N$. Overall, it's TBD whether stepdown will be our primary recommendation in the paper. Here, we should record and report its selected alpha, its failure rate, and any bias introduced by tempering the Fuller correction.

## Hypotheses

H1: Advantage of row-specific correction. At a given level of mean reliability, increasing reliability heterogeneity will not meaningfully increase the point-estimate bias of row-specific unshrinking plus Fuller or Lai 2S-PA.

H2: Effect of averaging measurement information. As the dispersion of cluster-specific \(\boldsymbol{\Lambda}_i\) and \(\boldsymbol{\Theta}_i\) increases, the difference between each row-specific estimator and its average-measurement counterpart will increase. (I think this only applies to the non-null structural conditions.)

H3: Naive OLS with unshrunk scores will not solve the problem. OLS regression using unshrunk predictor scores should be attenuated because the unshrunk scores contain additive sampling error. Fuller and 2S-PA should remove or substantially reduce this attenuation.

H4: BLUP predictor performance. Regression on the full BLUP vector will show little point-estimate bias as reliability heterogeneity increases. However, conventional homoskedastic standard errors may become miscalibrated as conditional score uncertainty becomes more heterogeneous.

​	(This could be strengthened by reporting an HC3 interval for the same BLUP regression. The primary concern under heterogeneous information is conditional heteroskedasticity. Comparing conventional and HC3 SEs would identify whether the inferential problem can be addressed without score unshrinking. But then we end up with the problem of when/whether to bring in HC3)

H5: Agreement when information is homogeneous. When cluster-specific measurement inputs are effectively constant, the mean difference between each row-specific estimator and its average-measurement counterpart will be negligible. Thus, row-specific and average-measurement Fuller should agree, as should 2S-PA and 2S-PAA. (In many respects, this should already be evident from the results of Studies 1 and 2.)

H6: Slope-information control condition. With zero intercept–slope correlation, the slope-information-matched control will produce smaller row-specific-versus-average differences than the ordinary severe-heterogeneity condition, despite having the same unequal cluster counts. That is, we should see much smaller PA–PAA and Fuller vs. Fuller-average differences than the ordinary severe profile.

H7: Numerical stability under sparse information. The largest 2S-PA and Fuller admissibility problems and extreme estimates should occur for $N=50$, mean reliability .25, and severe heterogeneity. For Fuller, alpha-stepdown should reduce failure rates and extreme RMSE, potentially at the cost of regularization bias.

- Perhaps we should also track the tails of the error distribution since RMSE tends to go wild with a few lousy estimates. Here, we could say alpha-stepdown will reduce the 95% upper error tail $|\hat{\beta} - \beta|$. 

H8: Direct-model benchmark. Under correct specification and adequate sample size, the row-specific two-stage methods should approach the MSEM and oracle estimator targets. 

## Monte Carlo Outcomes (stuff to store)

Report the common outcomes from Studies 1--3:

- bias and relative bias;
- RMSE;
- empirical SD of estimates;
- mean reported SE;
- mean-SE-to-empirical-SD ratio;
- 95% confidence-interval coverage;
- Type I error in null conditions;
- power in non-null conditions;
- first-stage singular-fit rate;
- Stage-2 nonconvergence or inadmissibility rate.

Add Study 4-specific outcomes:

- realized mean, SD, IQR, minimum, and maximum cluster reliability;
- dispersion of the loading and measurement-error matrices;
- within-replication 2S-PA minus 2S-PAA estimate difference;
- within-replication row-specific Fuller minus average-measurement Fuller
  estimate difference;
- comparison of the two row-specific-minus-average contrasts;
- within-replication Fuller minus 2S-PA estimate difference;
- Fuller correction and alpha-stepdown diagnostics;
- method runtime and convergence rate;
- cluster-size group-specific score bias and score RMSE.
- To bring in some robustness to RMSE's sensitivity to extreme estimates, add
  - median absolute error
  - 95th percentile of $|\hat{\beta} - \beta|$.

Side note: Bias, RMSE, and coverage among successful fits should always be displayed alongside failure rates. A method with many failed or inadmissible fits should not receive a favorable performance interpretation merely because its successful subset has low RMSE. This is part of what we've seen with Lai in some near-singular situations (it fails to converge, but the RMSE on converged samples looks good.)

----

# Supplemental Materials



## Information-Matched Falsification Control

### Why an information-matched control is needed

(This was mostly drafted by GPT-5.6. I've left the somewhat excessive detail here to let us evaluate it fully.)

The primary heterogeneous profiles deliberately make cluster size a source of
slope-information heterogeneity. With the ordinary centered and standardized
time design,

$$
\sum_j x_{ij}=0,
\qquad
\sum_j x_{ij}^2=m_i-1.
$$

Consequently, a cluster with 17 observations contributes substantially more
information about its random slope than a cluster with 3 observations. This is
the intended manipulation in the primary design, but it means that two features
change together: the number of observations and the reliability of the slope
estimate. If row-specific and average-measurement estimators separate in the
ordinary severe condition, that result alone does not establish that the
separation is caused by heterogeneous slope information. Unequal row counts
could instead affect numerical stability, leverage, or some other aspect of
the fitting algorithms.

The information-matched condition provides a focused falsification test of
that alternative explanation. It retains exactly the same severe cluster-size
distribution---half of the clusters have $m_i=3$ and half have $m_i=17$---but
removes the corresponding difference in slope information. The comparison is
therefore between two conditions with the same unequal counts but different
amounts of slope-reliability heterogeneity:

| Condition | Cluster sizes | Slope information |
|---|---:|---:|
| Ordinary severe profile | 3 and 17 | Heterogeneous |
| Information-matched severe profile | 3 and 17 | Equalized |

If the row-specific-versus-average contrasts are pronounced in the ordinary
severe profile but collapse in the information-matched profile, the evidence
more specifically implicates heterogeneous measurement information rather
than unequal counts themselves.

### Construction of the matched design

Use the severe $m_i=3/17$ profile, but multiply the centered, standardized
slope column in cluster $i$ by

$$
a_i
=
\sqrt{\frac{9}{m_i-1}}.
$$

Define the modified predictor as $x_{ij}^{*}=a_i x_{ij}$. Because multiplication
by a constant preserves centering,

$$
\sum_j x_{ij}^{*}=a_i\sum_jx_{ij}=0.
$$

Its total squared deviation is

$$
\sum_j(x_{ij}^{*})^2
=a_i^2\sum_jx_{ij}^2
=\frac{9}{m_i-1}(m_i-1)
=9.
$$

The value 9 is not arbitrary. It is the slope sum of squares in the homogeneous
$m=10$ design because its standardized predictor satisfies
$\sum_jx_{ij}^2=10-1=9$. Thus, the control gives both cluster-size groups the
same slope-information budget as the homogeneous reference condition rather
than introducing a new overall information level.

For iid residuals, $\mathbf R_i=\sigma^2\mathbf I$, and a centered random
intercept/random-slope design has

$$
\mathbf Z_i^\top\mathbf R_i^{-1}\mathbf Z_i
=
\frac{1}{\sigma^2}
\begin{bmatrix}
m_i & 0\\
0 & \sum_j(x_{ij}^{*})^2
\end{bmatrix}
=
\frac{1}{\sigma^2}
\begin{bmatrix}
m_i & 0\\
0 & 9
\end{bmatrix}.
$$

The slope-information entry is therefore exactly $9/\sigma^2$ for both the
$m_i=3$ and $m_i=17$ clusters. When the random intercept and slope are
uncorrelated, the slope block separates from the intercept block. Writing the
marginal slope variance as $\tau_1^2$, the posterior slope variance becomes

$$
V_{i,22}
=
\left(\frac{1}{\tau_1^2}+\frac{9}{\sigma^2}\right)^{-1},
$$

and the corresponding posterior reliability is

$$
\rho_{i,\mathrm{slope}}
=
1-\frac{V_{i,22}}{\tau_1^2}
=
\frac{9\tau_1^2}{\sigma^2+9\tau_1^2}.
$$

Neither expression depends on $m_i$. The sampling variance of the
cluster-specific OLS slope is likewise

$$
\operatorname{Var}(\widehat u_{1i}^{\mathrm{OLS}}\mid u_{1i})
=\frac{\sigma^2}{\sum_j(x_{ij}^{*})^2}
=\frac{\sigma^2}{9}.
$$

Thus, at the population parameter values, the control equalizes the principal
slope quantities used by both estimation frameworks: posterior slope
reliability, the slope element of the BLUP loading matrix, and the
sampling-error variance of the unshrunk slope score. This makes it a direct
check on whether the advantage of row-specific $\boldsymbol\Lambda_i$,
$\boldsymbol\Theta_i$, or score-error matrices is actually activated by
heterogeneous slope information.

### Why the rescaling is defensible and interpretable

The rescaling changes the predictor design before outcomes are generated; it
does not transform an estimated slope after seeing the data. Both the data
generator and the fitted mixed model use the same $x_{ij}^{*}$ values, so the
first-stage model remains correctly specified. The random coefficient $u_{1i}$
continues to mean the change in the outcome associated with a one-unit change
in the common predictor scale.

The design can be interpreted as trading the number of observations against
their leverage for estimating a linear trend. The three-observation clusters
are observed at predictor values farther from their cluster mean, whereas the
17-observation clusters are observed at more closely spaced predictor values.
The sparse clusters have fewer observations but broader predictor support; the
dense clusters have more observations but narrower support. Both schedules
have the same total squared distance from the predictor mean and therefore the
same information about a linear slope under homoscedastic errors. This is the
same principle used in regression and experimental design: slope precision
depends not only on sample size but also on the spread of the predictor.

This schedule need not be presented as the most typical applied sampling
pattern. Its value is diagnostic. A falsification control is useful precisely
because it changes the feature thought to be causal---slope information---while
retaining the most obvious alternative feature---unequal cluster counts. The
factor $a_i$ provides a transparent, deterministic way to make that
intervention without adding missingness, heteroscedastic residuals, or model
misspecification.

### Scope and interpretation of the control

The control equalizes slope information, not every aspect of the
random-intercept/random-slope measurement problem. The intercept-information
entry remains $m_i/\sigma^2$, so clusters with 17 observations still estimate
their random intercepts more precisely than clusters with 3 observations.
The full $\boldsymbol\Lambda_i$ and $\boldsymbol\Theta_i$ matrices can therefore
retain some cluster-size variation even when their focal slope components are
nearly identical. Restricting the control to zero intercept-slope correlation
makes the intercept and slope information blocks orthogonal and sharply limits
the extent to which the remaining intercept heterogeneity can contaminate the
focal slope comparison.

The control also does not remove the numerical consequences of having only
three observations in some clusters. That is a feature, not a defect: if
row-specific-versus-average differences persist after slope information is
matched, the result would point away from heterogeneous slope reliability and
toward sparse-cluster numerics, residual intercept heterogeneity, or another
mechanism. Accordingly, the result should be interpreted as evidence about
the slope-information mechanism, not as proof that all components of
reliability have been equalized.

The implementation should verify in every replication that:

- the $m_i=3/17$ count distribution is unchanged from the ordinary severe
  profile;
- $\sum_jx_{ij}^{*}=0$ and $\sum_j(x_{ij}^{*})^2=9$ within every cluster;
- the SD of cluster-specific slope reliability is numerically near zero;
- dispersion of `ols_var22` is numerically near zero;
- dispersion of the fitted `lambda22` values is much smaller than in the
  ordinary severe profile, recognizing that it need not be exactly zero when
  the Stage 1 intercept-slope covariance is estimated; and
- remaining dispersion in the full loading or error-covariance matrices is
  attributable primarily to the random-intercept block and finite-sample
  estimated intercept-slope coupling.

### Expected comparison

This is a falsification control rather than a primary design factor. If the
2S-PA minus 2S-PAA and row-specific Fuller minus average-measurement Fuller
contrasts are small here but large in the ordinary severe profile, that pattern
supports the claim that heterogeneous slope reliability, rather than unequal
counts by themselves, drives the primary method differences. If the contrasts
remain large, the proposed mechanism is not adequately isolated and the
remaining intercept-information and sparse-cluster explanations should be
investigated.

Run this control only for:

- $N\in\{50,150,300\}$;
- mean reliability .50;
- intercept-slope correlation 0;
- standardized focal effect .40.

## Average-measurement Fuller logic and implementation

The average-measurement Fuller estimator is included to make the test of row-specific measurement information internal to each estimation framework. Without it, row-specific Fuller can be compared with 2S-PAA, but any difference mixes two changes: the treatment of heterogeneous measurement information and the Stage-2 fitting criterion. Adding this estimator creates the paired design:

| Measurement treatment | Fuller | Lai/OpenMx |
|---|---|---|
| Row-specific | unshrinking plus Fuller | 2S-PA |
| Sample-average approximation | average-measurement Fuller | 2S-PAA |

For cluster $i$, write the BLUP measurement model as

$$
\tilde{\mathbf b}_i
=
\boldsymbol\Lambda_i\mathbf b_i+\boldsymbol\delta_i,
\qquad
\operatorname{Var}(\boldsymbol\delta_i)=\boldsymbol\Theta_i.
$$

The row-specific Fuller estimator uses

$$
\breve{\mathbf b}_i=\boldsymbol\Lambda_i^{-1}\tilde{\mathbf b}_i,
\qquad
\mathbf S_i=
\boldsymbol\Lambda_i^{-1}\boldsymbol\Theta_i
\boldsymbol\Lambda_i^{-\top}.
$$

In the implemented Gaussian mixed-model pathway, these are the full-vector closed-form corrected scores and their cluster-specific OLS/GLS sampling-error covariances. Average-measurement Fuller instead follows the same averaging order as the existing 2S-PAA implementation. First compute the elementwise sample averages:
$$
\bar{\boldsymbol\Lambda}
=\frac{1}{N}\sum_{i=1}^N\boldsymbol\Lambda_i,
\qquad
\bar{\boldsymbol\Theta}
=\frac{1}{N}\sum_{i=1}^N\boldsymbol\Theta_i.
$$

Then apply the common inverse loading to every original BLUP vector:

$$
\breve{\mathbf b}_i^{\mathrm{avg}}
=\bar{\boldsymbol\Lambda}^{-1}\tilde{\mathbf b}_i,
$$

and give every cluster the same additive measurement-error covariance

$$
\mathbf S^{\mathrm{avg}}
=
\bar{\boldsymbol\Lambda}^{-1}
\bar{\boldsymbol\Theta}
\bar{\boldsymbol\Lambda}^{-\top}.
$$

Standard dual-predictor Fuller estimation is then applied to \(\breve{\mathbf b}_i^{\mathrm{avg}}\), assigning every cluster the same three unique entries of the symmetric covariance matrix \(\mathbf S^{\mathrm{avg}}\). The raw slope coefficient and its conditional standard error are multiplied by the population slope SD, exactly as for row-specific Fuller, so both methods target the same standardized effect.

The implementation is split into two helpers in `R/stage2_estimators.R`:

1. `prepare_fuller_average_measurement()` constructs $\bar{\boldsymbol\Lambda}$, $\bar{\boldsymbol\Theta}$, the transformed BLUPs, and $\mathbf S^{\mathrm{avg}}$.
2. `fit_fuller_average_measurement()` calls the existing `fit_fuller_dual()` estimator with those transformed inputs and returns the standard Fuller result schema. A noninvertible average loading matrix is recorded as a method-specific failed fit rather than aborting the entire replication.

This estimator must not be implemented by first computing every $\boldsymbol\Lambda_i^{-1}\tilde{\mathbf b}_i$ and then replacing $\mathbf S_i$ by $N^{-1}\sum_i\mathbf S_i$. That alternative retains row-specific unshrinking and averages only the error covariance. In general,
$$
\bar{\boldsymbol\Lambda}^{-1}\bar{\boldsymbol\Theta}
\bar{\boldsymbol\Lambda}^{-\top}
\ne
\frac{1}{N}\sum_i
\boldsymbol\Lambda_i^{-1}\boldsymbol\Theta_i
\boldsymbol\Lambda_i^{-\top},
$$

so it does not reproduce the approximation made by 2S-PAA. Finally, average-measurement Fuller and 2S-PAA should not be expected to be numerically identical: they share the same average measurement approximation but use different finite-sample objectives and standard-error calculations.
