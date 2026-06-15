# Overview

## Four studies

1. BLUP as outcome
2. BLUP as predictor
3. Dual BLUPs: use dual-process growth as example, regressing slope from one growth model on intercept and slope from the other
4. Heterogeneous cluster sizes/reliabilities -- this is the advantage of Lai 2S-PA over the 2S-PAA. I still need to flesh this out...

## Notes

1. Now that I'm looking at everything, I think we can probably drop the diagonal-only correction approach. We know it's wrong in principle and we already have a lot of other estimators to test, including some other lousy ones. :)
2. Use slope reliability, not ICC, as the first-stage information factor. For each planned cluster-size and time/design matrix condition, calibrate the slope variance/residual variance ratio to achieve low, medium, and high expected posterior reliability:
$$
\rho_{\text{slope}} =
1 - \mathbb{E}_i \left[
\frac{\operatorname{Var}(b_{1i} \mid y_i)}
     {\operatorname{Var}(b_{1i})}
\right].
$$

# Study 1: BLUP as outcome

Data generating model
$$
y_{ij} = \beta_0 + \beta_z z_{ij} + b_{0i} + b_{1i} z_{ij} + e_{ij}
$$

$$
b_{1i} = \gamma x_i + u_{1i}
$$

The study focus is on $\gamma$, the effect of a covariate $x$ on the random slope $b_1$.

## Design matrix

- N: 30, 50, 100, 150, 300
- Cluster size (fixed within cell): 3, 5, 10, 25
- expected posterior slope reliability: low (.25), medium (.50), high (.80)
- intercept-slope correlation: -.5, 0, .5
- effect size (effect of x on on b): 0, .2, .4, .6

Total cells: 720

## Estimators

- Oracle regression of true `b_1i` on `x_i`
- Naive BLUP regression
- OLS Full-matrix-corrected BLUP
- Fuller Full-matrix-corrected BLUP
- Fuller stepdown full-matrix-corrected BLUP
- Single-subject OLS slope outcome
- Lai 2S-PA
- Direct MLM: `y ~ x + z + x:z + (1 + z | id)`

Should the direct MLE (GLS) scores be included as a separate estimator? Or just defer all of that -- which is more computational than statistical -- to a sensitivity study of MLE-generating method? I'm inclined to avoid the direct MLE scores too soon.

## Hypotheses

H1. Naive BLUP-as-outcome regressions should be attenuated, especially when
cluster size is small, slope reliability is low, and shrinkage is strong.

H2. Full-matrix unshrinking (and closed-form scores) should reduce or remove
point-estimate bias under correct first-stage specification.

H3. The direct MLM should be the best benchmark and
should agree with the oracle target asymptotically.

H4. The benefit of the full-matrix correction over diagonal correction should be
largest when random intercepts and slopes are correlated.

H5. Single-subject OLS slopes should be unbiased in principle but inefficient
and unstable when cluster sizes are low. Unshrinking should recover similar
likelihood-only information while using the first-stage model structure more
systematically.

# Study 2: BLUP as predictor

Data-generating model is same as above

Generate a cluster-level outcome that is related to the true latent random effects.
Include two structural targets:

1. Slope-only target:
$$
w_i = \alpha + \delta_1 b_{1i} + r_i
$$

2. Intercept-plus-slope target:
$$
w_i = \alpha + \delta_0 b_{0i} + \delta_1 b_{1i} + r_i
$$

The intercept-plus-slope target is essential because it tests whether slope effects are contaminated by intercept-slope covariance or by EB slope scores carrying information about the true intercept.

## Design matrix

Same as study 1.

- N: 30, 50, 100, 150, 300
- Cluster size (fixed within cell): 3, 5, 10, 25
- expected posterior slope reliability: low (.25), medium (.50), high (.80)
- intercept-slope correlation: -.5, 0, .5
- structural target: slope-only, intercept-plus-slope
- effect size (effect of b on on w): 0, .2, .4, .6

Total cells: 1440

## Estimators

- Oracle regression of `w_i` on true `b_1i` for the slope-only target, or on true `b_0i`, `b_1i` for the intercept-plus-slope target
- OLS Naive BLUP (predictor) regression with intercept + slope
- OLS unshrunk BLUP (predictor) regression with intercept + slope
- (Centered) Slope-only BLUP predictor (showing contamination)
- Fuller Full-matrix-corrected BLUP
- Fuller stepdown full-matrix-corrected BLUP
- Lai 2S-PA
- MSEM (as in Lai)

## Hypotheses

H1. Naive BLUPs as  predictors will have minimally biased point estimates  (esp. compared to Study 1), but their standard errors and coverage will be wrong.

​	H1a. Unshrunk scores used as predictors require Fuller or Lai 2S-PA to obtain 95% coverage.

H2. Full-vector methods should outperform slope-only methods when random
intercepts and slopes are correlated or when EB slope scores have cross-loadings
on the true intercept.

H3. The largest differences among methods will appear in low-reliability, small-cluster, and dual-predictor conditions.

H4. Null-effect conditions will show that naive OLS with BLUPs leads to inflated Type I errors.

# Study 3: Dual-ing BLUPs!

I think we should motivate this with the idea of dual-process growth models. This is a common application in psychology -- although it's also a case where this is easy enough to do in SEM, or with MLMs with process indicator variables. Still, I think that may be intuitive.

### Data-Generating Model

Generate two growth processes:

$$
\begin{align}
y_{ij} &= \beta_{y0} + \beta_{yz} z_{ij} + b_{y0i} + b_{y1i} z_{ij} + e_{yij} \\
q_{ik} &= \beta_{q0} + \beta_{qs} s_{ik} + b_{q0i} + b_{q1i} s_{ik} + e_{qik}
\end{align}
$$

Then define the second-process slope from a latent structural model. In the slope-only case,

$$
\eta_i = \theta b_{y1i}, \qquad
b_{q1i} = \eta_i + v_i.
$$

In the intercept-plus-slope case,

$$
\eta_i = \theta_0 b_{y0i} + \theta_1 b_{y1i}, \qquad
b_{q1i} = \eta_i + v_i.
$$

The target is $\theta$ in the slope-only case, or $\theta_1$ in the intercept-plus-slope case -- the relationship between latent slopes.

Generate the predictor-process random effects as

$$
\mathbf{b}_{yi} =
\begin{bmatrix}
b_{y0i} \\
b_{y1i}
\end{bmatrix}
\sim \mathcal{N}(\mathbf{0}, \mathbf{G}_y).
$$

Draw $(b_{q0i},v_i)$ independently of $\mathbf{b}_{yi}$, but allow their
covariance to preserve the planned Q intercept-slope covariance. Calibrate both process
slope variances from their posterior-reliability targets first. Define the focal
effect using the standardized slope coefficient

$$
\theta_{1,\mathrm{std}}
=
\theta_1
\frac{\operatorname{SD}(b_{y1i})}
     {\operatorname{SD}(b_{q1i})}.
$$

For the intercept-plus-slope condition, define $\theta_0$ on the same
standardized scale and hold it fixed as a nuisance effect. Convert both
coefficients to raw units, then set

$$
\sigma_v^2 =
\operatorname{Var}(b_{q1i})
-
\operatorname{Var}(\eta_i).
$$

The residual covariance between $b_{q0i}$ and $v_i$ preserves the
reliability-calibrated marginal covariance between the Q intercept and Q slope.
For null focal conditions, set $\theta_1=0$ while retaining $\theta_0$ in the
intercept-plus-slope condition.

This construction implies the full covariance matrix for
$(b_{y0i}, b_{y1i}, b_{q0i}, b_{q1i})^\top$. Let
$\mathbf{c} = (0, \theta)^\top$ in the slope-only case and
$\mathbf{c} = (\theta_0, \theta_1)^\top$ in the intercept-plus-slope case. Then

$$
\mathbf{G}_b =
\begin{bmatrix}
\mathbf{G}_y & \mathbf{0} & \mathbf{G}_y \mathbf{c} \\
\mathbf{0}^\top & G_{q,11} & G_{q,12} \\
\mathbf{c}^\top \mathbf{G}_y & G_{q,12} & G_{q,22}
\end{bmatrix}.
$$

Positive definiteness is enforced by requiring both the residual
$(b_{q0i},v_i)$ block and $\mathbf{G}_b$ to be positive definite before running
a condition.

## Design Matrix

Let's make this simpler in some respects since we now have cluster sizes for each process, which should presumably cross.

- N: 50, 100, 150, 300
- Cluster size for process y: 3, 5, 10
- Cluster size for process q: 3, 5, 10
- expected posterior slope reliability for Q and Y: low (.25) and high (.80) each. This leads to 4 combinations: high/high, low/low, low/high, high/low
- intercept-slope correlation: 0, .5
- effect size (effect of x on on b): 0, .2, .5

One way to cut down on the number of cells would be to have small symmetric (3 x 3), large symmetric (10 x 10), and asymmetric (3 x 10) cells for cluster sizes. We don't care if y vs. q is smaller, so we don't need to compare 5/3 to 3/5 and so on.

## Estimators

- Oracle regression using true latent slopes.
- OLS Naive BLUP-on-BLUP regression.
- OLS BLUP predictor with unshrunk outcome.
- OLS unshrunk BLUP outcome with BLUP predictor.
- OLS Unshrunk predictor with unshrunk outcome.
- Fuller correction on unshrunk outcome and predictor.
- Lai 2S-PA
- Direct joint latent model in SEM

## Hypotheses

H1. Naive BLUP-on-BLUP regressions should be unreliable because outcome-side
shrinkage, predictor-side measurement error, and heterogeneous reliability act
together.

H2. Correcting only one side should help only for the component of the problem it
targets. For example, unshrinking the outcome addresses outcome attenuation but
does not fix predictor measurement error.

H3. Full two-sided measurement correction, either through Lai 2S-PA or an
unshrunk-score Fuller moment correction, should best recover point estimates and coverage.

H4. The method differences should be largest when one first-stage process has
many observations and the other has few, because this creates asymmetric score
reliability.
