# Extending Zach's Score Correction from Computational Models to Gaussian MLMs

Michael Hallquist

17Apr2026

## Overview

These notes build on my review of MLM notation in Bates's [lme4 documentation](https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf), Hox's 2018 *Multilevel Analysis* book, and the technically weighty Demidenko 2013 book. Bates covers how to obtain posterior means and covariances under the frequentist Gaussian mixed-model setup. I also went through `mlm_random_slope_blup_correction_sim.R`, especially `unweight_random_effects` and `get_corrected_score`, and worked backward from the code to the underlying mixed-model algebra.

I also came across a paper from Wang and Merkle (2018) -- corresponding to the `merDeriv` package -- that I think may prove helpful in working toward standard errors since they've worked out key parts that would enter into a sandwich estimator in the second-stage/structural model. Here are two that seemed essential to have in hand, based on their paper:

`sandwich::estfun(mlm_obj, level = 2, ranpar = "var")` captures each subject’s contribution to the MLM (stage-1) score vector.

`sandwich::bread(mlm_obj, full = TRUE, information = "observed", ranpar = "var")` captures the MLM (stage-1) observed information structure.

Wang, T., & Merkle, E. C. (2018). merDeriv: Derivative Computations for Linear Mixed Effects Models with Application to Robust Standard Errors. *Journal of Statistical Software*, *87*, 1–16. https://doi.org/10.18637/jss.v087.c01

## MLM Notation

Here is the standard Gaussian mixed-model notation approximately matching Bates 2015. For cluster or subject $i$, let

$$
\mathbf{y}_i = \mathbf{X}_i \boldsymbol{\beta} + \mathbf{Z}_i \mathbf{b}_i + \boldsymbol{\varepsilon}_i,
$$

$$
\boldsymbol{\varepsilon}_i \sim \mathcal{N}(\mathbf{0}, \mathbf{R}_i),
$$

$$
\mathbf{b}_i \sim \mathcal{N}(\boldsymbol{\mu}_b, \mathbf{G}).
$$

Here:

- $\mathbf{y}_i$ is the vector of outcomes for subject $i$
- $\mathbf{X}_i$ is the fixed-effects design matrix
- $\mathbf{Z}_i$ is the random-effects design matrix
- $\boldsymbol{\beta}$ is the fixed-effects parameter vector
- $\mathbf{b}_i$ is the $q \times 1$ subject-specific random-effects vector
- $\mathbf{R}_i$ is the within-subject residual covariance matrix
- $\mathbf{G}$ is the between-subject covariance matrix of the random effects

It is also helpful to note the marginal covariance of $\mathbf{y}_i$ (cf. Wang and Merkle equations 3 and 5):

$$
\mathbf{V}_i = \mathbf{Z}_i \mathbf{G} \mathbf{Z}_i^\top + \mathbf{R}_i.
$$

That is the covariance of $\mathbf{y}_i$ after integrating over the random effects.

For the random intercept plus slope example already coded, $\mathbf{b}_i = (b_{0i}, b_{1i})^\top$ and $\mathbf{Z}_i$ has one column for the random intercept and one column for the random slope.

### Important point about centering of random effects

In standard MLM parameterizations, the random effects are centered and the population means live in $\boldsymbol{\beta}$. That is,

$\mathbb{E}\left[\mathbf{b}_i\right] = 0; \boldsymbol{\mu}_b = \mathbf{0}$

In other words, the model-implied prior mean of the random effects is 0. This simplifies the score correction.

## Thinking about the BLUPs in terms of the conditional distribution of the random effects

Following through Bates's section on conditional variances of random effects, we are interested in understanding the conditional distribution of the random effects, $\mathbf{b}_i$. Using plug-in estimates from the fitted model, the random effects have a Gaussian conditional distribution:

$$
\mathbf{b}_i \mid \mathbf{y}_i, \hat{\psi}
\sim
\mathcal{N}(\hat{\mathbf{m}}_i, \hat{\mathbf{V}}_i).
$$

where $\hat{\psi} = (\hat{\boldsymbol{\beta}}, \hat{\mathbf{G}}, \hat{\mathbf{R}}_i)$ collects the fitted MLM (stage-1) quantities being used as fixed plug-in values. $\hat{\mathbf{m}}_i$ is the posterior mean of $\mathbf{b}_i$ and $\hat{\mathbf{V}}_i$ is the posterior covariance of $\mathbf{b}_i$.

This follows most directly from Equations 56-59 in Bates. Demidenko's equations 3.51 and 3.52 are also useful because they give the BLUP/conditional-mean expression for the random effects in the conventional MLM form. The conditional distribution equation above is simply the subject-specific/blockwise version since we are focused on the subject-specific random effects/BLUPs.

In practical `lme4` terms:

- `ranef(fit, condVar = TRUE)[["id"]]` gives the posterior means $\hat{\mathbf{m}}_i$
- `attr(..., "postVar")` gives the posterior covariance blocks $\hat{\mathbf{V}}_i$
- `VarCorr(fit)$id` gives the prior covariance estimate $\hat{\mathbf{G}}$

Those are the three objects the current correction code is using.

## How the Existing Code Maps to the MLM Objects

The current correction function

$$
\texttt{unweight\_random\_effects(post\_mean, post\_vcov, prior\_mean, prior\_vcov)}
$$

is working with the following code-math mapping:

$$
\texttt{post\_mean} \leftrightarrow \hat{\mathbf{m}}_i,
$$

$$
\texttt{post\_vcov} \leftrightarrow \hat{\mathbf{V}}_i,
$$

$$
\texttt{prior\_mean} \leftrightarrow \hat{\boldsymbol{\mu}}_b,
$$

$$
\texttt{prior\_vcov} \leftrightarrow \hat{\mathbf{G}}.
$$

Now that we have these mappings, we have the ingredients to express the score correction for MLMs in the same algebraic form as Equation 32a! This is basically a plug-in empirical-Bayes case.

$$
\tilde{\mathbf{b}}_i
=
\left(
\hat{\mathbf{V}}_i^{-1} - \hat{\mathbf{G}}^{-1}
\right)^{-1}
\left(
\hat{\mathbf{V}}_i^{-1}\hat{\mathbf{m}}_i
- \hat{\mathbf{G}}^{-1}\hat{\boldsymbol{\mu}}_b
\right).
$$

This is the matrix form that matters for general/random-slope settings (i.e., not the scalar case of a random intercept-only model). It shows why the full intercept-slope covariance has to be carried along:

- the posterior covariance is a matrix, not two independent variances
- the prior covariance is a matrix, not two independent variances
- the corrected random slope is properly defined only after correcting the full vector

The equation above maps onto these two key pieces of code in the correction script:

```R
re_df <- ranef(fit_null, condVar = TRUE)[["id"]]
post_var_arr <- attr(re_df, "postVar")
prior_vcov <- as.matrix(VarCorr(fit_null)$id)
prior_mean <- c(0, 0)
re_names <- colnames(re_df)
...
solve(
      solve(post_vcov) - solve(prior_vcov),
      solve(post_vcov, post_mean) - solve(prior_vcov, prior_mean)
)
```



## Thinking about prior and posterior covariances in the frequentist MLM

The most important conceptual distinction is:

- **prior covariance:** the across-subject covariance of the 'G-side' random effects, $\hat{\mathbf{G}}$
- **posterior covariance:** the within-subject conditional uncertainty about each subject's random effects, $\hat{\mathbf{V}}_i$

### Prior Covariance

Echoing your treatment of hierarchical Bayesian models in the master's, the prior covariance tells you the distribution the BLUP was shrunk toward:

$$
\mathbf{b}_i \sim \mathcal{N}(\boldsymbol{\mu}_b, \mathbf{G}).
$$

In the current `lmer` setup, this is the matrix returned by `VarCorr(fit_null)$id`.

### Posterior Covariance

The posterior covariance tells you how uncertain the fitted model is about subject $i$ after seeing $\mathbf{y}_i$:

$$
\mathbf{b}_i \mid \mathbf{y}_i, \hat{\psi}
\sim
\mathcal{N}(\hat{\mathbf{m}}_i, \hat{\mathbf{V}}_i).
$$

In the current code, this is the subject-specific block from

$$
\texttt{attr(ranef(fit\_null, condVar = TRUE)[["id"]], "postVar")[,,i]}.
$$

### Intuition

The score correction relies on the interaction between these two covariance matrices. The correction is largely based on how the posterior precision differs from the prior precision, matching how you portrayed the problem in the master's.
