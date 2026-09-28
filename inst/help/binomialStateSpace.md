# Binomial State Space Models

Estimate how a success probability changes over time from the number of successes and the corresponding number of trials. The probability follows a local-level state space model on the log-odds scale. Optional continuous covariates can explain differences in probability with additive coefficients that are constant over time.

## Inputs and assumptions

- **Successes:** A numeric column of nonnegative whole counts. Each count must be no greater than the corresponding number of trials.
- **Number of trials:** Choose a numeric **Variable** or a **Constant** used for every row. Trials are the observed opportunities for success, not a population size. The number of trials may vary between rows. Binary data use successes of zero or one and a constant of one trial.
- **Time (optional):** A numeric column with unique, equally spaced values. Rows are sorted by this column before fitting. Without a time column, row order defines consecutive, equally spaced time points. The state-noise parameter refers to one row-to-row time step, not one unit of the numeric time variable.
- **Covariates:** Leave empty for the local-level model alone, or assign one or more numeric covariates. Each covariate must vary across time and is separately centered and divided by its sample standard deviation before fitting. Effects are additive; no interactions or additional regression intercept are included.

Assign each covariate only once. Covariates that are constant, non-finite, or linearly dependent after centering are rejected before sampling; variables are not silently dropped. Strongly correlated but linearly independent covariates are allowed. Their individual effects may be difficult to distinguish, and smoothly changing covariates may also be difficult to distinguish from the latent random-walk level. Inspect uncertainty and sampling diagnostics; the coefficients alone do not establish causal effects.

At least two time points are required. Assigned inputs must be finite and cannot contain missing values. Missing rows are not silently dropped. In variable-trials mode, keep an unobserved time step as a row with zero trials and zero successes; time and all covariates must still be supplied. At least one row must have positive trials.

Within a time point, the binomial likelihood assumes conditionally independent trials with a common success probability. Given the latent probabilities, observations at different time points are conditionally independent. Dependence over time is represented by the changing latent level. This model does not include a separate extra-binomial observation component, seasonality, a separate slope, categorical predictors, or time-varying regression coefficients.

## Model and priors

For time point `t`, the observed successes `y[t]` follow a binomial distribution with observed trials `n[t]` and probability `p[t]`:

```text
y[t] ~ Binomial(n[t], p[t])
level[1] = logit(initialProbability)
level[t] = level[t - 1] + innovation[t], for t > 1
innovation[t] ~ Normal(0, sigma^2)
```

Without covariates, `logit(p[t]) = level[t]`. With `K` covariates:

```text
logit(p[t]) = level[t] + sum(beta[j] * standardizedCovariate[t, j], j = 1, ..., K)
```

Each `beta[j]` is constant over time. Innovations are independent given the state-noise standard deviation `sigma`. Sampling uses Stan's Hamiltonian Monte Carlo with a non-centered representation of the innovations.

The **Model** section specifies:

- **Initial Probability:** A Beta prior with positive alpha and beta parameters, both one by default (a uniform prior). With covariates, this prior concerns the initial baseline probability when all covariates are at their respective means, not necessarily the probability at their first observed values.
- **State Noise:** A half-normal prior on `sigma`, with scale one by default. This is the standard deviation of changes in the latent log-odds level. Larger prior scales allow larger changes between consecutive time points.
- **Covariate Coefficients:** Choose a zero-centered Normal, Student t, or Cauchy prior, shown when covariates are assigned. The selected distribution and its settings apply independently to every coefficient; this is not a variable-selection or hierarchical prior. Each distribution has its own settings, enabled only when it is selected. The default is Student t with scale 0.5 and three degrees of freedom. Each coefficient is the log-odds change for a one-standard-deviation increase in its covariate, holding the other covariates and latent level fixed.
    - **Normal:** Set its positive **SD**, default 1.
    - **Student t:** Set its positive **Scale**, default 0.5, and positive **df**, default three. Fractional degrees of freedom are allowed. Smaller df give heavier tails; df = 1 gives a Cauchy distribution with the same scale. For df greater than two, the standard deviation is `scale * sqrt(df / (df - 2))`. For df of two or less, the prior has no finite variance; for df of one or less, its mean is undefined even though it is centered at zero.
    - **Cauchy:** Set its positive **Scale**, default 0.707. This is not a standard deviation: the Cauchy prior has no finite mean or variance.

All numeric prior settings must be strictly positive. Prior choices affect inference, particularly with sparse data, low trial counts, and unobserved time steps.

## Output

### Posterior Summary

The table reports the posterior mean, standard deviation, median, and 95% equal-tail credible interval for:

- the success probability at the final supplied time point;
- the state-noise standard deviation on the log-odds scale;
- each standardized covariate coefficient, labeled with its covariate name, when covariates are present.

The table reports probability on the zero-to-one scale. Probability plots use percentages. Estimates condition on **all supplied observations**: the trajectory contains smoothed probabilities, not sequential estimates based only on observations available at each time point. With covariates, probability summaries and plots include all their contributions.

The final probability refers to the final supplied row, including a final row with zero trials. It is not a forecast beyond the supplied series. Forecasting, future success-count prediction, acceptance decisions, and Bayes factors are not provided.

### Plots

- **Probability trajectory:** Posterior mean probability at each time point with a 95% equal-tail credible ribbon. Enabled by default.
- **Posterior distribution of final probability:** Density of the probability at the final supplied time point.
- **Posterior distributions of covariate coefficients:** A separate, named density plot for each time-constant standardized coefficient, with zero marked for reference.

### MCMC Diagnostics

Enable the **Overview table** to display posterior means, standard deviations, effective sample sizes (ESS), and R-hat for the selected parameters. State noise is preselected; probabilities at all time points and all covariate coefficients can also be selected. Each coefficient has a row labeled with its covariate name. Multiple chains are needed to assess between-chain convergence.

Review sampling diagnostics before interpreting results. R-hat values away from one, low effective sample sizes, or divergent transitions can indicate unreliable estimates. Summary footnotes flag divergent transitions and iterations that reach the maximum tree depth. More iterations alone do not necessarily resolve sampling problems.

## Advanced settings

- **Warmup iterations per chain:** Iterations used for adaptation and discarded; default 2000.
- **Sampling iterations per chain:** Additional post-warmup iterations before thinning; default 5000. At least two retained draws per chain are required.
- **Chains:** Independent chains; default three.
- **Thinning:** Retain every specified sampling iteration; default one keeps all draws.
- **Seed:** Default one. Reproducibility depends on the software environment as well as the seed.
- **Target acceptance probability:** Default 0.97. Increasing this can reduce divergent transitions but can slow sampling.
- **Maximum tree depth:** Default 15. Increasing this permits more computation per iteration and can substantially increase computation time. These settings do not guarantee divergence-free sampling; always review diagnostics.

## R packages

- `rstan`: posterior sampling and sampler summaries.
- `ggplot2` and `jaspGraphs`: graphical output.
