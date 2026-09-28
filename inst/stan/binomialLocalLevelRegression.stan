// Adapted from jaspAcceptanceSampling: inst/stan/complex.stan.
// Source commit: 7c7183f96578f26b148cff5affce1b3e7df6c986.
// Lot-size data and lot-completion predictions have been removed.
// Full density constants are retained for future marginal-likelihood estimation.

data {
  int<lower=1> T;                  // number of time points
  array[T] int<lower=0> n;         // number of observed trials
  array[T] int<lower=0> y;         // observed successes
  int<lower=1> K;                 // number of covariates
  matrix[T, K] predictor;         // each column centered/scaled before fitting

  // Hyperparameters for priors
  real<lower=0> prior_theta1_a;    // Beta(alpha) for theta1
  real<lower=0> prior_theta1_b;    // Beta(beta) for theta1
  real<lower=0> prior_sigma_sd;    // SD for (half-)normal prior on sigma_eta
  int<lower=1, upper=3> prior_beta_distribution; // 1: Normal, 2: Student t, 3: Cauchy
  real<lower=0> prior_beta_scale;  // scale for beta (SD for Normal)
  real<lower=0> prior_beta_df;     // degrees of freedom for Student t
}

parameters {
  // Initial baseline success probability when the predictor equals zero
  real<lower=0, upper=1> theta1;

  // Scale of time-step innovations (local level)
  real<lower=0> sigma_eta;

  // Non-centered innovations
  vector[T-1] z_raw;

  // Independent predictor effects on log odds, constant over time
  vector[K] beta;
}

transformed parameters {
  vector[T] nu;        // latent baseline level on logit scale (no predictor)
  vector[T] mu;        // logit including predictor effect
  vector[T-1] eta;     // actual innovations
  vector[T] theta;     // success probability

  // Random-walk innovations
  eta = sigma_eta * z_raw;

  // Baseline local-level process (no predictor)
  nu[1] = logit(theta1);
  nu[2:T] = nu[1] + cumulative_sum(eta);

  // Add all predictor effects; the latent level supplies the baseline.
  mu = nu + predictor * beta;

  // Map to probability scale
  theta = inv_logit(mu);
}

model {
  target += beta_lpdf(theta1 | prior_theta1_a, prior_theta1_b);

  // Normalize the half-normal prior on the positive state-noise scale.
  target += normal_lpdf(sigma_eta | 0, prior_sigma_sd)
            - normal_lccdf(0 | 0, prior_sigma_sd);

  target += normal_lpdf(z_raw | 0, 1);
  if (prior_beta_distribution == 1) {
    target += normal_lpdf(beta | 0, prior_beta_scale);
  } else if (prior_beta_distribution == 2) {
    target += student_t_lpdf(beta | prior_beta_df, 0, prior_beta_scale);
  } else {
    target += cauchy_lpdf(beta | 0, prior_beta_scale);
  }

  target += binomial_lpmf(y | n, theta);
}
