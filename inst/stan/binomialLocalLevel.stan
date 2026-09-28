// Adapted from jaspAcceptanceSampling: inst/stan/simple.stan.
// Source commit: 7c7183f96578f26b148cff5affce1b3e7df6c986.
// Lot-size data and lot-completion predictions have been removed.
// Full density constants are retained for future marginal-likelihood estimation.

data {
  int<lower=1> T;             // number of time points
  array[T] int<lower=0> n;    // number of observed trials
  array[T] int<lower=0> y;    // observed successes

  // Hyperparameters for priors
  real<lower=0> prior_theta1_a;   // Beta(alpha) for theta1
  real<lower=0> prior_theta1_b;   // Beta(beta) for theta1
  real<lower=0> prior_sigma_sd;   // SD for (half-)normal prior on sigma_eta
}

parameters {
  // Initial success probability
  real<lower=0, upper=1> theta1;

  // Scale of time-step innovations (constrained > 0)
  real<lower=0> sigma_eta;

  // Non-centered innovations
  vector[T-1] z_raw;
}

transformed parameters {
  vector[T] mu;
  vector[T] theta;
  vector[T-1] eta;

  eta = sigma_eta * z_raw;

  mu[1] = logit(theta1);
  mu[2:T] = mu[1] + cumulative_sum(eta);

  theta = inv_logit(mu);
}

model {
  target += beta_lpdf(theta1 | prior_theta1_a, prior_theta1_b);

  // Normalize the half-normal prior on the positive state-noise scale.
  target += normal_lpdf(sigma_eta | 0, prior_sigma_sd)
            - normal_lccdf(0 | 0, prior_sigma_sd);

  target += normal_lpdf(z_raw | 0, 1);

  target += binomial_lpmf(y | n, theta);
}
