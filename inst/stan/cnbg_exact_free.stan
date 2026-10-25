// cnbg_exact_free.stan -- Phase 4, Task 4.1 (exact-posterior, FREE accident-year levels)
// Same model as clean_run/cnbg_exact.stan except that the expected row count mu_i is a free
// parameter with a vague prior; nu_ij = mu_i * w_j * Ybar with the development weights w_j
// supplied as data at their true values and Ybar = A alpha / (alpha + 1).
data {
  int<lower=1> N;
  vector<lower=0>[N] x;
  vector<lower=0>[N] w;                 // true development weight of each observed cell
  int<lower=1> I;
  array[N] int<lower=1, upper=I> row;
  real<lower=0> alpha;
  int<lower=1> n_max;
  real<lower=0> prior_shape;
  real<lower=0> prior_rate;
  real<lower=0> A_prior_shape;
  real<lower=0> A_prior_rate;
  vector<lower=0>[I] mu0;               // rough row-count scale for the vague lognormal prior
  real<lower=0> mu_prior_sdlog;
}
transformed data {
  vector[N] is_zero;
  for (k in 1:N) is_zero[k] = (x[k] <= 0) ? 1 : 0;
}
parameters {
  real<lower=0> kappa;
  real<lower=0> A;
  vector<lower=0>[I] U;
  vector<lower=0>[I] mu;
}
transformed parameters {
  real Ybar = A * alpha / (alpha + 1);
  real sev_rate = (alpha + 1) / A;
}
model {
  kappa ~ gamma(prior_shape, prior_rate);
  A ~ gamma(A_prior_shape, A_prior_rate);
  U ~ gamma(kappa, kappa);
  mu ~ lognormal(log(mu0), mu_prior_sdlog);
  for (k in 1:N) {
    real c = mu[row[k]] * w[k] * U[row[k]];          // expected latent count
    if (is_zero[k] > 0.5) target += -c;
    else {
      vector[n_max] lp;
      for (n in 1:n_max) lp[n] = poisson_lpmf(n | c) + gamma_lpdf(x[k] | n * alpha, sev_rate);
      target += log_sum_exp(lp);
    }
  }
}
