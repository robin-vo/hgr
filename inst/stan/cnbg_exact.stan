data {
  int<lower=1> N;
  vector<lower=0>[N] x;
  vector<lower=0>[N] nu;
  int<lower=1> I;
  array[N] int<lower=1, upper=I> row;
  real<lower=0> alpha;
  int<lower=1> n_max;
  real<lower=0> prior_shape;
  real<lower=0> prior_rate;
  real<lower=0> A_prior_shape;
  real<lower=0> A_prior_rate;
}
transformed data {
  vector[N] is_zero;
  for (k in 1:N)
    is_zero[k] = (x[k] <= 0) ? 1 : 0;
}
parameters {
  real<lower=0> kappa;
  real<lower=0> A;
  vector<lower=0>[I] U;
}
transformed parameters {
  real Ybar     = A * alpha / (alpha + 1);
  real sev_rate = (alpha + 1) / A;
}
model {
  kappa ~ gamma(prior_shape, prior_rate);
  A     ~ gamma(A_prior_shape, A_prior_rate);
  U     ~ gamma(kappa, kappa);
  for (k in 1:N) {
    real c = nu[k] * U[row[k]] / Ybar;
    if (is_zero[k] > 0.5) {
      target += -c;
    } else {
      vector[n_max] lp;
      for (n in 1:n_max)
        lp[n] = poisson_lpmf(n | c)
              + gamma_lpdf(x[k] | n * alpha, sev_rate);
      target += log_sum_exp(lp);
    }
  }
}
generated quantities {
  real Ybar_out = Ybar;
}
