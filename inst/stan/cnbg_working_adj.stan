// cnbg_working_adj.stan -- Phase 2 (Task 2.1 and 2.3)
// Working likelihood z_ij ~ Gamma(s_ij, s_ij / mu_z_ij), mu_z = A nu + nu^2 / kappa,
// with a per-cell shape vector s_cell (legacy specification: s_cell = rep(2, n))
// and the Ribatet et al. (2012) curvature adjustment evaluated on the
// unconstrained scale phi = (log kappa, log A):
//   phi* = phi_hat + C (phi - phi_hat),  likelihood evaluated at exp(phi*).
// adjust = 0 reproduces the unadjusted working posterior (C ignored).
data {
  int<lower=1> n_obs;
  vector<lower=0>[n_obs] nu;
  vector<lower=0>[n_obs] z;
  vector<lower=0>[n_obs] s_cell;
  real<lower=0> kappa_a;
  real<lower=0> kappa_b;
  real<lower=0> A_a;
  real<lower=0> A_b;
  real<lower=0> kappa_lb;          // legacy lower bound 0.1
  int<lower=0, upper=1> adjust;
  vector[2] phi_hat;               // (log kappa_hat, log A_hat) at the posterior mode
  matrix[2, 2] C;                  // M^{-1} M_A
  int<lower=0, upper=1> absorbed;  // 1: mu_z = A nu only (free levels absorb the frailty; kappa = prior)
}
parameters {
  real<lower=kappa_lb> kappa;
  real<lower=0> A;
}
transformed parameters {
  vector[2] phi_star;
  {
    vector[2] phi = [log(kappa), log(A)]';
    if (adjust == 1) phi_star = phi_hat + C * (phi - phi_hat);
    else             phi_star = phi;
  }
}
model {
  vector[n_obs] mu_z = exp(phi_star[2]) * nu + (1 - absorbed) * square(nu) / exp(phi_star[1]);
  kappa ~ gamma(kappa_a, kappa_b);
  A     ~ gamma(A_a, A_b);
  target += gamma_lpdf(z | s_cell, s_cell ./ mu_z);
}
