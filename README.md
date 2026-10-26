# hgr

An R package implementing the Grammar of Reserving, a generative construction for macro-level non-life claims reserving, currently including:

- the Negative Binomial Chain-Ladder (NB-CL) model
- a Unified Credibility Reserving (UCR) framework based on the insight that **classical reserving methods are credibility estimators**
- **model-agnostic conditional predictive intervals** via Dirichlet-Multinomial allocation, with a portfolio-specific operability bound for when the framework can report at all
- the Compound Negative Binomial-Gamma (CNBG) **amount-triangle** model of *What Paid Triangles Identify about Reserve Uncertainty*, as a port of that paper's replication library: the information boundaries of a paid triangle, the modular predictive CNBG-Mod, and the diagnostic for when paid amounts alone can be trusted

## Multinomial Parametric Bootstrap

Conditional predictive intervals for claims reserving. The observed triangle is held fixed and only the allocation proportion is simulated, so the interval answers "how uncertain is the reserve on **this** triangle" rather than "how would the estimator vary across triangles I did not observe".

- **Method-agnostic**: Supply cumulative development proportions from Chain-Ladder, Bornhuetter-Ferguson, Cape Cod, GLMs or expert judgment — the framework adds the uncertainty quantification. For BF and Cape Cod, analytic prediction errors exist in the literature but no conditional predictive bootstrap did
- **One dispersion parameter**: Given the development pattern, allocation uncertainty is governed by a single Dirichlet concentration `c`, estimated from the triangle by a partial-column moment method
- **Two diagnostics, not one threshold**: `diagnose_c()` reports the portfolio-specific operability bound `c†(τ) = τ/π₀` and, separately, the heterogeneity screen with its measured error rates. They answer different questions and a portfolio can pass one while failing the other
- **Tiered moment reporting**: `(1-W)/W` is Beta-prime, whose k-th moment is finite only when `c·F > k`. The package enforces this and returns `NA` for the mean and standard error where they do not exist, rather than a number with no population counterpart
- **Joint count-amount**: The same structure handles claim counts (Multinomial) and claim amounts (Dirichlet)
- **Fully parametric**: Drawn from the Gamma-Dirichlet generative model. Note that a residual bootstrap with free accident-year effects reproduces the diffuse-limit posterior of the row level and is a correctly targeted benchmark under that hierarchy; the distinction that survives is the informative anchor, which a residual bootstrap has no channel to carry

## NB-CL

The NB-CL model provides a full likelihood framework for claim count triangles, addressing limitations of classical Chain-Ladder methods:

- **Full likelihood**: Enables likelihood-based inference, AIC/BIC comparison, and LR tests
- **Overdispersion**: Models variance through dispersion parameter κ with structural interpretation
- **Bias correction**: REML-like correction for finite-sample bias in κ estimation
- **Prediction intervals**: Parametric bootstrap incorporating both process and parameter uncertainty

## UCR

Unified Credibility Reserving (UCR) provides a data-adaptive framework that nests all classical reserving methods as special cases:

- **Unification**: Chain-Ladder, Bornhuetter-Ferguson and Cape Cod are exact posterior means under a single Poisson-Gamma-Multinomial-Dirichlet model, differing only in hyperparameter choice; Mack's distribution-free variance formula is recovered in a restricted special case
- **Adaptive weights**: Estimates between-year heterogeneity τ² from data, by marginal maximum likelihood on the exposure-anchored Negative Binomial marginal of the row totals (`method = "mml"`, which removes the negative-estimate truncation of Bühlmann-Straub) or by Bühlmann-Straub (`method = "bs"`), selecting the blend of individual vs pooled information
- **Method selection**: When τ² is large → UCR ≈ Chain-Ladder; when τ² ≈ 0 → UCR ≈ Cape Cod; with external prior → UCR generalises Bornhuetter-Ferguson
- **Accuracy, honestly stated**: UCR is not offered as more accurate. In simulation it tracks Chain-Ladder on the median across the heterogeneity range; its durable advantage is the damping of Chain-Ladder's extrapolation tail (the mean-MSE reduction reaches ~38% only at exact homogeneity, where UCR leans toward Cape Cod). Its value is diagnostic and protective: the fitted weights make the degree of pooling auditable
- **Diagnostics**: Credibility weights Z and estimated τ² provide interpretable diagnostics for method appropriateness

## CNBG

The Compound Negative Binomial-Gamma (CNBG) model is the distributional model for **incremental paid-loss** triangles of "What Paid Triangles Identify about Reserve Uncertainty: Information Boundaries in Loss Reserving" (Van Oirbeek, 2026). Its variance function `Var(X) = A nu + nu^2 / kappa_B` adds to the linear severity term the quadratic heterogeneity term that the ODP and Tweedie bootstraps omit, with `1 + 1/kappa_B = (1 + 1/kappa_r)(1 + 1/kappa_c)` for a row frailty and a cell shock. The module is a port of the library behind the paper's tables (`working_fit.R`, `predictive.R`, `cell_shapes.R`, `curvature_adjustment.R`, `oracle.R`, `legacy_methods.R` of the replication material), with the same conventions: with `engine = "nuts"` and the same seed, a fit here is the fit the paper reports.

- **What the triangle identifies**: with free accident-year levels the row frailty `kappa_r` is absorbed by the levels and is not identified from paid amounts (Theorem 1 of the paper, a parameterisation fact); the cell shock `kappa_c` survives in the residual curvature only below a threshold of about 10 to 20 (Remark 4.5). The second-moment estimator `fit_cnbg_wls()` therefore sits at its cap on most correctly specified triangles; `two_stage_target()` gives its population target `B* + c_q / kappa_c`
- **The working-likelihood posterior** `fit_cnbg_working()`: squared Stage-1 residuals, leverage corrected, with Gamma working shapes set cell by cell from the fourth moment, the posterior mode, and the curvature adjustment of Ribatet, Cooley and Davison (2012) in matrix form. `mean = "absorbed"` is the free-level fit that supplies `A`; `mean = "full"` carries `nu^2 / kappa_B` and feeds the diagnostic. `engine = "grid"` evaluates the posterior deterministically (no Stan); `engine = "nuts"` is the sampler the paper ran (`rstan`)
- **The modular predictive CNBG-Mod** `fit_cnbg()` + `reserve_cnbg()`: `A` from the absorbed fit, `kappa_B` from the prior (reference `Gamma(2, 0.4)` or the sensitivity prior `Gamma(1.78, 0.22)`), the mean parameters drawn with the retransformation correction, the row frailty restored by the conjugate update given the observed row, and moment-matched Gamma cells. The observed triangle is held fixed throughout
- **Point-estimate procedures**: `bootstrap_cnbg_cond()` is the paper's WLS-cond (one unconditioned frailty per accident year and replication at the second-moment estimate); `bootstrap_cnbg_refit()` is the Frequentist CNBG refit bootstrap, which under free levels reproduces to first order the conditional predictive variance under a diffuse prior on levels and pattern
- **Trust diagnostic, not a method chooser**: `diagnose_kappa()` returns the contraction `c` (posterior to prior 95% interval width), the shift `s`, the moment-matched KL divergence and its tier (`prior-driven`, `weak`, `moderate`, `strong`) from the full-mean fit, with the regime `informative` (`c < 1`), `not_informative` or `no_bayes` (no usable posterior). On most paid triangles the verdict is not informative, which is what the identification results predict; the next step is then the count triangle (NB-CL), not another amount-only method
- **Information on `kappa_B`**: `cnbg_information_kappa()` gives the Fisher and Godambe information of the squared-residual equation at a calibration (Proposition 6.1), the coefficient of variation of the `nu^-2` WLS estimator and the heterogeneity-limit information; at the package defaults it reproduces the `G_avg` column of the paper's Table godambe. The paper states no van Trees bound from `G`
- **The exact posterior** `fit_cnbg_exact()` (requires `rstan`): the frailty-augmented posterior by HMC with the cell means anchored at the Stage-1 fit (the paper's Schedule P comparison) or with free accident-year levels under a vague lognormal prior (the paper's Section S5 and its Phase 17 sensitivity). It costs roughly 200 times the working fit and is the per-portfolio reference, not the production estimator
- **Data and design**: `prepare_triangle()` applies the paper's Schedule P rules (negative increments to `NA`); `simulate_cnbg()` draws from the model at the paper's calibration with `kappa_r`, `kappa_c` and the calendar drift of Phase 16, in the same RNG order as the paper's scripts

The model is appropriate for paid (claim-emergence) triangles only, not cell-level incurred-loss triangles, which mix payments with case-reserve and IBNR revisions.

## Installation

```r
devtools::install_github("robin-vo/hgr")
```

The exact CNBG engine `fit_cnbg_exact()` and the `engine = "nuts"` option of the working fit require `rstan` (in Suggests); the default grid engine of `fit_cnbg()` does not. `bootstrap_cnbg_refit()` requires `tweedie`.

```r
install.packages("rstan")   # only for fit_cnbg_exact() and engine = "nuts"
```

## Quick Start — Multinomial Bootstrap

```r
library(hgr)
library(ChainLadder)

incr <- cum2incr(GenIns)

# 1. Diagnostics FIRST: operability bound and heterogeneity screen.
#    The bound is portfolio-specific; there is no universal threshold.
diagnose_c(incr)

# 2. Development proportions (shown for transparency; estimated
#    automatically by the bootstrap if not supplied)
dev <- estimate_dev_proportions(incr)
dev$pi_hat
dev$F_hat

# 3. Concentration parameter
c_hat <- estimate_c(incr, pi_hat = dev$pi_hat)
c_hat

# 4. Conditional predictive intervals
boot <- multinomial_bootstrap(incr, B = 10000)
print(boot)

#    Per-accident-year detail, including the moment tier of each year.
#    reserve_mean and reserve_se are NA where the moment does not exist.
boot$by_origin

# 5. Any development proportions: Bornhuetter-Ferguson, Cape Cod, a GLM,
#    or expert judgment. Note that calibration under an informative anchor
#    depends on the anchor being right (see "What this does not do").
earned_premium <- c(10e6, 11e6, 12e6, 13e6, 14e6, 15e6, 16e6, 17e6, 18e6, 19e6)
elr <- 0.65
bf_ultimate <- earned_premium * elr
bf_pi <- colSums(incr, na.rm = TRUE) / sum(bf_ultimate)
bf_pi <- bf_pi / sum(bf_pi)

boot_bf <- multinomial_bootstrap(incr, pi_hat = bf_pi, B = 10000)
print(boot_bf)

# 6. Regularised version: integrates over the posterior of c rather than
#    plugging it in. Use where diagnose_c()$operable is FALSE, or I < 10.
#    Returns quantiles only -- no moments exist under the posterior mixture.
reg <- bayesian_bootstrap(incr, B = 5000)
reg$quantiles
reg$c_posterior_mean

# 7. Closed-form variance. Retained for speed; not recommended over the
#    bootstrap at any c_hat.
delta_method_var(incr)
```
## Quick Start — NB-CL

```r
library(hgr)
library(ChainLadder)

# Convert cumulative triangle to incremental (NB-CL works on incremental data)
incr <- cum2incr(GenIns)

# Fit Negative Binomial Chain-Ladder
fit <- fit_nbcl(incr)
print(fit)          # model summary, kappa correction, AIC/BIC
summary(fit)        # coefficients and diagnostics

# Deterministic reserve estimate (mean prediction)
reserve_nbcl(fit)

# Predict lower-triangle means
pred <- predict_nbcl(fit)
head(pred)

# Parametric bootstrap for full predictive distribution
boot <- bootstrap_nbcl(fit, B = 5000, correct_kappa = TRUE)

# Prediction interval for the total reserve
predict_interval(boot, level = 0.95)

# Distribution of accident-year reserves
boot$reserves_by_ay[1:5, ]

# Inspect bootstrap dispersion estimates
hist(boot$kappas, main = "Bootstrap kappa", xlab = "kappa")

# Diagnostics: residuals, AY/DY patterns, profile likelihood for kappa
plot_diagnostics(fit, which = 1:4)

# Profile likelihood for kappa (standalone)
prof <- profile_kappa(fit)
plot(prof)
```

## Quick Start — UCR

```r
library(hgr)
library(ChainLadder)

# Load a standard cumulative triangle
triangle <- GenIns   # UCR works on cumulative data

# Optional: choose an exposure measure (default = first development period)
exposure <- triangle[, 1]

# 1. Fit the Unified Credibility Reserving model
ucr_fit <- ucr(triangle, exposure = exposure)

# 1b. Dispersion-corrected UCR for real amount data (Sec. 4.8):
#     required on premium-scaled losses, where the Poisson convention
#     collapses UCR onto Chain-Ladder
ucr_fit <- ucr(triangle, exposure = exposure, dispersion = "odp")

# 1c. MML estimator of tau^2 instead of Bühlmann-Straub (Sec. 4.6):
ucr_fit <- ucr(triangle, exposure = exposure, method = "mml")

# 2. Inspect credibility structure and reserve estimates
print(ucr_fit)

# 3. Visual diagnostics (rates, credibility weights, reserve comparison)
plot(ucr_fit, which = 1:3)

# 4. Compare classical methods and UCR
comparison <- compare_reserves(triangle, exposure = exposure)
print(comparison)

# 5. Extract credibility weights directly
ucr_fit$Z1        # accident-year credibility
ucr_fit$Lambda_CL # individual CL rates
ucr_fit$Lambda_UCR# credibility-adjusted rates

# 6. Access variance components
ucr_fit$sigma_sq   # process variance
ucr_fit$tau_sq_hat # between-year variance
ucr_fit$k          # sigma^2 / tau^2

# 7. Mack comparison (distribution-free benchmark)
mack_fit <- fit_mack(triangle)
print(mack_fit)
```

## Quick Start — CNBG

```r
library(hgr)

# Simulate an incremental paid-loss triangle at the paper's calibration
# (Y_bar = 5000, mu_i = 30 e^{0.05(i-1)}, explicit decreasing development
# weights; kappa_c = Inf is the single-frailty main design)
sim <- simulate_cnbg(I = 10, J = 10, kappa_r = 5, kappa_c = Inf, alpha = 3)
tri <- sim$triangle

# 1. Stage 1 (free accident-year levels) and the second-moment estimator,
#    which sits at the cap on most correctly specified triangles
ms  <- cnbg_stage1(tri)
fit_cnbg_wls(ms)
two_stage_target(kappa_r = c(3, 5, 10, 20, 50), kappa_c = 5)   # its free-level target

# 2. The full pipeline: absorbed working fit for A, full-mean fit for the
#    diagnostic. engine = "grid" needs no Stan; engine = "nuts" is the
#    sampler the paper ran.
fit <- fit_cnbg(tri, prior = "ref", alpha = 3, engine = "grid")
print(fit)

# 3. Trust diagnostic: informative | not_informative | no_bayes, with the
#    contraction c, the shift s and the KL tier
diagnose_kappa(fit)

# 4. CNBG-Mod predictive reserve (A from the absorbed fit, kappa_B from the
#    prior, conjugate frailty update, moment-matched Gamma cells)
reserve_cnbg(fit, B = 4000)
reserve_cnbg(fit, B = 4000, prior = "empb")   # the sensitivity prior

# 5. The point-estimate procedures of the paper
bootstrap_cnbg_cond(tri, B = 4000)             # WLS-cond
## bootstrap_cnbg_refit(tri, B = 500)          # Frequentist CNBG (needs tweedie)

# 6. Information on kappa_B at the calibration (the paper's Table godambe)
cnbg_information_kappa(kappa_grid = c(3, 5, 10, 20, 50), alpha = 3)

# 7. Real data: apply the paper's data rules first
## ta  <- prepare_triangle(taylor_ashe)         # negatives to NA
## fit <- fit_cnbg(ta$triangle)

# 8. Exact frailty-augmented posterior (requires rstan; minutes per fit)
## fx <- fit_cnbg_exact(tri, levels = "anchored")   # or levels = "free"
## reserve_cnbg(fx); diagnose_kappa(fx)
```

## Key Functions

### Multinomial Bootstrap

| Function                     | Description                                              |
| ---------------------------- | -------------------------------------------------------- |
| `multinomial_bootstrap()`    | Conditional predictive intervals, tiered moment reporting |
| `bayesian_bootstrap()`       | Regularised variant: `c` integrated over its posterior    |
| `diagnose_c()`               | Operability bound `c†(τ)` and heterogeneity screen         |
| `estimate_c()`               | Dirichlet concentration (subcomposition-corrected)        |
| `estimate_dev_proportions()` | Chain-Ladder development proportions                      |
| `delta_method_var()`         | Closed-form variance (not recommended over the bootstrap)  |

### Reserving Methods

| Function                 | Description                    |
| ------------------------ | ------------------------------ |
| `ucr()`                  | Unified Credibility Reserving  |
| `chain_ladder()`         | Chain-Ladder method            |
| `cape_cod()`             | Cape Cod method                |
| `bornhuetter_ferguson()` | Bornhuetter-Ferguson method    |
| `fit_mack()`             | Mack's distribution-free model |
| `compare_reserves()`     | Compare all methods            |

### NB-CL Model

| Function           | Description                               |
| ------------------ | ----------------------------------------- |
| `fit_nbcl()`       | Fit Negative Binomial Chain-Ladder        |
| `bootstrap_nbcl()` | Parametric bootstrap with bias correction |
| `profile_kappa()`  | Profile likelihood for dispersion         |

### CNBG Model

| Function                   | Description                                                                          |
| -------------------------- | ------------------------------------------------------------------------------------ |
| `simulate_cnbg()`          | CNBG triangle at the paper's calibration (`kappa_r`, `kappa_c`, calendar drift)        |
| `prepare_triangle()`       | The paper's data rules: incremental amounts, negatives to `NA`                        |
| `cnbg_stage1()`            | Stage-1 quasi-Poisson mean structure, free or anchored levels                         |
| `fit_cnbg_wls()`           | Two-stage second-moment estimator (`nu^-2` NNLS, cap 1000)                            |
| `two_stage_target()`       | Population target of that estimator under free levels (`B* + c_q / kappa_c`)          |
| `fit_cnbg_working()`       | Working-likelihood posterior: cell shapes, curvature adjustment; `grid` or `nuts`     |
| `fit_cnbg()`               | Pipeline: Stage 1, second-moment estimate, absorbed and full-mean working fits         |
| `reserve_cnbg()`           | CNBG-Mod predictive (`A` from the fit, `kappa_B` from the prior); exact predictive     |
| `predict_cnbg()`           | The CNBG-Mod predictive from a Stage-1 object and posterior draws of `A`              |
| `bootstrap_cnbg_cond()`    | WLS-cond                                                                               |
| `bootstrap_cnbg_refit()`   | Frequentist CNBG refit bootstrap                                                       |
| `diagnose_kappa()`         | Contraction, shift, KL tier and regime of the `kappa_B` posterior                      |
| `cnbg_information_kappa()` | Fisher and Godambe information on `kappa_B` (the paper's Table godambe)                |
| `fit_cnbg_exact()`         | Exact frailty-augmented HMC posterior, anchored or free levels (requires `rstan`)      |

## Replication material

The R scripts, seeds, stored result files and Schedule P inputs behind the tables of *What Paid Triangles Identify about Reserve Uncertainty* (about 70 MB) are not part of this repository; they are available from the author on request. The CNBG module of the package is a port of those scripts, so the procedures themselves can be run from the package.

## References

Van Oirbeek, R. and Verdonck, T. (2026). Conditional Prediction for Macro-Level Claims Reserving: What a Single Paid Triangle Identifies. *arXiv:2605.15896*.

Van Oirbeek, R. (2026). The Negative Binomial Chain-Ladder: A Full Likelihood Model for Claim Count Reserving. *CAS Forum*, 2026(1); CAS Bornhuetter Loss Reserving Award. Preprint *arXiv:2605.15811*.

Van Oirbeek, R. (2026). Classical Reserving Methods as Credibility Estimators: A Unified Bayesian Framework. *Working Paper*.

Van Oirbeek, R. (2026). What Paid Triangles Identify about Reserve Uncertainty: Information Boundaries in Loss Reserving. *Working Paper*.
