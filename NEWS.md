# hgr 0.6.0

## CNBG module rebuilt as a port of the paper's replication library

The CNBG module now implements "What Paid Triangles Identify about Reserve
Uncertainty: Information Boundaries in Loss Reserving" (Van Oirbeek, 2026)
as that paper's scripts do. Fits obtained with `engine = "nuts"` and the
same seed are identical to the scripts' output; `engine = "grid"` (the
default) evaluates the same posterior deterministically.

* `simulate_cnbg()`: two frailties (`kappa_r`, `kappa_c`), calendar drift
  (`delta`, `onset`), latent counts returned; draws in the scripts' RNG order.
* `cnbg_stage1()`: the Stage-1 quasi-Poisson GLM with `quasi(link = "log",
  variance = "mu")`, the paper's convergence control, leverages, and anchored
  levels through an exposure offset.
* `fit_cnbg_wls()`: the second-moment estimator is now the single-pass
  non-negative least squares with weights `nu^-2`, the cap `kappa_B <=
  1000`, the Pearson fallback and the `(n - q) / n` factor; it replaces the
  iterated weighted regression of 0.5.0.
* `two_stage_target()`: the population target of that estimator under free
  levels, `B* + c_q / kappa_c` (Proposition 4.2 and Remark 4.5), and the
  target `A (1 - mean s_ij)` of the constrained intercept at the cap.
* `fit_cnbg_working()`: the working-likelihood posterior with the Gamma
  shapes set cell by cell from the fourth moment, leverage-corrected squared
  residuals, dropped zero-fit cells, the posterior mode, and the matrix
  curvature adjustment of Ribatet, Cooley and Davison (2012); `mean =
  "absorbed"` is the free-level fit of the paper. The joint `(kappa, A)` grid
  with Gamma shape 2 and the scalar "adjusted" engine of 0.5.0 are gone.
* `fit_cnbg()` and `reserve_cnbg()`: the modular predictive CNBG-Mod. `A`
  comes from the absorbed fit and `kappa_B` from the prior (floor 0.5), the
  mean parameters are drawn with the retransformation correction, the row
  frailty is restored by the conjugate update, and the future cells are
  moment-matched Gamma (`cells = "compound"` gives the compound law). In
  0.5.0 `reserve_cnbg()` sampled `(kappa, A)` jointly from the working
  posterior and drew compound cells; that construction is no longer offered.
* `bootstrap_cnbg_cond()` is now the paper's WLS-cond (one unconditioned
  frailty per accident year and replication, cells at the fitted means);
  the 0.5.0 version conditioned the frailty on the observed row.
* `bootstrap_cnbg_refit()`: the Frequentist CNBG refit bootstrap (new).
* `diagnose_kappa()`: contraction, shift, moment-matched KL and tier as in
  the paper's Schedule P diagnostics; regimes are `informative`,
  `not_informative` and `no_bayes`. The `incipient_flight` and
  `strong_flight` labels and the grid-exact KL are gone.
* `cnbg_information_kappa()`: Fisher and Godambe information, the
  coefficient of variation of the `nu^-2` WLS estimator and the
  heterogeneity-limit information. The van Trees columns are gone: the paper
  states no van Trees bound from the Godambe information.
* `fit_cnbg_exact()`: `levels = "free"` fits the free-level exact model
  (`inst/stan/cnbg_exact_free.stan`) with the lognormal level prior of the
  paper's Phase 17; `reserve_cnbg()` on an exact fit draws the compound
  law at the posterior draws.
* `prepare_triangle()`: the paper's data rules (negative increments to `NA`).
* `nnls` moved to Imports; `tweedie`, `statmod` and `R.utils` to Suggests.

## Other

* The object-less roxygen headers of `R/ucr.R` and `R/nbcl.R` are plain
  comments, so `roxygen2::roxygenise()` regenerates `NAMESPACE` cleanly.
* `tests/testthat` added for the CNBG module.
