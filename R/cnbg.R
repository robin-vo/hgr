#' Compound Negative Binomial-Gamma (CNBG) reserving
#'
#' Implementation of the amount-triangle methodology of "What Paid Triangles
#' Can and Cannot Identify" (the CNBG paper). Provides the CNBG
#' data-generating process at the paper's calibration, two-stage estimation
#' of the structural parameters (kappa, A), the conditioning-respecting
#' posterior predictive reserve distribution, the kappa posterior diagnostic
#' with the paper's regime and KL-tier classification, and the
#' Fisher/Godambe/van Trees information bound (which reproduces the paper's
#' van Trees table exactly at its defaults).
#'
#' Two Bayesian engines are provided. `engine = "paper"` (the default) is a
#' deterministic joint grid quadrature of the paper's working-likelihood
#' posterior over (kappa, A): the naive independent-cell composite
#' likelihood, prior-anchored, with no curvature adjustment -- the estimator
#' behind the paper's simulation and Schedule P results (there fitted by
#' NUTS; the grid is the same target integrated deterministically).
#' `engine = "adjusted"` additionally tempers the kappa likelihood by the
#' cluster-robust Godambe magnitude adjustment of Ribatet, Cooley and
#' Davison (2012); the paper discusses this correction but does not apply
#' it, so the adjusted engine is a variant, not a reproduction. The exact
#' frailty-augmented posterior is available separately via
#' [fit_cnbg_exact()] (requires rstan, which is in Suggests -- the core
#' package carries no Stan dependency). The variance function is
#' V(nu) = A nu + nu^2 / kappa throughout.
#'
#' @keywords internal
#' @name cnbg
#' @importFrom stats glm quasipoisson predict coef vcov fitted lm qgamma dgamma
#'   rgamma rpois rnorm quantile integrate digamma lgamma sd median pgamma
#'   model.matrix approx
#' @importFrom MASS mvrnorm
NULL
## NOTE (K8): rstan is deliberately NOT imported. All rstan calls are
## namespace-qualified and gated by requireNamespace() at call time, so
## rstan belongs in Suggests in DESCRIPTION, matching the header claim.

## Default kappa priors (shape, rate) -----------------------------------------
.cnbg_prior_ref  <- c(shape = 2.00, rate = 0.40)   # weakly informative reference
.cnbg_prior_empb <- c(shape = 1.78, rate = 0.22)   # empirical Bayes (Sched. P)
.cnbg_kappa_cap  <- 1000                            # WLS boundary cap
## Saturation convention: with n = 55 observed cells and q = 19 mean
## parameters on a 10x10 triangle, the df adjustment maps the cap to
## 1000 * 36/55 = 654.5; "at cap" counts estimates at this boundary, and
## posterior edge mass is reported above kappa = 600, both matching the
## paper's convention.
.cnbg_edge_thresh <- 600

## Paper calibration (Section "Design" of the simulation study) ---------------
.cnbg_paper_w <- c(0.35, 0.25, 0.15, 0.10, 0.07, 0.04, 0.02, 0.01,
                   0.005, 0.005)


## ===========================================================================
## DATA-GENERATING PROCESS
## ===========================================================================

#' Simulate a CNBG paid-loss triangle
#'
#' Draws an incremental paid-loss triangle from the CNBG generative model:
#' accident-year frailty `U_i ~ Gamma(kappa, kappa)`, cell counts
#' `N_ij ~ Poisson(mu_i U_i w_j)`, Gamma severities, and `X_ij` their sum.
#'
#' Defaults reproduce the paper's simulation calibration exactly:
#' accident-year expected claim counts `mu_i = mu_base * exp(mu_growth (i-1))`
#' with `mu_base = 30` (roughly 30--47 expected claims per row), the explicit
#' decreasing development-weight vector
#' `(.35, .25, .15, .10, .07, .04, .02, .01, .005, .005)`, and mean severity
#' `Y_bar = 5000`, for an expected outstanding reserve of about \$355k per
#' triangle. Note the late development cells are sparse by design (~0.15
#' expected claims, ~86\% zero probability): the paper's information bounds
#' are derived for this regime. For a geometric weight pattern pass, e.g.,
#' `w = 0.7^(0:(J-1))`.
#'
#' @param I,J Number of accident and development years.
#' @param kappa Frailty dispersion (small = heterogeneous accident years).
#' @param alpha Gamma severity shape.
#' @param Y_bar Mean severity.
#' @param mu_base,mu_growth Accident-year expected claim count is
#'   `mu_base * exp(mu_growth (i - 1))`.
#' @param w Development weights, length `J` (normalised internally).
#' @return A list with `triangle` (incremental amounts, lower triangle `NA`),
#'   `full` (complete I x J matrix), `R_true` (realised outstanding reserve),
#'   and `U` (the drawn frailties).
#' @examples
#' set.seed(1)
#' tri <- simulate_cnbg(I = 10, J = 10, kappa = 5, alpha = 3)
#' tri$R_true
#' @export
simulate_cnbg <- function(I = 10L, J = 10L, kappa = 5, alpha = 3,
                          Y_bar = 5000, mu_base = 30, mu_growth = 0.05,
                          w = .cnbg_paper_w) {
  stopifnot(length(w) == J, all(w > 0))
  mu_i <- mu_base * exp(mu_growth * (seq_len(I) - 1L))
  w    <- w / sum(w)
  g_rate <- alpha / Y_bar

  U    <- rgamma(I, shape = kappa, rate = kappa)      # E[U]=1, Var=1/kappa
  full <- matrix(0, I, J)
  for (i in seq_len(I)) for (j in seq_len(J)) {
    n <- rpois(1L, mu_i[i] * U[i] * w[j])
    if (n > 0L) full[i, j] <- rgamma(1L, shape = n * alpha, rate = g_rate)
  }
  upper <- outer(seq_len(I), seq_len(J), function(i, j) (i + j - 1L) <= I)
  tri <- full; tri[!upper] <- NA_real_
  list(triangle = tri, full = full,
       R_true = sum(full[!upper]), U = U)
}


## ===========================================================================
## STAGE 1: MEAN STRUCTURE  (free-row Chain-Ladder quasi-Poisson GLM)
## ===========================================================================

.cnbg_long <- function(tri) {
  I <- nrow(tri); J <- ncol(tri)
  idx <- which(!is.na(tri), arr.ind = TRUE)
  data.frame(x  = tri[!is.na(tri)],
             ay = factor(idx[, 1], levels = seq_len(I)),
             dy = factor(idx[, 2], levels = seq_len(J)))
}

.cnbg_stage1 <- function(tri) {
  I <- nrow(tri); J <- ncol(tri)
  d <- .cnbg_long(tri)
  fit1 <- glm(x ~ ay + dy, data = d,
              family = quasipoisson(link = "log"))
  nu_obs <- fitted(fit1)
  nd <- expand.grid(ay = factor(seq_len(I), levels = seq_len(I)),
                    dy = factor(seq_len(J), levels = seq_len(J)))
  nu_mat <- matrix(predict(fit1, nd, type = "response"), I, J)
  Xd     <- model.matrix(~ ay + dy, nd)

  obs    <- !is.na(tri)
  X_obs_by_ay  <- tapply(d$x, d$ay, sum)
  nu_obs_by_ay <- vapply(seq_len(I), function(i) sum(nu_mat[i, obs[i, ]]),
                         numeric(1))
  list(fit1 = fit1, nu_mat = nu_mat, tri = tri, obs = obs,
       I = I, J = J, n_obs = nrow(d), q = length(coef(fit1)),
       X_design = Xd,
       sum_X_by_ay  = as.numeric(X_obs_by_ay),
       sum_nu_by_ay = nu_obs_by_ay)
}


## ===========================================================================
## STAGE 2a: WLS point estimate of (kappa, A)   (paper's weighted estimator)
## ===========================================================================
## The paper's estimator (Sec. "Estimation", Davidian--Carroll tradition):
## through-origin regression of the squared residual z = (x - nu)^2 on
## (nu, nu^2), with working-variance weights 1/mu_z^2 where
## mu_z = A nu + B nu^2, iterated to convergence of the weights. This is
## the estimator that returned identical (kappa, A) from two independent
## Schedule P implementations; the previous package version used an
## unweighted transformed regression and would not have matched them (K5).

.cnbg_wls <- function(ms, kappa_cap = .cnbg_kappa_cap, n_iter = 4L) {
  nu <- ms$nu_mat[ms$obs]; x <- ms$tri[ms$obs]
  z  <- (x - nu)^2
  nu2 <- nu^2

  ## iteration 0: unweighted through-origin fit for starting values
  cf <- coef(lm(z ~ 0 + nu + nu2))
  A_hat <- max(cf[[1]], .Machine$double.eps)
  B_hat <- max(cf[[2]], 1 / kappa_cap)

  for (it in seq_len(n_iter)) {
    mu_z <- A_hat * nu + B_hat * nu2
    wgt  <- 1 / pmax(mu_z, .Machine$double.eps)^2
    cf   <- coef(lm(z ~ 0 + nu + nu2, weights = wgt))
    A_new <- max(cf[[1]], .Machine$double.eps)
    B_new <- max(cf[[2]], 1 / kappa_cap)
    if (abs(log(A_new / A_hat)) < 1e-8 &&
        abs(log(B_new / B_hat)) < 1e-8) { A_hat <- A_new; B_hat <- B_new; break }
    A_hat <- A_new; B_hat <- B_new
  }
  if (A_hat <= .Machine$double.eps) A_hat <- mean(z / nu)  # severity fallback

  kappa_raw <- min(1 / B_hat, kappa_cap)
  kappa_adj <- kappa_raw * (ms$n_obs - ms$q) / ms$n_obs
  list(kappa = kappa_adj, A = A_hat,
       at_cap = (1 / B_hat) >= kappa_cap)
}


## ===========================================================================
## STAGE 2b: Bayesian posterior of (kappa, A)
## ===========================================================================

## Exact KL of a discrete grid posterior against the (same-grid) prior (K4).
.cnbg_grid_kl <- function(kg, p, prior_shape, prior_rate) {
  pr <- dgamma(kg, prior_shape, prior_rate) * kg      # log-grid measure
  pr <- pr / sum(pr)
  ok <- p > 0 & pr > 0
  sum(p[ok] * (log(p[ok]) - log(pr[ok])))
}

## ---- engine = "paper": joint (kappa, A) grid quadrature --------------------
## Target: the paper's working likelihood -- each observed cell contributes
## z = (x - nu)^2 ~ Gamma(2, 2 / mu_z), mu_z = A nu + nu^2 / kappa, cells
## treated as independent (the composite likelihood of the paper), anchored
## by the kappa prior and A ~ Gamma(2, 2 / A0) (the same weakly-informative
## A-prior convention as the exact engine). No curvature adjustment: this is
## the estimator whose posterior the paper's published fits integrate.
.cnbg_posterior_paper <- function(ms, prior_shape, prior_rate,
                                  kappa_cap = .cnbg_kappa_cap,
                                  n_k = 400L, n_A = 120L, n_draw = 4000L) {
  nu <- ms$nu_mat[ms$obs]; x <- ms$tri[ms$obs]
  z  <- (x - nu)^2
  A0 <- max(stats::median(z / nu), .Machine$double.eps)

  kg <- exp(seq(log(0.1), log(kappa_cap), length.out = n_k))
  Ag <- A0 * exp(seq(log(1 / 8), log(8), length.out = n_A))

  lp <- matrix(-Inf, n_k, n_A)
  nu2 <- nu^2
  lprior_k <- dgamma(kg, prior_shape, prior_rate, log = TRUE) + log(kg)
  lprior_A <- dgamma(Ag, 2, rate = 2 / A0, log = TRUE) + log(Ag)
  for (a in seq_len(n_A)) {
    ## vectorised over the kappa grid: cells x n_k matrix of mu_z
    mu_z <- outer(nu, rep(Ag[a], n_k)) + outer(nu2, 1 / kg)
    ll   <- colSums(dgamma(z, shape = 2, rate = 2 / mu_z, log = TRUE))
    lp[, a] <- ll + lprior_k + lprior_A[a]
  }
  lp[!is.finite(lp)] <- -Inf
  mx <- max(lp)
  P  <- if (is.finite(mx)) exp(lp - mx) else matrix(1, n_k, n_A)
  P[!is.finite(P)] <- 0
  if (sum(P) <= 0) P <- matrix(1, n_k, n_A)
  P <- P / sum(P)

  pk <- rowSums(P)                                    # kappa marginal
  cdf <- cumsum(pk)
  qk  <- function(q) approx(cdf, kg, xout = q, rule = 2, ties = min)$y
  kmean <- sum(pk * kg)

  ## joint draws preserve the (kappa, A) posterior dependence
  cell <- sample.int(n_k * n_A, n_draw, replace = TRUE, prob = as.numeric(P))
  ki <- ((cell - 1L) %% n_k) + 1L
  ai <- ((cell - 1L) %/% n_k) + 1L

  list(kappa_draws = kg[ki], A_draws = Ag[ai],
       kappa_mean = kmean, kappa_median = qk(0.5),
       kappa_q025 = qk(0.025), kappa_q975 = qk(0.975),
       kappa_sd = sqrt(sum(pk * (kg - kmean)^2)),
       A_emp = A0, A_hat = sum(colSums(P) * Ag),
       adjusted = FALSE, godambe_r = NA_real_, widen = NA_real_,
       KL_grid = .cnbg_grid_kl(kg, pk, prior_shape, prior_rate),
       edge_mass = sum(pk[kg >= .cnbg_edge_thresh]),
       engine = "paper")
}

## ---- engine = "adjusted": 1-D tempered variant (kept, clearly labelled) ----
## Fixes A at a robust estimate and tempers the 1-D kappa working likelihood
## by the cluster-robust Godambe magnitude adjustment r = I_naive/J_cluster
## (Ribatet, Cooley and Davison, 2012), computed from row-summed scores at
## the prior mean. The paper DISCUSSES this correction but does not apply
## it; use this engine when a deliberately deflated, cluster-honest kappa
## posterior is wanted, and engine = "paper" to reproduce the paper.
.cnbg_posterior_adjusted <- function(ms, prior_shape, prior_rate,
                                     kappa_cap = .cnbg_kappa_cap,
                                     n_grid = 600L, n_draw = 4000L) {
  nu <- ms$nu_mat[ms$obs]; x <- ms$tri[ms$obs]
  z  <- (x - nu)^2
  rows  <- which(ms$obs, arr.ind = TRUE)[, 1]         # accident-year cluster
  A_emp <- max(mean(z) / mean(nu), .Machine$double.eps)
  ## robust severity dispersion; the /50 floor guards degenerate medians
  ## on near-deterministic triangles (documented, m3)
  A_hat <- max(stats::median(z / nu), A_emp / 50)

  llk <- function(k) {
    mu_z <- A_hat * nu + nu^2 / k
    sum(dgamma(z, shape = 2, rate = 2 / mu_z, log = TRUE))
  }

  k0 <- prior_shape / prior_rate
  hh <- 1e-4 * k0
  sc <- (dgamma(z, 2, rate = 2 / (A_hat * nu + nu^2 / (k0 + hh)), log = TRUE) -
         dgamma(z, 2, rate = 2 / (A_hat * nu + nu^2 / (k0 - hh)), log = TRUE)) /
        (2 * hh)
  I_naive   <- sum(sc^2)
  J_cluster <- sum(tapply(sc, rows, sum)^2)
  r <- 1
  if (is.finite(I_naive) && is.finite(J_cluster) && J_cluster > 0)
    r <- min(max(I_naive / J_cluster, 0.02), 1)

  kg <- exp(seq(log(0.1), log(kappa_cap), length.out = n_grid))
  ll <- vapply(kg, llk, numeric(1))
  ll[!is.finite(ll)] <- -Inf
  lp <- r * (ll - max(ll[is.finite(ll)])) +
        dgamma(kg, shape = prior_shape, rate = prior_rate, log = TRUE) +
        log(kg)
  lp[!is.finite(lp)] <- -Inf
  mx <- suppressWarnings(max(lp))
  p  <- if (!is.finite(mx)) rep(1, n_grid) else exp(lp - mx)
  p[!is.finite(p)] <- 0; if (sum(p) <= 0) p <- rep(1, n_grid); p <- p / sum(p)

  cdf <- cumsum(p)
  qk  <- function(q) approx(cdf, kg, xout = q, rule = 2, ties = min)$y
  kmean <- sum(p * kg)
  di <- sample(kg, n_draw, replace = TRUE, prob = p)

  ## A uncertainty proxy: lognormal jitter with sd tied to the effective
  ## cell count (replaces the previous hardcoded 0.08, m3/K2): the
  ## working-model CV of a mean of n_obs Gamma(2) cells.
  A_sd <- sqrt(1 / (2 * ms$n_obs))
  A_draws <- A_hat * exp(rnorm(n_draw, -A_sd^2 / 2, A_sd))

  list(kappa_draws = di, A_draws = A_draws,
       kappa_mean = kmean, kappa_median = qk(0.5),
       kappa_q025 = qk(0.025), kappa_q975 = qk(0.975),
       kappa_sd = sqrt(sum(p * (kg - kmean)^2)),
       A_emp = A_emp, A_hat = A_hat,
       adjusted = TRUE, godambe_r = r, widen = sqrt(1 / r),
       KL_grid = .cnbg_grid_kl(kg, p, prior_shape, prior_rate),
       edge_mass = sum(p[kg >= .cnbg_edge_thresh]),
       engine = "adjusted")
}


## ===========================================================================
## FIT
## ===========================================================================

#' Fit the CNBG model to a paid-loss triangle
#'
#' Estimates the mean structure by a free-row quasi-Poisson GLM (Chain-Ladder
#' means) and the structural pair `(kappa, A)` by the Bayesian
#' working-likelihood posterior (`method = "bayes"`, recommended) or the
#' two-stage weighted least-squares point estimate (`method = "wls"`, the
#' paper's weighted Davidian--Carroll regression with df-adjusted cap).
#'
#' @param triangle Incremental paid-loss triangle (matrix); observed cells
#'   non-`NA`, future cells `NA`.
#' @param method `"bayes"` or `"wls"`.
#' @param engine For `method = "bayes"`: `"paper"` (default) is the joint
#'   `(kappa, A)` grid quadrature of the paper's working-likelihood
#'   posterior -- naive composite likelihood, prior-anchored, no curvature
#'   adjustment -- i.e. the estimator behind the paper's published fits.
#'   `"adjusted"` fixes `A` and tempers the kappa likelihood by the
#'   cluster-robust Godambe magnitude adjustment (Ribatet, Cooley and
#'   Davison, 2012), a deliberately deflated variant the paper discusses
#'   but does not apply.
#' @param prior `"ref"` (Gamma(2, 0.4)) or `"empb"` (Gamma(1.78, 0.22)), or a
#'   length-2 numeric `c(shape, rate)`.
#' @param kappa_cap WLS boundary cap on `kappa`.
#' @return An object of class `cnbg_fit`.
#'
#' @note The working likelihood treats the cells of an accident-year row as
#'   independent although they share a row frailty; the paper's information
#'   bound (see [cnbg_information_kappa()]) quantifies what this can and
#'   cannot identify. Amount triangles identify `kappa` only weakly, so the
#'   posterior is largely prior-governed under either engine -- treat point
#'   summaries with caution and lean on [diagnose_kappa()].
#' @examples
#' set.seed(1)
#' tri <- simulate_cnbg(kappa = 5, alpha = 3)$triangle
#' fit <- fit_cnbg(tri, method = "bayes")
#' fit
#' @export
fit_cnbg <- function(triangle, method = c("bayes", "wls"),
                     engine = c("paper", "adjusted"),
                     prior = "ref", kappa_cap = .cnbg_kappa_cap) {
  method <- match.arg(method)
  engine <- match.arg(engine)
  pr <- if (is.numeric(prior)) prior
        else switch(prior, ref = .cnbg_prior_ref, empb = .cnbg_prior_empb,
                    stop("prior must be 'ref', 'empb', or c(shape, rate)"))
  names(pr) <- c("shape", "rate")

  ms <- tryCatch(.cnbg_stage1(triangle), error = function(e) NULL)
  if (is.null(ms))
    return(structure(list(ok = FALSE, method = method,
                          message = "Stage-1 GLM failed"),
                     class = "cnbg_fit"))

  out <- list(ok = TRUE, method = method, ms = ms,
              prior = pr,
              prior_label = if (is.numeric(prior)) "custom" else prior,
              kappa_cap = kappa_cap)

  if (method == "wls") {
    w <- .cnbg_wls(ms, kappa_cap)
    out$kappa <- w$kappa; out$A <- w$A; out$at_cap <- w$at_cap
  } else {
    g <- if (engine == "paper")
      .cnbg_posterior_paper(ms, pr[["shape"]], pr[["rate"]], kappa_cap)
    else
      .cnbg_posterior_adjusted(ms, pr[["shape"]], pr[["rate"]], kappa_cap)
    out <- c(out, g)
    out$kappa <- g$kappa_median; out$A <- mean(g$A_draws)
  }
  structure(out, class = "cnbg_fit")
}

#' @export
print.cnbg_fit <- function(x, ...) {
  if (!isTRUE(x$ok)) { cat("CNBG fit FAILED:", x$message, "\n"); return(invisible(x)) }
  cat(sprintf("CNBG fit (%s)\n", x$method))
  if (x$method %in% c("bayes", "exact")) {
    cat(sprintf("  kappa: median %.2f, 95%% CI [%.2f, %.2f]   A (mean): %.4g\n",
                x$kappa_median, x$kappa_q025, x$kappa_q975, x$A))
    curv <- if (identical(x$engine, "exact-rstan"))
      sprintf("exact frailty-augmented (Rhat %.3f, n_div %d)",
              x$rhat_max, x$n_divergent)
    else if (identical(x$engine, "adjusted"))
      sprintf("Godambe-adjusted variant (kappa SD x%.1f)", x$widen)
    else "paper working-likelihood (joint grid)"
    cat(sprintf("  prior: %s   engine: %s", x$prior_label, curv))
    cat(sprintf("   edge mass (kappa >= %d): %.1f%%\n",
                .cnbg_edge_thresh, 100 * x$edge_mass))
  } else {
    cat(sprintf("  kappa (adj): %.2f%s   A: %.4g\n",
                x$kappa, if (isTRUE(x$at_cap)) " [at cap]" else "", x$A))
  }
  invisible(x)
}


## ===========================================================================
## PREDICTIVE ENGINE  (conditioning-respecting posterior predictive reserve)
## ===========================================================================
## Single engine used by both reserve_cnbg() (parameter draws from the
## posterior) and bootstrap_cnbg_cond() (point estimate replicated). They
## therefore differ ONLY in whether (kappa, A) carry posterior uncertainty,
## which is exactly the comparison the paper makes.
##
## Future cells are drawn from the EXACT conditional compound Poisson-Gamma
## (K6): given the frailty draw U_i, the expected claim count of cell (i,j)
## is nu_ij U_i / Y_bar with Y_bar = A alpha / (alpha + 1), the count is
## Poisson, and the amount is Gamma(n alpha, (alpha+1)/A). This carries the
## zero atom the previous Gamma moment-match lacked -- material at the
## paper's calibration, where late development cells are ~86% zero.

.cnbg_predict <- function(ms, kappa_v, A_v, B, alpha = 3,
                          propagate_mean = FALSE, level = 0.95) {
  I <- ms$I; J <- ms$J
  fut_mask <- outer(seq_len(I), seq_len(J), function(i, j) (i + j - 1L) > I)
  future <- which(fut_mask, arr.ind = TRUE)

  coef_hat <- coef(ms$fit1)
  V        <- vcov(ms$fit1)
  if (propagate_mean &&
      min(eigen(V, symmetric = TRUE, only.values = TRUE)$values) < 1e-10)
    V <- V + diag(1e-8, nrow(V))

  R <- numeric(B)
  for (b in seq_len(B)) {
    kap <- max(kappa_v[b], 0.1)                # aligned with grid floor (K6)
    A   <- max(A_v[b], 1e-8)
    Yb  <- A * alpha / (alpha + 1)             # implied mean severity
    sev_rate <- (alpha + 1) / A

    if (propagate_mean) {
      th <- MASS::mvrnorm(1, coef_hat, V)
      nu_b <- matrix(exp(ms$X_design %*% th), I, J)
      nu_b <- pmin(nu_b, 10 * pmax(ms$nu_mat, 1))
      snu  <- vapply(seq_len(I), function(i) sum(nu_b[i, ms$obs[i, ]]),
                     numeric(1))
    } else {
      nu_b <- ms$nu_mat; snu <- ms$sum_nu_by_ay
    }

    ## conjugate frailty posterior given the observed row totals
    U <- rgamma(I, shape = kap + ms$sum_X_by_ay / A,
                rate  = pmax(kap + snu / A, 1e-10))
    Rb <- 0
    for (r in seq_len(nrow(future))) {
      i <- future[r, 1]; j <- future[r, 2]; nu0 <- nu_b[i, j]
      if (nu0 > 0) {
        n <- rpois(1L, nu0 * U[i] / Yb)
        if (n > 0L) Rb <- Rb + rgamma(1L, shape = n * alpha, rate = sev_rate)
      }
    }
    R[b] <- Rb
  }
  a <- (1 - level) / 2
  list(reserve = R, mean = mean(R), sd = sd(R),
       pi = quantile(R, c(a, 1 - a)), level = level, alpha = alpha)
}


#' Posterior predictive reserve distribution (CNBG-Bayes)
#'
#' Conditioning-respecting posterior predictive reserve: the observed triangle
#' is held fixed, accident-year frailties are drawn from their conjugate
#' posterior, future cells are drawn from the exact conditional compound
#' Poisson--Gamma (with its zero atom), and parameter uncertainty in
#' `(kappa, A)` is integrated over the posterior. The severity shape `alpha`
#' is held fixed (paid triangles do not identify it), as in the exact engine.
#'
#' Setting `propagate_mean = TRUE` additionally propagates estimation
#' uncertainty in the Stage-1 mean coefficients by drawing them from their
#' asymptotic distribution; no data are resampled, so this remains distinct
#' from the refit-on-pseudo-triangle bootstrap. The default `FALSE` is the
#' strictly conditional predictive of the paper's Algorithm.
#'
#' @param fit A `cnbg_fit` from [fit_cnbg()] with `method = "bayes"` or from
#'   [fit_cnbg_exact()].
#' @param B Number of predictive draws.
#' @param alpha Severity shape, held fixed.
#' @param propagate_mean Propagate Stage-1 mean-coefficient uncertainty.
#' @param level Predictive interval level.
#' @return An object of class `cnbg_reserve`.
#' @examples
#' set.seed(1)
#' tri <- simulate_cnbg(kappa = 5, alpha = 3)$triangle
#' reserve_cnbg(fit_cnbg(tri, "bayes"))
#' @export
reserve_cnbg <- function(fit, B = 4000L, alpha = 3,
                         propagate_mean = FALSE, level = 0.95) {
  if (!isTRUE(fit$ok)) stop("fit failed; cannot predict")
  if (!fit$method %in% c("bayes", "exact"))
    stop("reserve_cnbg requires a Bayesian or exact fit; use bootstrap_cnbg_cond for WLS")
  if (!is.null(fit$alpha)) alpha <- fit$alpha   # exact fits carry alpha
  idx <- sample.int(length(fit$kappa_draws), B, replace = TRUE)
  res <- .cnbg_predict(fit$ms, fit$kappa_draws[idx], fit$A_draws[idx],
                       B, alpha, propagate_mean, level)
  structure(c(res, list(method = "CNBG-Bayes")), class = "cnbg_reserve")
}


#' Conditioning-respecting point-estimate bootstrap (WLS-cond)
#'
#' Point-estimate counterpart to [reserve_cnbg()]: identical predictive engine
#' (shared accident-year frailty, conjugate update, exact conditional
#' compound Poisson--Gamma cell draws) but with `(kappa, A)` fixed at the WLS
#' point estimate instead of integrated over the posterior. Because the two
#' share the engine, the gap between them isolates the value of posterior
#' integration over `(kappa, A)`.
#'
#' @param triangle Incremental paid-loss triangle.
#' @param B Number of bootstrap draws.
#' @param alpha Severity shape, held fixed.
#' @param propagate_mean Propagate Stage-1 mean-coefficient uncertainty
#'   (default `FALSE`, matching the point-estimate design).
#' @param level Predictive interval level.
#' @param kappa_cap WLS boundary cap.
#' @return An object of class `cnbg_reserve`.
#' @examples
#' set.seed(1)
#' tri <- simulate_cnbg(kappa = 5, alpha = 3)$triangle
#' bootstrap_cnbg_cond(tri)
#' @export
bootstrap_cnbg_cond <- function(triangle, B = 4000L, alpha = 3,
                                propagate_mean = FALSE, level = 0.95,
                                kappa_cap = .cnbg_kappa_cap) {
  fit <- fit_cnbg(triangle, method = "wls", kappa_cap = kappa_cap)
  if (!isTRUE(fit$ok)) stop("fit failed; cannot bootstrap")
  res <- .cnbg_predict(fit$ms, rep(fit$kappa, B), rep(fit$A, B),
                       B, alpha, propagate_mean, level)
  structure(c(res, list(method = "WLS-cond",
                        kappa = fit$kappa, A = fit$A)),
            class = "cnbg_reserve")
}

#' @export
print.cnbg_reserve <- function(x, ...) {
  cat(sprintf("%s reserve\n", x$method))
  cat(sprintf("  mean %.0f   %.0f%% PI [%.0f, %.0f]   width %.0f\n",
              x$mean, 100 * x$level, x$pi[1], x$pi[2], diff(x$pi)))
  invisible(x)
}


## ===========================================================================
## KAPPA POSTERIOR DIAGNOSTIC   (paper's regime + KL-tier classification)
## ===========================================================================

.kl_gamma <- function(a1, b1, a2, b2)
  (a1 - a2) * digamma(a1) - lgamma(a1) + lgamma(a2) +
    a2 * (log(b1) - log(b2)) + a1 * (b2 - b1) / b1

#' Kappa posterior diagnostic: regime and KL tier
#'
#' Compares the kappa posterior to its prior via three statistics --
#' credible-interval contraction `c`, standardised mean shift `s`, and
#' KL divergence -- and returns both classifications the paper uses:
#'
#' * `regime`: `"informative"` (`c < 1`; note this INCLUDES prior-driven
#'   portfolios, exactly as in the paper's binary split),
#'   `"incipient_flight"` (`c > 1`, `KL <= 5`), `"strong_flight"`
#'   (`c > 1`, `KL > 5`), or `"no_bayes"` (no usable posterior).
#' * `kl_tier`: `"prior-driven"` (`KL <= 1`), `"weak"` (`<= 5`),
#'   `"moderate"` (`<= 20`), `"strong"` (`> 20`) -- the four-level grading
#'   whose coverage gradient is the paper's validated instrument.
#'
#' KL is computed exactly on the posterior grid where the fit provides one
#' (`KL_grid`); for engines without a grid (the exact rstan engine) a
#' Gamma moment-match fallback is used and flagged in the output.
#'
#' Note the `"no_bayes"` regime of the paper is a sampler-failure
#' phenomenon; the deterministic grid engines cannot exhibit it, so on grid
#' fits it arises only from Stage-1 failure. Calling this on a WLS fit
#' returns `"not_applicable"`.
#'
#' @param fit A `cnbg_fit` from [fit_cnbg()] (`method = "bayes"`) or
#'   [fit_cnbg_exact()].
#' @return An object of class `cnbg_kappa_diag`.
#' @examples
#' set.seed(1)
#' tri <- simulate_cnbg(kappa = 5, alpha = 3)$triangle
#' diagnose_kappa(fit_cnbg(tri, "bayes"))
#' @export
diagnose_kappa <- function(fit) {
  if (!isTRUE(fit$ok))
    return(structure(list(regime = "no_bayes", kl_tier = NA_character_,
                          c = NA_real_, s = NA_real_, KL = NA_real_,
                          informative = FALSE, kl_exact = FALSE),
                     class = "cnbg_kappa_diag"))
  if (!fit$method %in% c("bayes", "exact"))
    return(structure(list(regime = "not_applicable", kl_tier = NA_character_,
                          c = NA_real_, s = NA_real_, KL = NA_real_,
                          informative = NA, kl_exact = FALSE),
                     class = "cnbg_kappa_diag"))

  pr <- fit$prior
  prior_mean <- pr[["shape"]] / pr[["rate"]]
  prior_sd   <- sqrt(pr[["shape"]]) / pr[["rate"]]
  prior_w    <- diff(qgamma(c(0.025, 0.975), pr[["shape"]], pr[["rate"]]))

  cc <- (fit$kappa_q975 - fit$kappa_q025) / prior_w
  ss <- (fit$kappa_mean - prior_mean) / prior_sd

  kl_exact <- !is.null(fit$KL_grid) && is.finite(fit$KL_grid)
  KL <- if (kl_exact) fit$KL_grid else {
    post_sd <- (fit$kappa_q975 - fit$kappa_q025) / 3.92
    if (post_sd > 0)
      .kl_gamma(fit$kappa_mean^2 / post_sd^2, fit$kappa_mean / post_sd^2,
                pr[["shape"]], pr[["rate"]])
    else NA_real_
  }

  kl_tier <- if (is.na(KL)) NA_character_
             else if (KL <= 1)  "prior-driven"
             else if (KL <= 5)  "weak"
             else if (KL <= 20) "moderate"
             else               "strong"

  regime <- if (is.na(cc)) "no_bayes"
            else if (cc < 1) "informative"
            else if (!is.na(KL) && KL > 5) "strong_flight"
            else "incipient_flight"

  structure(list(regime = regime, kl_tier = kl_tier,
                 c = cc, s = ss, KL = KL,
                 informative = isTRUE(cc < 1), kl_exact = kl_exact),
            class = "cnbg_kappa_diag")
}

#' @export
print.cnbg_kappa_diag <- function(x, ...) {
  cat(sprintf("kappa diagnostic: regime = %s", x$regime))
  if (!is.na(x$kl_tier)) cat(sprintf("   KL tier = %s", x$kl_tier))
  cat("\n")
  if (!is.na(x$c))
    cat(sprintf("  contraction c = %.3f   mean shift s = %.2f   KL = %.2f%s\n",
                x$c, x$s, x$KL, if (isTRUE(x$kl_exact)) "" else " (moment-match)"))
  invisible(x)
}


## ===========================================================================
## INFORMATION BOUND  (Fisher diagonal, Godambe, van Trees)
## ===========================================================================

.cnbg_m4 <- function(nu, alpha, kappa, Y_bar) {
  A  <- (alpha + 1) * Y_bar / alpha
  c1 <- (alpha + 1) * (alpha + 2) * (alpha + 3) / alpha^3 * Y_bar^3
  c2 <- 3 * A^2 * (1 + 1 / kappa) +
        4 * (alpha + 1) * (alpha + 2) * Y_bar^2 / (alpha^2 * kappa)
  c3 <- 6 * A * (1 + 2 / kappa) / kappa
  c4 <- 3 * (kappa + 2) / kappa^3
  c1 * nu + c2 * nu^2 + c3 * nu^3 + c4 * nu^4
}

.cnbg_cov_z <- function(n1, n2, A, kappa)
  n1^2 * n2^2 * (6 / kappa^3 + 2 / kappa^2) +
    A^2 * n1 * n2 / kappa +
    2 * A * n1 * n2 * (n1 + n2) / kappa^2

#' Information bound on kappa from the amount triangle
#'
#' Computes, at a given calibration, the working-independence (diagonal)
#' Fisher information for `kappa`, the within-row score-covariance
#' correction, the Godambe information `G`, the prior information `I_pi`
#' (van Trees), the prior fraction and the van Trees standard-deviation
#' floor on `kappa`, across a grid of true `kappa`.
#'
#' At its defaults (the paper's calibration) this function reproduces the
#' paper's van Trees table exactly; averaging `G` over
#' `alpha` in `c(1, 3, 10)` gives the published `G_avg` column. `I_pi`
#' depends only on the (truncated) prior and is independent of the data
#' calibration.
#'
#' @param kappa_grid True `kappa` values.
#' @param alpha Severity shape.
#' @param I,J,Y_bar,mu_base,mu_growth,w Calibration (see [simulate_cnbg()]).
#' @param prior Length-2 `c(shape, rate)` kappa prior; default reference.
#' @param kappa_trunc Lower truncation used by the prior (`>= 0.1`).
#' @return A data frame: `kappa`, `I_diag`, `Delta`, `G`, `I_pi`,
#'   `prior_frac`, `vantrees_sd`.
#' @examples
#' ## the published table (G averaged over alpha):
#' tabs <- lapply(c(1, 3, 10), function(a)
#'   cnbg_information_kappa(alpha = a)$G)
#' round(rowMeans(do.call(cbind, tabs)), 4)   # 0.1644 0.0610 0.0120 ...
#' @export
cnbg_information_kappa <- function(kappa_grid = c(3, 5, 10, 20, 50), alpha = 3,
                                   I = 10L, J = 10L, Y_bar = 5000,
                                   mu_base = 30, mu_growth = 0.05,
                                   w = .cnbg_paper_w,
                                   prior = .cnbg_prior_ref,
                                   kappa_trunc = 0.1) {
  stopifnot(length(w) == J)
  mu_i <- mu_base * exp(mu_growth * (seq_len(I) - 1L))
  w    <- w / sum(w)
  nu   <- outer(mu_i, w) * Y_bar
  obs  <- outer(seq_len(I), seq_len(J), function(i, j) (i + j - 1L) <= I)
  A    <- (alpha + 1) * Y_bar / alpha

  a <- prior[["shape"]]; b <- prior[["rate"]]
  integ <- function(k) ((a - 1) / k - b)^2 * dgamma(k, a, b)
  I_pi <- integrate(integ, kappa_trunc, Inf)$value /
          (1 - pgamma(kappa_trunc, a, b))

  rows <- lapply(kappa_grid, function(kappa) {
    mz <- A * nu + nu^2 / kappa
    cv <- .cnbg_m4(nu, alpha, kappa, Y_bar) / mz^2 - 1
    I_diag <- sum(((nu^2 / kappa^2)^2 / (cv * mz^2))[obs])
    Delta <- 0
    for (i in seq_len(I)) {
      js <- which(obs[i, ]); if (length(js) < 2L) next
      for (p in seq_along(js)[-length(js)]) for (q in js[(p + 1L):length(js)]) {
        j1 <- js[p]; n1 <- nu[i, j1]; n2 <- nu[i, q]
        Delta <- Delta + 2 * (n1^2 / kappa^2) * (n2^2 / kappa^2) *
                 .cnbg_cov_z(n1, n2, A, kappa) /
                 (cv[i, j1] * mz[i, j1]^2 * cv[i, q] * mz[i, q]^2)
      }
    }
    G <- I_diag^2 / (I_diag + Delta)
    data.frame(kappa = kappa, I_diag = I_diag, Delta = Delta, G = G,
               I_pi = I_pi, prior_frac = I_pi / (G + I_pi),
               vantrees_sd = 1 / sqrt(G + I_pi))
  })
  do.call(rbind, rows)
}

## --- regression test making the calibration error unrepeatable (K1) ---------
## Internal; call from testthat. Asserts the information bound at package
## defaults reproduces the PUBLISHED van Trees table (G averaged over
## alpha in {1, 3, 10}; prop33_constants.csv, 2026-07-23).
.cnbg_selftest_information <- function(tol_G = 0.01, tol_sd = 5e-3) {
  target_G  <- c(0.1644, 0.06096, 0.01199, 0.001653, 7.538e-5)
  target_sd <- c(1.4972, 1.7084, 1.8454, 1.8788, 1.8840)
  Gs <- sapply(c(1, 3, 10), function(a) cnbg_information_kappa(alpha = a)$G)
  G_avg <- rowMeans(Gs)
  I_pi  <- cnbg_information_kappa(alpha = 3)$I_pi[1]
  vt    <- 1 / sqrt(G_avg + I_pi)
  ok <- all(abs(G_avg / target_G - 1) < tol_G) &&
        all(abs(vt - target_sd) < tol_sd)
  if (!ok) stop("cnbg information bound does not reproduce the published table",
                "\n  G_avg = ", paste(signif(G_avg, 4), collapse = ", "))
  invisible(TRUE)
}


## ===========================================================================
## EXACT FRAILTY-AUGMENTED POSTERIOR ENGINE  (rstan; theoretical reference)
## ===========================================================================
## Returns the SAME list shape as the grid engines (kappa_draws, A_draws,
## kappa_median, CIs, ...), so reserve_cnbg() and diagnose_kappa() consume
## it unchanged. rstan is gated at call time via requireNamespace and lives
## in Suggests (K8).
##
## This is the exact posterior of the single-frailty (kappa_c -> infinity)
## CNBG model: one row frailty U_i per accident year, no cell shock,
## matching simulate_cnbg() and the paper's coverage experiment. It
## converges cleanly on individual 10x10 triangles (the paper's benchmark
## runs: all Rhat <= 1.001, zero divergences) at minutes per fit against
## seconds for the working likelihood -- priced out of population-scale
## deployment by the ~200-400x runtime, not broken (K7). The working-
## likelihood fit_cnbg() remains the recommended default.
##
## SINGLE-SOURCE NOTE (K7): inst/stan/cnbg_exact.stan is the canonical
## model; the embedded string below is the fallback for source()-level use
## only. A package test should assert the two are byte-identical so they
## cannot drift (the four-copy lesson).

.cnbg_exact_stan_code <- '
data {
  int<lower=1> N;
  vector<lower=0>[N] x;
  vector<lower=0>[N] nu;
  int<lower=1> I;
  int<lower=1, upper=I> row[N];
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
'
## NOTE: pre-2.26 array syntax `int row[N]`. On rstan >= 2.26 change to
## `array[N] int<lower=1,upper=I> row;` (here and in inst/stan/cnbg_exact.stan).

.cnbg_stan_env <- new.env(parent = emptyenv())

.cnbg_exact_model <- function() {
  if (!is.null(.cnbg_stan_env$mod)) return(.cnbg_stan_env$mod)
  if (!requireNamespace("rstan", quietly = TRUE))
    stop("fit_cnbg_exact() requires the 'rstan' package. install.packages('rstan').")

  f <- system.file("stan", "cnbg_exact.stan", package = "hgr")
  mod <- if (nzchar(f) && file.exists(f))
    rstan::stan_model(file = f, model_name = "cnbg_exact")
  else
    rstan::stan_model(model_code = .cnbg_exact_stan_code,
                      model_name = "cnbg_exact")

  .cnbg_stan_env$mod <- mod
  mod
}

.cnbg_exact <- function(ms, prior_shape, prior_rate, alpha = 3,
                        A_init = NULL, n_max = NULL,
                        chains = 4L, iter = 2000L, warmup = 1000L,
                        seed = 2026L, kappa_cap = .cnbg_kappa_cap,
                        adapt_delta = 0.95, max_treedepth = 12L,
                        verbose = FALSE, ...) {
  obs  <- ms$obs
  nu   <- as.numeric(ms$nu_mat[obs])
  x    <- as.numeric(ms$tri[obs])
  rows <- which(obs, arr.ind = TRUE)[, 1]

  if (is.null(A_init)) {
    z <- (x - nu)^2
    A_init <- max(stats::median(z / nu), .Machine$double.eps)
  }
  Ybar0 <- A_init * alpha / (alpha + 1)

  if (is.null(n_max)) {
    approx_count <- x / max(Ybar0, .Machine$double.eps)
    n_max <- as.integer(ceiling(3 * max(approx_count, 1) + 50))
  }

  standata <- list(
    N = length(x), x = x, nu = nu, I = ms$I,
    row = as.integer(rows), alpha = alpha, n_max = n_max,
    prior_shape = prior_shape, prior_rate = prior_rate,
    A_prior_shape = 2, A_prior_rate = 2 / A_init       # A ~ Gamma(2, 2/A_init)
  )

  init_fun <- function()
    list(kappa = prior_shape / prior_rate, A = A_init, U = rep(1, ms$I))

  fit <- rstan::sampling(
    .cnbg_exact_model(), data = standata,
    chains = chains, iter = iter, warmup = warmup, seed = seed,
    init = init_fun, refresh = if (verbose) max(iter %/% 10L, 1L) else 0L,
    control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
    ...)

  ## convergence guard: surface problems as findings, never hide them
  sm <- rstan::summary(fit, pars = c("kappa", "A"))$summary
  rhat_max <- suppressWarnings(max(sm[, "Rhat"], na.rm = TRUE))
  ess_min  <- suppressWarnings(min(sm[, "n_eff"], na.rm = TRUE))
  ndiv <- sum(vapply(rstan::get_sampler_params(fit, inc_warmup = FALSE),
                     function(p) sum(p[, "divergent__"]), numeric(1)))
  if (is.finite(rhat_max) && rhat_max > 1.01)
    warning(sprintf("exact fit: max Rhat = %.3f (> 1.01); treat draws with care.",
                    rhat_max))
  if (ndiv > 0)
    warning(sprintf("exact fit: %d divergent transition(s); consider adapt_delta > %.2f.",
                    ndiv, adapt_delta))

  post <- rstan::extract(fit, pars = c("kappa", "A"))
  kd <- pmin(post$kappa, kappa_cap)
  Ad <- post$A
  q  <- function(p) unname(stats::quantile(kd, p))

  list(
    kappa_draws  = kd,            A_draws    = Ad,
    kappa_mean   = mean(kd),      kappa_median = stats::median(kd),
    kappa_q025   = q(0.025),      kappa_q975 = q(0.975),
    kappa_sd     = stats::sd(kd),
    A_emp = A_init, A_hat = mean(Ad),
    adjusted = NA, godambe_r = NA, widen = NA,   # n/a to exact engine
    KL_grid = NULL,                              # diagnose falls back (doc'd)
    edge_mass = mean(kd >= .cnbg_edge_thresh),
    stanfit = fit, n_max = n_max, alpha = alpha,
    rhat_max = rhat_max, n_eff_min = ess_min, n_divergent = ndiv,
    engine = "exact-rstan"
  )
}

#' Fit the exact frailty-augmented CNBG posterior (rstan)
#'
#' Fits the exact, frailty-augmented Bayesian posterior by Hamiltonian Monte
#' Carlo. Conditional on the accident-year frailty, each cell is compound
#' Poisson-Gamma (Tweedie, \eqn{1<p<2}); the density is evaluated by a
#' truncated latent-count summation, numerically equivalent to the Dunn and
#' Smyth (2005) series. The mean structure is the same Stage-1 free-row GLM
#' used by [fit_cnbg()]; the severity shape `alpha` is held fixed (paid
#' triangles do not identify it). This is the single-frailty
#' (\eqn{\kappa_c\to\infty}) model: one row frailty per accident year, no
#' cell shock.
#'
#' The returned object is a `cnbg_fit` with `method = "exact"`, carrying the
#' same `(kappa, A)` posterior summaries and `kappa_draws`/`A_draws` as the
#' grid engines, so [reserve_cnbg()] and [diagnose_kappa()] apply unchanged.
#' The exact engine is the theoretical reference: it converges cleanly on
#' individual \eqn{10\times10} triangles (the paper's benchmark runs report
#' all \eqn{\hat R \le 1.001} with zero divergences) at minutes per fit
#' against seconds for the working likelihood, so it is priced out of
#' population-scale deployment by runtime, and [fit_cnbg()] remains the
#' recommended default.
#'
#' @param triangle Incremental paid-loss triangle (matrix); future cells `NA`.
#' @param prior `"ref"`, `"empb"`, or `c(shape, rate)`.
#' @param alpha Severity shape, held fixed.
#' @param kappa_cap Upper cap applied to `kappa` draws.
#' @param ... Passed to the sampler (`chains`, `iter`, `warmup`, `seed`,
#'   `adapt_delta`, ...).
#' @return A `cnbg_fit` (`method = "exact"`).
#' @examples
#' \dontrun{
#' tri <- simulate_cnbg(kappa = 5, alpha = 3)$triangle
#' fit <- fit_cnbg_exact(tri, prior = "ref", alpha = 3)
#' fit
#' reserve_cnbg(fit)
#' diagnose_kappa(fit)
#' }
#' @export
fit_cnbg_exact <- function(triangle, prior = "ref", alpha = 3,
                           kappa_cap = .cnbg_kappa_cap, ...) {
  pr <- if (is.numeric(prior)) prior
  else switch(prior, ref = .cnbg_prior_ref, empb = .cnbg_prior_empb,
              stop("prior must be 'ref', 'empb', or c(shape, rate)"))
  names(pr) <- c("shape", "rate")

  ms <- tryCatch(.cnbg_stage1(triangle), error = function(e) NULL)
  if (is.null(ms))
    return(structure(list(ok = FALSE, method = "exact",
                          message = "Stage-1 GLM failed"), class = "cnbg_fit"))

  g <- .cnbg_exact(ms, pr[["shape"]], pr[["rate"]], alpha = alpha,
                   kappa_cap = kappa_cap, ...)

  out <- c(list(ok = TRUE, method = "exact", ms = ms, prior = pr,
                prior_label = if (is.numeric(prior)) "custom" else prior,
                kappa_cap = kappa_cap), g)
  out$kappa <- g$kappa_median
  out$A     <- mean(g$A_draws)
  structure(out, class = "cnbg_fit")
}
