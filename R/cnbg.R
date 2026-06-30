#' Compound Negative Binomial-Gamma (CNBG) reserving
#'
#' Implementation of the amount-triangle methodology of "What Paid Triangles
#' Can and Cannot Identify: An Information Ceiling for Loss-Reserve
#' Uncertainty". Provides the CNBG data-generating process, two-stage
#' estimation of the structural parameters (kappa, A), the
#' conditioning-respecting posterior predictive reserve distribution, the
#' kappa posterior diagnostic with regime classification, and the
#' Fisher/Godambe/van Trees information bound used in the paper.
#'
#' The Bayesian working-likelihood posterior is computed by dense grid
#' quadrature over the two parameters (kappa, A) rather than MCMC: the target
#' is the Gamma(2, 2/mu_z) composite likelihood of the paper, the integration
#' is deterministic, and the package carries no Stan dependency. The variance
#' function is V(nu) = A nu + nu^2 / kappa throughout.
#'
#' @keywords internal
#' @name cnbg
#' @importFrom stats glm quasipoisson predict coef vcov fitted lm qgamma dgamma
#'   rgamma rpois quantile integrate digamma lgamma sd median pgamma model.matrix
#' @importFrom MASS mvrnorm
#' @importFrom rstan stan_model sampling extract summary get_sampler_params
NULL

## Default kappa priors (shape, rate) -----------------------------------------
.cnbg_prior_ref  <- c(shape = 2.00, rate = 0.40)   # weakly informative reference
.cnbg_prior_empb <- c(shape = 1.78, rate = 0.22)   # empirical Bayes (Sched. P)
.cnbg_kappa_cap  <- 1000                            # WLS boundary cap


## ===========================================================================
## DATA-GENERATING PROCESS
## ===========================================================================

#' Simulate a CNBG paid-loss triangle
#'
#' Draws an incremental paid-loss triangle from the CNBG generative model:
#' accident-year frailty `U_i ~ Gamma(kappa, kappa)`, cell counts
#' `N_ij ~ Poisson(mu_i U_i w_j)`, Gamma severities, and `X_ij` their sum.
#'
#' Defaults reproduce the paper's stated calibration: `mu_i = exp(mu_log_base +
#' mu_growth (i-1))`, geometric development weights, and mean severity `Y_bar`.
#' This places the cells in the heterogeneity-capable regime the information
#' bounds are derived for; lowering `mu_log_base` moves the triangle into the
#' sparse, severity-dominated regime.
#'
#' @param I,J Number of accident and development years.
#' @param kappa Frailty dispersion (small = heterogeneous accident years).
#' @param alpha Gamma severity shape.
#' @param Y_bar Mean severity.
#' @param mu_log_base,mu_growth Log accident-year intensity is
#'   `mu_log_base + mu_growth (i - 1)`.
#' @param w_ratio Geometric ratio for development weights (decreasing).
#' @return A list with `triangle` (incremental amounts, lower triangle `NA`),
#'   `full` (complete I x J matrix), `R_true` (realised outstanding reserve),
#'   and `U` (the drawn frailties).
#' @examples
#' set.seed(1)
#' tri <- simulate_cnbg(I = 10, J = 10, kappa = 5, alpha = 3)
#' tri$R_true
#' @export
simulate_cnbg <- function(I = 10L, J = 10L, kappa = 5, alpha = 3,
                          Y_bar = 1e4, mu_log_base = 7, mu_growth = 0.05,
                          w_ratio = 0.7) {
  mu_i  <- exp(mu_log_base + mu_growth * (seq_len(I) - 1L))
  w     <- w_ratio^(seq_len(J) - 1L); w <- w / sum(w)
  g_rate <- alpha / Y_bar

  U   <- rgamma(I, shape = kappa, rate = kappa)        # E[U]=1, Var=1/kappa
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
## STAGE 2a: WLS point estimate of (kappa, A)
## ===========================================================================

.cnbg_wls <- function(ms, kappa_cap = .cnbg_kappa_cap) {
  nu <- ms$nu_mat[ms$obs]; x <- ms$tri[ms$obs]
  z  <- (x - nu)^2 / nu^2          # z = B + A * (1/nu)
  u  <- 1 / nu
  cf <- coef(lm(z ~ u))
  B_hat <- max(cf[[1]], 0); A_hat <- max(cf[[2]], 0)
  if (B_hat < 1 / kappa_cap) B_hat <- 1 / kappa_cap
  if (A_hat <= 0) A_hat <- mean(z * nu)            # fallback: severity term
  kappa <- min(1 / B_hat, kappa_cap)
  kappa_adj <- kappa * (ms$n_obs - ms$q) / ms$n_obs
  list(kappa = kappa_adj, A = A_hat,
       at_cap = (1 / B_hat) >= kappa_cap)
}


## ===========================================================================
## STAGE 2b: Bayesian posterior of (kappa, A)
## ===========================================================================
## The CNBG cell likelihood is a working (composite) likelihood: cells of an
## accident-year row share a frailty, so they are correlated and the naive
## independent likelihood over-counts information. Amount triangles barely
## identify kappa at all, so the raw working likelihood is near-flat in kappa
## and drifts toward the boundary; a naive Bayesian fit is therefore either
## falsely sharp (if it locks onto a spurious mode) or boundary-seeking. We
## form an honest kappa posterior by (i) fixing A at its well-identified
## estimate, (ii) tempering the 1-D kappa working likelihood by the
## cluster-robust Godambe deflation r = I_naive / J_cluster -- computed from
## row-summed scores at the prior mean, where the likelihood is not flat -- and
## (iii) anchoring with the prior. With r < 1 the over-counted information is
## removed and the posterior collapses toward the prior, reflecting that the
## triangle cannot pin kappa. Set adjust = FALSE for the naive (overconfident)
## likelihood, for comparison only.

.cnbg_posterior <- function(ms, prior_shape, prior_rate,
                            kappa_cap = .cnbg_kappa_cap,
                            adjust = TRUE, n_grid = 600L, n_draw = 4000L) {
  nu <- ms$nu_mat[ms$obs]; x <- ms$tri[ms$obs]
  z  <- (x - nu)^2
  rows  <- which(ms$obs, arr.ind = TRUE)[, 1]        # accident-year cluster id
  A_emp <- max(mean(z) / mean(nu), .Machine$double.eps)
  A_hat <- max(median(z / nu), A_emp / 50)               # robust severity dispersion

  ## 1-D kappa working log-likelihood with A fixed
  llk <- function(k) {
    mu_z <- A_hat * nu + nu^2 / k
    sum(dgamma(z, shape = 2, rate = 2 / mu_z, log = TRUE))
  }

  ## cluster-robust Godambe deflation r, evaluated at the prior mean (non-flat)
  r <- 1
  if (adjust) {
    k0 <- prior_shape / prior_rate
    hh <- 1e-4 * k0
    sc <- (dgamma(z, 2, rate = 2 / (A_hat * nu + nu^2 / (k0 + hh)), log = TRUE) -
           dgamma(z, 2, rate = 2 / (A_hat * nu + nu^2 / (k0 - hh)), log = TRUE)) / (2 * hh)
    I_naive   <- sum(sc^2)                              # independent-cell info
    J_cluster <- sum(tapply(sc, rows, sum)^2)           # row-clustered info
    if (is.finite(I_naive) && is.finite(J_cluster) && J_cluster > 0)
      r <- min(max(I_naive / J_cluster, 0.02), 1)
  }

  kg <- exp(seq(log(0.1), log(kappa_cap), length.out = n_grid))
  ll <- vapply(kg, llk, numeric(1))
  ll[!is.finite(ll)] <- -Inf
  lp <- r * (ll - max(ll[is.finite(ll)])) +
        dgamma(kg, shape = prior_shape, rate = prior_rate, log = TRUE) +
        log(kg)                                         # log-grid Jacobian
  lp[!is.finite(lp)] <- -Inf
  mx <- suppressWarnings(max(lp))
  p  <- if (!is.finite(mx)) rep(1, n_grid) else exp(lp - mx)
  p[!is.finite(p)] <- 0; if (sum(p) <= 0) p <- rep(1, n_grid); p <- p / sum(p)

  cdf <- cumsum(p)
  qk  <- function(q) stats::approx(cdf, kg, xout = q, rule = 2, ties = min)$y
  kmean <- sum(p * kg)
  di <- sample(kg, n_draw, replace = TRUE, prob = p)

  ## A draws: well identified; mild log-normal jitter around A_hat
  A_sd <- 0.08
  A_draws <- A_hat * exp(rnorm(n_draw, -A_sd^2 / 2, A_sd))

  list(kappa_draws = di, A_draws = A_draws,
       kappa_mean = kmean, kappa_median = qk(0.5),
       kappa_q025 = qk(0.025), kappa_q975 = qk(0.975),
       kappa_sd = sqrt(sum(p * (kg - kmean)^2)),
       A_emp = A_emp, A_hat = A_hat, adjusted = adjust,
       godambe_r = r, widen = sqrt(1 / r),
       edge_mass = sum(p[kg >= 0.95 * kappa_cap]))
}


## ===========================================================================
## FIT
## ===========================================================================

#' Fit the CNBG model to a paid-loss triangle
#'
#' Estimates the mean structure by a free-row quasi-Poisson GLM (Chain-Ladder
#' means) and the structural pair `(kappa, A)` either by the Bayesian
#' working-likelihood posterior (`method = "bayes"`, the recommended default)
#' or by the two-stage weighted least-squares point estimate
#' (`method = "wls"`).
#'
#' @param triangle Incremental paid-loss triangle (matrix); observed cells
#'   non-`NA`, future cells `NA`.
#' @param method `"bayes"` (posterior over `(kappa, A)`) or `"wls"` (point
#'   estimate).
#' @param prior `"ref"` (Gamma(2, 0.4)) or `"empb"` (Gamma(1.78, 0.22)), or a
#'   length-2 numeric `c(shape, rate)`.
#' @param kappa_cap WLS boundary cap on `kappa`.
#' @param adjust For `method = "bayes"`, apply the cluster-robust Godambe
#'   curvature adjustment (default `TRUE`). `FALSE` returns the naive
#'   working-likelihood posterior, which is overconfident and should be used
#'   only for comparison.
#' @return An object of class `cnbg_fit`.
#'
#' @note The Bayesian engine is a deterministic grid quadrature of the
#'   working-likelihood posterior of `(kappa, A)` (target `z ~ Gamma(2, 2/mu_z)`,
#'   `mu_z = A nu + nu^2 / kappa`). That working likelihood treats the cells of
#'   an accident-year row as independent, but they share a row frailty, so the
#'   naive posterior manufactures information and its credible intervals are
#'   several-fold too narrow (they badly undercover the true `kappa`). With
#'   `adjust = TRUE` the curvature is rescaled to the cluster-robust Godambe
#'   information (rows as independent clusters) via the magnitude adjustment of
#'   Ribatet, Cooley and Davison (2012); the reported `widen` factor is how much
#'   the `kappa` posterior SD is inflated relative to the naive fit. Even
#'   adjusted, amount triangles identify `kappa` only weakly, so the posterior
#'   is largely prior-governed -- treat point summaries with caution and lean on
#'   [diagnose_kappa()].
#' @examples
#' set.seed(1)
#' tri <- simulate_cnbg(kappa = 5, alpha = 3)$triangle
#' fit <- fit_cnbg(tri, method = "bayes")
#' fit
#' @export
fit_cnbg <- function(triangle, method = c("bayes", "wls"),
                     prior = "ref", kappa_cap = .cnbg_kappa_cap,
                     adjust = TRUE) {
  method <- match.arg(method)
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
              prior = pr, prior_label = if (is.numeric(prior)) "custom" else prior,
              kappa_cap = kappa_cap)

  if (method == "wls") {
    w <- .cnbg_wls(ms, kappa_cap)
    out$kappa <- w$kappa; out$A <- w$A; out$at_cap <- w$at_cap
  } else {
    g <- .cnbg_posterior(ms, pr[["shape"]], pr[["rate"]], kappa_cap, adjust = adjust)
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
    else if (isTRUE(x$adjusted))
      sprintf("Godambe-adjusted (kappa SD x%.1f)", x$widen)
    else "NAIVE (overconfident)"
    cat(sprintf("  prior: %s   curvature: %s", x$prior_label, curv))
    cat(sprintf("   edge mass: %.1f%%\n", 100 * x$edge_mass))
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

.cnbg_predict <- function(ms, kappa_v, A_v, B, propagate_mean = FALSE,
                          level = 0.95) {
  I <- ms$I; J <- ms$J
  fut_mask <- outer(seq_len(I), seq_len(J), function(i, j) (i + j - 1L) > I)
  future <- which(fut_mask, arr.ind = TRUE)

  coef_hat <- coef(ms$fit1)
  V        <- vcov(ms$fit1)
  if (propagate_mean && min(eigen(V, symmetric = TRUE, only.values = TRUE)$values) < 1e-10)
    V <- V + diag(1e-8, nrow(V))

  R <- numeric(B)
  for (b in seq_len(B)) {
    kap <- max(kappa_v[b], 0.5); A <- max(A_v[b], 1e-8)

    if (propagate_mean) {
      th <- MASS::mvrnorm(1, coef_hat, V)
      nu_b <- matrix(exp(ms$X_design %*% th), I, J)
      nu_b <- pmin(nu_b, 10 * pmax(ms$nu_mat, 1))
      snu  <- vapply(seq_len(I), function(i) sum(nu_b[i, ms$obs[i, ]]), numeric(1))
    } else {
      nu_b <- ms$nu_mat; snu <- ms$sum_nu_by_ay
    }

    ## conjugate frailty posterior given the observed row totals
    U <- rgamma(I, shape = kap + ms$sum_X_by_ay / A,
                rate  = pmax(kap + snu / A, 1e-10))
    Rb <- 0
    for (r in seq_len(nrow(future))) {
      i <- future[r, 1]; j <- future[r, 2]; nu0 <- nu_b[i, j]
      if (nu0 > 0) Rb <- Rb + rgamma(1, shape = nu0 * U[i] / A, rate = 1 / A)
    }
    R[b] <- Rb
  }
  a <- (1 - level) / 2
  list(reserve = R, mean = mean(R), sd = sd(R),
       pi = quantile(R, c(a, 1 - a)), level = level)
}


#' Posterior predictive reserve distribution (CNBG-Bayes)
#'
#' Conditioning-respecting posterior predictive reserve: the observed triangle
#' is held fixed, accident-year frailties are drawn from their conjugate
#' posterior, future cells are drawn from the CNBG conditional moments, and
#' parameter uncertainty in `(kappa, A)` is integrated over the posterior.
#'
#' Setting `propagate_mean = TRUE` additionally propagates estimation
#' uncertainty in the Stage-1 mean coefficients by drawing them from their
#' asymptotic distribution; no data are resampled, so this remains distinct
#' from the refit-on-pseudo-triangle bootstrap. The default `FALSE` is the
#' strictly conditional predictive of the paper's Algorithm.
#'
#' @param fit A `cnbg_fit` from [fit_cnbg()] with `method = "bayes"`.
#' @param B Number of predictive draws.
#' @param propagate_mean Propagate Stage-1 mean-coefficient uncertainty.
#' @param level Predictive interval level.
#' @return An object of class `cnbg_reserve`.
#' @examples
#' set.seed(1)
#' tri <- simulate_cnbg(kappa = 5, alpha = 3)$triangle
#' reserve_cnbg(fit_cnbg(tri, "bayes"))
#' @export
reserve_cnbg <- function(fit, B = 4000L, propagate_mean = FALSE, level = 0.95) {
  if (!isTRUE(fit$ok)) stop("fit failed; cannot predict")
  if (!fit$method %in% c("bayes", "exact"))
    stop("reserve_cnbg requires a Bayesian or exact fit; use bootstrap_cnbg_cond for WLS")
  idx <- sample.int(length(fit$kappa_draws), B, replace = TRUE)
  res <- .cnbg_predict(fit$ms, fit$kappa_draws[idx], fit$A_draws[idx],
                       B, propagate_mean, level)
  structure(c(res, list(method = "CNBG-Bayes")), class = "cnbg_reserve")
}


#' Conditioning-respecting point-estimate bootstrap (WLS-cond)
#'
#' Point-estimate counterpart to [reserve_cnbg()]: identical predictive engine
#' (shared accident-year frailty, conjugate update, CNBG conditional cell
#' draws) but with `(kappa, A)` fixed at the WLS point estimate instead of
#' integrated over the posterior. Because the two share the engine, the gap
#' between them isolates the value of posterior integration over `(kappa, A)`.
#'
#' @param triangle Incremental paid-loss triangle.
#' @param B Number of bootstrap draws.
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
bootstrap_cnbg_cond <- function(triangle, B = 4000L, propagate_mean = FALSE,
                                level = 0.95, kappa_cap = .cnbg_kappa_cap) {
  fit <- fit_cnbg(triangle, method = "wls", kappa_cap = kappa_cap)
  if (!isTRUE(fit$ok)) stop("fit failed; cannot bootstrap")
  res <- .cnbg_predict(fit$ms, rep(fit$kappa, B), rep(fit$A, B),
                       B, propagate_mean, level)
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
## KAPPA POSTERIOR DIAGNOSTIC
## ===========================================================================

.kl_gamma <- function(a1, b1, a2, b2)
  (a1 - a2) * digamma(a1) - lgamma(a1) + lgamma(a2) +
    a2 * (log(b1) - log(b2)) + a1 * (b2 - b1) / b1

#' Kappa posterior diagnostic and regime classification
#'
#' Compares the kappa posterior to its prior via three statistics --
#' credible-interval contraction `c`, standardised mean shift `s`, and
#' KL divergence -- and assigns the portfolio to one of four regimes:
#' `"no_bayes"`, `"informative"`, `"prior_driven"`, or `"flight"`. The binary
#' informative/not-informative split is driven by `c < 1`.
#'
#' @param fit A `cnbg_fit` from [fit_cnbg()] with `method = "bayes"`.
#' @return An object of class `cnbg_kappa_diag`.
#' @examples
#' set.seed(1)
#' tri <- simulate_cnbg(kappa = 5, alpha = 3)$triangle
#' diagnose_kappa(fit_cnbg(tri, "bayes"))
#' @export
diagnose_kappa <- function(fit) {
  if (!isTRUE(fit$ok) || !fit$method %in% c("bayes", "exact"))
    return(structure(list(regime = "no_bayes", c = NA, s = NA, KL = NA),
                     class = "cnbg_kappa_diag"))

  pr <- fit$prior
  prior_mean <- pr[["shape"]] / pr[["rate"]]
  prior_sd   <- sqrt(pr[["shape"]]) / pr[["rate"]]
  prior_w    <- diff(qgamma(c(0.025, 0.975), pr[["shape"]], pr[["rate"]]))

  cc <- (fit$kappa_q975 - fit$kappa_q025) / prior_w
  ss <- (fit$kappa_mean - prior_mean) / prior_sd
  post_sd <- (fit$kappa_q975 - fit$kappa_q025) / 3.92
  KL <- if (post_sd > 0)
          .kl_gamma(fit$kappa_mean^2 / post_sd^2, fit$kappa_mean / post_sd^2,
                    pr[["shape"]], pr[["rate"]])
        else NA_real_

  regime <- if (is.na(cc)) "no_bayes"
            else if (cc > 1) "flight"
            else if (!is.na(KL) && KL > 1) "informative"
            else "prior_driven"

  structure(list(regime = regime, c = cc, s = ss, KL = KL,
                 informative = isTRUE(cc < 1)),
            class = "cnbg_kappa_diag")
}

#' @export
print.cnbg_kappa_diag <- function(x, ...) {
  cat(sprintf("kappa diagnostic: regime = %s\n", x$regime))
  if (!is.na(x$c))
    cat(sprintf("  contraction c = %.3f   mean shift s = %.2f   KL = %.2f\n",
                x$c, x$s, x$KL))
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
#' Computes, at a given calibration, the working-independence (diagonal) Fisher
#' information for `kappa`, the within-row score-covariance correction, the
#' Godambe information `G`, the prior information `I_pi` (van Trees), the prior
#' fraction and the van Trees standard-deviation floor on `kappa`, across a
#' grid of true `kappa`. Reproduces the paper's van Trees table.
#'
#' `I_pi` depends only on the (truncated) prior and is independent of the
#' data calibration.
#'
#' @param kappa_grid True `kappa` values.
#' @param alpha Severity shape.
#' @param I,J,Y_bar,mu_log_base,mu_growth,w_ratio Calibration (see
#'   [simulate_cnbg()]).
#' @param prior Length-2 `c(shape, rate)` kappa prior; default reference.
#' @param kappa_trunc Lower truncation used by the prior (Stan-style `>= 0.1`).
#' @return A data frame: `kappa`, `I_diag`, `Delta`, `G`, `I_pi`,
#'   `prior_frac`, `vantrees_sd`.
#' @examples
#' cnbg_information_kappa(kappa_grid = c(3, 10, 50), alpha = 3)
#' @export
cnbg_information_kappa <- function(kappa_grid = c(3, 5, 10, 20, 50), alpha = 3,
                                   I = 10L, J = 10L, Y_bar = 1e4,
                                   mu_log_base = 7, mu_growth = 0.05,
                                   w_ratio = 0.7,
                                   prior = .cnbg_prior_ref,
                                   kappa_trunc = 0.1) {
  mu_i <- exp(mu_log_base + mu_growth * (seq_len(I) - 1L))
  w    <- w_ratio^(seq_len(J) - 1L); w <- w / sum(w)
  nu   <- outer(mu_i, w) * Y_bar
  obs  <- outer(seq_len(I), seq_len(J), function(i, j) (i + j - 1L) <= I)
  A    <- (alpha + 1) * Y_bar / alpha

  ## prior (van Trees) information, truncated and renormalised
  a <- prior[["shape"]]; b <- prior[["rate"]]
  integ <- function(k) ((a - 1) / k - b)^2 * dgamma(k, a, b)
  I_pi <- integrate(integ, kappa_trunc, Inf)$value / (1 - pgamma(kappa_trunc, a, b))

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

## ===========================================================================
## EXACT FRAILTY-AUGMENTED POSTERIOR ENGINE  (rstan; Remark 5.12 / Table 1)
## ===========================================================================
## Returns the SAME list shape as .cnbg_posterior() (kappa_draws, A_draws,
## kappa_median, CIs, ...), so reserve_cnbg() and diagnose_kappa() consume it
## unchanged once their method guards admit "exact" (see part B / the guards
## above). rstan is gated at call time via requireNamespace, so it can sit in
## Suggests rather than Imports.
##
## This is the exact posterior of the single-frailty (kappa_c -> infinity)
## CNBG model: one row frailty U_i per accident year, no cell shock, matching
## simulate_cnbg() and the coverage experiment of the paper. The working-
## likelihood engine (.cnbg_posterior) is the recommended default; this engine
## is the theoretical reference, fragile on 10x10 triangles (kappa-lambda
## funnel) and not for population-scale deployment.
## ===========================================================================

## --- embedded Stan source (kept in sync with inst/stan/cnbg_exact.stan) ------
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

## --- locate / compile the Stan model (cached, lazy) -------------------------
.cnbg_stan_env <- new.env(parent = emptyenv())

.cnbg_exact_model <- function() {
  if (!is.null(.cnbg_stan_env$mod)) return(.cnbg_stan_env$mod)
  if (!requireNamespace("rstan", quietly = TRUE))
    stop("fit_cnbg_exact() requires the 'rstan' package. install.packages('rstan').")
  
  ## prefer the packaged model file if the package is built; else compile the
  ## embedded string. Both are read lazily, never at source time.
  f <- system.file("stan", "cnbg_exact.stan", package = "hgr")
  mod <- if (nzchar(f) && file.exists(f))
    rstan::stan_model(file = f, model_name = "cnbg_exact")
  else
    rstan::stan_model(model_code = .cnbg_exact_stan_code, model_name = "cnbg_exact")
  
  .cnbg_stan_env$mod <- mod
  mod
}

## --- core fit ---------------------------------------------------------------
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
    warning(sprintf("exact fit: max Rhat = %.3f (> 1.01); treat draws with care.", rhat_max))
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
    adjusted = NA, godambe_r = NA, widen = NA,   # not applicable to exact engine
    edge_mass = mean(kd >= 0.95 * kappa_cap),
    stanfit = fit, n_max = n_max, alpha = alpha,
    rhat_max = rhat_max, n_eff_min = ess_min, n_divergent = ndiv,
    engine = "exact-rstan"
  )
}

## --- user-facing fit --------------------------------------------------------
#' Fit the exact frailty-augmented CNBG posterior (rstan)
#'
#' Fits the exact, frailty-augmented Bayesian posterior of Remark 5.12 by
#' Hamiltonian Monte Carlo. Conditional on the accident-year frailty, each cell
#' is compound Poisson-Gamma (Tweedie, \eqn{1<p<2}); the density is evaluated by
#' a truncated latent-count summation, numerically equivalent to the Dunn and
#' Smyth (2005) series. The mean structure is the same Stage-1 free-row GLM used
#' by [fit_cnbg()]; the severity shape `alpha` is held fixed (paid triangles do
#' not identify it). This is the single-frailty (\eqn{\kappa_c\to\infty}) model:
#' one row frailty \eqn{\kappa_r} per accident year, no cell shock.
#'
#' The returned object is a `cnbg_fit` with `method = "exact"`, carrying the
#' same `(kappa, A)` posterior summaries and `kappa_draws`/`A_draws` as the
#' working-likelihood engine, so [reserve_cnbg()] and [diagnose_kappa()] apply
#' unchanged. The exact engine is a theoretical reference: on \eqn{10\times10}
#' triangles the \eqn{\kappa}-\eqn{\lambda} funnel makes HMC fragile, and the
#' working-likelihood [fit_cnbg()] is the recommended default for deployment.
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

