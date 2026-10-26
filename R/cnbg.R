#' Compound Negative Binomial-Gamma (CNBG) reserving
#'
#' Implementation of the amount-triangle methodology of "What Paid Triangles
#' Identify about Reserve Uncertainty: Information Boundaries in Loss
#' Reserving" (Van Oirbeek, 2026; the CNBG paper). The module is a port of
#' the library behind the paper's tables (`R/working_fit.R`, `R/predictive.R`,
#' `R/cell_shapes.R`, `R/curvature_adjustment.R`, `R/oracle.R` and
#' `R/legacy_methods.R` of the replication material), with the same
#' conventions, so that a fit obtained here is the fit the paper reports.
#'
#' The pieces, in the order the paper uses them:
#' \itemize{
#'   \item [simulate_cnbg()]: the data-generating process at the paper's
#'     calibration, with the row frailty `kappa_r`, the optional cell shock
#'     `kappa_c` (Section 7 of the paper, Phase 15) and the optional calendar
#'     drift (Phase 16).
#'   \item [cnbg_stage1()]: the Stage-1 mean structure, a quasi-Poisson GLM
#'     with free accident-year levels (Chain-Ladder means) or with the levels
#'     anchored to a known exposure.
#'   \item [fit_cnbg_wls()]: the two-stage second-moment estimator of
#'     `(A, B = 1/kappa_B)`, a non-negative least-squares fit of the squared
#'     residuals on `(nu, nu^2)` with weights `nu^-2` and the cap
#'     `kappa_B <= 1000`; [two_stage_target()] gives its population target
#'     under free levels (Proposition 4.2 and Remark 4.5 of the paper).
#'   \item [fit_cnbg_working()]: the working-likelihood posterior, with the
#'     Gamma shapes set cell by cell from the fourth moment, the squared
#'     residuals corrected for leverage, and the curvature adjustment of
#'     Ribatet, Cooley and Davison (2012). `mean = "absorbed"` is the fit the
#'     predictive uses under free levels; `mean = "full"` carries the
#'     `nu^2 / kappa_B` term and feeds [diagnose_kappa()].
#'   \item [fit_cnbg()] and [reserve_cnbg()]: the modular predictive CNBG-Mod
#'     (`A` from the absorbed fit, `kappa_B` from the prior, the row frailty
#'     restored by the conjugate update, moment-matched Gamma cells).
#'   \item [bootstrap_cnbg_cond()] (WLS-cond) and [bootstrap_cnbg_refit()]
#'     (Frequentist CNBG), the two point-estimate procedures of the paper.
#'   \item [diagnose_kappa()]: contraction, shift, Kullback--Leibler
#'     divergence and tier of the `kappa_B` posterior against its prior.
#'   \item [cnbg_information_kappa()]: the Fisher and Godambe information on
#'     `kappa_B` from squared residuals at a calibration (Section 6).
#'   \item [fit_cnbg_exact()]: the exact frailty-augmented posterior by
#'     Hamiltonian Monte Carlo, with anchored or free accident-year levels
#'     (requires `rstan`).
#' }
#'
#' The variance function is `Var(X) = A nu + nu^2 / kappa_B` throughout, with
#' `1 + 1/kappa_B = (1 + 1/kappa_r)(1 + 1/kappa_c)`. Under free accident-year
#' levels the row frailty is absorbed by the levels (Theorem 1 of the paper):
#' the amount triangle then supplies `A` and the predictive takes `kappa_B`
#' from its prior.
#'
#' @keywords internal
#' @name cnbg
#' @importFrom stats glm quasi glm.control predict coef vcov fitted lm qgamma
#'   dgamma rgamma rpois rnorm quantile integrate digamma lgamma sd median
#'   pgamma model.matrix approx optim influence runif
#' @importFrom stats AIC BIC aggregate dnorm logLik lowess model.frame pchisq
#'   poisson qchisq quasipoisson rbeta residuals rmultinom rnbinom setNames var
#' @importFrom graphics abline barplot boxplot legend lines par points polygon
#' @importFrom grDevices rgb
#' @importFrom MASS mvrnorm
NULL

`%||%` <- function(a, b) if (is.null(a)) b else a

## Default kappa priors (shape, rate) -----------------------------------------
.cnbg_prior_ref  <- c(shape = 2.00, rate = 0.40)   # reference prior: median 4.2, 95% mass on (0.6, 13.9)
.cnbg_prior_empb <- c(shape = 1.78, rate = 0.22)   # sensitivity prior: mean 8, sd about 6
.cnbg_kappa_cap  <- 1000                            # cap of the second-moment estimator

## Paper calibration (Section 5 of the paper) ---------------------------------
.cnbg_paper_w <- c(0.35, 0.25, 0.15, 0.10, 0.07, 0.04, 0.02, 0.01,
                   0.005, 0.005)

.cnbg_prior <- function(prior) {
  pr <- if (is.numeric(prior)) prior
        else switch(prior, ref = .cnbg_prior_ref, empb = .cnbg_prior_empb,
                    stop("prior must be 'ref', 'empb', or c(shape, rate)"))
  pr <- as.numeric(pr); names(pr) <- c("shape", "rate"); pr
}

.cnbg_obs_mask <- function(I, J) outer(seq_len(I), seq_len(J), "+") <= I + 1L

## ===========================================================================
## DATA-GENERATING PROCESS
## ===========================================================================

#' Simulate a CNBG paid-loss triangle
#'
#' Draws an incremental paid-loss triangle from the CNBG generative model:
#' accident-year frailty `U_i ~ Gamma(kappa_r, kappa_r)`, optional cell shock
#' `eps_ij ~ Gamma(kappa_c, kappa_c)`, cell counts
#' `N_ij ~ Poisson(mu_i U_i eps_ij w_j g_ij)`, Gamma severities with shape
#' `alpha` and mean `Y_bar`, and `X_ij` their sum. `g_ij = exp(-delta
#' max(0, t - onset))` with `t = i + j - 1` is the calendar drift of the
#' paper's Phase 16; `delta = 0` (the default) switches it off.
#'
#' Defaults reproduce the paper's simulation calibration: accident-year
#' expected claim counts `mu_i = mu_base exp(mu_growth (i - 1))` with
#' `mu_base = 30` (30 to 47 expected claims per row), the decreasing
#' development-weight vector `(.35, .25, .15, .10, .07, .04, .02, .01, .005,
#' .005)`, mean severity `Y_bar = 5000` and severity shape `alpha = 3`. The
#' random draws are made in the same order as the paper's scripts
#' (`oracle_dgp()`, `dgp_cellshock()`, `dgp_drift()`), so that a seed
#' reproduces the paper's triangles.
#'
#' @param I,J Number of accident and development years.
#' @param kappa_r Row-frailty dispersion (small = heterogeneous accident
#'   years). `kappa` is accepted as an alias.
#' @param kappa_c Cell-shock dispersion; `Inf` (the default) is the
#'   single-frailty model of the paper's main design.
#' @param alpha Gamma severity shape.
#' @param Y_bar Mean severity.
#' @param mu_base,mu_growth Accident-year expected claim count is
#'   `mu_base * exp(mu_growth (i - 1))`.
#' @param w Development weights, length `J` (normalised internally). For
#'   `J < 10` the first `J` paper weights are used.
#' @param delta,onset Calendar drift rate and onset (Phase 16 of the paper).
#' @param kappa Alias of `kappa_r`.
#' @return A list with `triangle` (incremental amounts, lower triangle `NA`),
#'   `full` (complete `I x J` matrix), `R_true` (realised outstanding
#'   reserve), `U` (row frailties), `eps` (cell shocks, all `1` when
#'   `kappa_c = Inf`), `N` (latent counts) and the calibration (`mu`, `w`,
#'   `nu`, `A`).
#' @examples
#' set.seed(1)
#' tri <- simulate_cnbg(kappa_r = 5, alpha = 3)
#' tri$R_true
#' @export
simulate_cnbg <- function(I = 10L, J = 10L, kappa_r = 5, kappa_c = Inf, alpha = 3,
                          Y_bar = 5000, mu_base = 30, mu_growth = 0.05,
                          w = NULL, delta = 0, onset = 7L, kappa = NULL) {
  if (!is.null(kappa)) kappa_r <- kappa
  if (is.null(w)) {
    stopifnot(J <= length(.cnbg_paper_w))
    w <- .cnbg_paper_w[seq_len(J)]
  }
  stopifnot(length(w) == J, all(w > 0), kappa_r > 0, kappa_c > 0, delta >= 0)
  mu_i <- mu_base * exp(mu_growth * (seq_len(I) - 1L))
  w    <- w / sum(w)
  sev_rate <- alpha / Y_bar

  lambdas <- rgamma(I, shape = kappa_r, rate = kappa_r / mu_i)   # lambda_i = mu_i U_i
  eps <- if (is.finite(kappa_c)) matrix(rgamma(I * J, shape = kappa_c, rate = kappa_c), I, J)
         else matrix(1, I, J)
  N <- matrix(0L, I, J); full <- matrix(0, I, J)
  for (i in seq_len(I)) for (j in seq_len(J)) {
    g <- if (delta > 0) exp(-delta * max(0, i + j - 1L - onset)) else 1
    n <- rpois(1L, lambdas[i] * eps[i, j] * w[j] * g)
    if (n > 0L) full[i, j] <- rgamma(1L, shape = n * alpha, rate = sev_rate)
    N[i, j] <- n
  }
  obs <- .cnbg_obs_mask(I, J)
  tri <- full; tri[!obs] <- NA_real_
  list(triangle = tri, full = full, R_true = sum(full[!obs]),
       U = lambdas / mu_i, eps = eps, N = N,
       mu = mu_i, w = w, nu = outer(mu_i, w) * Y_bar,
       A = (1 + 1 / alpha) * Y_bar, kappa_r = kappa_r, kappa_c = kappa_c,
       alpha = alpha, Y_bar = Y_bar)
}

#' Prepare a paid triangle for the CNBG fits
#'
#' Applies the data rules of the paper's Schedule P analysis: an optional
#' conversion from cumulative to incremental amounts, the lower triangle set to
#' `NA`, and negative increments either set to `NA` (the rule of the paper's
#' outcome-free panel, `negatives = "na"`) or floored at a small positive
#' value (`negatives = "floor"`, the earlier rule of the paper's Section
#' S8.3 sensitivity). Zero increments are kept.
#'
#' @param x Square matrix of amounts, accident years in rows and development
#'   years in columns.
#' @param cumulative `TRUE` if `x` holds cumulative paid amounts.
#' @param negatives `"na"` or `"floor"`.
#' @param floor_value Value used by `negatives = "floor"`.
#' @return A list with `triangle` (incremental, lower triangle `NA`),
#'   `n_negative`, `n_zero` and, when the lower triangle of `x` was
#'   observed, `R_actual` (the realised outstanding amount).
#' @examples
#' x <- simulate_cnbg()$full
#' prepare_triangle(x)$n_zero
#' @export
prepare_triangle <- function(x, cumulative = FALSE, negatives = c("na", "floor"),
                             floor_value = 1) {
  negatives <- match.arg(negatives)
  x <- as.matrix(x); I <- nrow(x); J <- ncol(x)
  incr <- if (cumulative) cbind(x[, 1, drop = FALSE], x[, -1, drop = FALSE] - x[, -J, drop = FALSE]) else x
  obs <- .cnbg_obs_mask(I, J)
  R_actual <- if (all(is.finite(incr[!obs]))) sum(incr[!obs]) else NA_real_
  tri <- incr; tri[!obs] <- NA_real_
  neg <- !is.na(tri) & tri < 0
  if (negatives == "na") tri[neg] <- NA_real_ else tri[neg] <- floor_value
  list(triangle = tri, n_negative = sum(neg), n_zero = sum(!is.na(tri) & tri == 0),
       R_actual = R_actual)
}

## ===========================================================================
## STAGE 1: MEAN STRUCTURE
## ===========================================================================

.cnbg_long <- function(tri, mu_anchor = NULL) {
  I <- nrow(tri); J <- ncol(tri)
  idx <- which(!is.na(tri), arr.ind = TRUE)
  d <- data.frame(ay = factor(idx[, 1], levels = seq_len(I)),
                  dy = factor(idx[, 2], levels = seq_len(J)),
                  x = tri[idx], i = idx[, 1], j = idx[, 2])
  if (!is.null(mu_anchor)) d$off <- log(mu_anchor[idx[, 1]])
  d
}

#' Stage 1: the mean structure of a paid triangle
#'
#' Fits the multiplicative mean `nu_ij = exp(a_i + b_j)` by a quasi-Poisson
#' GLM (`quasi(link = "log", variance = "mu")`, `maxit = 50`, `epsilon =
#' 1e-7`), the Stage-1 fit of the paper. `levels = "free"` estimates one
#' level per accident year (the Chain-Ladder means); `levels = "anchored"`
#' fixes the levels at a known exposure `mu_anchor` through an offset and
#' estimates one common scale, the anchored design of the paper's Section 4.
#'
#' @param triangle Incremental paid-loss triangle (matrix), future cells
#'   `NA`; negative increments must have been handled by
#'   [prepare_triangle()].
#' @param levels `"free"` or `"anchored"`.
#' @param mu_anchor Exposure (expected level) per accident year, used when
#'   `levels = "anchored"`.
#' @return A list of class `cnbg_stage1` with the GLM (`fit1`), the long data
#'   (`dat`), fitted means on the observed cells (`nu_hat`) and on the full
#'   square (`nu_mat`), leverages (`hat`), the prediction design
#'   (`X_pred`, `off_pred`), `q` (number of mean parameters), `n_obs`,
#'   `converged`, `I`, `J`, `tri`, `obs`.
#' @examples
#' set.seed(1)
#' ms <- cnbg_stage1(simulate_cnbg()$triangle)
#' ms$converged
#' @export
cnbg_stage1 <- function(triangle, levels = c("free", "anchored"), mu_anchor = NULL) {
  levels <- match.arg(levels)
  tri <- as.matrix(triangle); I <- nrow(tri); J <- ncol(tri)
  if (levels == "anchored") stopifnot(length(mu_anchor) == I, all(mu_anchor > 0))
  dat <- .cnbg_long(tri, if (levels == "anchored") mu_anchor else NULL)
  ctrl <- glm.control(maxit = 50, epsilon = 1e-7)
  fit1 <- if (levels == "free")
    glm(x ~ ay + dy, data = dat, family = quasi(link = "log", variance = "mu"), control = ctrl)
  else
    glm(x ~ dy + offset(off), data = dat, family = quasi(link = "log", variance = "mu"), control = ctrl)
  newdat <- expand.grid(ay = factor(seq_len(I), levels = seq_len(I)),
                        dy = factor(seq_len(J), levels = seq_len(J)))
  if (levels == "anchored") newdat$off <- log(mu_anchor[as.integer(newdat$ay)])
  nu_mat <- matrix(predict(fit1, newdata = newdat, type = "response"), I, J)
  hat <- tryCatch(influence(fit1)$hat, error = function(e) rep(0, nrow(dat)))
  X_pred <- if (levels == "free") model.matrix(~ ay + dy, data = newdat) else model.matrix(~ dy, data = newdat)
  off_pred <- if (levels == "free") rep(0, nrow(newdat)) else newdat$off
  structure(list(fit1 = fit1, dat = dat, nu_hat = fitted(fit1), hat = hat, nu_mat = nu_mat,
                 X_pred = X_pred, off_pred = off_pred, q = length(coef(fit1)), n_obs = nrow(dat),
                 levels = levels, converged = fit1$converged, I = I, J = J,
                 tri = tri, obs = !is.na(tri)),
            class = "cnbg_stage1")
}

#' @export
print.cnbg_stage1 <- function(x, ...) {
  fut <- !.cnbg_obs_mask(x$I, x$J)
  cat(sprintf("CNBG Stage 1 (%s levels): %d observed cells, %d mean parameters, %s\n",
              x$levels, x$n_obs, x$q, if (x$converged) "converged" else "NOT converged"))
  cat(sprintf("  Chain-Ladder-type reserve estimate: %.0f\n", sum(x$nu_mat[fut])))
  invisible(x)
}

## ===========================================================================
## STAGE 2: THE SECOND-MOMENT ESTIMATOR
## ===========================================================================

#' Two-stage second-moment estimator of (A, kappa_B)
#'
#' The second-stage estimator of the paper (eq. wls): non-negative least
#' squares of the squared Stage-1 residuals `z = (x - nu_hat)^2` on
#' `(nu_hat, nu_hat^2)` with weights `nu_hat^-2`, on all observed cells; the
#' quadratic coefficient is floored at `1 / kappa_cap` (the cap `kappa_B <=
#' 1000`); a non-positive linear coefficient falls back to the Pearson-type
#' estimate; the dispersion is scaled by `(n - q) / n`.
#'
#' Under free accident-year levels the population target of this estimator
#' is not the true `kappa_B`: the first-stage error of the row total gives
#' the quadratic coefficient the negative target `B*` of Proposition 4.2,
#' and with a finite cell shock the target is `B* + c_q / kappa_c` (Remark
#' 4.5). See [two_stage_target()]. The estimator therefore sits at the cap on
#' most correctly specified triangles, which is why the predictive takes
#' `kappa_B` from the prior. With the slope at the cap the intercept is the
#' mean of `z / nu_hat`, whose target is `A (1 - mean s_ij)`, about `0.82 A`
#' at the paper's design: the second-moment fit also understates `A`. The
#' working fit ([fit_cnbg_working()]) removes that deflation by the
#' leverage factor.
#'
#' @param ms A `cnbg_stage1` object.
#' @param kappa_cap Cap on `kappa_B`.
#' @return A list with `kappa` (df-adjusted), `A`, `B` (raw quadratic
#'   coefficient) and `at_cap`.
#' @examples
#' set.seed(1)
#' ms <- cnbg_stage1(simulate_cnbg()$triangle)
#' fit_cnbg_wls(ms)
#' @export
fit_cnbg_wls <- function(ms, kappa_cap = .cnbg_kappa_cap) {
  nu_hat <- as.numeric(ms$nu_hat); sq_res <- (ms$dat$x - nu_hat)^2
  wts <- 1 / pmax(nu_hat^2, 1e-8); sw <- sqrt(wts)
  nn <- tryCatch(nnls::nnls(cbind(sw * nu_hat, sw * nu_hat^2), sw * sq_res), error = function(e) NULL)
  if (!is.null(nn)) { A_hat <- nn$x[1]; B_hat <- nn$x[2] } else {
    fit2 <- lm(sq_res ~ 0 + nu_hat + I(nu_hat^2), weights = wts)
    A_hat <- max(coef(fit2)[1], 0); B_hat <- max(coef(fit2)[2], 0) }
  at_cap <- B_hat <= 1 / kappa_cap
  if (B_hat < 1 / kappa_cap) B_hat <- 1 / kappa_cap
  if (A_hat < 1e-8) A_hat <- max(sum(sq_res / pmax(nu_hat, 1e-8)) / (ms$n_obs - ms$q), 1e-6)
  list(kappa = unname((1 / B_hat) * (ms$n_obs - ms$q) / ms$n_obs), A = unname(A_hat),
       B = unname(B_hat), at_cap = unname(at_cap))
}

#' Population target of the free-level second-moment estimator
#'
#' With the pattern known and free accident-year levels, the fitted level of
#' a row is its observed total, so `nu_hat_ij = nu_ij U_i (1 + e_i)` and
#' `E[(X_ij - nu_hat_ij)^2 | U_i] = A nu_ij U_i (1 - s_ij) + (nu_ij U_i)^2
#' q_ij / kappa_c`, with `s_ij = nu_ij / nu_i+` the cell's share of its row
#' and `q_ij = 1 - 2 s_ij + sum_l s_il^2`. The population slope of the
#' second-stage regression is then `B_target = B* + c_q / kappa_c` with
#' `B* = -A Cov(s, nu U) / Var(nu U) < 0` on a decreasing pattern and
#' `c_q = Cov(nu U q, nu U) / Var(nu U)` (Proposition 4.2 and Remark 4.5 of
#' the paper; `make_table_two_stage_target.R` and
#' `make_table_two_stage_cellshock.R` of the replication material). Where
#' `B_target <= 0` the constrained estimator converges to the cap.
#'
#' @param kappa_r Row-frailty dispersion (vector).
#' @param kappa_c Cell-shock dispersion (vector; `Inf` for the single-frailty
#'   design).
#' @inheritParams simulate_cnbg
#' @return A data frame with one row per `(kappa_r, kappa_c)` pair:
#'   `B_star`, `c_q`, `B_target`, `kappa_target` (`Inf` where the target is
#'   the cap), `A_star_over_A` (the unconstrained intercept relative to the
#'   truth), `A_cap_over_A` (the target of the implemented constrained
#'   estimator of `A` when the slope sits at the cap: the mean of `1 - s_ij`
#'   over the cells, `0.82` at the paper's design, since the NNLS intercept
#'   is then the mean of `z / nu_hat`) and `kappa_c_threshold = c_q / |B*|`,
#'   the cell dispersion below which the second stage reads it.
#' @examples
#' two_stage_target(kappa_r = c(3, 5, 10, 20, 50), kappa_c = 5)
#' @export
two_stage_target <- function(kappa_r = c(3, 5, 10, 20, 50), kappa_c = Inf,
                             I = 10L, J = 10L, alpha = 3, Y_bar = 5000,
                             mu_base = 30, mu_growth = 0.05, w = NULL) {
  if (is.null(w)) w <- .cnbg_paper_w[seq_len(J)]
  w <- w / sum(w); mu <- mu_base * exp(mu_growth * (seq_len(I) - 1L))
  A <- (alpha + 1) * Y_bar / alpha; nu <- outer(mu, w) * Y_bar
  obs <- .cnbg_obs_mask(I, J)
  nup <- rowSums(nu * obs); s <- nu / nup
  q <- 1 - 2 * s + matrix(rowSums((s * obs)^2), I, J)
  x_cells <- nu[obs]; s_cells <- s[obs]; q_cells <- q[obs]
  Ex <- mean(x_cells); Es <- mean(s_cells); Exs <- mean(x_cells * s_cells)
  Exq <- mean(x_cells * q_cells); Ex2q <- mean(x_cells^2 * q_cells)
  grid <- expand.grid(kappa_r = kappa_r, kappa_c = kappa_c)
  out <- do.call(rbind, lapply(seq_len(nrow(grid)), function(r) {
    kr <- grid$kappa_r[r]; kc <- grid$kappa_c[r]
    Vx <- mean(x_cells^2) * (1 + 1 / kr) - Ex^2
    Bs <- -A * (Exs - Ex * Es) / Vx
    cq <- (Ex2q * (1 + 1 / kr) - Exq * Ex) / Vx
    Bt <- Bs + if (is.finite(kc)) cq / kc else 0
    Ast <- A * (1 - Es) - Bs * Ex
    data.frame(kappa_r = kr, kappa_c = kc, B_star = Bs, c_q = cq, B_target = Bt,
               kappa_target = if (Bt > 0) 1 / Bt else Inf,
               A_star_over_A = Ast / A, A_cap_over_A = 1 - Es, kappa_c_threshold = cq / abs(Bs))
  }))
  rownames(out) <- NULL
  out
}
