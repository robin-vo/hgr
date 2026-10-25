## ===========================================================================
## PREDICTIVE PROCEDURES  (port of R/predictive.R and R/legacy_methods.R of the
## paper's replication material)
## ===========================================================================

.cnbg_stage1_design <- function(ms) {
  I <- ms$I %||% max(ms$dat$i); J <- ms$J %||% max(ms$dat$j)
  newdat <- expand.grid(ay = factor(seq_len(I), levels = seq_len(I)), dy = factor(seq_len(J), levels = seq_len(J)))
  list(X = model.matrix(~ ay + dy, data = newdat), I = I, J = J)
}

#' The modular predictive CNBG-Mod
#'
#' The conditioning-respecting predictive of the paper's Section 5
#' ("Draws: the modular benchmark"). Each of the `B` draws proceeds as
#' follows: `A` is drawn from the posterior draws supplied (the absorbed
#' working fit under free levels) and `kappa_B` from the prior, floored at
#' `kappa_floor`; the mean parameters are drawn from `N(theta_hat, V)` with
#' the retransformation correction `-v/2` on the log scale (`theta_mode =
#' "biascorr"`), the drawn cell means capped at `nu_cap_mult` times the fitted
#' ones; the row frailty is restored by the conjugate update `U_i ~
#' Gamma(kappa + X_i+/A, kappa + nu_i+/A)`; the future cells are drawn from
#' the moment-matched `Gamma(nu U / A, 1/A)` (`cells = "moment"`) or from
#' the compound Poisson--Gamma with its zero atom at the severity convention
#' `alpha` (`cells = "compound"`, the law of the paper's oracle and exact
#' arms). The reserve is the sum over the future cells.
#'
#' `kappa_draws` replaces the prior by posterior draws of `kappa_B` paired
#' with `A_draws` by index (the full-mean variant of the paper's Section
#' S8); `frailty = "fixed"` holds `U_i = 1` (the step-3 ablation).
#'
#' @param ms A `cnbg_stage1` object (free or anchored levels).
#' @param A_draws Posterior draws of `A`.
#' @param prior Prior on `kappa_B`: `"ref"`, `"empb"` or `c(shape, rate)`.
#' @param B Number of predictive draws.
#' @param kappa_floor Lower bound on the drawn `kappa_B`.
#' @param nu_cap_mult Cap on the drawn cell means relative to the fitted ones.
#' @param theta_mode `"biascorr"` or `"fixed"`.
#' @param kappa_draws Optional posterior draws of `kappa_B`.
#' @param frailty `"conjugate"` or `"fixed"`.
#' @param cells `"moment"` or `"compound"`.
#' @param alpha Severity shape used by `cells = "compound"`.
#' @param level Interval level.
#' @return An object of class `cnbg_reserve` with the predictive draws
#'   (`reserve`), `mean`, `median`, `sd`, the central interval `pi`,
#'   `width`, the Chain-Ladder-type point estimate `cl_est` and the
#'   posterior standard deviation of the row frailties (`U_sd`).
#' @examples
#' set.seed(1)
#' ms <- cnbg_stage1(simulate_cnbg()$triangle)
#' f  <- fit_cnbg_working(ms, mean = "absorbed")
#' predict_cnbg(ms, f$post$A, B = 500)
#' @export
predict_cnbg <- function(ms, A_draws, prior = "ref", B = 4000L, kappa_floor = 0.5,
                         nu_cap_mult = 10, theta_mode = c("biascorr", "fixed"), kappa_draws = NULL,
                         frailty = c("conjugate", "fixed"), cells = c("moment", "compound"),
                         alpha = 3, level = 0.95) {
  theta_mode <- match.arg(theta_mode); frailty <- match.arg(frailty); cells <- match.arg(cells)
  pr <- .cnbg_prior(prior)
  I <- nrow(ms$nu_mat); J <- ncol(ms$nu_mat)
  X_design <- if (!is.null(ms$X_pred)) ms$X_pred else .cnbg_stage1_design(ms)$X
  off <- if (!is.null(ms$off_pred)) ms$off_pred else rep(0, nrow(X_design))
  fut <- outer(seq_len(I), seq_len(J), "+") > I + 1L
  coef_hat <- coef(ms$fit1); V <- vcov(ms$fit1)
  if (min(eigen(V, symmetric = TRUE, only.values = TRUE)$values) < 1e-10) V <- V + diag(1e-8, nrow(V))
  v_eta <- rowSums((X_design %*% V) * X_design)
  sum_X_by_ay <- tapply(ms$dat$x, factor(ms$dat$i, levels = seq_len(I)), sum); sum_X_by_ay[is.na(sum_X_by_ay)] <- 0
  obs <- outer(seq_len(I), seq_len(J), "+") <= I + 1L
  fut_idx <- which(fut); fut_i <- ((fut_idx - 1L) %% I) + 1L
  idx <- sample.int(length(A_draws), B, replace = TRUE); A_v <- pmax(A_draws[idx], 1e-6)
  kappa_v <- if (is.null(kappa_draws)) pmax(rgamma(B, shape = pr[["shape"]], rate = pr[["rate"]]), kappa_floor)
             else pmax(kappa_draws[idx], kappa_floor)
  reserves <- numeric(B); U_sd <- matrix(NA_real_, B, I)
  for (b in seq_len(B)) {
    theta_b <- if (theta_mode == "fixed") coef_hat else MASS::mvrnorm(1, mu = coef_hat, Sigma = V)
    eta_b <- as.numeric(X_design %*% theta_b) + off
    if (theta_mode == "biascorr") eta_b <- eta_b - 0.5 * v_eta
    nu_b <- matrix(exp(eta_b), I, J)
    if (is.finite(nu_cap_mult)) nu_b <- pmin(nu_b, nu_cap_mult * pmax(ms$nu_mat, 1))
    sum_nu_b <- rowSums(nu_b * obs)
    U <- if (frailty == "fixed") rep(1, I)
         else rgamma(I, shape = kappa_v[b] + sum_X_by_ay / A_v[b], rate = pmax(kappa_v[b] + sum_nu_b / A_v[b], 1e-10))
    U_sd[b, ] <- U
    m <- nu_b[fut_idx] * U[fut_i]
    if (cells == "moment") {
      sh <- m / A_v[b]; ok <- sh > 0
      reserves[b] <- sum(rgamma(sum(ok), shape = sh[ok], rate = 1 / A_v[b]))
    } else {
      Yb <- A_v[b] * alpha / (alpha + 1)
      N <- rpois(length(m), m / Yb); pos <- N > 0L
      reserves[b] <- if (any(pos)) sum(rgamma(sum(pos), shape = N[pos] * alpha, rate = alpha / Yb)) else 0
    }
  }
  a <- (1 - level) / 2
  q <- unname(quantile(reserves, c(a, 1 - a)))
  structure(list(reserve = reserves, mean = mean(reserves), median = median(reserves), sd = sd(reserves),
                 pi = q, width = q[2] - q[1], level = level, cl_est = sum(ms$nu_mat[fut]),
                 U_sd = apply(U_sd, 2, sd), method = "CNBG-Mod", prior = pr, cells = cells),
            class = "cnbg_reserve")
}

#' @export
print.cnbg_reserve <- function(x, ...) {
  cat(sprintf("%s reserve\n", x$method))
  cat(sprintf("  mean %.0f   median %.0f   %.0f%% interval [%.0f, %.0f]   width %.0f\n",
              x$mean, x$median, 100 * x$level, x$pi[1], x$pi[2], x$width))
  if (!is.null(x$cl_est)) cat(sprintf("  Chain-Ladder-type point estimate %.0f\n", x$cl_est))
  if (!is.null(x$kappa)) cat(sprintf("  at kappa_B = %.2f, A = %.4g\n", x$kappa, x$A))
  invisible(x)
}

#' Fit the CNBG model to a paid triangle (CNBG-Mod)
#'
#' Convenience wrapper that runs the pipeline of the paper's Schedule P and
#' simulation studies: [cnbg_stage1()], the second-moment estimator
#' [fit_cnbg_wls()], the absorbed working fit [fit_cnbg_working()] that
#' supplies the posterior of `A`, and, when `diagnose = TRUE`, the full-mean
#' working fit that feeds [diagnose_kappa()]. The predictive is obtained by
#' [reserve_cnbg()].
#'
#' @param triangle Incremental paid-loss triangle (matrix), future cells
#'   `NA`; see [prepare_triangle()] for the data rules.
#' @param prior Prior on `kappa_B`: `"ref"` (Gamma(2, 0.4), the reference
#'   prior), `"empb"` (Gamma(1.78, 0.22), the sensitivity prior) or
#'   `c(shape, rate)`.
#' @param alpha Severity-shape convention.
#' @param engine `"grid"` or `"nuts"` (see [fit_cnbg_working()]).
#' @param levels,mu_anchor Passed to [cnbg_stage1()].
#' @param diagnose Also fit the full-mean working posterior.
#' @param seed Seed passed to the working fits.
#' @param ... Further arguments to [fit_cnbg_working()].
#' @return An object of class `cnbg_fit` with the Stage-1 object (`ms`),
#'   the second-moment estimate (`wls`), the absorbed fit (`abs`), the
#'   full-mean fit (`full`, or `NULL`), the posterior draws of `A`
#'   (`A_draws`), `prior`, `alpha`, `engine` and `ok`.
#' @examples
#' set.seed(1)
#' tri <- simulate_cnbg(kappa_r = 5, alpha = 3)$triangle
#' fit <- fit_cnbg(tri)
#' fit
#' @export
fit_cnbg <- function(triangle, prior = "ref", alpha = 3, engine = c("grid", "nuts"),
                     levels = c("free", "anchored"), mu_anchor = NULL, diagnose = TRUE,
                     seed = 1L, ...) {
  engine <- match.arg(engine); levels <- match.arg(levels)
  pr <- .cnbg_prior(prior)
  ms <- tryCatch(cnbg_stage1(triangle, levels, mu_anchor), error = function(e) NULL)
  if (is.null(ms) || !ms$converged)
    return(structure(list(ok = FALSE, message = "Stage-1 GLM failed or did not converge", ms = ms),
                     class = "cnbg_fit"))
  wls <- tryCatch(fit_cnbg_wls(ms), error = function(e) NULL)
  f_abs <- tryCatch(fit_cnbg_working(ms, alpha, "absorbed", pr, engine, seed = seed + 2L, ...), error = function(e) NULL)
  f_full <- if (diagnose) tryCatch(fit_cnbg_working(ms, alpha, "full", pr, engine, seed = seed, ...), error = function(e) NULL) else NULL
  ok <- .cnbg_working_usable(f_abs)
  structure(list(ok = ok, message = if (ok) NULL else "absorbed working fit not usable",
                 ms = ms, wls = wls, abs = f_abs, full = f_full,
                 A_draws = if (ok) f_abs$post$A else NULL,
                 prior = pr, prior_label = if (is.numeric(prior)) "custom" else prior,
                 alpha = alpha, engine = engine, levels = levels),
            class = "cnbg_fit")
}

#' @export
print.cnbg_fit <- function(x, ...) {
  if (is.null(x$ms)) { cat("CNBG fit FAILED:", x$message, "\n"); return(invisible(x)) }
  if (identical(x$method, "exact")) {
    cat(sprintf("CNBG exact frailty-augmented posterior (%s levels, prior %s, alpha = %g)\n",
                x$levels, x$prior_label, x$alpha))
    cat(sprintf("  kappa_B: median %.2f   95%% [%.2f, %.2f]   A: median %.4g\n",
                x$kappa_median, x$kappa_q025, x$kappa_q975, median(x$A_draws)))
    cat(sprintf("  sampler: Rhat %.3f, min n_eff %.0f, divergences %d\n", x$rhat_max, x$n_eff_min, x$n_divergent))
    return(invisible(x))
  }
  cat(sprintf("CNBG-Mod fit (%s levels, %s engine, prior %s = Gamma(%g, %g), alpha = %g)\n",
              x$levels, x$engine, x$prior_label, x$prior[["shape"]], x$prior[["rate"]], x$alpha))
  if (!isTRUE(x$ok)) cat("  NOT usable:", x$message, "\n")
  fut <- !.cnbg_obs_mask(x$ms$I, x$ms$J)
  cat(sprintf("  Stage 1: %d cells, Chain-Ladder-type reserve %.0f\n", x$ms$n_obs, sum(x$ms$nu_mat[fut])))
  if (!is.null(x$wls))
    cat(sprintf("  second-moment estimate: kappa_B %.1f%s, A %.4g\n", x$wls$kappa, if (x$wls$at_cap) " (at the cap)" else "", x$wls$A))
  if (!is.null(x$abs) && is.null(x$abs$failed))
    cat(sprintf("  absorbed fit: A median %.4g, 95%% [%.4g, %.4g] (%d cells)\n",
                x$abs$adj$A[2], x$abs$adj$A[1], x$abs$adj$A[3], x$abs$n_cells))
  if (!is.null(x$full) && is.null(x$full$failed))
    cat(sprintf("  full-mean fit: kappa_B median %.2f, 95%% [%.2f, %.2f]; see diagnose_kappa()\n",
                x$full$adj$kappa[2], x$full$adj$kappa[1], x$full$adj$kappa[3]))
  cat("  predictive: reserve_cnbg() draws A from the absorbed fit and kappa_B from the prior\n")
  invisible(x)
}

#' Predictive reserve distribution of a CNBG fit
#'
#' For a [fit_cnbg()] object, the CNBG-Mod predictive of [predict_cnbg()]
#' with `A` from the absorbed working fit and `kappa_B` from the prior (the
#' fit's prior unless `prior` is given, so that the sensitivity prior can be
#' applied to the same fit, as in the paper's CNBG-Mod-EmpB). For a
#' [fit_cnbg_exact()] object, the exact predictive: at each posterior draw of
#' `(kappa_B, A, U, mu)` the future cells are drawn from the compound
#' Poisson--Gamma law.
#'
#' @param fit A `cnbg_fit`.
#' @param B Number of predictive draws.
#' @param prior Prior on `kappa_B` for the CNBG-Mod predictive (default: the
#'   fit's prior).
#' @param ... Further arguments to [predict_cnbg()] (`kappa_floor`,
#'   `cells`, `theta_mode`, `level`, ...).
#' @return A `cnbg_reserve` object.
#' @examples
#' set.seed(1)
#' tri <- simulate_cnbg(kappa_r = 5, alpha = 3)$triangle
#' reserve_cnbg(fit_cnbg(tri, diagnose = FALSE), B = 1000)
#' @export
reserve_cnbg <- function(fit, B = 4000L, prior = NULL, ...) {
  if (identical(fit$method, "exact")) return(.cnbg_exact_predict(fit, B, ...))
  if (!isTRUE(fit$ok)) stop("fit not usable; cannot predict (", fit$message, ")")
  pr <- if (is.null(prior)) fit$prior else .cnbg_prior(prior)
  res <- predict_cnbg(fit$ms, fit$A_draws, pr, B = B, alpha = fit$alpha, ...)
  res$method <- if (identical(fit$prior_label, "empb") || (is.character(prior) && prior == "empb")) "CNBG-Mod-EmpB"
                else if (is.null(prior) && identical(fit$prior_label, "ref")) "CNBG-Mod-Ref" else "CNBG-Mod"
  res
}

#' The conditional WLS bootstrap (WLS-cond)
#'
#' The point-estimate procedure of the paper's Section 5.3 and Algorithm 1:
#' the mean structure held at the Stage-1 fit, `(kappa_B, A)` at the
#' second-moment estimate, one frailty `U_i ~ Gamma(kappa_hat, kappa_hat)`
#' drawn per accident year and per replication (not conditioned on the
#' observed row), and the future cells drawn from `Gamma(nu / A, 1 / (A
#' U_i))`. It carries neither parameter nor mean-structure uncertainty and
#' is the lower edge of the paper's design space.
#'
#' @param x A triangle (matrix) or a `cnbg_stage1` object.
#' @param B Number of replications.
#' @param level Interval level.
#' @param kappa_cap Cap of the second-moment estimator.
#' @return A `cnbg_reserve` object carrying `kappa` and `A`.
#' @examples
#' set.seed(1)
#' bootstrap_cnbg_cond(simulate_cnbg()$triangle, B = 1000)
#' @export
bootstrap_cnbg_cond <- function(x, B = 4000L, level = 0.95, kappa_cap = .cnbg_kappa_cap) {
  ms <- if (inherits(x, "cnbg_stage1")) x else cnbg_stage1(x, "free")
  wls <- fit_cnbg_wls(ms, kappa_cap)
  I <- nrow(ms$nu_mat); J <- ncol(ms$nu_mat)
  fut <- outer(seq_len(I), seq_len(J), "+") > I + 1L
  fut_idx <- which(fut); fut_i <- ((fut_idx - 1L) %% I) + 1L
  nu_f <- pmax(ms$nu_mat[fut], 1e-10); sh <- nu_f / wls$A; ok <- sh > 1e-10
  k <- pmax(wls$kappa, 1e-6)
  reserves <- vapply(seq_len(B), function(b) {
    U <- rgamma(I, shape = k, rate = k)
    sum(rgamma(sum(ok), shape = sh[ok], rate = 1 / (wls$A * pmax(U[fut_i][ok], 1e-10)))) }, numeric(1))
  a <- (1 - level) / 2; q <- unname(quantile(reserves, c(a, 1 - a)))
  structure(list(reserve = reserves, mean = mean(reserves), median = median(reserves), sd = sd(reserves),
                 pi = q, width = q[2] - q[1], level = level, cl_est = sum(ms$nu_mat[fut]),
                 kappa = wls$kappa, A = wls$A, at_cap = wls$at_cap, method = "WLS-cond"),
            class = "cnbg_reserve")
}

#' The refit bootstrap (Frequentist CNBG)
#'
#' The paper's Frequentist CNBG (F-CNBG): pseudo-triangles are drawn cell by
#' cell from a Tweedie law matched to the CNBG mean and variance at the
#' second-moment estimate, the Stage-1 GLM and the second-moment estimator
#' are refitted on each, and the future cells are drawn from the refitted
#' law. Under free accident-year levels this reproduces, to first order, the
#' conditional predictive variance under a diffuse prior on the levels and
#' the pattern (Section 5 of the paper). Requires the `tweedie` package.
#'
#' @param triangle Incremental paid-loss triangle.
#' @param B Number of replications.
#' @param level Interval level.
#' @return A `cnbg_reserve` object carrying the second-moment estimate.
#' @examples
#' \donttest{
#' set.seed(1)
#' if (requireNamespace("tweedie", quietly = TRUE))
#'   bootstrap_cnbg_refit(simulate_cnbg()$triangle, B = 200)
#' }
#' @export
bootstrap_cnbg_refit <- function(triangle, B = 4000L, level = 0.95) {
  if (!requireNamespace("tweedie", quietly = TRUE))
    stop("bootstrap_cnbg_refit() requires the 'tweedie' package")
  ms <- cnbg_stage1(triangle, "free"); fc <- fit_cnbg_wls(ms)
  I <- nrow(triangle); J <- ncol(triangle); obs <- .cnbg_obs_mask(I, J); fut <- !obs
  tp <- function(nu0, kap, Ahat) {
    p <- 1 + (nu0 / kap) / (Ahat + nu0 / kap); p <- min(max(p, 1.01), 1.99)
    list(p = p, phi = max((Ahat * nu0 + nu0^2 / kap) / nu0^p, 1e-6)) }
  rc <- function(nu0, kap, Ahat) {
    if (nu0 <= 0) return(0); t <- tp(nu0, kap, Ahat)
    tryCatch(tweedie::rtweedie(1, power = t$p, mu = nu0, phi = t$phi), error = function(e) nu0) }
  quiet <- function(expr) withCallingHandlers(tryCatch(expr, error = function(e) NULL),
                                              warning = function(w) invokeRestart("muffleWarning"))
  cl_est <- sum(ms$nu_mat[fut]); reserves <- numeric(B)
  for (b in seq_len(B)) {
    X_b <- matrix(NA_real_, I, J)
    for (i in seq_len(I)) for (j in seq_len(J)) if (obs[i, j]) X_b[i, j] <- rc(ms$nu_mat[i, j], fc$kappa, fc$A)
    ms_b <- quiet(cnbg_stage1(X_b, "free"))
    fc_b <- if (!is.null(ms_b) && ms_b$converged) quiet(fit_cnbg_wls(ms_b)) else NULL
    if (is.null(fc_b)) { reserves[b] <- cl_est; next }
    R_b <- 0
    for (i in seq_len(I)) for (j in seq_len(J)) if (fut[i, j]) R_b <- R_b + rc(ms_b$nu_mat[i, j], fc_b$kappa, fc_b$A)
    reserves[b] <- R_b
  }
  a <- (1 - level) / 2; q <- unname(quantile(reserves, c(a, 1 - a)))
  structure(list(reserve = reserves, mean = mean(reserves), median = median(reserves), sd = sd(reserves),
                 pi = q, width = q[2] - q[1], level = level, cl_est = cl_est,
                 kappa = fc$kappa, A = fc$A, at_cap = fc$at_cap, method = "Frequentist CNBG"),
            class = "cnbg_reserve")
}

## ===========================================================================
## KAPPA POSTERIOR DIAGNOSTIC
## ===========================================================================

.cnbg_kl_gamma <- function(a1, b1, a2, b2)
  (a1 - a2) * digamma(a1) - lgamma(a1) + lgamma(a2) + a2 * (log(b1) - log(b2)) + a1 * (b2 - b1) / b1

.cnbg_prior_moments <- function(pr) {
  q <- qgamma(c(0.025, 0.975), shape = pr[["shape"]], rate = pr[["rate"]])
  list(shape = pr[["shape"]], rate = pr[["rate"]], mean = pr[["shape"]] / pr[["rate"]],
       sd = sqrt(pr[["shape"]]) / pr[["rate"]], q025 = q[1], q975 = q[2], ci_width = diff(q))
}

#' Kappa posterior diagnostic
#'
#' Compares the posterior of `kappa_B` from the full-mean working fit with
#' its prior, as in the paper's Schedule P diagnostics: the contraction
#' `c` (ratio of the posterior to the prior 95 percent interval width), the
#' standardised mean shift `s`, the Kullback--Leibler divergence from the
#' prior to a Gamma moment-matched to the posterior, and its tier
#' (`prior-driven` below 1, then `weak`, `moderate` above 5, `strong` above
#' 20). The regime is `informative` when `c < 1`, `not_informative` when
#' `c >= 1`, and `no_bayes` when the working fit returned no usable posterior
#' (sampler failure, or Rhat above 1.05 or bulk ESS below 100 for the
#' `nuts` engine).
#'
#' The paper reads this as a trust diagnostic for amount-only inference, not
#' as a method chooser: under free levels `kappa_B` is identified only
#' through the cell shock and only below a threshold (Remark 4.5), so an
#' uninformative verdict is the expected outcome on most paid triangles.
#'
#' @param fit A `cnbg_fit` from [fit_cnbg()] with `diagnose = TRUE`, a
#'   `cnbg_working` object with `mean = "full"`, or a `cnbg_fit` from
#'   [fit_cnbg_exact()].
#' @return An object of class `cnbg_kappa_diag` with `regime`, `tier`,
#'   `c`, `s`, `KL`, `informative` and the posterior summaries used.
#' @examples
#' set.seed(1)
#' tri <- simulate_cnbg(kappa_r = 5, alpha = 3)$triangle
#' diagnose_kappa(fit_cnbg(tri))
#' @export
diagnose_kappa <- function(fit) {
  none <- function(msg) structure(list(regime = "no_bayes", tier = NA_character_, c = NA_real_, s = NA_real_,
                                       KL = NA_real_, informative = FALSE, message = msg), class = "cnbg_kappa_diag")
  if (inherits(fit, "cnbg_working")) {
    if (fit$mean != "full") return(none("needs the full-mean working fit"))
    if (!.cnbg_working_usable(fit)) return(none("no usable posterior"))
    pr <- fit$prior; pm <- fit$adj$kappa_mean; q025 <- fit$adj$kappa[1]; q975 <- fit$adj$kappa[3]
  } else if (identical(fit$method, "exact")) {
    pr <- fit$prior; pm <- fit$kappa_mean; q025 <- fit$kappa_q025; q975 <- fit$kappa_q975
  } else {
    if (is.null(fit$full)) return(none("fit_cnbg() was called with diagnose = FALSE"))
    if (!.cnbg_working_usable(fit$full)) return(none("no usable posterior"))
    pr <- fit$prior; pm <- fit$full$adj$kappa_mean; q025 <- fit$full$adj$kappa[1]; q975 <- fit$full$adj$kappa[3]
  }
  if (any(!is.finite(c(pm, q025, q975))) || pm <= 0) return(none("non-finite posterior summary"))
  prior <- .cnbg_prior_moments(pr)
  w <- q975 - q025; sd_ <- w / 3.92; v <- sd_^2
  cc <- w / prior$ci_width; s <- (pm - prior$mean) / prior$sd
  KL <- if (v > 0) tryCatch(.cnbg_kl_gamma(pm^2 / v, pm / v, prior$shape, prior$rate),
                            error = function(e) NA_real_, warning = function(w) NA_real_) else NA_real_
  tier <- if (is.na(KL)) NA_character_ else if (KL > 20) "strong" else if (KL > 5) "moderate" else if (KL > 1) "weak" else "prior-driven"
  structure(list(regime = if (cc < 1) "informative" else "not_informative", tier = tier,
                 c = cc, s = s, KL = KL, informative = cc < 1,
                 kappa_mean = pm, kappa_q025 = q025, kappa_q975 = q975, prior = pr, message = NULL),
            class = "cnbg_kappa_diag")
}

#' @export
print.cnbg_kappa_diag <- function(x, ...) {
  cat(sprintf("kappa_B diagnostic: regime = %s", x$regime))
  if (!is.na(x$tier)) cat(sprintf("   KL tier = %s", x$tier))
  cat("\n")
  if (!is.null(x$message)) cat("  ", x$message, "\n")
  if (!is.na(x$c))
    cat(sprintf("  contraction c = %.3f   mean shift s = %.2f   KL = %.2f   posterior median-range [%.2f, %.2f]\n",
                x$c, x$s, x$KL, x$kappa_q025, x$kappa_q975))
  invisible(x)
}
