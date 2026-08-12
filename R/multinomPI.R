# ===========================================================================
# multinomPI.R
# Model-agnostic predictive intervals for claims reserving via
# Dirichlet-Multinomial allocation.
#
# Reference:
# Van Oirbeek, R. and Verdonck, T. (2026). Model-Agnostic Predictive
# Intervals for Claims Reserving via Dirichlet-Multinomial Allocation.
# ===========================================================================


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Sample from a Dirichlet distribution
#' @keywords internal
rdirichlet_ <- function(n, alpha) {
  k <- length(alpha)
  out <- matrix(0, n, k)
  for (i in seq_len(n)) {
    y <- rgamma(k, shape = alpha, rate = 1)
    out[i, ] <- y / sum(y)
  }
  out
}

#' Dirichlet log-density of one observation under concentration `conc`
#'
#' @param w Observed proportions, summing to 1.
#' @param conc Scalar TOTAL concentration of this composition. For a
#'   subcomposition through horizon k this is c*F_k, not c (Lemma A.3).
#' @param p Base proportions of this composition, summing to 1.
#' @keywords internal
ddirichlet_log_ <- function(w, conc, p) {
  a <- conc * p
  sum((a - 1) * log(w + 1e-300)) + lgamma(sum(a)) - sum(lgamma(a))
}

#' Moment-reporting tier of Remark 4.5
#'
#' For (1-W)/W ~ BetaPrime(c(1-F), cF) the k-th moment is finite iff cF > k.
#' Returns the highest finite moment order (0 = none).
#' @keywords internal
moment_tier_ <- function(c_param, F_obs) {
  cf <- c_param * F_obs
  ifelse(cf > 4, 4L, ifelse(cf > 2, 2L, ifelse(cf > 1, 1L, 0L)))
}

# ===========================================================================
# DEVELOPMENT PROPORTIONS
# ===========================================================================

#' Estimate Development Proportions from a Run-Off Triangle
#'
#' @param triangle Incremental run-off triangle, NA for unobserved cells.
#' @return List with `pi_hat` (incremental, summing to 1) and `F_hat`
#'   (cumulative, with `F_hat[J] = 1`).
#' @details Any method producing development proportions may be used instead.
#'   The bootstrap and diagnostic functions accept a user-supplied `pi_hat`,
#'   which is what makes the framework method-agnostic.
#' @export
estimate_dev_proportions <- function(triangle) {
  J <- ncol(triangle)
  cum <- t(apply(triangle, 1, cumsum))
  cum[is.na(triangle)] <- NA
  f <- numeric(J - 1)
  for (j in seq_len(J - 1)) {
    ok <- which(!is.na(cum[, j]) & !is.na(cum[, j + 1]))
    f[j] <- if (length(ok) > 0) sum(cum[ok, j + 1]) / sum(cum[ok, j]) else 1
  }
  F_hat <- numeric(J); F_hat[J] <- 1
  for (j in (J - 1):1) F_hat[j] <- F_hat[j + 1] / f[j]
  list(pi_hat = c(F_hat[1], diff(F_hat)), F_hat = F_hat)
}

# ===========================================================================
# CONCENTRATION PARAMETER ESTIMATION
# ===========================================================================

#' Estimate the Dirichlet Concentration Parameter
#'
#' Partial-column moment estimator of the concentration `c`, corrected for the
#' subcomposition property.
#'
#' @param triangle Incremental run-off triangle.
#' @param pi_hat Development proportions; estimated from the triangle if NULL.
#' @return Scalar estimate of `c`. Returns 50 with a warning when the estimator
#'   is undefined (fewer than three rows at every horizon, i.e. I < 5).
#'
#' @details
#' For each horizon k the rows observed through k supply the subcomposition
#' W_ij = X_ij / sum_{l<=k} X_il, which by Lemma A.3 is Dirichlet with TOTAL
#' concentration c*F_k. The moment estimator of that quantity is therefore
#' divided by F_k to obtain an estimator of c. Omitting that division leaves
#' an estimator inconsistent for c.
#'
#' The estimator is consistent but converges slowly: E[c_hat]/c is about 1.89
#' at I = 7, 1.31 at I = 10 and 1.15 at I = 15, with the ratio scale-free in
#' c. Where the decision is material, prefer the posterior of
#' `bayesian_bootstrap()`, which is defined at every triangle size.
#'
#' @export
estimate_c <- function(triangle, pi_hat = NULL) {
  J <- ncol(triangle)
  if (is.null(pi_hat)) pi_hat <- estimate_dev_proportions(triangle)$pi_hat

  c_est <- numeric(0)
  for (k in 1:(J - 1)) {                       # k_min = 1, not 2
    rows_k <- which(apply(triangle[, 1:(k + 1), drop = FALSE], 1,
                          function(x) all(!is.na(x))))
    if (length(rows_k) < 3) next
    S  <- rowSums(triangle[rows_k, 1:(k + 1), drop = FALSE])
    ok <- rows_k[S > 0]; S <- S[S > 0]
    if (length(ok) < 3) next
    W  <- triangle[ok, 1:(k + 1), drop = FALSE] / S
    Fk <- sum(pi_hat[1:(k + 1)])
    if (Fk <= 0 || Fk > 1.01) next
    pi_k <- pi_hat[1:(k + 1)] / Fk
    v    <- apply(W, 2, var)
    # bracket = moment estimator of the SUBCOMPOSITION concentration c*F_k;
    # divide by F_k to get c (Lemma A.3, Prop. 4.4)
    cj <- (pi_k * (1 - pi_k) / (v + 1e-10) - 1) / Fk
    c_est <- c(c_est, cj[cj > 0 & is.finite(cj)])
  }

  if (length(c_est) == 0) {
    warning("Insufficient data to estimate c (needs three rows at some ",
            "horizon, i.e. I >= 5); returning 50. Intervals will not ",
            "reflect uncertainty in c. Consider bayesian_bootstrap().",
            call. = FALSE)
    return(50)
  }
  c_med <- median(c_est)
  if (c_med < 1) {
    warning(sprintf("Estimated c = %.2f truncated to 1: allocation is highly ",
                    c_med),
            "heterogeneous. See diagnose_c() and consider frailty models.",
            call. = FALSE)
    return(1)
  }
  c_med
}

# ===========================================================================
# PARAMETRIC BOOTSTRAP
# ===========================================================================

#' Multinomial Parametric Bootstrap for Predictive Intervals
#'
#' Conditional predictive distribution of the total outstanding reserve. The
#' observed triangle is held fixed and only the allocation proportion
#' W_i ~ Beta(c F_{I-i}, c(1-F_{I-i})) is sampled per accident year per
#' replication. Any vector of cumulative development proportions may be
#' supplied, from Chain-Ladder, Bornhuetter-Ferguson, Cape Cod or elsewhere.
#'
#' @param triangle Incremental run-off triangle, NA for unobserved cells.
#' @param pi_hat Development proportions summing to 1; estimated if NULL.
#' @param c_param Concentration; estimated if NULL.
#' @param B Bootstrap replications (default 10000).
#' @param level Confidence level (default 0.95).
#' @param tier_min Minimum moment tier for reporting the mean and standard
#'   error of an accident year (default 4, i.e. c*F > 4; the package uses
#'   c*F >= 5 as that boundary with a margin). Years below the tier are
#'   summarised by quantiles only.
#'
#' @return Object of class `multinomial_boot`. `reserve_mean` and
#'   `reserve_se` are NA when any reported accident year falls below
#'   `tier_min`; quantiles are always available.
#'
#' @details
#' This is Algorithm 1 of the paper: it satisfies the conditioning requirement
#' exactly and substitutes (c_hat, pi_hat) rather than integrating over them.
#' The cost of that substitution is measured at 6.5 coverage points at I = 10
#' and 2.1 at I = 30; where c_hat is not well determined, use
#' `bayesian_bootstrap()`.
#'
#' Under compound Poisson misspecification the construction is conservative:
#' the predictive standard deviation exceeds the truth by the factor
#' 1/sqrt(F_{I-i}) in the compound Poisson limit p -> 1, and by more for
#' p in (1,2).
#'
#' @export
multinomial_bootstrap <- function(triangle, pi_hat = NULL, c_param = NULL,
                                  B = 10000, level = 0.95, tier_min = 4L) {
  I <- nrow(triangle); J <- ncol(triangle)
  if (is.null(pi_hat)) {
    dev <- estimate_dev_proportions(triangle)
    pi_hat <- dev$pi_hat; F_hat <- dev$F_hat
  } else F_hat <- cumsum(pi_hat)
  if (is.null(c_param)) c_param <- estimate_c(triangle, pi_hat)

  X_obs <- rowSums(triangle, na.rm = TRUE)
  k_obs <- pmin(I - seq_len(I) + 1, J)
  F_obs <- F_hat[k_obs]
  open  <- which(k_obs < J & F_obs > 0 & F_obs < 1)

  reserve_cl <- sum(ifelse(k_obs < J, X_obs * (1 - F_obs) / F_obs, 0))
  tier <- rep(NA_integer_, I)
  tier[open] <- moment_tier_(c_param, F_obs[open])

  origin_res <- matrix(0, B, I)
  for (i in open) {
    W <- rbeta(B, c_param * F_obs[i], c_param * (1 - F_obs[i]))
    # clamp and KEEP: dropping the draw truncates the upper tail, which is
    # where a 95% interval is determined
    W <- pmin(pmax(W, 1e-12), 1 - 1e-12)
    origin_res[, i] <- X_obs[i] * (1 - W) / W
  }
  total_res <- rowSums(origin_res)

  a  <- (1 - level) / 2
  ci <- quantile(total_res, c(a, 1 - a))
  q5 <- quantile(total_res, c(.05, .25, .50, .75, .95), names = FALSE)

  ok_all <- length(open) == 0 || all(tier[open] >= tier_min)
  if (!ok_all)
    warning(sprintf(paste0("Accident year(s) %s have c*F below the reporting ",
                           "tier; mean and standard error are not defined and ",
                           "are returned as NA. Quantiles are unaffected. ",
                           "See diagnose_c()."),
                    paste(open[tier[open] < tier_min], collapse = ", ")),
            call. = FALSE)

  by_origin <- data.frame(
    origin = seq_len(I), observed = X_obs, F_obs = F_obs,
    cF = c_param * F_obs, tier = tier,
    reserve_cl = ifelse(k_obs < J, X_obs * (1 - F_obs) / F_obs, 0),
    reserve_med = apply(origin_res, 2, median),
    reserve_mean = ifelse(!is.na(tier) & tier >= tier_min,
                          colMeans(origin_res), NA_real_),
    reserve_se = ifelse(!is.na(tier) & tier >= tier_min,
                        apply(origin_res, 2, sd), NA_real_))

  structure(list(
    reserve_cl = reserve_cl,
    reserve_mean = if (ok_all) mean(total_res) else NA_real_,
    reserve_se   = if (ok_all) sd(total_res)   else NA_real_,
    reserve_median = q5[3], quantiles = q5,
    ci_lower = unname(ci[1]), ci_upper = unname(ci[2]),
    cv = if (ok_all) sd(total_res) / mean(total_res) else NA_real_,
    c_hat = c_param, pi_hat = pi_hat, F_hat = F_hat,
    tier = tier, tier_min = tier_min, all_reportable = ok_all,
    level = level, B = B, samples = total_res, by_origin = by_origin
  ), class = "multinomial_boot")
}

#' @export
print.multinomial_boot <- function(x, ...) {
  cat("Multinomial Parametric Bootstrap (conditional predictive)\n")
  cat(sprintf("  B = %d, level = %.0f%%\n", x$B, 100 * x$level))
  cat(sprintf("  Concentration c = %.1f\n", x$c_hat))
  d <- diagnose_c(pi_hat = x$pi_hat, c_hat = x$c_hat)
  cat(sprintf("  Operability bound c+(5) = %.1f  ->  %s\n",
              d$c_dagger5, if (d$operable) "all years reportable"
                           else "some years below tier"))
  cat("\n")
  cat(sprintf("  CL reserve:      %12.0f\n", x$reserve_cl))
  cat(sprintf("  Predictive med:  %12.0f\n", x$reserve_median))
  if (x$all_reportable) {
    cat(sprintf("  Bootstrap mean:  %12.0f\n", x$reserve_mean))
    cat(sprintf("  Bootstrap SE:    %12.0f\n", x$reserve_se))
    cat(sprintf("  CV:              %11.1f%%\n", 100 * x$cv))
  } else {
    cat("  Bootstrap mean/SE: not defined (see $by_origin$tier)\n")
  }
  cat(sprintf("  5/25/50/75/95:   %s\n",
              paste(format(round(x$quantiles), big.mark = ","),
                    collapse = " / ")))
  cat(sprintf("  %.0f%% PI:         [%.0f, %.0f]\n",
              100 * x$level, x$ci_lower, x$ci_upper))
  invisible(x)
}

# ===========================================================================
# DELTA METHOD APPROXIMATION
# ===========================================================================

#' Delta Method Variance for IBNP Reserves
#'
#' Closed-form process variance under the Dirichlet-Gamma model.
#'
#' @inheritParams multinomial_bootstrap
#' @return Data frame with per-origin `origin`, `X_obs`, `F_k`, `cF`,
#'   `reserve`, `process_se`, `cv`. `process_se` is NA where c*F <= 2, the
#'   variance not existing there.
#'
#' @details Var(S | X) ~ X^2 (1-F) / (F^3 (c+1)). We do not recommend this in
#'   preference to `multinomial_bootstrap()` at any c_hat: the moment
#'   expansion is least reliable in exactly the small-c*F regime where a
#'   shortcut is most tempting, and its coverage has not been evaluated.
#' @export
delta_method_var <- function(triangle, pi_hat = NULL, c_param = NULL) {
  I <- nrow(triangle); J <- ncol(triangle)
  if (is.null(pi_hat)) {
    dev <- estimate_dev_proportions(triangle)
    pi_hat <- dev$pi_hat; F_hat <- dev$F_hat
  } else F_hat <- cumsum(pi_hat)
  if (is.null(c_param)) c_param <- estimate_c(triangle, pi_hat)

  X_obs <- rowSums(triangle, na.rm = TRUE)
  k_obs <- pmin(I - seq_len(I) + 1, J)
  Fk    <- F_hat[k_obs]
  R     <- ifelse(k_obs < J & Fk > 0 & Fk < 1, X_obs * (1 - Fk) / Fk, 0)
  se    <- ifelse(k_obs < J & Fk > 0 & Fk < 1 & c_param * Fk > 2,
                  sqrt(X_obs^2 * (1 - Fk) / (Fk^3 * (c_param + 1))), NA_real_)
  data.frame(origin = seq_len(I), X_obs = X_obs, F_k = Fk,
             cF = c_param * Fk, reserve = R, process_se = se, cv = se / R)
}

# ===========================================================================
# REGULARISED (BAYESIAN) PREDICTIVE
# ===========================================================================

#' Regularised Conditional Predictive
#'
#' Algorithm 2 of the paper: the conditional predictive with `c` integrated
#' over its posterior rather than substituted. `pi` may optionally be sampled
#' as well.
#'
#' @inheritParams multinomial_bootstrap
#' @param B Posterior predictive samples (default 5000).
#' @param mu_c,sigma_c Prior mean and SD of log c (default log(50), 1).
#' @param n_mcmc,burnin MCMC length after burn-in, and burn-in.
#' @param lambda_pi Dirichlet random-walk concentration for the pi step.
#' @param bayesian_c_only If TRUE, only c is sampled and pi is plugged in
#'   (default TRUE, matching Algorithm 2 as stated in the paper).
#'
#' @return List with predictive quantiles, posterior summaries of c, and
#'   acceptance rates.
#'
#' @details
#' The likelihood is the exact subcomposition likelihood: a row observed
#' through development year k contributes a Dirichlet with TOTAL concentration
#' c*F_k (Lemma A.3), so partially developed rows enter without approximation.
#'
#' No conjugate update exists for pi. The Dirichlet is a full exponential
#' family in alpha = c*pi and not in pi, since the log-partition contributes
#' sum_j lgamma(c*pi_j), which is not a Dirichlet kernel in pi. pi is therefore
#' sampled by a Dirichlet random-walk Metropolis step with the Hastings
#' correction, under a uniform prior on the simplex.
#'
#' NO MOMENTS EXIST. Under the posterior mixture the predictive mean is the
#' integral of c(1-F)/(cF-1) against the posterior of c, which has a simple
#' pole at c = 1/F; a log-normal prior puts positive posterior density there,
#' so the integral diverges, as does every higher moment. Only quantiles are
#' returned, and this is a property of the construction rather than a
#' reporting convention.
#'
#' @export
bayesian_bootstrap <- function(triangle, B = 5000,
                               mu_c = log(50), sigma_c = 1,
                               n_mcmc = 5000, burnin = 1000,
                               prop_sd = 0.4, lambda_pi = 600,
                               bayesian_c_only = TRUE, level = 0.95) {
  I <- nrow(triangle); J <- ncol(triangle)
  dev <- estimate_dev_proportions(triangle)
  pi_hat <- dev$pi_hat; F_hat <- dev$F_hat
  X_obs <- rowSums(triangle, na.rm = TRUE)
  k_obs <- pmin(I - seq_len(I) + 1, J)

  # subcompositions: one per row with at least two observed, positive cells
  W_list <- list(); k_list <- integer(0)
  for (i in seq_len(I)) {
    oc <- which(!is.na(triangle[i, ]))
    if (length(oc) < 2) next
    if (!identical(oc, seq_along(oc))) next          # leading run only
    x <- triangle[i, oc]
    if (any(x <= 0)) next
    W_list[[length(W_list) + 1L]] <- as.numeric(x / sum(x))
    k_list <- c(k_list, length(oc))
  }
  if (length(W_list) < 2) {
    message("Insufficient rows for MCMC; falling back to the plug-in bootstrap.")
    return(multinomial_bootstrap(triangle, B = B, level = level))
  }

  loglik <- function(cc, pv) {
    if (!is.finite(cc) || cc <= 0 || any(pv <= 0)) return(-Inf)
    s <- 0
    for (r in seq_along(W_list)) {
      k  <- k_list[r]
      Fk <- sum(pv[1:k])                              # subcomposition total
      s  <- s + ddirichlet_log_(W_list[[r]], cc * Fk, pv[1:k] / Fk)
    }
    if (is.finite(s)) s else -Inf
  }
  ldir <- function(x, a) sum((a - 1) * log(x)) + lgamma(sum(a)) - sum(lgamma(a))

  logc <- log(max(estimate_c(triangle, pi_hat), 1))
  if (!is.finite(logc)) logc <- mu_c
  pv <- pi_hat
  ll <- loglik(exp(logc), pv)
  if (!is.finite(ll)) {
    message("Likelihood not finite at the starting value; falling back.")
    return(multinomial_bootstrap(triangle, B = B, level = level))
  }

  keep <- max(n_mcmc, 1L)
  c_draws <- numeric(keep); pi_draws <- matrix(NA_real_, keep, J)
  acc_c <- 0L; acc_pi <- 0L; kk <- 0L
  for (it in seq_len(n_mcmc + burnin)) {
    lc2 <- logc + rnorm(1, 0, prop_sd)
    ll2 <- loglik(exp(lc2), pv)
    lp  <- dnorm(lc2, mu_c, sigma_c, log = TRUE) -
           dnorm(logc, mu_c, sigma_c, log = TRUE)
    if (is.finite(ll2) && log(runif(1)) < (ll2 - ll + lp)) {
      logc <- lc2; ll <- ll2; if (it > burnin) acc_c <- acc_c + 1L
    }
    if (!bayesian_c_only) {
      pv2 <- as.vector(rdirichlet_(1, pmax(lambda_pi * pv, 1e-8)))
      if (all(pv2 > 0)) {
        ll2 <- loglik(exp(logc), pv2)
        if (is.finite(ll2)) {
          hast <- ldir(pv,  pmax(lambda_pi * pv2, 1e-8)) -
                  ldir(pv2, pmax(lambda_pi * pv,  1e-8))
          if (log(runif(1)) < (ll2 - ll + hast)) {
            pv <- pv2; ll <- ll2; if (it > burnin) acc_pi <- acc_pi + 1L
          }
        }
      }
    }
    if (it > burnin) { kk <- kk + 1L; c_draws[kk] <- exp(logc); pi_draws[kk, ] <- pv }
  }
  ok <- is.finite(c_draws) & c_draws > 0
  c_draws <- c_draws[ok]; pi_draws <- pi_draws[ok, , drop = FALSE]

  idx <- sample.int(nrow(pi_draws), B, replace = TRUE)
  total_res <- numeric(B)
  for (b in seq_len(B)) {
    cb <- c_draws[idx[b]]; Fb <- cumsum(pi_draws[idx[b], ])
    tot <- 0
    for (i in seq_len(I)) {
      if (k_obs[i] >= J) next
      Fi <- min(max(Fb[k_obs[i]], 1e-6), 1 - 1e-6)
      W  <- rbeta(1, cb * Fi, cb * (1 - Fi))
      W  <- min(max(W, 1e-12), 1 - 1e-12)
      tot <- tot + X_obs[i] * (1 - W) / W
    }
    total_res[b] <- tot
  }
  total_res <- total_res[is.finite(total_res) & total_res > 0]

  a  <- (1 - level) / 2
  ci <- quantile(total_res, c(a, 1 - a))
  list(quantiles = quantile(total_res, c(.05, .25, .50, .75, .95), names = FALSE),
       reserve_median = median(total_res),
       ci_lower = unname(ci[1]), ci_upper = unname(ci[2]),
       c_posterior_mean = mean(c_draws), c_posterior_sd = sd(c_draws),
       accept_c = acc_c / n_mcmc, accept_pi = acc_pi / n_mcmc,
       n_rows_used = length(W_list), samples = total_res,
       note = paste("No moments exist under the posterior mixture;",
                    "quantiles only. See ?bayesian_bootstrap."))
}

# ===========================================================================
# DIAGNOSTICS
# ===========================================================================

#' Concentration Diagnostics: Operability Bound and Heterogeneity Screen
#'
#' Two distinct questions, with different answers.
#'
#' @param triangle Incremental run-off triangle (or NULL if `pi_hat` and
#'   `c_hat` are supplied).
#' @param pi_hat,c_hat Optionally supplied directly.
#' @param tau Moment tier required (default 5, i.e. the c*F > 4 boundary with
#'   a margin).
#' @param screen_t Heterogeneity screen threshold (default 30).
#'
#' @return Object of class `c_diagnostic`.
#'
#' @details
#' \strong{Operability.} The k-th predictive moment exists iff c*F_{I-i} > k,
#' so the requirement on c is portfolio-specific: c >= c_dagger(tau) =
#' tau / pi_0, with pi_0 the share of the ultimate reported in the first
#' development year. On standard benchmarks c_dagger(5) is 44.6 for RAA, 72.2
#' for Taylor-Ashe and 632.5 for Mortgage. A fixed threshold such as
#' c_hat < 30 implicitly assumes pi_0 ~ 1/6 and is not used here.
#'
#' \strong{Heterogeneity.} Separately, a small c_hat indicates that allocation
#' varies across accident years, favouring a model with explicit frailty. This
#' is a screen and not a test: because c_hat is upward-biased at small I, a
#' portfolio sitting exactly at the threshold is flagged with probability
#' 0.32-0.45 rather than 0.50, and at I = 7 the rule detects genuinely
#' heterogeneous development on barely half of triangles.
#'
#' The two questions are independent: a portfolio can clear one and fail the
#' other.
#'
#' @export
diagnose_c <- function(triangle = NULL, pi_hat = NULL, c_hat = NULL,
                       tau = 5, screen_t = 30) {
  if (is.null(pi_hat)) pi_hat <- estimate_dev_proportions(triangle)$pi_hat
  if (is.null(c_hat))  c_hat  <- estimate_c(triangle, pi_hat)
  pi0 <- pi_hat[1]
  cd5 <- tau / pi0; cd2 <- 2 / pi0
  operable <- c_hat >= cd5
  flagged  <- c_hat < screen_t
  structure(list(
    c_hat = c_hat, pi_0 = pi0, c_dagger2 = cd2, c_dagger5 = cd5,
    operable = operable, screen_flag = flagged, screen_t = screen_t,
    message = paste0(
      sprintf("c_hat = %.1f;  pi_0 = %.4f;  c+(2) = %.1f, c+(5) = %.1f.\n",
              c_hat, pi0, cd2, cd5),
      sprintf("Operability: %s.\n", if (operable)
        "all accident years support full moment reporting" else
        "at least one accident year falls below the tier; report quantiles or exclude"),
      sprintf("Screen (c_hat < %g): %s. Note the screen's error rates: at I = 7 it\n",
              screen_t, if (flagged) "FLAGGED" else "not flagged"),
      "detects c < 30 on about half of triangles, and a portfolio at the\n",
      "threshold is flagged with probability 0.32-0.45 rather than 0.50.")
  ), class = "c_diagnostic")
}

#' @export
print.c_diagnostic <- function(x, ...) { cat(x$message, "\n"); invisible(x) }
