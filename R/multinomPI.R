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
#'
#' @param n Number of samples.
#' @param alpha Numeric vector of concentration parameters.
#' @return Matrix with \code{n} rows and \code{length(alpha)} columns.
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

#' Dirichlet log-likelihood for a single observation
#'
#' @param w Numeric vector of observed proportions (sums to 1).
#' @param c_val Scalar concentration parameter.
#' @param pi_vec Numeric vector of base proportions.
#' @return Scalar log-likelihood.
#' @keywords internal
dirichlet_loglik <- function(w, c_val, pi_vec) {
  a <- c_val * pi_vec
  sum((a - 1) * log(w + 1e-300)) + lgamma(sum(a)) - sum(lgamma(a))
}


# ===========================================================================
# DEVELOPMENT PROPORTIONS
# ===========================================================================

#' Estimate Development Proportions from a Run-Off Triangle
#'
#' Computes cumulative development proportions \eqn{F} and incremental
#' proportions \eqn{\pi} from an incremental run-off triangle using
#' volume-weighted Chain-Ladder link ratios.
#'
#' @param triangle Numeric matrix. Incremental run-off triangle with
#'   \code{NA} for unobserved cells. Rows are accident years, columns
#'   are development periods.
#'
#' @return A list with:
#' \describe{
#'   \item{pi_hat}{Numeric vector of incremental development proportions
#'     summing to 1.}
#'   \item{F_hat}{Numeric vector of cumulative development proportions
#'     with \code{F_hat[J] = 1}.}
#' }
#'
#' @details
#' Any method producing development proportions can be used instead of
#' Chain-Ladder. The bootstrap and diagnostic functions accept
#' user-supplied \code{pi_hat} vectors, making the framework
#' method-agnostic.
#'
#' @export
estimate_dev_proportions <- function(triangle) {
  I <- nrow(triangle)
  J <- ncol(triangle)

  cum <- t(apply(triangle, 1, cumsum))
  cum[is.na(triangle)] <- NA

  f <- numeric(J - 1)
  for (j in seq_len(J - 1)) {
    ok <- which(!is.na(cum[, j]) & !is.na(cum[, j + 1]))
    f[j] <- if (length(ok) > 0) sum(cum[ok, j + 1]) / sum(cum[ok, j]) else 1
  }

  F_hat <- numeric(J)
  F_hat[J] <- 1
  for (j in (J - 1):1) F_hat[j] <- F_hat[j + 1] / f[j]

  pi_hat <- c(F_hat[1], diff(F_hat))

  list(pi_hat = pi_hat, F_hat = F_hat)
}


# ===========================================================================
# CONCENTRATION PARAMETER ESTIMATION
# ===========================================================================

#' Estimate the Dirichlet Concentration Parameter
#'
#' Estimates the Dirichlet concentration parameter \eqn{c} from an
#' incremental run-off triangle using a partial-column moment method.
#' Works for any triangle shape including \eqn{I = J}.
#'
#' @param triangle Numeric matrix. Incremental run-off triangle.
#' @param pi_hat Numeric vector. Development proportions. If \code{NULL},
#'   estimated from the triangle.
#'
#' @return Scalar estimate of \eqn{c}. Returns 50 (diffuse default) if
#'   insufficient data.
#'
#' @details
#' For each development horizon \eqn{k = 2, \ldots, J-1}, partial
#' proportions are computed from the rows observed through period \eqn{k},
#' their variance is matched to the Dirichlet formula, and the median
#' across columns and horizons is returned.
#'
#' The diagnostic threshold \eqn{\hat{c} < 30} signals substantial
#' accident-year heterogeneity.
#'
#' @export
estimate_c <- function(triangle, pi_hat = NULL) {
  I <- nrow(triangle)
  J <- ncol(triangle)

  if (is.null(pi_hat)) {
    pi_hat <- estimate_dev_proportions(triangle)$pi_hat
  }

  c_est <- numeric(0)

  for (k in 2:(J - 1)) {
    rows_k <- which(apply(triangle[, 1:(k + 1), drop = FALSE], 1,
                           function(x) all(!is.na(x))))
    if (length(rows_k) < 3) next

    S <- rowSums(triangle[rows_k, 1:(k + 1), drop = FALSE])
    ok <- rows_k[S > 0]
    S <- S[S > 0]
    if (length(ok) < 3) next

    W <- triangle[ok, 1:(k + 1), drop = FALSE] / S

    Fk <- sum(pi_hat[1:(k + 1)])
    if (Fk <= 0 || Fk > 1.01) next
    pi_k <- pi_hat[1:(k + 1)] / Fk

    v <- apply(W, 2, var)
    cj <- pi_k * (1 - pi_k) / (v + 1e-10) - 1
    cj <- cj[cj > 0 & is.finite(cj)]
    c_est <- c(c_est, cj)
  }

  if (length(c_est) == 0) return(50)
  max(median(c_est), 1)
}


# ===========================================================================
# PARAMETRIC BOOTSTRAP
# ===========================================================================

#' Multinomial Parametric Bootstrap for Predictive Intervals
#'
#' Model-agnostic predictive intervals for total outstanding reserves.
#' Any development proportions can be supplied; the framework adds a
#' single Dirichlet concentration parameter \eqn{c} to produce full
#' predictive distributions.
#'
#' @param triangle Numeric matrix. Incremental run-off triangle.
#' @param pi_hat Numeric vector. Development proportions. If \code{NULL},
#'   estimated via Chain-Ladder.
#' @param c_param Scalar. Dirichlet concentration. If \code{NULL},
#'   estimated from the triangle.
#' @param B Integer. Bootstrap replications (default 10000).
#' @param level Numeric. Confidence level (default 0.95).
#'
#' @return A list with:
#' \describe{
#'   \item{reserve_cl}{Chain-Ladder point estimate.}
#'   \item{reserve_mean}{Bootstrap mean.}
#'   \item{reserve_se}{Bootstrap standard error.}
#'   \item{ci_lower, ci_upper}{Predictive interval bounds.}
#'   \item{cv}{Coefficient of variation.}
#'   \item{c_hat}{Concentration parameter used.}
#'   \item{pi_hat, F_hat}{Development proportions used.}
#'   \item{samples}{Vector of bootstrap reserve samples.}
#'   \item{by_origin}{Per-accident-year breakdown.}
#' }
#'
#' @details
#' The bootstrap samples ultimates from
#' \eqn{S_i \sim \mathrm{Gamma}(\alpha, \alpha / \hat{S}_i)} and
#' allocations from \eqn{W_i \sim \mathrm{Dir}(c \hat{\pi})}, truncates
#' to the run-off shape, re-estimates development factors, and computes
#' reserves. Process risk is captured through Gamma-Dirichlet sampling;
#' parameter risk through re-estimation.
#'
#' @export
multinomial_bootstrap <- function(triangle, pi_hat = NULL, c_param = NULL,
                                   B = 10000, level = 0.95) {
  I <- nrow(triangle)
  J <- ncol(triangle)

  if (is.null(pi_hat)) {
    dev <- estimate_dev_proportions(triangle)
    pi_hat <- dev$pi_hat
    F_hat  <- dev$F_hat
  } else {
    F_hat <- cumsum(pi_hat)
  }

  if (is.null(c_param)) c_param <- estimate_c(triangle, pi_hat)

  # Ultimates
  X_obs <- rowSums(triangle, na.rm = TRUE)
  k_obs <- pmin(I - seq_len(I) + 1, J)
  S_ult <- X_obs / F_hat[k_obs]

  # Gamma shape
  cv_ult <- sd(S_ult) / mean(S_ult)
  if (cv_ult < 0.01) cv_ult <- 0.1
  alpha_tot <- 1 / cv_ult^2

  # CL point estimate
  reserve_cl <- sum(ifelse(k_obs < J, X_obs * (1 - F_hat[k_obs]) / F_hat[k_obs], 0))

  # Bootstrap
  total_res <- numeric(B)
  origin_res <- matrix(0, B, I)
  alpha_dir <- c_param * pi_hat

  for (b in seq_len(B)) {
    X_star <- matrix(0, I, J)
    for (i in seq_len(I)) {
      S <- rgamma(1, shape = alpha_tot, rate = alpha_tot / S_ult[i])
      W <- rdirichlet_(1, alpha_dir)[1, ]
      X_star[i, ] <- S * W
    }
    for (i in seq_len(I)) {
      k <- I - i + 1
      if (k < J) X_star[i, (k + 1):J] <- NA
    }

    dev_s <- estimate_dev_proportions(X_star)
    Fs <- dev_s$F_hat

    for (i in seq_len(I)) {
      k <- I - i + 1
      if (k < J) {
        Xo <- sum(X_star[i, 1:k], na.rm = TRUE)
        if (Xo > 0 && Fs[k] > 0 && Fs[k] < 1) {
          origin_res[b, i] <- Xo * (1 - Fs[k]) / Fs[k]
        }
      }
    }
    total_res[b] <- sum(origin_res[b, ])
  }

  a <- (1 - level) / 2
  ci <- quantile(total_res, c(a, 1 - a))

  by_origin <- data.frame(
    origin       = seq_len(I),
    observed     = X_obs,
    ultimate     = S_ult,
    reserve_cl   = ifelse(k_obs < J, X_obs * (1 - F_hat[k_obs]) / F_hat[k_obs], 0),
    reserve_mean = colMeans(origin_res),
    reserve_se   = apply(origin_res, 2, sd)
  )

  structure(list(
    reserve_cl   = reserve_cl,
    reserve_mean = mean(total_res),
    reserve_se   = sd(total_res),
    ci_lower     = unname(ci[1]),
    ci_upper     = unname(ci[2]),
    cv           = sd(total_res) / mean(total_res),
    c_hat        = c_param,
    pi_hat       = pi_hat,
    F_hat        = F_hat,
    level        = level,
    B            = B,
    samples      = total_res,
    by_origin    = by_origin
  ), class = "multinomial_boot")
}


#' @export
print.multinomial_boot <- function(x, ...) {
  cat("Multinomial Parametric Bootstrap\n")
  cat(sprintf("  B = %d, level = %.0f%%\n", x$B, 100 * x$level))
  cat(sprintf("  Concentration c = %.1f", x$c_hat))
  if (x$c_hat < 30) cat("  [WARNING: c < 30, consider richer models]")
  cat("\n\n")
  cat(sprintf("  CL reserve:       %12.0f\n", x$reserve_cl))
  cat(sprintf("  Bootstrap mean:   %12.0f\n", x$reserve_mean))
  cat(sprintf("  Bootstrap SE:     %12.0f\n", x$reserve_se))
  cat(sprintf("  CV:               %11.1f%%\n", 100 * x$cv))
  cat(sprintf("  %.0f%% PI:       [%12.0f, %12.0f]\n",
              100 * x$level, x$ci_lower, x$ci_upper))
  invisible(x)
}


# ===========================================================================
# DELTA METHOD APPROXIMATION
# ===========================================================================

#' Delta Method Variance for IBNP Reserves
#'
#' Closed-form process variance under the Dirichlet-Gamma model.
#' Avoids bootstrapping when speed is needed.
#'
#' @param triangle Numeric matrix. Incremental run-off triangle.
#' @param pi_hat Numeric vector. Development proportions (or \code{NULL}).
#' @param c_param Scalar. Concentration parameter (or \code{NULL}).
#'
#' @return Data frame with per-origin columns: \code{origin},
#'   \code{X_obs}, \code{F_k}, \code{reserve}, \code{process_se},
#'   \code{cv}.
#'
#' @details
#' \deqn{\mathrm{Var}(S_i^{\mathrm{IBNP}} \mid X_i^{\mathrm{obs}})
#'   \approx (X_i^{\mathrm{obs}})^2 \cdot \frac{1 - F}{F^3(c+1)}}
#'
#' @export
delta_method_var <- function(triangle, pi_hat = NULL, c_param = NULL) {
  I <- nrow(triangle)
  J <- ncol(triangle)

  if (is.null(pi_hat)) {
    dev <- estimate_dev_proportions(triangle)
    pi_hat <- dev$pi_hat
    F_hat  <- dev$F_hat
  } else {
    F_hat <- cumsum(pi_hat)
  }

  if (is.null(c_param)) c_param <- estimate_c(triangle, pi_hat)

  X_obs <- rowSums(triangle, na.rm = TRUE)
  k_obs <- pmin(I - seq_len(I) + 1, J)

  out <- data.frame(origin = seq_len(I))
  out$X_obs <- X_obs
  out$F_k <- F_hat[k_obs]
  out$reserve <- out$process_se <- out$cv <- 0

  for (i in seq_len(I)) {
    Fk <- F_hat[k_obs[i]]
    if (k_obs[i] < J && Fk > 0 && Fk < 1) {
      R <- X_obs[i] * (1 - Fk) / Fk
      V <- X_obs[i]^2 * (1 - Fk) / (Fk^3 * (c_param + 1))
      out$reserve[i] <- R
      out$process_se[i] <- sqrt(V)
      out$cv[i] <- sqrt(V) / R
    }
  }

  out
}


# ===========================================================================
# BAYESIAN PREDICTIVE BOOTSTRAP
# ===========================================================================

#' Bayesian Predictive Bootstrap for Claims Reserves
#'
#' Places priors on \eqn{c} and \eqn{\pi} and uses MCMC to produce
#' posterior predictive reserve distributions.
#'
#' @param triangle Numeric matrix. Incremental run-off triangle.
#' @param B Integer. Predictive samples (default 5000).
#' @param mu_c Scalar. Prior mean of \eqn{\log c} (default \code{log(50)}).
#' @param sigma_c Scalar. Prior SD of \eqn{\log c} (default 1).
#' @param a0 Scalar. Dirichlet prior for \eqn{\pi} (default 1).
#' @param n_mcmc Integer. MCMC iterations after burn-in (default 2000).
#' @param burnin Integer. Burn-in (default 500).
#' @param bayesian_c_only Logical. If \code{TRUE}, only \eqn{c} is
#'   sampled; \eqn{\pi} fixed at plug-in.
#' @param level Numeric. Confidence level (default 0.95).
#'
#' @return A list with: \code{reserve_mean}, \code{reserve_se},
#'   \code{ci_lower}, \code{ci_upper}, \code{c_posterior_mean},
#'   \code{c_posterior_sd}, \code{accept_rate}, \code{samples}.
#'
#' @export
bayesian_bootstrap <- function(triangle, B = 5000,
                                mu_c = log(50), sigma_c = 1, a0 = 1,
                                n_mcmc = 2000, burnin = 500,
                                bayesian_c_only = FALSE, level = 0.95) {

  I <- nrow(triangle)
  J <- ncol(triangle)

  dev <- estimate_dev_proportions(triangle)
  pi_hat <- dev$pi_hat
  F_hat  <- dev$F_hat

  # Usable rows (>= 3 observed columns)
  min_cols <- min(3, J)
  usable <- which(apply(triangle, 1, function(x) sum(!is.na(x)) >= min_cols))

  if (length(usable) < 2) {
    message("Insufficient rows for MCMC; falling back to plug-in bootstrap.")
    return(multinomial_bootstrap(triangle, B = B, level = level))
  }

  # Partial-row data
  W_list <- list()
  cols_list <- list()
  for (i in usable) {
    oc <- which(!is.na(triangle[i, ]))
    s <- sum(triangle[i, oc])
    if (s <= 0) next
    Fo <- sum(pi_hat[oc])
    if (Fo <= 0) next
    W_list[[length(W_list) + 1]] <- triangle[i, oc] / s
    cols_list[[length(cols_list) + 1]] <- oc
  }

  if (length(W_list) < 2) {
    message("Insufficient valid rows for MCMC; falling back to plug-in bootstrap.")
    return(multinomial_bootstrap(triangle, B = B, level = level))
  }

  # Ultimates
  X_obs <- rowSums(triangle, na.rm = TRUE)
  k_obs <- pmin(I - seq_len(I) + 1, J)
  S_ult <- X_obs / F_hat[k_obs]
  cv_ult <- sd(S_ult) / mean(S_ult)
  if (cv_ult < 0.01) cv_ult <- 0.1
  alpha_tot <- 1 / cv_ult^2

  # MCMC
  c_cur <- estimate_c(triangle, pi_hat)
  if (c_cur < 1) c_cur <- 10
  pi_cur <- pi_hat

  c_samp  <- numeric(n_mcmc)
  pi_samp <- matrix(0, n_mcmc, J)
  prop_sd <- 0.5
  n_acc   <- 0L

  for (iter in seq_len(n_mcmc + burnin)) {

    # Update pi | c
    if (!bayesian_c_only) {
      col_w <- numeric(J)
      for (idx in seq_along(W_list)) {
        oc <- cols_list[[idx]]
        col_w[oc] <- col_w[oc] + c_cur * W_list[[idx]]
      }
      pi_cur <- rdirichlet_(1, pmax(a0 + col_w, 0.01))[1, ]
    }

    # Update c | pi (Metropolis)
    lc_prop <- rnorm(1, log(c_cur), prop_sd)
    c_prop  <- exp(lc_prop)

    ll_cur <- ll_prop <- 0
    for (idx in seq_along(W_list)) {
      oc <- cols_list[[idx]]
      Fo <- sum(pi_cur[oc])
      if (Fo <= 0) next
      pp <- pi_cur[oc] / Fo
      ll_cur  <- ll_cur  + dirichlet_loglik(W_list[[idx]], c_cur, pp)
      ll_prop <- ll_prop + dirichlet_loglik(W_list[[idx]], c_prop, pp)
    }

    lp_cur  <- dnorm(log(c_cur),  mu_c, sigma_c, log = TRUE)
    lp_prop <- dnorm(log(c_prop), mu_c, sigma_c, log = TRUE)

    if (log(runif(1)) < (ll_prop + lp_prop) - (ll_cur + lp_cur)) {
      c_cur <- c_prop
      if (iter > burnin) n_acc <- n_acc + 1L
    }

    if (iter > burnin) {
      idx_out <- iter - burnin
      c_samp[idx_out]    <- c_cur
      pi_samp[idx_out, ] <- pi_cur
    }
  }

  # Predictive bootstrap from posterior
  thin <- round(seq(1, n_mcmc, length.out = B))
  total_res <- numeric(B)

  for (b in seq_len(B)) {
    cb  <- c_samp[thin[b]]
    pib <- if (bayesian_c_only) pi_hat else pi_samp[thin[b], ]

    X_star <- matrix(0, I, J)
    for (i in seq_len(I)) {
      S <- rgamma(1, shape = alpha_tot, rate = alpha_tot / S_ult[i])
      W <- rdirichlet_(1, cb * pib)[1, ]
      X_star[i, ] <- S * W
    }
    for (i in seq_len(I)) {
      k <- I - i + 1
      if (k < J) X_star[i, (k + 1):J] <- NA
    }

    dev_s <- estimate_dev_proportions(X_star)
    Fs <- dev_s$F_hat
    res <- 0
    for (i in seq_len(I)) {
      k <- I - i + 1
      if (k < J) {
        Xo <- sum(X_star[i, 1:k], na.rm = TRUE)
        if (Xo > 0 && Fs[k] > 0 && Fs[k] < 1) {
          res <- res + Xo * (1 - Fs[k]) / Fs[k]
        }
      }
    }
    total_res[b] <- res
  }

  a <- (1 - level) / 2
  ci <- quantile(total_res, c(a, 1 - a))

  list(
    reserve_mean     = mean(total_res),
    reserve_se       = sd(total_res),
    ci_lower         = unname(ci[1]),
    ci_upper         = unname(ci[2]),
    c_posterior_mean = mean(c_samp),
    c_posterior_sd   = sd(c_samp),
    accept_rate      = n_acc / n_mcmc,
    samples          = total_res
  )
}


# ===========================================================================
# DIAGNOSTIC
# ===========================================================================

#' Concentration Parameter Diagnostic
#'
#' Estimates \eqn{c} and assesses whether the homogeneous development
#' assumption is appropriate.
#'
#' @param triangle Numeric matrix. Incremental run-off triangle.
#' @param threshold Scalar. Diagnostic threshold (default 30).
#'
#' @return A list with: \code{c_hat}, \code{threshold}, \code{adequate}
#'   (logical), \code{message}.
#'
#' @export
diagnose_c <- function(triangle, threshold = 30) {
  c_hat <- estimate_c(triangle)
  adequate <- c_hat >= threshold

  msg <- if (adequate) {
    sprintf("c_hat = %.1f >= %d: stable development. Multinomial framework appropriate.",
            c_hat, threshold)
  } else {
    sprintf("c_hat = %.1f < %d: heterogeneous development. Consider frailty models.",
            c_hat, threshold)
  }

  structure(list(
    c_hat     = c_hat,
    threshold = threshold,
    adequate  = adequate,
    message   = msg
  ), class = "c_diagnostic")
}

#' @export
print.c_diagnostic <- function(x, ...) {
  cat(x$message, "\n")
  invisible(x)
}
