#' @title Multinomial Parametric Inference for Claims Reserving
#' @description
#' Model-agnostic predictive intervals for claims reserving via
#' Dirichlet-Multinomial allocation. Any macro-level reserving method
#' producing development proportions can be embedded in this framework.
#'
#' @references
#' Van Oirbeek, R. and Verdonck, T. (2026). Model-Agnostic Predictive
#' Intervals for Claims Reserving via Dirichlet-Multinomial Allocation.
#' Working paper.
#' @name multinomPI
NULL

# =============================================================================
# DEVELOPMENT PROPORTIONS
# =============================================================================

#' Estimate Development Proportions from a Run-Off Triangle
#'
#' Computes cumulative development proportions F and incremental
#' proportions pi from an incremental run-off triangle using
#' volume-weighted Chain-Ladder link ratios.
#'
#' @param triangle Numeric matrix. Incremental run-off triangle with
#'   \code{NA} for unobserved cells. Rows are accident years, columns
#'   are development periods.
#'
#' @return A list with components:
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
estimate_development_proportions <- function(triangle) {
  I <- nrow(triangle)
  J <- ncol(triangle)

  cum_triangle <- t(apply(triangle, 1, cumsum))
  cum_triangle[is.na(triangle)] <- NA

  link_ratios <- numeric(J - 1)
  for (j in 1:(J - 1)) {
    valid <- which(!is.na(cum_triangle[, j]) & !is.na(cum_triangle[, j + 1]))
    if (length(valid) > 0) {
      link_ratios[j] <- sum(cum_triangle[valid, j + 1]) /
        sum(cum_triangle[valid, j])
    } else {
      link_ratios[j] <- 1
    }
  }

  F_hat <- numeric(J)
  F_hat[J] <- 1
  for (j in (J - 1):1) {
    F_hat[j] <- F_hat[j + 1] / link_ratios[j]
  }

  pi_hat <- c(F_hat[1], diff(F_hat))

  list(pi_hat = pi_hat, F_hat = F_hat)
}


# =============================================================================
# CONCENTRATION PARAMETER ESTIMATION
# =============================================================================

#' Estimate the Dirichlet Concentration Parameter
#'
#' Estimates the Dirichlet concentration parameter \eqn{c} from an
#' incremental run-off triangle using a partial-column moment method.
#' This estimator works for any triangle shape, including square
#' triangles (\eqn{I = J}) where only one row is fully observed.
#'
#' @param triangle Numeric matrix. Incremental run-off triangle.
#' @param pi_hat Numeric vector. Development proportions. If \code{NULL},
#'   estimated via \code{\link{estimate_development_proportions}}.
#'
#' @return Scalar estimate of \eqn{c}. Returns 50 (diffuse default) if
#'   insufficient data for estimation.
#'
#' @details
#' For each development horizon \eqn{k = 2, \ldots, J-1}, the estimator
#' computes partial proportions \eqn{W_{ij}^{(k)} = X_{ij} / \sum_{l=0}^{k} X_{il}}
#' from the \eqn{I - k} rows observed through period \eqn{k}, matches
#' their variance to the Dirichlet formula
#' \eqn{\mathrm{Var}(W_j^{(k)}) = \pi_j^{(k)}(1 - \pi_j^{(k)}) / (c + 1)},
#' and takes the median across columns and horizons.
#'
#' The diagnostic threshold \eqn{\hat{c} < 30} signals substantial
#' accident-year heterogeneity requiring richer models.
#'
#' @export
estimate_dirichlet_c <- function(triangle, pi_hat = NULL) {
  I <- nrow(triangle)
  J <- ncol(triangle)

  if (is.null(pi_hat)) {
    pi_hat <- estimate_development_proportions(triangle)$pi_hat
  }

  c_estimates <- numeric(0)

  for (k in 2:(J - 1)) {
    rows_with_k <- which(apply(
      triangle[, 1:(k + 1), drop = FALSE], 1,
      function(x) all(!is.na(x))
    ))
    if (length(rows_with_k) < 3) next

    partial_sums <- rowSums(triangle[rows_with_k, 1:(k + 1), drop = FALSE])
    valid <- rows_with_k[partial_sums > 0]
    partial_sums <- partial_sums[partial_sums > 0]
    if (length(valid) < 3) next

    W_partial <- triangle[valid, 1:(k + 1), drop = FALSE] / partial_sums

    F_k <- sum(pi_hat[1:(k + 1)])
    if (F_k <= 0 || F_k > 1.01) next
    pi_partial <- pi_hat[1:(k + 1)] / F_k

    var_W <- apply(W_partial, 2, var)
    c_j <- pi_partial * (1 - pi_partial) / (var_W + 1e-10) - 1
    c_j <- c_j[c_j > 0 & is.finite(c_j)]
    c_estimates <- c(c_estimates, c_j)
  }

  if (length(c_estimates) == 0) return(50)
  max(median(c_estimates), 1)
}


# =============================================================================
# PARAMETRIC BOOTSTRAP
# =============================================================================

#' Multinomial Parametric Bootstrap for Predictive Intervals
#'
#' Generates predictive intervals for total outstanding reserves using
#' the Dirichlet-Gamma parametric bootstrap. The method is
#' model-agnostic: any development proportions can be supplied.
#'
#' @param triangle Numeric matrix. Incremental run-off triangle.
#' @param pi_hat Numeric vector. Development proportions. If \code{NULL},
#'   estimated from the triangle via Chain-Ladder.
#' @param c_param Scalar. Dirichlet concentration parameter. If \code{NULL},
#'   estimated from the triangle via \code{\link{estimate_dirichlet_c}}.
#' @param B Integer. Number of bootstrap replications (default 10000).
#' @param level Numeric. Confidence level for predictive interval
#'   (default 0.95).
#'
#' @return A list with components:
#' \describe{
#'   \item{reserve_cl}{Chain-Ladder point estimate of total reserve.}
#'   \item{reserve_mean}{Bootstrap mean of total reserve.}
#'   \item{reserve_se}{Bootstrap standard error.}
#'   \item{ci_lower}{Lower bound of predictive interval.}
#'   \item{ci_upper}{Upper bound of predictive interval.}
#'   \item{cv}{Coefficient of variation (se / mean).}
#'   \item{c_hat}{Estimated (or supplied) concentration parameter.}
#'   \item{pi_hat}{Development proportions used.}
#'   \item{F_hat}{Cumulative development proportions used.}
#'   \item{samples}{Numeric vector of length \code{B} with bootstrap
#'     reserve samples.}
#'   \item{by_origin}{Data frame with per-accident-year reserve
#'     estimates.}
#' }
#'
#' @details
#' The algorithm:
#' \enumerate{
#'   \item Estimate ultimates \eqn{S_i = X_{i,\mathrm{obs}} / F_{I-i}}.
#'   \item For each bootstrap replication:
#'     \enumerate{
#'       \item Sample \eqn{S_i^* \sim \mathrm{Gamma}(\alpha, \alpha / S_i)}.
#'       \item Sample \eqn{W_i^* \sim \mathrm{Dir}(c \hat{\pi})}.
#'       \item Set \eqn{X_{ij}^* = S_i^* W_{ij}^*}, truncate to
#'         run-off shape.
#'       \item Re-estimate \eqn{\hat{F}^*} from truncated triangle.
#'       \item Compute reserve \eqn{R^* = \sum_i X_{i,\mathrm{obs}}^*
#'         (1 - F^*_{I-i}) / F^*_{I-i}}.
#'     }
#'   \item Return quantiles of \eqn{\{R^{*(b)}\}_{b=1}^B}.
#' }
#'
#' Process risk is captured through the Gamma-Dirichlet sampling;
#' parameter risk through re-estimation. No residual resampling is
#' required.
#'
#' @importFrom MCMCpack rdirichlet
#' @export
multinomial_bootstrap <- function(triangle, pi_hat = NULL, c_param = NULL,
                                   B = 10000, level = 0.95) {

  if (!requireNamespace("MCMCpack", quietly = TRUE)) {
    stop("Package 'MCMCpack' is required. Install with install.packages('MCMCpack').")
  }

  I <- nrow(triangle)
  J <- ncol(triangle)

  # Development proportions
  if (is.null(pi_hat)) {
    dev <- estimate_development_proportions(triangle)
    pi_hat <- dev$pi_hat
    F_hat  <- dev$F_hat
  } else {
    F_hat <- cumsum(pi_hat)
  }

  # Concentration parameter
  if (is.null(c_param)) {
    c_param <- estimate_dirichlet_c(triangle, pi_hat)
  }

  # Ultimate estimates
  X_obs_total <- rowSums(triangle, na.rm = TRUE)
  S_ult <- numeric(I)
  k_obs <- integer(I)
  for (i in 1:I) {
    k_obs[i] <- min(I - i + 1, J)
    S_ult[i] <- X_obs_total[i] / F_hat[k_obs[i]]
  }

  # Gamma shape from CV of ultimates
  cv_ult <- sd(S_ult) / mean(S_ult)
  if (cv_ult < 0.01) cv_ult <- 0.1
  alpha_total <- 1 / cv_ult^2

  # Chain-Ladder point estimate
  reserve_cl <- 0
  for (i in 1:I) {
    if (k_obs[i] < J) {
      reserve_cl <- reserve_cl + X_obs_total[i] * (1 - F_hat[k_obs[i]]) / F_hat[k_obs[i]]
    }
  }

  # Bootstrap
  total_reserve <- numeric(B)
  origin_reserve <- matrix(0, B, I)

  alpha_dir <- c_param * pi_hat

  for (b in 1:B) {
    X_star <- matrix(0, I, J)
    for (i in 1:I) {
      S <- rgamma(1, shape = alpha_total, rate = alpha_total / S_ult[i])
      W <- as.vector(MCMCpack::rdirichlet(1, alpha_dir))
      X_star[i, ] <- S * W
    }

    # Truncate to run-off shape
    for (i in 1:I) {
      k <- I - i + 1
      if (k < J) X_star[i, (k + 1):J] <- NA
    }

    # Re-estimate development
    dev_star <- estimate_development_proportions(X_star)
    F_star <- dev_star$F_hat

    # Compute reserves
    for (i in 1:I) {
      k <- I - i + 1
      if (k < J) {
        X_obs_i <- sum(X_star[i, 1:k], na.rm = TRUE)
        if (X_obs_i > 0 && F_star[k] > 0 && F_star[k] < 1) {
          origin_reserve[b, i] <- X_obs_i * (1 - F_star[k]) / F_star[k]
        }
      }
    }
    total_reserve[b] <- sum(origin_reserve[b, ])
  }

  # Quantiles
  alpha_ci <- (1 - level) / 2
  ci <- quantile(total_reserve, c(alpha_ci, 1 - alpha_ci))

  # Per-origin summary
  by_origin <- data.frame(
    origin     = 1:I,
    observed   = X_obs_total,
    ultimate   = S_ult,
    reserve_cl = ifelse(k_obs < J,
                        X_obs_total * (1 - F_hat[k_obs]) / F_hat[k_obs], 0),
    reserve_mean = colMeans(origin_reserve),
    reserve_se   = apply(origin_reserve, 2, sd)
  )

  list(
    reserve_cl   = reserve_cl,
    reserve_mean = mean(total_reserve),
    reserve_se   = sd(total_reserve),
    ci_lower     = unname(ci[1]),
    ci_upper     = unname(ci[2]),
    cv           = sd(total_reserve) / mean(total_reserve),
    c_hat        = c_param,
    pi_hat       = pi_hat,
    F_hat        = F_hat,
    samples      = total_reserve,
    by_origin    = by_origin
  )
}


# =============================================================================
# DELTA METHOD APPROXIMATION
# =============================================================================

#' Delta Method Variance Approximation for IBNP Reserves
#'
#' Closed-form approximation of the process variance of outstanding
#' reserves under the Dirichlet-Gamma model, avoiding the need for
#' bootstrapping.
#'
#' @param triangle Numeric matrix. Incremental run-off triangle.
#' @param pi_hat Numeric vector. Development proportions. If \code{NULL},
#'   estimated from the triangle.
#' @param c_param Scalar. Concentration parameter. If \code{NULL},
#'   estimated from the triangle.
#'
#' @return A data frame with per-origin-year columns:
#' \describe{
#'   \item{origin}{Accident year index.}
#'   \item{X_obs}{Observed cumulative amount.}
#'   \item{F_k}{Cumulative development proportion at observed horizon.}
#'   \item{reserve}{Chain-Ladder reserve estimate.}
#'   \item{process_var}{Delta method process variance.}
#'   \item{process_se}{Process standard error.}
#'   \item{cv}{Coefficient of variation.}
#' }
#'
#' @details
#' Under the Dirichlet-Gamma model:
#' \deqn{\mathrm{Var}(S_i^{\mathrm{IBNP}} \mid X_i^{\mathrm{obs}})
#'   \approx (X_i^{\mathrm{obs}})^2 \cdot \frac{1 - F}{F^3 (c + 1)}}
#'
#' @export
delta_method_variance <- function(triangle, pi_hat = NULL, c_param = NULL) {

  I <- nrow(triangle)
  J <- ncol(triangle)

  if (is.null(pi_hat)) {
    dev <- estimate_development_proportions(triangle)
    pi_hat <- dev$pi_hat
    F_hat  <- dev$F_hat
  } else {
    F_hat <- cumsum(pi_hat)
  }

  if (is.null(c_param)) {
    c_param <- estimate_dirichlet_c(triangle, pi_hat)
  }

  X_obs <- rowSums(triangle, na.rm = TRUE)

  result <- data.frame(origin = 1:I)
  result$X_obs <- X_obs
  result$F_k <- numeric(I)
  result$reserve <- numeric(I)
  result$process_var <- numeric(I)
  result$process_se <- numeric(I)
  result$cv <- numeric(I)

  for (i in 1:I) {
    k <- min(I - i + 1, J)
    Fk <- F_hat[k]
    result$F_k[i] <- Fk

    if (k < J && Fk > 0 && Fk < 1) {
      R_i <- X_obs[i] * (1 - Fk) / Fk
      V_i <- X_obs[i]^2 * (1 - Fk) / (Fk^3 * (c_param + 1))
      result$reserve[i] <- R_i
      result$process_var[i] <- V_i
      result$process_se[i] <- sqrt(V_i)
      result$cv[i] <- sqrt(V_i) / R_i
    }
  }

  result
}


# =============================================================================
# BAYESIAN PREDICTIVE BOOTSTRAP
# =============================================================================

#' Dirichlet Log-Likelihood
#'
#' Evaluates the Dirichlet log-likelihood for a single observation.
#'
#' @param w Numeric vector. Observed proportions (sums to 1).
#' @param c_val Scalar. Concentration parameter.
#' @param pi_vec Numeric vector. Base proportions.
#'
#' @return Scalar log-likelihood value.
#'
#' @keywords internal
dirichlet_loglik <- function(w, c_val, pi_vec) {
  a <- c_val * pi_vec
  sum((a - 1) * log(w + 1e-300)) + lgamma(sum(a)) - sum(lgamma(a))
}


#' Bayesian Predictive Bootstrap for Claims Reserves
#'
#' Full Bayesian extension of the multinomial bootstrap that places
#' priors on the concentration parameter \eqn{c} and development
#' proportions \eqn{\pi}, using MCMC for posterior computation.
#'
#' @param triangle Numeric matrix. Incremental run-off triangle.
#' @param B Integer. Number of predictive samples (default 5000).
#' @param mu_c Scalar. Prior mean of \eqn{\log c} (default \code{log(50)}).
#' @param sigma_c Scalar. Prior SD of \eqn{\log c} (default 1).
#' @param a0 Scalar. Dirichlet prior concentration for \eqn{\pi}
#'   (default 1, i.e., uniform).
#' @param n_mcmc Integer. Number of MCMC iterations after burn-in
#'   (default 2000).
#' @param burnin Integer. Number of burn-in iterations (default 500).
#' @param bayesian_c_only Logical. If \code{TRUE}, only \eqn{c} is
#'   sampled; \eqn{\pi} is fixed at the plug-in estimate.
#' @param level Numeric. Confidence level (default 0.95).
#'
#' @return A list with components:
#' \describe{
#'   \item{reserve_mean}{Predictive mean of total reserve.}
#'   \item{reserve_se}{Predictive standard error.}
#'   \item{ci_lower}{Lower bound of predictive interval.}
#'   \item{ci_upper}{Upper bound of predictive interval.}
#'   \item{c_posterior_mean}{Posterior mean of \eqn{c}.}
#'   \item{c_posterior_sd}{Posterior SD of \eqn{c}.}
#'   \item{accept_rate}{Metropolis acceptance rate for \eqn{c}.}
#'   \item{samples}{Numeric vector of predictive reserve samples.}
#' }
#'
#' @details
#' The prior on \eqn{c} is log-normal with default median 50 and 95\%
#' interval approximately (7, 370). The MCMC uses Metropolis-within-Gibbs:
#' \eqn{\pi \mid c} is conjugate Dirichlet; \eqn{c \mid \pi} uses a
#' random-walk Metropolis step on \eqn{\log c}. Partial rows (with
#' \eqn{\geq 3} observed columns) contribute to the likelihood with
#' \eqn{\pi} rescaled to the observed subset.
#'
#' @importFrom MCMCpack rdirichlet
#' @export
bayesian_predictive_bootstrap <- function(triangle, B = 5000,
                                           mu_c = log(50), sigma_c = 1,
                                           a0 = 1,
                                           n_mcmc = 2000, burnin = 500,
                                           bayesian_c_only = FALSE,
                                           level = 0.95) {

  if (!requireNamespace("MCMCpack", quietly = TRUE)) {
    stop("Package 'MCMCpack' is required. Install with install.packages('MCMCpack').")
  }

  I <- nrow(triangle)
  J <- ncol(triangle)

  dev <- estimate_development_proportions(triangle)
  pi_hat <- dev$pi_hat
  F_hat  <- dev$F_hat

  # Identify usable rows (at least 3 observed columns)
  min_cols <- min(3, J)
  usable_rows <- which(apply(triangle, 1, function(x) sum(!is.na(x)) >= min_cols))

  if (length(usable_rows) < 2) {
    message("Insufficient rows for MCMC; falling back to plug-in bootstrap.")
    return(multinomial_bootstrap(triangle, B = B, level = level))
  }

  # Prepare partial-row data for likelihood
  W_obs_list <- list()
  obs_cols_list <- list()
  for (idx in seq_along(usable_rows)) {
    i <- usable_rows[idx]
    obs_cols <- which(!is.na(triangle[i, ]))
    row_total <- sum(triangle[i, obs_cols])
    if (row_total <= 0) next
    F_obs <- sum(pi_hat[obs_cols])
    if (F_obs <= 0) next
    W_obs_list[[length(W_obs_list) + 1]] <- triangle[i, obs_cols] / row_total
    obs_cols_list[[length(obs_cols_list) + 1]] <- obs_cols
  }

  if (length(W_obs_list) < 2) {
    message("Insufficient valid rows for MCMC; falling back to plug-in bootstrap.")
    return(multinomial_bootstrap(triangle, B = B, level = level))
  }

  # Ultimates for bootstrap generation
  X_obs_total <- rowSums(triangle, na.rm = TRUE)
  S_ult <- numeric(I)
  for (i in 1:I) {
    k <- min(I - i + 1, J)
    S_ult[i] <- X_obs_total[i] / F_hat[k]
  }
  cv_ult <- sd(S_ult) / mean(S_ult)
  if (cv_ult < 0.01) cv_ult <- 0.1
  alpha_total <- 1 / cv_ult^2

  # --- MCMC ---
  c_cur <- estimate_dirichlet_c(triangle, pi_hat)
  if (c_cur < 1) c_cur <- 10
  pi_cur <- pi_hat

  c_samples  <- numeric(n_mcmc)
  pi_samples <- matrix(0, n_mcmc, J)
  prop_sd    <- 0.5
  n_accept   <- 0L

  for (iter in 1:(n_mcmc + burnin)) {

    # Update pi | c (approximate conjugate using partial rows)
    if (!bayesian_c_only) {
      col_contrib <- numeric(J)
      for (idx in seq_along(W_obs_list)) {
        obs_cols <- obs_cols_list[[idx]]
        w <- W_obs_list[[idx]]
        col_contrib[obs_cols] <- col_contrib[obs_cols] + c_cur * w
      }
      alpha_post <- a0 + col_contrib
      alpha_post <- pmax(alpha_post, 0.01)
      pi_cur <- as.vector(MCMCpack::rdirichlet(1, alpha_post))
    }

    # Update c | pi (Metropolis on log c)
    log_c_prop <- rnorm(1, log(c_cur), prop_sd)
    c_prop <- exp(log_c_prop)

    ll_cur <- 0
    ll_prop <- 0
    for (idx in seq_along(W_obs_list)) {
      w <- W_obs_list[[idx]]
      obs_cols <- obs_cols_list[[idx]]
      F_obs <- sum(pi_cur[obs_cols])
      if (F_obs <= 0) next
      pi_partial <- pi_cur[obs_cols] / F_obs
      ll_cur  <- ll_cur  + dirichlet_loglik(w, c_cur, pi_partial)
      ll_prop <- ll_prop + dirichlet_loglik(w, c_prop, pi_partial)
    }

    lp_cur  <- dnorm(log(c_cur),  mu_c, sigma_c, log = TRUE)
    lp_prop <- dnorm(log(c_prop), mu_c, sigma_c, log = TRUE)
    log_alpha <- (ll_prop + lp_prop) - (ll_cur + lp_cur)

    if (log(runif(1)) < log_alpha) {
      c_cur <- c_prop
      if (iter > burnin) n_accept <- n_accept + 1L
    }

    if (iter > burnin) {
      idx_out <- iter - burnin
      c_samples[idx_out]    <- c_cur
      pi_samples[idx_out, ] <- pi_cur
    }
  }

  # --- Predictive bootstrap from posterior samples ---
  thin_idx <- round(seq(1, n_mcmc, length.out = B))
  total_reserve <- numeric(B)

  for (b in 1:B) {
    c_b  <- c_samples[thin_idx[b]]
    pi_b <- if (bayesian_c_only) pi_hat else pi_samples[thin_idx[b], ]

    X_star <- matrix(0, I, J)
    for (i in 1:I) {
      S <- rgamma(1, shape = alpha_total, rate = alpha_total / S_ult[i])
      W <- as.vector(MCMCpack::rdirichlet(1, c_b * pi_b))
      X_star[i, ] <- S * W
    }
    for (i in 1:I) {
      k <- I - i + 1
      if (k < J) X_star[i, (k + 1):J] <- NA
    }
    dev_star <- estimate_development_proportions(X_star)
    F_star <- dev_star$F_hat
    reserve <- 0
    for (i in 1:I) {
      k <- I - i + 1
      if (k < J) {
        X_obs_i <- sum(X_star[i, 1:k], na.rm = TRUE)
        if (X_obs_i > 0 && F_star[k] > 0 && F_star[k] < 1) {
          reserve <- reserve + X_obs_i * (1 - F_star[k]) / F_star[k]
        }
      }
    }
    total_reserve[b] <- reserve
  }

  alpha_ci <- (1 - level) / 2
  ci <- quantile(total_reserve, c(alpha_ci, 1 - alpha_ci))

  list(
    reserve_mean     = mean(total_reserve),
    reserve_se       = sd(total_reserve),
    ci_lower         = unname(ci[1]),
    ci_upper         = unname(ci[2]),
    c_posterior_mean = mean(c_samples),
    c_posterior_sd   = sd(c_samples),
    accept_rate      = n_accept / n_mcmc,
    samples          = total_reserve
  )
}


# =============================================================================
# DIAGNOSTIC SUMMARY
# =============================================================================

#' Concentration Parameter Diagnostic
#'
#' Estimates the Dirichlet concentration parameter and provides a
#' diagnostic assessment of whether the homogeneous development
#' assumption is appropriate.
#'
#' @param triangle Numeric matrix. Incremental run-off triangle.
#' @param threshold Scalar. Diagnostic threshold (default 30).
#'
#' @return A list with components:
#' \describe{
#'   \item{c_hat}{Estimated concentration parameter.}
#'   \item{threshold}{Diagnostic threshold used.}
#'   \item{adequate}{Logical. \code{TRUE} if \eqn{\hat{c} \geq} threshold.}
#'   \item{message}{Interpretive message.}
#' }
#'
#' @export
diagnose_concentration <- function(triangle, threshold = 30) {
  dev <- estimate_development_proportions(triangle)
  c_hat <- estimate_dirichlet_c(triangle, dev$pi_hat)
  adequate <- c_hat >= threshold

  msg <- if (adequate) {
    sprintf("c_hat = %.1f >= %d: development patterns are stable. Multinomial framework is appropriate.",
            c_hat, threshold)
  } else {
    sprintf("c_hat = %.1f < %d: substantial accident-year heterogeneity detected. Consider models with explicit frailty.",
            c_hat, threshold)
  }

  list(
    c_hat     = c_hat,
    threshold = threshold,
    adequate  = adequate,
    message   = msg
  )
}
