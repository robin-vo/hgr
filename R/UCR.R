#' @title Classical Reserving Methods as Credibility Estimators
#' @description 
#' Implementation of the unified credibility framework for claims reserving.
#' This file contains all functions from "Classical Reserving Methods as 
#' Credibility Estimators: A Unified Bayesian Framework" (Van Oirbeek, 2026).
#' 
#' The key insight is that Chain-Ladder, Cape Cod, Bornhuetter-Ferguson, and
#' Mack are all credibility estimators under a Poisson-Gamma-Multinomial model.
#' 
#' @author Robin Van Oirbeek
#' @import stats

# =============================================================================
# PART 1: UTILITY FUNCTIONS
# =============================================================================

#' Estimate development pattern from triangle
#'
#' @param triangle A matrix representing the run-off triangle (cumulative claims)
#' @return A list with age-to-age factors and cumulative proportions
#' @export
estimate_development <- function(triangle) {
  I <- nrow(triangle)
  J <- ncol(triangle)
  
  # Age-to-age factors
  dev_factors <- numeric(J - 1)
  for (j in 1:(J - 1)) {
    num <- 0
    den <- 0
    for (i in 1:(I - j)) {
      if (!is.na(triangle[i, j]) && !is.na(triangle[i, j + 1])) {
        num <- num + triangle[i, j + 1]
        den <- den + triangle[i, j]
      }
    }
    dev_factors[j] <- if (den > 0) num / den else NA
  }
  
  # Cumulative development factors (to ultimate)
  cum_factors <- rev(cumprod(rev(dev_factors[!is.na(dev_factors)])))
  cum_factors <- c(cum_factors, 1)
  
  # Pad if needed
  if (length(cum_factors) < J) {
    cum_factors <- c(rep(NA, J - length(cum_factors)), cum_factors)
  }
  
  # Proportions emerged
  F_cum <- 1 / cum_factors
  
  list(
    dev_factors = dev_factors,
    cum_factors = cum_factors,
    F_cum = F_cum,
    I = I,
    J = J
  )
}

# =============================================================================
# PART 2: CLASSICAL RESERVING METHODS
# =============================================================================

#' Chain-Ladder reserving method
#'
#' The Chain-Ladder method corresponds to the credibility limit Z -> 1
#' (full weight on individual experience, diffuse prior).
#'
#' @param triangle A matrix representing the run-off triangle (cumulative claims)
#' @param F_cum Optional vector of cumulative development proportions
#' @return An object of class "reserve_estimate" with reserve estimates
#' @export
chain_ladder <- function(triangle, F_cum = NULL) {
  I <- nrow(triangle)
  J <- ncol(triangle)
  
  if (is.null(F_cum)) {
    dev <- estimate_development(triangle)
    F_cum <- dev$F_cum
  }
  
  # Get observed cumulative for each accident year
  C_obs <- numeric(I)
  j_obs <- numeric(I)
  for (i in 1:I) {
    j_obs[i] <- max(which(!is.na(triangle[i, ])))
    C_obs[i] <- triangle[i, j_obs[i]]
  }
  
  F_obs <- F_cum[j_obs]
  ultimate <- C_obs / F_obs
  reserve <- ultimate - C_obs
  
  result <- list(
    ultimate = ultimate,
    reserve = reserve,
    total_reserve = sum(reserve),
    C_obs = C_obs,
    j_obs = j_obs,
    F_obs = F_obs,
    F_cum = F_cum,
    triangle = triangle,
    method = "Chain-Ladder"
  )
  
  class(result) <- "reserve_estimate"
  return(result)
}

#' Cape Cod reserving method
#'
#' The Cape Cod method corresponds to the credibility limit Z -> 0
#' (full weight on pooled experience, homogeneous rates assumption).
#' This is the Bühlmann-Straub credibility estimator.
#'
#' @param triangle A matrix representing the run-off triangle (cumulative claims)
#' @param exposure Vector of exposures (e.g., premium, first period claims)
#' @param F_cum Optional vector of cumulative development proportions
#' @return An object of class "reserve_estimate" with reserve estimates
#' @export
cape_cod <- function(triangle, exposure = NULL, F_cum = NULL) {
  I <- nrow(triangle)
  J <- ncol(triangle)
  
  if (is.null(F_cum)) {
    dev <- estimate_development(triangle)
    F_cum <- dev$F_cum
  }
  
  # Default exposure: first development period claims
  if (is.null(exposure)) {
    exposure <- triangle[, 1]
  }
  
  # Get observed cumulative for each accident year
  C_obs <- numeric(I)
  j_obs <- numeric(I)
  for (i in 1:I) {
    j_obs[i] <- max(which(!is.na(triangle[i, ])))
    C_obs[i] <- triangle[i, j_obs[i]]
  }
  
  F_obs <- F_cum[j_obs]
  
  # Pooled claim rate
  q_pool <- sum(C_obs) / sum(exposure * F_obs)
  
  # Cape Cod ultimate and reserve
  ultimate <- exposure * q_pool
  reserve <- ultimate - C_obs
  
  result <- list(
    ultimate = ultimate,
    reserve = reserve,
    total_reserve = sum(reserve),
    q_pool = q_pool,
    C_obs = C_obs,
    j_obs = j_obs,
    F_obs = F_obs,
    F_cum = F_cum,
    exposure = exposure,
    triangle = triangle,
    method = "Cape Cod"
  )
  
  class(result) <- "reserve_estimate"
  return(result)
}

#' Bornhuetter-Ferguson reserving method
#'
#' The Bornhuetter-Ferguson method uses separation credibility:
#' observed claims get full credibility, unreported claims use prior.
#' Ultimate = Observed + (1-F) * Prior_Ultimate
#'
#' @param triangle A matrix representing the run-off triangle (cumulative claims)
#' @param exposure Vector of exposures
#' @param prior_elr Prior expected loss ratio
#' @param ultimate_prior Optional vector of prior ultimates (overrides prior_elr)
#' @param F_cum Optional vector of cumulative development proportions
#' @return An object of class "reserve_estimate" with reserve estimates
#' @export
bornhuetter_ferguson <- function(triangle, exposure = NULL, prior_elr = NULL, 
                                  ultimate_prior = NULL, F_cum = NULL) {
  I <- nrow(triangle)
  J <- ncol(triangle)
  
  if (is.null(F_cum)) {
    dev <- estimate_development(triangle)
    F_cum <- dev$F_cum
  }
  
  # Default exposure: first development period claims
  if (is.null(exposure)) {
    exposure <- triangle[, 1]
  }
  
  # Get observed cumulative for each accident year
  C_obs <- numeric(I)
  j_obs <- numeric(I)
  for (i in 1:I) {
    j_obs[i] <- max(which(!is.na(triangle[i, ])))
    C_obs[i] <- triangle[i, j_obs[i]]
  }
  
  F_obs <- F_cum[j_obs]
  
  # Prior ultimate
  if (is.null(ultimate_prior)) {
    if (is.null(prior_elr)) {
      # Default: use Cape Cod pooled rate as prior
      prior_elr <- sum(C_obs) / sum(exposure * F_obs)
    }
    ultimate_prior <- exposure * prior_elr
  }
  
  # BF formula: Reserve = (1 - F) * Prior Ultimate
  reserve <- (1 - F_obs) * ultimate_prior
  ultimate <- C_obs + reserve
  
  result <- list(
    ultimate = ultimate,
    reserve = reserve,
    total_reserve = sum(reserve),
    ultimate_prior = ultimate_prior,
    C_obs = C_obs,
    j_obs = j_obs,
    F_obs = F_obs,
    F_cum = F_cum,
    exposure = exposure,
    triangle = triangle,
    method = "Bornhuetter-Ferguson"
  )
  
  class(result) <- "reserve_estimate"
  return(result)
}

#' Fit Mack's Chain-Ladder model
#'
#' Mack's distribution-free model. The variance formula equals the
#' Bayesian posterior predictive variance under diffuse priors.
#'
#' @param triangle A matrix representing the run-off triangle (cumulative claims)
#' @param tail_factor Tail factor for development beyond the triangle (default 1.0)
#' @param alpha Vector of variance weights (default: all 1)
#' @return An object of class "mack" with estimates and standard errors
#' @export
fit_mack <- function(triangle, tail_factor = 1.0, alpha = NULL) {
  I <- nrow(triangle)
  J <- ncol(triangle)
  
  if (is.null(alpha)) {
    alpha <- rep(1, J - 1)
  }
  
  # Age-to-age factors
  f <- numeric(J - 1)
  for (j in 1:(J - 1)) {
    num <- 0
    den <- 0
    for (i in 1:(I - j)) {
      if (!is.na(triangle[i, j]) && !is.na(triangle[i, j + 1])) {
        num <- num + triangle[i, j + 1]
        den <- den + triangle[i, j]
      }
    }
    f[j] <- if (den > 0) num / den else NA
  }
  
  # Cumulative development factors
  F_cum <- c(rev(cumprod(rev(f))), 1) * tail_factor
  
  # Estimate sigma_j^2 (Mack's formula)
  sigma_sq <- numeric(J - 1)
  for (j in 1:(J - 1)) {
    n_obs <- I - j
    if (n_obs > 1) {
      ss <- 0
      for (i in 1:n_obs) {
        if (!is.na(triangle[i, j]) && !is.na(triangle[i, j + 1])) {
          residual <- triangle[i, j + 1] / triangle[i, j] - f[j]
          ss <- ss + triangle[i, j]^alpha[j] * residual^2
        }
      }
      sigma_sq[j] <- ss / (n_obs - 1)
    } else {
      if (j >= 2 && !is.na(sigma_sq[j - 1])) {
        sigma_sq[j] <- min(sigma_sq[j - 1]^2 / sigma_sq[max(1, j - 2)], 
                          sigma_sq[j - 1])
      } else {
        sigma_sq[j] <- sigma_sq[max(1, j - 1)]
      }
    }
  }
  
  # Get observed cumulative for each accident year
  C_obs <- numeric(I)
  j_obs <- numeric(I)
  for (i in 1:I) {
    j_obs[i] <- max(which(!is.na(triangle[i, ])))
    C_obs[i] <- triangle[i, j_obs[i]]
  }
  
  # Ultimate estimates
  ultimate <- C_obs * F_cum[j_obs]
  reserve <- ultimate - C_obs
  
  # Mean squared error of prediction (MSEP)
  msep <- numeric(I)
  
  for (i in 2:I) {
    se_sq <- 0
    C_hat <- C_obs[i]
    
    for (j in j_obs[i]:(J - 1)) {
      if (!is.na(f[j]) && !is.na(sigma_sq[j])) {
        S_j <- sum(triangle[1:(I - j), j], na.rm = TRUE)
        
        if (S_j > 0 && C_hat > 0) {
          se_sq <- se_sq + (sigma_sq[j] / f[j]^2) * (1/C_hat + 1/S_j)
        }
        
        C_hat <- C_hat * f[j]
      }
    }
    
    msep[i] <- ultimate[i]^2 * se_sq
  }
  
  se <- sqrt(msep)
  total_se <- sqrt(sum(msep))
  
  result <- list(
    ultimate = ultimate,
    reserve = reserve,
    total_reserve = sum(reserve),
    se = se,
    total_se = total_se,
    msep = msep,
    f = f,
    F_cum = F_cum,
    sigma_sq = sigma_sq,
    alpha = alpha,
    C_obs = C_obs,
    j_obs = j_obs,
    tail_factor = tail_factor,
    triangle = triangle,
    I = I,
    J = J
  )
  
  class(result) <- "mack"
  return(result)
}

# =============================================================================
# PART 3: CREDIBILITY THEORY
# =============================================================================

#' Bühlmann credibility estimator
#'
#' @param x Vector of observations for the risk
#' @param mu Collective mean (prior mean)
#' @param sigma_sq Within-risk variance (process variance)
#' @param tau_sq Between-risk variance (structural variance)
#' @return A list with credibility estimate, weight, and components
#' @export
buhlmann <- function(x, mu, sigma_sq, tau_sq) {
  n <- length(x)
  x_bar <- mean(x)
  
  k <- sigma_sq / tau_sq
  Z <- n / (n + k)
  estimate <- Z * x_bar + (1 - Z) * mu
  
  result <- list(
    estimate = estimate,
    Z = Z,
    k = k,
    x_bar = x_bar,
    mu = mu,
    n = n,
    sigma_sq = sigma_sq,
    tau_sq = tau_sq
  )
  
  class(result) <- "buhlmann"
  return(result)
}

#' Bühlmann-Straub credibility estimator
#'
#' @param x Vector of observations (claim rates)
#' @param w Vector of exposure weights
#' @param mu Collective mean
#' @param sigma_sq Within-risk variance per unit weight
#' @param tau_sq Between-risk variance
#' @return A list with credibility estimates and weights
#' @export
buhlmann_straub <- function(x, w, mu, sigma_sq, tau_sq) {
  I <- length(x)
  
  if (tau_sq > 0) {
    k <- sigma_sq / tau_sq
    Z <- w / (w + k)
  } else {
    k <- Inf
    Z <- rep(0, I)
  }
  
  estimate <- Z * x + (1 - Z) * mu
  
  result <- list(
    estimate = estimate,
    Z = Z,
    k = k,
    x = x,
    w = w,
    mu = mu,
    I = I,
    sigma_sq = sigma_sq,
    tau_sq = tau_sq
  )
  
  class(result) <- "buhlmann_straub"
  return(result)
}

#' Estimate Bühlmann-Straub variance components
#'
#' @param x Vector of observations (claim rates)
#' @param w Vector of exposure weights
#' @param sigma_sq Optional fixed process variance
#' @return A list with estimated parameters
#' @export
estimate_bs_params <- function(x, w, sigma_sq = NULL) {
  I <- length(x)
  
  mu_hat <- sum(w * x) / sum(w)
  
  if (is.null(sigma_sq)) {
    sigma_sq_hat <- mu_hat  # Poisson assumption
  } else {
    sigma_sq_hat <- sigma_sq
  }
  
  SS_between <- sum(w * (x - mu_hat)^2)
  w_total <- sum(w)
  w_factor <- w_total - sum(w^2) / w_total
  
  tau_sq_hat <- (SS_between - (I - 1) * sigma_sq_hat) / w_factor
  tau_sq_hat <- max(tau_sq_hat, 0)
  
  list(
    mu = mu_hat,
    sigma_sq = sigma_sq_hat,
    tau_sq = tau_sq_hat,
    k = if (tau_sq_hat > 0) sigma_sq_hat / tau_sq_hat else Inf,
    I = I,
    w = w
  )
}

#' Three-source credibility estimator
#'
#' @param x Vector of individual observations/rates
#' @param w Vector of exposure weights
#' @param mu_pool Pooled (collective) mean
#' @param mu_prior External prior mean
#' @param sigma_sq Process variance per unit weight
#' @param tau_sq Between-risk variance
#' @param tau_prior_sq Prior variance
#' @return A list with credibility estimates and weights
#' @export
three_source_credibility <- function(x, w, mu_pool, mu_prior, 
                                      sigma_sq, tau_sq, tau_prior_sq) {
  I <- length(x)
  
  prec_ind <- w / sigma_sq
  prec_pool <- if (tau_sq > 0) 1 / tau_sq else 0
  prec_prior <- if (tau_prior_sq > 0) 1 / tau_prior_sq else 0
  prec_total <- prec_ind + prec_pool + prec_prior
  
  Z1 <- prec_ind / prec_total
  Z2 <- prec_pool / prec_total
  Z0 <- prec_prior / prec_total
  
  estimate <- Z1 * x + Z2 * mu_pool + Z0 * mu_prior
  
  list(
    estimate = estimate,
    Z1 = Z1,
    Z2 = Z2,
    Z0 = Z0,
    x = x,
    mu_pool = mu_pool,
    mu_prior = mu_prior
  )
}

# =============================================================================
# PART 4: UNIFIED CREDIBILITY RESERVING (UCR)
# =============================================================================

#' Unified Credibility Reserving (UCR)
#'
#' UCR is a general credibility-based reserving method that nests Chain-Ladder,
#' Cape Cod, and Bornhuetter-Ferguson as special cases. It adaptively estimates
#' credibility weights from the data using the Bühlmann-Straub framework.
#'
#' @param triangle A matrix representing the run-off triangle (cumulative claims)
#' @param exposure Vector of exposures (default: first period claims)
#' @param F_cum Optional vector of cumulative development proportions
#' @param prior_mean Optional external prior mean (mu_0)
#' @param prior_var Optional external prior variance (tau_0^2)
#' @param sigma_sq Optional process variance (default: estimated as pooled rate)
#' @return An object of class "ucr" with reserve estimates and diagnostics
#' @export
ucr <- function(triangle, exposure = NULL, F_cum = NULL, 
                prior_mean = NULL, prior_var = NULL, sigma_sq = NULL) {
  I <- nrow(triangle)
  J <- ncol(triangle)
  
  if (is.null(F_cum)) {
    dev <- estimate_development(triangle)
    F_cum <- dev$F_cum
  }
  
  if (is.null(exposure)) {
    exposure <- triangle[, 1]
  }
  
  # Get observed cumulative for each accident year
  C_obs <- numeric(I)
  j_obs <- numeric(I)
  for (i in 1:I) {
    j_obs[i] <- max(which(!is.na(triangle[i, ])))
    C_obs[i] <- triangle[i, j_obs[i]]
  }
  
  F_obs <- F_cum[j_obs]
  w <- exposure * F_obs
  
  # Individual Chain-Ladder rates
  Lambda_CL <- C_obs / w
  
  # Pooled rate (Cape Cod)
  mu_pool <- sum(w * Lambda_CL) / sum(w)
  
  # Process variance
  if (is.null(sigma_sq)) {
    sigma_sq <- mu_pool
  }
  
  # Estimate between-year variance (Bühlmann-Straub)
  SS_between <- sum(w * (Lambda_CL - mu_pool)^2)
  w_total <- sum(w)
  w_factor <- w_total - sum(w^2) / w_total
  tau_sq_hat <- max((SS_between - (I - 1) * sigma_sq) / w_factor, 0)
  
  # Compute credibility weights
  if (tau_sq_hat > 0) {
    k <- sigma_sq / tau_sq_hat
    
    if (!is.null(prior_mean) && !is.null(prior_var)) {
      # Three-source credibility
      prec_ind <- w / sigma_sq
      prec_pool <- 1 / tau_sq_hat
      prec_prior <- 1 / prior_var
      prec_total <- prec_ind + prec_pool + prec_prior
      
      Z1 <- prec_ind / prec_total
      Z2 <- prec_pool / prec_total
      Z0 <- prec_prior / prec_total
      
      Lambda_UCR <- Z1 * Lambda_CL + Z2 * mu_pool + Z0 * prior_mean
    } else {
      # Two-source credibility
      Z1 <- w / (w + k)
      Z2 <- 1 - Z1
      Z0 <- rep(0, I)
      
      Lambda_UCR <- Z1 * Lambda_CL + Z2 * mu_pool
    }
  } else {
    # tau^2 = 0: full pooling (Cape Cod)
    Z1 <- rep(0, I)
    Z2 <- rep(1, I)
    Z0 <- rep(0, I)
    k <- Inf
    
    Lambda_UCR <- rep(mu_pool, I)
  }
  
  ultimate <- Lambda_UCR * exposure
  reserve <- ultimate - C_obs
  
  result <- list(
    ultimate = ultimate,
    reserve = reserve,
    total_reserve = sum(reserve),
    Lambda_CL = Lambda_CL,
    Lambda_UCR = Lambda_UCR,
    mu_pool = mu_pool,
    sigma_sq = sigma_sq,
    tau_sq_hat = tau_sq_hat,
    k = k,
    Z1 = Z1,
    Z2 = Z2,
    Z0 = Z0,
    C_obs = C_obs,
    j_obs = j_obs,
    F_obs = F_obs,
    F_cum = F_cum,
    w = w,
    exposure = exposure,
    triangle = triangle,
    prior_mean = prior_mean,
    prior_var = prior_var,
    method = "UCR"
  )
  
  class(result) <- c("ucr", "reserve_estimate")
  return(result)
}

# =============================================================================
# PART 5: COMPARISON AND DIAGNOSTICS
# =============================================================================

#' Compare all reserving methods
#'
#' @param triangle A matrix representing the run-off triangle (cumulative claims)
#' @param exposure Vector of exposures (optional)
#' @param prior_elr Prior expected loss ratio for BF (optional)
#' @return A data frame comparing all methods
#' @export
compare_reserves <- function(triangle, exposure = NULL, prior_elr = NULL) {
  cl <- chain_ladder(triangle)
  cc <- cape_cod(triangle, exposure)
  ucr_fit <- ucr(triangle, exposure)
  
  if (!is.null(prior_elr) || !is.null(exposure)) {
    bf <- bornhuetter_ferguson(triangle, exposure, prior_elr)
  } else {
    bf <- bornhuetter_ferguson(triangle, exposure, prior_elr = cc$q_pool)
  }
  
  I <- nrow(triangle)
  
  comparison <- data.frame(
    AY = 1:I,
    Observed = cl$C_obs,
    CL = round(cl$reserve, 0),
    CC = round(cc$reserve, 0),
    BF = round(bf$reserve, 0),
    UCR = round(ucr_fit$reserve, 0),
    Z1 = round(ucr_fit$Z1, 3)
  )
  
  totals <- data.frame(
    AY = "Total",
    Observed = sum(cl$C_obs),
    CL = round(cl$total_reserve, 0),
    CC = round(cc$total_reserve, 0),
    BF = round(bf$total_reserve, 0),
    UCR = round(ucr_fit$total_reserve, 0),
    Z1 = round(mean(ucr_fit$Z1), 3)
  )
  
  result <- list(
    by_ay = comparison,
    totals = totals,
    cl = cl,
    cc = cc,
    bf = bf,
    ucr = ucr_fit
  )
  
  class(result) <- "reserve_comparison"
  return(result)
}

# =============================================================================
# PART 6: PRINT AND PLOT METHODS
# =============================================================================

#' @export
print.reserve_estimate <- function(x, ...) {
  cat(x$method, "Reserve Estimates\n")
  cat(paste(rep("=", nchar(x$method) + 18), collapse = ""), "\n\n")
  
  I <- length(x$ultimate)
  df <- data.frame(
    AY = 1:I,
    Observed = x$C_obs,
    Ultimate = round(x$ultimate, 0),
    Reserve = round(x$reserve, 0)
  )
  print(df, row.names = FALSE)
  cat("\nTotal Reserve:", format(round(x$total_reserve, 0), big.mark = ","), "\n")
}

#' @export
print.ucr <- function(x, ...) {
  cat("Unified Credibility Reserving (UCR)\n")
  cat("===================================\n\n")
  
  cat("Variance Estimation:\n")
  cat("  sigma^2 (process):", round(x$sigma_sq, 4), "\n")
  cat("  tau^2 (between-year):", round(x$tau_sq_hat, 6), "\n")
  if (is.finite(x$k)) {
    cat("  k = sigma^2/tau^2:", round(x$k, 4), "\n")
  } else {
    cat("  k = sigma^2/tau^2: Inf (homogeneous rates)\n")
  }
  cat("\n")
  
  I <- length(x$ultimate)
  df <- data.frame(
    AY = 1:I,
    Observed = format(round(x$C_obs, 0), big.mark = ","),
    Lambda_CL = round(x$Lambda_CL, 4),
    Z1 = round(x$Z1, 3),
    Lambda_UCR = round(x$Lambda_UCR, 4),
    Reserve = format(round(x$reserve, 0), big.mark = ",")
  )
  print(df, row.names = FALSE)
  
  cat("\nTotal Reserve:", format(round(x$total_reserve, 0), big.mark = ","), "\n")
  
  mean_Z1 <- mean(x$Z1)
  if (mean_Z1 > 0.9) {
    cat("\nInterpretation: High Z1 -> UCR ≈ Chain-Ladder\n")
  } else if (mean_Z1 < 0.1) {
    cat("\nInterpretation: Low Z1 -> UCR ≈ Cape Cod\n")
  } else {
    cat("\nInterpretation: Moderate Z1 -> UCR interpolates between CL and CC\n")
  }
}

#' @export
print.mack <- function(x, ...) {
  cat("Mack's Chain-Ladder Model\n")
  cat("=========================\n\n")
  
  df <- data.frame(
    AY = 1:x$I,
    Observed = format(round(x$C_obs, 0), big.mark = ","),
    Ultimate = format(round(x$ultimate, 0), big.mark = ","),
    Reserve = format(round(x$reserve, 0), big.mark = ","),
    SE = format(round(x$se, 0), big.mark = ","),
    CV = paste0(round(100 * x$se / pmax(x$reserve, 1), 1), "%")
  )
  print(df, row.names = FALSE)
  
  cat("\nTotal Reserve:", format(round(x$total_reserve, 0), big.mark = ","), "\n")
  cat("Total S.E.:   ", format(round(x$total_se, 0), big.mark = ","), "\n")
}

#' @export
print.reserve_comparison <- function(x, ...) {
  cat("Reserve Comparison\n")
  cat("==================\n\n")
  print(x$by_ay, row.names = FALSE)
  cat("\nTotals:\n")
  cat("  Chain-Ladder:", format(x$totals$CL, big.mark = ","), "\n")
  cat("  Cape Cod:    ", format(x$totals$CC, big.mark = ","), "\n")
  cat("  BF:          ", format(x$totals$BF, big.mark = ","), "\n")
  cat("  UCR:         ", format(x$totals$UCR, big.mark = ","), "\n")
  cat("\nUCR: tau^2 =", round(x$ucr$tau_sq_hat, 6), ", mean Z1 =", round(mean(x$ucr$Z1), 3), "\n")
}

#' @export
plot.ucr <- function(x, which = 1:3, ...) {
  n_plots <- length(which)
  if (n_plots > 1) {
    old_par <- par(mfrow = c(1, min(n_plots, 3)), mar = c(4.5, 4.5, 3, 1))
    on.exit(par(old_par))
  }
  
  I <- length(x$ultimate)
  
  if (1 %in% which) {
    plot(1:I, x$Lambda_CL, type = "b", pch = 19, col = "blue",
         xlab = "Accident Year", ylab = "Claim Rate",
         main = "Individual vs Pooled Rates", ylim = range(c(x$Lambda_CL, x$mu_pool)))
    abline(h = x$mu_pool, lty = 2, col = "red", lwd = 2)
    points(1:I, x$Lambda_UCR, pch = 17, col = "darkgreen")
    legend("topright", legend = c("Chain-Ladder", "Cape Cod", "UCR"),
           pch = c(19, NA, 17), lty = c(1, 2, NA), col = c("blue", "red", "darkgreen"), bty = "n")
  }
  
  if (2 %in% which) {
    barplot(x$Z1, names.arg = 1:I, col = "steelblue",
            xlab = "Accident Year", ylab = expression(Z^{(1)}),
            main = "Credibility Weights", ylim = c(0, 1))
    abline(h = mean(x$Z1), lty = 2, col = "red")
  }
  
  if (3 %in% which) {
    cl <- chain_ladder(x$triangle, x$F_cum)
    cc <- cape_cod(x$triangle, x$exposure, x$F_cum)
    reserves <- rbind(cl$reserve, cc$reserve, x$reserve)
    barplot(reserves, beside = TRUE, names.arg = 1:I,
            col = c("blue", "red", "darkgreen"),
            xlab = "Accident Year", ylab = "Reserve", main = "Reserve Comparison")
    legend("topleft", legend = c("CL", "CC", "UCR"),
           fill = c("blue", "red", "darkgreen"), bty = "n")
  }
}

# =============================================================================
# PART 7: SIMULATION
# =============================================================================

#' Generate claims triangle from Poisson-Gamma-Multinomial model
#'
#' @param I Number of accident years
#' @param J Number of development periods
#' @param exposure Vector of exposures (or single value)
#' @param mu Mean frailty (claim rate)
#' @param tau_sq Between-year variance of frailty
#' @param pi Development pattern (sums to 1)
#' @param seed Random seed (optional)
#' @return A list with triangle, ultimates, and frailties
#' @export
simulate_triangle <- function(I, J, exposure = 1000, mu = 0.1, tau_sq = 0.01, 
                               pi = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  if (length(exposure) == 1) exposure <- rep(exposure, I)
  
  if (is.null(pi)) {
    pi <- c(0.15, 0.25, 0.20, 0.15, 0.10, 0.06, 0.04, 0.03, 0.015, 0.005)
    if (J != 10) {
      pi <- pi[1:min(J, length(pi))]
      if (length(pi) < J) pi <- c(pi, rep(0.01, J - length(pi)))
      pi <- pi / sum(pi)
    }
  }
  
  # Gamma frailties
  if (tau_sq > 0) {
    alpha <- mu^2 / tau_sq
    beta <- mu / tau_sq
    Lambda <- rgamma(I, shape = alpha, rate = beta)
  } else {
    Lambda <- rep(mu, I)
  }
  
  # Ultimate counts
  N_ult <- rpois(I, lambda = Lambda * exposure)
  
  # Multinomial allocation
  N_incr <- matrix(0, nrow = I, ncol = J)
  for (i in 1:I) {
    if (N_ult[i] > 0) {
      N_incr[i, ] <- rmultinom(1, size = N_ult[i], prob = pi)
    }
  }
  
  # Cumulative
  triangle <- t(apply(N_incr, 1, cumsum))
  
  # Lower triangle mask
  for (i in 1:I) {
    if (i < I) {
      triangle[i, (I - i + 2):J] <- NA
    }
  }
  
  rownames(triangle) <- paste0("AY", 1:I)
  colnames(triangle) <- paste0("Dev", 1:J)
  
  list(
    triangle = triangle,
    N_ult = N_ult,
    Lambda = Lambda,
    exposure = exposure,
    mu = mu,
    tau_sq = tau_sq,
    pi = pi,
    F_cum = cumsum(pi)
  )
}

#' Run UCR simulation study
#'
#' @param n_reps Number of replications
#' @param tau_sq_grid Grid of tau^2 values to test
#' @param I Number of accident years
#' @param J Number of development periods
#' @param mu Mean frailty
#' @param exposure Exposure level
#' @param pi Development pattern
#' @param verbose Print progress
#' @return A data frame with simulation results
#' @export
run_ucr_simulation <- function(n_reps = 500, 
                                tau_sq_grid = c(0, 0.0001, 0.0005, 0.001, 0.005, 0.01),
                                I = 10, J = 10, mu = 0.1, exposure = 1000,
                                pi = NULL, verbose = TRUE) {
  
  results <- data.frame()
  
  for (tau_sq in tau_sq_grid) {
    if (verbose) cat("  tau^2 =", tau_sq, "...")
    
    for (r in 1:n_reps) {
      sim <- simulate_triangle(I, J, exposure, mu, tau_sq, pi)
      
      cl <- chain_ladder(sim$triangle)
      cc <- cape_cod(sim$triangle, sim$exposure)
      ucr_fit <- ucr(sim$triangle, sim$exposure)
      
      incomplete <- which(ucr_fit$j_obs < J)
      if (length(incomplete) > 0) {
        MSE_CL <- mean((cl$ultimate[incomplete] - sim$N_ult[incomplete])^2)
        MSE_CC <- mean((cc$ultimate[incomplete] - sim$N_ult[incomplete])^2)
        MSE_UCR <- mean((ucr_fit$ultimate[incomplete] - sim$N_ult[incomplete])^2)
      } else {
        MSE_CL <- MSE_CC <- MSE_UCR <- NA
      }
      
      results <- rbind(results, data.frame(
        tau_sq_true = tau_sq,
        tau_sq_hat = ucr_fit$tau_sq_hat,
        Z1_mean = mean(ucr_fit$Z1),
        MSE_CL = MSE_CL,
        MSE_CC = MSE_CC,
        MSE_UCR = MSE_UCR
      ))
    }
    
    if (verbose) cat(" done\n")
  }
  
  return(results)
}

# =============================================================================
# PART 8: CREDIBILITY-BAYES EQUIVALENCE
# =============================================================================

#' Posterior distribution under Poisson-Gamma model
#'
#' @param n_obs Observed count
#' @param exposure Exposure
#' @param alpha Prior shape parameter
#' @param beta Prior rate parameter
#' @return A list with posterior parameters and summary statistics
#' @export
poisson_gamma_posterior <- function(n_obs, exposure, alpha, beta) {
  alpha_post <- alpha + n_obs
  beta_post <- beta + exposure
  
  post_mean <- alpha_post / beta_post
  post_var <- alpha_post / beta_post^2
  Z <- exposure / (exposure + beta)
  prior_mean <- alpha / beta
  x_bar <- n_obs / exposure
  
  list(
    alpha_post = alpha_post,
    beta_post = beta_post,
    mean = post_mean,
    var = post_var,
    sd = sqrt(post_var),
    Z = Z,
    prior_mean = prior_mean,
    x_bar = x_bar,
    credibility_check = Z * x_bar + (1 - Z) * prior_mean
  )
}

#' Demonstrate credibility = posterior mean
#'
#' @param n_obs Observed count
#' @param exposure Exposure
#' @param mu Prior mean of Lambda
#' @param tau_sq Prior variance of Lambda
#' @return A list comparing credibility and Bayesian estimates
#' @export
credibility_bayes_equivalence <- function(n_obs, exposure, mu, tau_sq) {
  alpha <- mu^2 / tau_sq
  beta <- mu / tau_sq
  
  bayes <- poisson_gamma_posterior(n_obs, exposure, alpha, beta)
  
  sigma_sq <- mu
  k <- sigma_sq / tau_sq
  Z <- exposure / (exposure + k)
  x_bar <- n_obs / exposure
  cred_estimate <- Z * x_bar + (1 - Z) * mu
  
  list(
    bayesian_mean = bayes$mean,
    credibility_estimate = cred_estimate,
    difference = bayes$mean - cred_estimate,
    Z_bayesian = bayes$Z,
    Z_credibility = Z,
    note = "Difference should be ~0: credibility = posterior mean"
  )
}

# MML estimator of (mu, tau^2) under PGMD via NegBin marginal
# Drop-in alongside hgr's BS estimator.

mml_estimator <- function(N_obs, w) {
  # N_obs: vector of observed cumulative counts per AY
  # w:     vector of effective exposures, w_i = E_i * F_{j_i}
  # Returns list with mu_hat, tau2_hat.
  #
  # Uses MASS::glm.nb with log-link and offset log(w).
  # MASS parameterisation: var(N) = mu + mu^2/theta.
  # PGMD parameterisation: var(N|mu,tau2) = m + (tau2/mu) * m^2 where m = mu*w.
  # Matching: 1/theta = tau2/mu, i.e. tau2 = mu/theta. (Note: NOT mu^2/theta.)
  # Equivalently, theta = mu^2/tau2 (the NB size parameter from Lemma 2.4).

  fit <- MASS::glm.nb(N_obs ~ 1 + offset(log(w)),
                      link = log,
                      control = glm.control(maxit = 200))
  mu_hat   <- as.numeric(exp(coef(fit)))
  theta    <- fit$theta                # NB size; theta = mu^2 / tau2
  tau2_hat <- mu_hat^2 / theta

  list(mu = mu_hat, tau2 = tau2_hat, theta = theta, fit = fit)
}

# Sanity check on Taylor-Ashe (should give numbers close to but not identical to BS):
# library(ChainLadder); data(GenIns)
# tri <- as.triangle(GenIns)  # adapt to hgr conventions
# ... compute N_obs and w from triangle ...
# bs  <- hgr::ucr_bs(N_obs, w)
# mml <- mml_estimator(N_obs, w)
# print(c(bs_tau2 = bs$tau2, mml_tau2 = mml$tau2))

