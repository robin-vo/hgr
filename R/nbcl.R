#' @title Negative Binomial Chain-Ladder Model
#' @description Functions for fitting and predicting with the NB-CL model
#' @author Robin Van Oirbeek
#' @import MASS
#' @import stats
#' @importFrom R.utils withTimeout

# =============================================================================
# Utility functions
# =============================================================================

#' Convert run-off triangle to long format
#'
#' @param triangle A matrix representing the run-off triangle (incremental claims)
#' @return A data frame with columns AY, DY, A
#' @export
triangle_to_long <- function(triangle) {
  I <- nrow(triangle)
  J <- ncol(triangle)

  data <- data.frame(
    AY = integer(),
    DY = integer(),
    A = numeric()
  )

  for (i in 1:I) {
    for (j in 1:J) {
      if (!is.na(triangle[i, j])) {
        data <- rbind(data, data.frame(AY = i, DY = j - 1, A = triangle[i, j]))
      }
    }
  }

  data$AY <- factor(data$AY)
  data$DY <- factor(data$DY)

  return(data)
}

#' Convert long format back to triangle
#'
#' @param df A data frame with columns AY, DY, A
#' @param I Number of accident years
#' @param J Number of development years
#' @return A matrix representing the run-off triangle
#' @export
long_to_triangle <- function(df, I = NULL, J = NULL) {
  if (is.null(I)) I <- max(as.integer(as.character(df$AY)))
  if (is.null(J)) J <- max(as.integer(as.character(df$DY))) + 1

  triangle <- matrix(NA, nrow = I, ncol = J)

  for (k in 1:nrow(df)) {
    i <- as.integer(as.character(df$AY[k]))
    j <- as.integer(as.character(df$DY[k])) + 1
    triangle[i, j] <- df$A[k]
  }

  rownames(triangle) <- 1:I
  colnames(triangle) <- 0:(J - 1)

  return(triangle)
}

#' Reparameterise glm.nb output to the simplex convention
#'
#' Standard glm.nb output uses treatment contrasts (beta_0 = 0). This function
#' converts to the simplex convention sum_j exp(beta_j) = 1, under which
#' mu_i = exp(alpha_i) is the expected total claim count for accident year i
#' and w_j = exp(beta_j) is the proportion of ultimate claims reported in
#' development year j.
#'
#' The cell means mu_{ij} = exp(alpha_i + beta_j) are invariant under this
#' transformation; only the interpretation of individual coefficients changes.
#'
#' @param fit A fitted glm.nb object with AY and DY factors
#' @return A list with components alpha, beta, mu, w (see Details).
#' @export
reparam_simplex <- function(fit) {
  cf <- coef(fit)
  intercept <- unname(cf["(Intercept)"])

  ay_levels <- levels(model.frame(fit)$AY)
  dy_levels <- levels(model.frame(fit)$DY)

  # Build full alpha and beta vectors (reference level = 0)
  alpha_tilde <- setNames(rep(0, length(ay_levels)), ay_levels)
  beta_tilde <- setNames(rep(0, length(dy_levels)), dy_levels)

  alpha_tilde[-1] <- unname(cf[grep("^AY", names(cf))])
  beta_tilde[-1] <- unname(cf[grep("^DY", names(cf))])

  # Absorb intercept into alpha
  alpha_tilde <- alpha_tilde + intercept

  # Simplex transformation
  S <- sum(exp(beta_tilde))
  alpha <- alpha_tilde + log(S)
  beta <- beta_tilde - log(S)

  list(
    alpha = alpha,  # log AY totals; mu_i = exp(alpha_i)
    beta = beta,    # log development weights; sum_j exp(beta_j) = 1
    mu = exp(alpha),  # AY totals
    w = exp(beta)     # development weights, sum to 1
  )
}

# =============================================================================
# Core model fitting
# =============================================================================

#' Fit Negative Binomial Chain-Ladder model
#'
#' @param triangle A matrix representing the run-off triangle (incremental claims)
#' @return An object of class "nbcl" containing the fitted model and metadata
#' @export
fit_nbcl <- function(triangle) {
  df <- triangle_to_long(triangle)
  n <- nrow(df)

  fit <- MASS::glm.nb(A ~ AY + DY, data = df, link = log)

  p <- length(coef(fit))
  kappa_mle <- fit$theta
  kappa_adj <- kappa_mle * (n - p) / n

  # Simplex parameterisation: mu_i = exp(alpha_i), sum_j exp(beta_j) = 1
  simplex <- reparam_simplex(fit)

  result <- list(
    fit = fit,
    triangle = triangle,
    data = df,
    n = n,
    p = p,
    kappa_mle = kappa_mle,
    kappa_corrected = kappa_adj,
    correction_factor = (n - p) / n,
    alpha = simplex$alpha,  # log AY totals
    beta = simplex$beta,    # log development weights
    mu = simplex$mu,        # AY totals
    w = simplex$w,          # development weights
    aic = AIC(fit),
    bic = BIC(fit),
    loglik = as.numeric(logLik(fit))
  )

  class(result) <- "nbcl"
  return(result)
}

#' Coefficients in the simplex parameterisation
#'
#' Returns alpha (log AY totals) and beta (log development weights satisfying
#' sum_j exp(beta_j) = 1). This is the parameterisation used throughout the
#' NB-CL paper and differs from the treatment-contrast coefficients returned
#' by the underlying glm.nb fit.
#'
#' @param object An object of class "nbcl"
#' @param ... Additional arguments (ignored)
#' @return A named numeric vector of simplex-parameterised coefficients.
#' @export
coef.nbcl <- function(object, ...) {
  c(setNames(object$alpha, paste0("alpha[", names(object$alpha), "]")),
    setNames(object$beta, paste0("beta[", names(object$beta), "]")))
}

#' Extract bias-corrected kappa
#'
#' @param object An object of class "nbcl"
#' @return The bias-corrected dispersion parameter
#' @export
kappa_corrected <- function(object) {
  if (inherits(object, "nbcl")) {
    return(object$kappa_corrected)
  } else if (inherits(object, "negbin")) {
    n <- length(fitted(object))
    p <- length(coef(object))
    return(object$theta * (n - p) / n)
  } else {
    stop("Object must be of class 'nbcl' or 'negbin'")
  }
}

#' Print method for nbcl objects
#'
#' @param x An object of class "nbcl"
#' @param ... Additional arguments (ignored)
#' @export
print.nbcl <- function(x, ...) {
  cat("Negative Binomial Chain-Ladder Model\n")
  cat("=====================================\n\n")
  cat("Triangle size:", nrow(x$triangle), "x", ncol(x$triangle), "\n")
  cat("Observations:", x$n, "\n")
  cat("Parameters:", x$p, "\n\n")
  cat("Dispersion parameter (kappa):\n")
  cat("  MLE:       ", round(x$kappa_mle, 3), "\n")
  cat("  Corrected: ", round(x$kappa_corrected, 3), "\n")
  cat("  Correction factor (n-p)/n:", round(x$correction_factor, 3), "\n\n")
  cat("Accident-year totals (mu_i = exp(alpha_i)):\n")
  print(round(x$mu, 1))
  cat("\nDevelopment weights (w_j = exp(beta_j), sum to 1):\n")
  print(round(x$w, 4))
  cat("  sum =", round(sum(x$w), 6), "\n\n")
  cat("Model fit:\n")
  cat("  Log-likelihood:", round(x$loglik, 2), "\n")
  cat("  AIC:", round(x$aic, 2), "\n")
  cat("  BIC:", round(x$bic, 2), "\n")
}

#' Summary method for nbcl objects
#'
#' @param object An object of class "nbcl"
#' @param ... Additional arguments (ignored)
#' @export
summary.nbcl <- function(object, ...) {
  cat("Negative Binomial Chain-Ladder Model Summary\n")
  cat("=============================================\n\n")
  print(object)
  cat("\nUnderlying GLM coefficients (treatment contrasts, beta_0 = 0):\n")
  print(summary(object$fit)$coefficients)
  cat("\nNote: paper convention uses simplex parameterisation\n")
  cat("      (sum_j exp(beta_j) = 1; see fields $alpha, $beta, $mu, $w).\n")
}

# =============================================================================
# Prediction
# =============================================================================

#' Predict lower triangle
#'
#' @param object An object of class "nbcl"
#' @param newdata Optional new data for prediction
#' @return A data frame with predicted values for the lower triangle
#' @export
predict_nbcl <- function(object, newdata = NULL) {
  if (!inherits(object, "nbcl")) {
    stop("Object must be of class 'nbcl'")
  }

  I <- nrow(object$triangle)
  J <- ncol(object$triangle)

  # Create prediction grid for lower triangle
  future <- data.frame(AY = integer(), DY = integer())
  for (i in 2:I) {
    for (j in (I - i + 1):(J - 1)) {
      future <- rbind(future, data.frame(AY = i, DY = j))
    }
  }

  future$AY <- factor(future$AY, levels = levels(object$data$AY))
  future$DY <- factor(future$DY, levels = levels(object$data$DY))

  future$fitted <- predict(object$fit, newdata = future, type = "response")
  future$kappa_mle <- object$kappa_mle
  future$kappa_corrected <- object$kappa_corrected

  return(future)
}

reserve_nbcl <- function(object) {
  if (!inherits(object, "nbcl")) {
    stop("Object must be of class 'nbcl'")
  }

  triangle <- object$triangle
  I <- nrow(triangle)
  J <- ncol(triangle)

  # Compute cumulative triangle (treating NAs as missing future cells)
  cum_tri <- t(apply(triangle, 1, function(row) {
    cumsum(replace(row, is.na(row), 0))
  }))

  # Restore NAs in the lower triangle
  for (i in 1:I) {
    for (j in 1:J) {
      if (is.na(triangle[i, j])) cum_tri[i, j] <- NA
    }
  }

  # Chain-Ladder development factors
  dev_factors <- numeric(J - 1)
  for (j in 1:(J - 1)) {
    num <- sum(cum_tri[1:(I - j), j + 1], na.rm = TRUE)
    den <- sum(cum_tri[1:(I - j), j], na.rm = TRUE)
    dev_factors[j] <- num / den
  }

  # Project ultimates
  ultimates <- numeric(I)
  ultimates[1] <- cum_tri[1, J]
  for (i in 2:I) {
    last_obs <- J - i + 1
    ultimates[i] <- cum_tri[i, last_obs]
    if (last_obs < J) {
      for (j in last_obs:(J - 1)) {
        ultimates[i] <- ultimates[i] * dev_factors[j]
      }
    }
  }

  # Reserves = ultimate - latest observed cumulative
  reserves <- numeric(I)
  for (i in 1:I) {
    last_obs <- J - i + 1
    reserves[i] <- ultimates[i] - cum_tri[i, last_obs]
  }

  by_ay <- data.frame(
    AY = 1:I,
    Reserve = reserves
  )
  by_ay <- by_ay[by_ay$AY > 1, ]  # AY 1 has no reserve
  rownames(by_ay) <- NULL

  list(
    total = sum(reserves[-1]),
    by_accident_year = by_ay,
    ultimates = ultimates,
    dev_factors = dev_factors
  )
}

# =============================================================================
# Bootstrap
# =============================================================================

#' Parametric bootstrap for NB-CL model
#'
#' Implements the unconditional parametric bootstrap (Algorithm 1 of the
#' NB-CL paper). Each expensive refit is wrapped in a hard per-fit timeout;
#' timeouts and failed refits are recorded as NA and reported via
#' \code{n_failed}, and never hang the run. If \code{checkpoint_file} is
#' supplied, results are checkpointed incrementally via atomic
#' write-then-rename, and an interrupted run resumes from the last
#' checkpoint on re-invocation with the same arguments.
#'
#' @param object An object of class "nbcl"
#' @param B Number of bootstrap samples (default 5000)
#' @param correct_kappa Whether to apply bias correction (default TRUE)
#' @param seed Random seed for reproducibility
#' @param verbose Print progress (default TRUE)
#' @param timeout Hard per-refit timeout in seconds (default 30). On timeout
#'   the replicate is recorded as failed (NA) and the run continues.
#' @param checkpoint_file Optional path to an RDS checkpoint file. If supplied,
#'   partial results are saved every \code{checkpoint_every} replicates and a
#'   re-run resumes from the last checkpoint.
#' @param checkpoint_every Checkpoint frequency in replicates (default 500).
#' @return An object of class "nbcl_boot" containing bootstrap results
#' @export
bootstrap_nbcl <- function(object, B = 5000, correct_kappa = TRUE,
                           seed = NULL, verbose = TRUE,
                           timeout = 30, checkpoint_file = NULL,
                           checkpoint_every = 500) {
  if (!inherits(object, "nbcl")) {
    stop("Object must be of class 'nbcl'")
  }

  if (!is.null(seed)) set.seed(seed)

  I <- nrow(object$triangle)
  n <- object$n
  p <- object$p

  # Use corrected or MLE kappa for initial sampling
  kappa_init <- if (correct_kappa) object$kappa_corrected else object$kappa_mle

  reserves_total <- rep(NA_real_, B)
  reserves_by_ay <- matrix(NA, nrow = B, ncol = I - 1)
  kappas <- rep(NA_real_, B)
  b_start <- 1L

  # Resume from checkpoint if present and compatible
  if (!is.null(checkpoint_file) && file.exists(checkpoint_file)) {
    cp <- readRDS(checkpoint_file)
    if (identical(cp$B, B)) {
      reserves_total <- cp$reserves_total
      reserves_by_ay <- cp$reserves_by_ay
      kappas <- cp$kappas
      b_start <- cp$b_done + 1L
      if (verbose) cat("Resuming from checkpoint at b =", b_start, "\n")
    }
  }

  save_ckpt <- function(b_done) {
    if (is.null(checkpoint_file)) return(invisible(NULL))
    tmp <- paste0(checkpoint_file, ".tmp")
    saveRDS(list(B = B, b_done = b_done, reserves_total = reserves_total,
                 reserves_by_ay = reserves_by_ay, kappas = kappas), tmp)
    file.rename(tmp, checkpoint_file)  # atomic write-then-rename
  }

  future_grid <- predict_nbcl(object)[, c("AY", "DY")]

  if (b_start > B && verbose) cat("All", B, "samples already complete.\n")
  for (b in seq.int(b_start, length.out = max(0L, B - b_start + 1L))) {

    # Step 1: Simulate observed triangle
    object$data$A_boot <- rnbinom(n, size = kappa_init, mu = fitted(object$fit))

    # Step 2: Re-fit model (hard per-fit timeout; timeout or error -> NA, continue)
    fit_boot <- tryCatch({
      R.utils::withTimeout(
        suppressWarnings(
          MASS::glm.nb(A_boot ~ AY + DY, data = object$data, link = log,
                       control = glm.control(maxit = 50))
        ),
        timeout = timeout, onTimeout = "silent"
      )
    }, error = function(e) NULL)

    if (is.null(fit_boot)) {
      reserves_total[b] <- NA
      next
    }

    # Step 3: Apply correction if requested
    kappa_boot <- if (correct_kappa) {
      fit_boot$theta * (n - p) / n
    } else {
      fit_boot$theta
    }
    kappas[b] <- kappa_boot

    # Step 4: Predict future and simulate
    mu_boot <- predict(fit_boot, newdata = future_grid, type = "response")
    future_claims <- rnbinom(length(mu_boot), size = kappa_boot, mu = mu_boot)

    # Store results (index per-AY row explicitly by accident year, so a
    # non-square or partially observed triangle cannot desync the columns)
    reserves_total[b] <- sum(future_claims)
    future_grid$claims <- future_claims
    by_ay <- aggregate(claims ~ AY, data = future_grid, FUN = sum)
    by_ay <- by_ay[order(as.integer(as.character(by_ay$AY))), ]
    reserves_by_ay[b, as.integer(as.character(by_ay$AY)) - 1] <- by_ay$claims

    if (!is.null(checkpoint_file) && b %% checkpoint_every == 0) save_ckpt(b)

    if (verbose && b %% 1000 == 0) {
      cat("Completed", b, "of", B, "bootstrap samples\n")
    }
  }
  save_ckpt(B)

  # Remove failed iterations
  valid <- !is.na(reserves_total)

  result <- list(
    reserves_total = reserves_total[valid],
    reserves_by_ay = reserves_by_ay[valid, , drop = FALSE],
    kappas = kappas[valid],
    B = sum(valid),
    B_requested = B,
    n_failed = sum(!valid),
    correct_kappa = correct_kappa,
    point_estimate = reserve_nbcl(object)$total
  )

  class(result) <- "nbcl_boot"
  return(result)
}

#' Extract prediction intervals from bootstrap
#'
#' @param object An object of class "nbcl_boot"
#' @param level Confidence level (default 0.95)
#' @return A data frame with point estimate and prediction intervals
#' @export
predict_interval <- function(object, level = 0.95) {
  if (!inherits(object, "nbcl_boot")) {
    stop("Object must be of class 'nbcl_boot'")
  }

  alpha <- 1 - level
  probs <- c(alpha / 2, 1 - alpha / 2)

  total_quantiles <- quantile(object$reserves_total, probs = probs)

  result <- data.frame(
    point_estimate = object$point_estimate,
    mean = mean(object$reserves_total),
    sd = sd(object$reserves_total),
    lower = total_quantiles[1],
    upper = total_quantiles[2],
    level = level
  )
  rownames(result) <- NULL

  return(result)
}

#' Print method for nbcl_boot objects
#'
#' @param x An object of class "nbcl_boot"
#' @param ... Additional arguments (ignored)
#' @export
print.nbcl_boot <- function(x, ...) {
  cat("NB-CL Bootstrap Results\n")
  cat("=======================\n\n")
  cat("Bootstrap samples:", x$B, "of", x$B_requested, "requested",
      if (x$n_failed > 0) paste0("(", x$n_failed, " failed/timed out)") else "",
      "\n")
  cat("Kappa correction:", if (x$correct_kappa) "Yes (REML-like)" else "No (MLE)", "\n\n")
  cat("Total Reserve:\n")
  cat("  Point estimate:", format(round(x$point_estimate), big.mark = ","), "\n")
  cat("  Bootstrap mean:", format(round(mean(x$reserves_total)), big.mark = ","), "\n")
  cat("  Bootstrap S.E.:", format(round(sd(x$reserves_total)), big.mark = ","), "\n")
  q95 <- quantile(x$reserves_total, probs = c(0.025, 0.975))
  cat("  95% PI: [", format(round(q95[1]), big.mark = ","), ", ",
      format(round(q95[2]), big.mark = ","), "]\n", sep = "")
}

# =============================================================================
# Diagnostics
# =============================================================================

#' Profile likelihood for kappa
#'
#' @param object An object of class "nbcl"
#' @param kappa_range Range of kappa values to evaluate (default: auto)
#' @param n_points Number of points in the grid (default 50)
#' @param level Confidence level for CI (default 0.95)
#' @return An object of class "nbcl_profile" with profile likelihood results
#' @export
profile_kappa <- function(object, kappa_range = NULL, n_points = 50, level = 0.95) {
  if (!inherits(object, "nbcl")) {
    stop("Object must be of class 'nbcl'")
  }

  if (is.null(kappa_range)) {
    kappa_range <- c(
      max(1, object$kappa_mle / 4),
      object$kappa_mle * 4
    )
  }

  kappa_grid <- seq(kappa_range[1], kappa_range[2], length.out = n_points)
  profile_ll <- numeric(n_points)

  for (k in seq_along(kappa_grid)) {
    fit_fixed <- glm(A ~ AY + DY, data = object$data,
                     family = MASS::negative.binomial(theta = kappa_grid[k]))
    profile_ll[k] <- as.numeric(logLik(fit_fixed))
  }

  profile_ll_norm <- profile_ll - max(profile_ll)

  # Find CI using chi-squared threshold
  threshold <- -qchisq(level, df = 1) / 2
  ci_idx <- which(profile_ll_norm >= threshold)
  ci_lower <- min(kappa_grid[ci_idx])
  ci_upper <- max(kappa_grid[ci_idx])

  result <- list(
    kappa_grid = kappa_grid,
    profile_ll = profile_ll,
    profile_ll_norm = profile_ll_norm,
    kappa_mle = object$kappa_mle,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    level = level
  )

  class(result) <- "nbcl_profile"
  return(result)
}

#' Plot profile likelihood
#'
#' @param x An object of class "nbcl_profile"
#' @param ... Additional arguments passed to plot
#' @export
plot.nbcl_profile <- function(x, ...) {
  threshold <- -qchisq(x$level, df = 1) / 2

  plot(x$kappa_grid, x$profile_ll_norm,
       type = "l", lwd = 2,
       xlab = expression(kappa),
       ylab = "Profile log-likelihood (normalized)",
       main = expression(paste("Profile likelihood for ", kappa)),
       ...)
  abline(v = x$kappa_mle, lty = 2, col = "blue", lwd = 1.5)
  abline(h = threshold, lty = 3, col = "darkgray")

  ci_idx <- which(x$profile_ll_norm >= threshold)
  polygon(c(x$ci_lower, x$kappa_grid[ci_idx], x$ci_upper),
          c(-10, x$profile_ll_norm[ci_idx], -10),
          col = rgb(0, 0, 1, 0.1), border = NA)

  legend("topright",
         legend = c(paste("MLE =", round(x$kappa_mle, 1)),
                    paste0(x$level * 100, "% CI: [", round(x$ci_lower, 1),
                           ", ", round(x$ci_upper, 1), "]")),
         lty = c(2, NA), col = c("blue", NA), lwd = c(1.5, NA),
         bty = "n")
}

#' Extract residuals from NB-CL model
#'
#' @param object An object of class "nbcl"
#' @param type Type of residuals: "pearson", "deviance", or "response"
#' @param ... Additional arguments (ignored)
#' @return A numeric vector of residuals
#' @export
residuals.nbcl <- function(object, type = "pearson", ...) {
  residuals(object$fit, type = type)
}

#' Diagnostic plots for NB-CL model
#'
#' @param object An object of class "nbcl"
#' @param which Which plots to produce (1 = residuals vs fitted,
#'   2 = by AY, 3 = by DY, 4 = profile likelihood)
#' @export
plot_diagnostics <- function(object, which = 1:4) {
  if (!inherits(object, "nbcl")) {
    stop("Object must be of class 'nbcl'")
  }

  pearson_resid <- residuals(object$fit, type = "pearson")
  fitted_vals <- fitted(object$fit)

  n_plots <- length(which)
  if (n_plots > 1) {
    old_par <- par(mfrow = c(ceiling(n_plots / 2), min(n_plots, 2)),
                   mar = c(4.5, 4.5, 2, 1))
    on.exit(par(old_par))
  }

  if (1 %in% which) {
    plot(log(fitted_vals), pearson_resid,
         xlab = "log(Fitted values)",
         ylab = "Pearson residuals",
         main = "Residuals vs. Fitted",
         pch = 19, col = rgb(0, 0, 0, 0.6))
    abline(h = 0, lty = 2, col = "darkgray")
    abline(h = c(-2, 2), lty = 3, col = "darkgray")
    lines(lowess(log(fitted_vals), pearson_resid), col = "red", lwd = 2)
  }

  if (2 %in% which) {
    boxplot(pearson_resid ~ object$data$AY,
            xlab = "Accident Year",
            ylab = "Pearson residuals",
            main = "Residuals by Accident Year",
            col = "lightblue")
    abline(h = 0, lty = 2, col = "darkgray")
  }

  if (3 %in% which) {
    boxplot(pearson_resid ~ object$data$DY,
            xlab = "Development Year",
            ylab = "Pearson residuals",
            main = "Residuals by Development Year",
            col = "lightgreen")
    abline(h = 0, lty = 2, col = "darkgray")
  }

  if (4 %in% which) {
    profile <- profile_kappa(object)
    plot(profile)
  }
}

# =============================================================================
# Model comparison
# =============================================================================

#' Compare Poisson, ODP, and NB-CL models
#'
#' @param triangle A matrix representing the run-off triangle
#' @return A data frame comparing model fits
#' @export
compare_models <- function(triangle) {
  df <- triangle_to_long(triangle)
  n <- nrow(df)

  # Poisson
  fit_pois <- glm(A ~ AY + DY, data = df, family = poisson(link = "log"))

  # Quasi-Poisson (ODP)
  fit_odp <- glm(A ~ AY + DY, data = df, family = quasipoisson(link = "log"))
  phi_odp <- summary(fit_odp)$dispersion

  # Negative Binomial
  fit_nb <- MASS::glm.nb(A ~ AY + DY, data = df, link = log)
  p <- length(coef(fit_nb))
  kappa_mle <- fit_nb$theta
  kappa_adj <- kappa_mle * (n - p) / n

  result <- data.frame(
    Model = c("Poisson", "ODP (quasi-Poisson)", "NB-CL (MLE)", "NB-CL (corrected)"),
    Dispersion = c(1, phi_odp, kappa_mle, kappa_adj),
    LogLik = c(as.numeric(logLik(fit_pois)), NA, as.numeric(logLik(fit_nb)), NA),
    AIC = c(AIC(fit_pois), NA, AIC(fit_nb), NA),
    BIC = c(BIC(fit_pois), NA, BIC(fit_nb), NA)
  )

  return(result)
}

#' Likelihood ratio test for overdispersion
#'
#' @param object An object of class "nbcl"
#' @return A list with test statistic, df, and p-value
#' @export
lr_test_overdispersion <- function(object) {
  if (!inherits(object, "nbcl")) {
    stop("Object must be of class 'nbcl'")
  }

  # Fit Poisson model
  fit_pois <- glm(A ~ AY + DY, data = object$data, family = poisson(link = "log"))

  # LR statistic
  lr_stat <- 2 * (object$loglik - as.numeric(logLik(fit_pois)))

  # P-value (mixture of chi-squared, conservative approximation)
  p_value <- 0.5 * pchisq(lr_stat, df = 1, lower.tail = FALSE)

  result <- list(
    statistic = lr_stat,
    df = 1,
    p_value = p_value,
    kappa_mle = object$kappa_mle
  )

  cat("Likelihood Ratio Test for Overdispersion\n")
  cat("=========================================\n")
  cat("H0: Poisson (kappa = Inf) vs H1: Negative Binomial\n\n")
  cat("LR statistic:", round(lr_stat, 2), "\n")
  cat("P-value:", format.pval(p_value), "\n")
  cat("Estimated kappa:", round(object$kappa_mle, 2), "\n")

  if (p_value < 0.05) {
    cat("\nConclusion: Reject H0, significant overdispersion detected.\n")
  } else {
    cat("\nConclusion: Fail to reject H0, no significant overdispersion.\n")
  }

  invisible(result)
}
