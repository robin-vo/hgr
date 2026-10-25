## ===========================================================================
## EXACT FRAILTY-AUGMENTED POSTERIOR  (rstan; the paper's exact arms:
## clean_run/cnbg_exact.R, R/run_phase4_1_exact.R, R/cnbg_exact_free.stan)
## ===========================================================================

#' Exact frailty-augmented CNBG posterior (rstan)
#'
#' Fits the exact posterior of the single-frailty CNBG model by Hamiltonian
#' Monte Carlo: `kappa_B ~ Gamma(prior)`, `A ~ Gamma(2, 2 / A_init)`,
#' `U_i ~ Gamma(kappa_B, kappa_B)`, and each cell compound Poisson--Gamma
#' given its mean, the density evaluated by a truncated latent-count sum.
#' The severity shape `alpha` is held fixed.
#'
#' `levels = "anchored"` holds the cell means at the Stage-1 free-level
#' fitted values (the paper's Schedule P comparison of the exact and working
#' posteriors, Table exact_panel). `levels = "free"` makes the accident-year
#' levels `mu_i` parameters with a vague lognormal prior `lognormal(log mu0_i,
#' mu_prior_sdlog)` around a row-total start, the free-level exact arm of the
#' paper's Section S5 (and its Phase 17 sensitivity to `mu_prior_sdlog`),
#' with the development pattern `w` held fixed (the Stage-1 pattern unless
#' supplied).
#'
#' The exact posterior costs roughly 200 times the working fit and is the
#' paper's per-portfolio reference, not its production estimator.
#'
#' @param triangle Incremental paid-loss triangle.
#' @param prior `"ref"`, `"empb"` or `c(shape, rate)`.
#' @param alpha Severity shape, held fixed.
#' @param levels `"anchored"` or `"free"`.
#' @param w Development weights for `levels = "free"` (default: Stage-1).
#' @param mu_prior_sdlog Spread of the lognormal level prior (`2` in the
#'   paper; `0.5` and `4` in its sensitivity).
#' @param kappa_cap Cap applied to the `kappa_B` draws.
#' @param chains,iter,warmup,seed,adapt_delta,max_treedepth Sampler settings
#'   (the paper's Phase 4 used two chains of 1,000 iterations with 500
#'   warm-up; the defaults here are more conservative).
#' @param ... Passed to `rstan::sampling()`.
#' @return A `cnbg_fit` with `method = "exact"`, the posterior draws of
#'   `kappa_B`, `A`, `U` (and `mu` under free levels), their summaries, the
#'   sampler diagnostics and the `stanfit`.
#' @examples
#' \dontrun{
#' tri <- simulate_cnbg(kappa_r = 5, alpha = 3)$triangle
#' fx <- fit_cnbg_exact(tri, levels = "anchored", chains = 2, iter = 1000, warmup = 500)
#' fx
#' reserve_cnbg(fx)
#' diagnose_kappa(fx)
#' }
#' @export
fit_cnbg_exact <- function(triangle, prior = "ref", alpha = 3, levels = c("anchored", "free"),
                           w = NULL, mu_prior_sdlog = 2, kappa_cap = .cnbg_kappa_cap,
                           chains = 4L, iter = 2000L, warmup = 1000L, seed = 2026L,
                           adapt_delta = 0.95, max_treedepth = 12L, ...) {
  levels <- match.arg(levels)
  pr <- .cnbg_prior(prior)
  ms <- tryCatch(cnbg_stage1(triangle, "free"), error = function(e) NULL)
  if (is.null(ms) || !ms$converged)
    return(structure(list(ok = FALSE, method = "exact", message = "Stage-1 GLM failed", ms = ms), class = "cnbg_fit"))
  obs <- ms$obs; I <- ms$I; J <- ms$J
  nu   <- as.numeric(ms$nu_mat[obs]); x <- as.numeric(ms$tri[obs])
  ij   <- which(obs, arr.ind = TRUE); rows <- ij[, 1]; cols <- ij[, 2]
  A_init <- max(median((x - nu)^2 / nu), .Machine$double.eps)
  Ybar0 <- A_init * alpha / (alpha + 1)
  n_max <- as.integer(ceiling(3 * max(x / Ybar0, 1) + 50))
  standata <- list(N = length(x), x = x, I = I, row = as.integer(rows), alpha = alpha, n_max = n_max,
                   prior_shape = pr[["shape"]], prior_rate = pr[["rate"]], A_prior_shape = 2, A_prior_rate = 2 / A_init)
  if (levels == "free") {
    if (is.null(w)) w <- ms$nu_mat[1, ] / sum(ms$nu_mat[1, ])
    w <- w / sum(w)
    rt <- tapply(x, rows, sum); ws <- tapply(w[cols], rows, sum)
    standata$w <- w[cols]; standata$mu0 <- as.numeric(pmax(rt, 1) / (Ybar0 * ws)); standata$mu_prior_sdlog <- mu_prior_sdlog
    init_fun <- function() list(kappa = pr[["shape"]] / pr[["rate"]], A = A_init, U = rep(1, I), mu = standata$mu0)
    model <- "cnbg_exact_free"
  } else {
    standata$nu <- nu
    init_fun <- function() list(kappa = pr[["shape"]] / pr[["rate"]], A = A_init, U = rep(1, I))
    model <- "cnbg_exact"
  }
  fit <- rstan::sampling(.cnbg_stan_model(model), data = standata, chains = chains, iter = iter, warmup = warmup,
                         seed = seed, init = init_fun, refresh = 0,
                         control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth), ...)
  sm <- rstan::summary(fit, pars = c("kappa", "A"))$summary
  rhat_max <- suppressWarnings(max(sm[, "Rhat"], na.rm = TRUE))
  ess_min  <- suppressWarnings(min(sm[, "n_eff"], na.rm = TRUE))
  ndiv <- sum(vapply(rstan::get_sampler_params(fit, inc_warmup = FALSE),
                     function(p) sum(p[, "divergent__"]), numeric(1)))
  if (is.finite(rhat_max) && rhat_max > 1.01)
    warning(sprintf("exact fit: max Rhat = %.3f (> 1.01); treat the draws with care", rhat_max))
  if (ndiv > 0) warning(sprintf("exact fit: %d divergent transition(s); consider a larger adapt_delta", ndiv))
  post <- rstan::extract(fit, pars = c("kappa", "A", "U"))
  mu_post <- if (levels == "free") rstan::extract(fit, pars = "mu")$mu else NULL
  kd <- pmin(post$kappa, kappa_cap)
  q  <- function(p) unname(quantile(kd, p))
  structure(list(ok = TRUE, method = "exact", ms = ms, prior = pr,
                 prior_label = if (is.numeric(prior)) "custom" else prior, alpha = alpha, levels = levels,
                 w = if (levels == "free") w else NULL,
                 kappa_draws = kd, A_draws = post$A, U_draws = post$U, mu_draws = mu_post,
                 kappa_mean = mean(kd), kappa_median = median(kd), kappa_q025 = q(0.025), kappa_q975 = q(0.975),
                 kappa_sd = sd(kd), A_init = A_init, n_max = n_max,
                 edge_mass = mean(kd >= 0.95 * kappa_cap),
                 rhat_max = rhat_max, n_eff_min = ess_min, n_divergent = ndiv, stanfit = fit),
            class = "cnbg_fit")
}

## exact predictive: compound law at the posterior draws (R/run_phase4_1_exact.R)
.cnbg_exact_predict <- function(fit, B = 4000L, level = 0.95, ...) {
  ms <- fit$ms; I <- ms$I; J <- ms$J; a <- fit$alpha
  fut <- !.cnbg_obs_mask(I, J); fut_idx <- which(fut)
  fut_i <- ((fut_idx - 1L) %% I) + 1L; fut_j <- ((fut_idx - 1L) %/% I) + 1L
  nd <- length(fit$kappa_draws); idx <- sample.int(nd, B, replace = TRUE); R <- numeric(B)
  for (b in seq_len(B)) {
    d <- idx[b]; Yb <- fit$A_draws[d] * a / (a + 1)
    nf <- if (fit$levels == "free") fit$mu_draws[d, fut_i] * fit$w[fut_j] * Yb else ms$nu_mat[fut_idx]
    Nn <- rpois(length(fut_idx), nf * fit$U_draws[d, fut_i] / Yb); pos <- Nn > 0L
    R[b] <- if (any(pos)) sum(rgamma(sum(pos), shape = Nn[pos] * a, rate = a / Yb)) else 0
  }
  al <- (1 - level) / 2; q <- unname(quantile(R, c(al, 1 - al)))
  structure(list(reserve = R, mean = mean(R), median = median(R), sd = sd(R), pi = q, width = q[2] - q[1],
                 level = level, cl_est = sum(ms$nu_mat[fut]),
                 method = if (fit$levels == "free") "CNBG exact (free levels)" else "CNBG exact (anchored levels)"),
            class = "cnbg_reserve")
}
