## ===========================================================================
## WORKING-LIKELIHOOD POSTERIOR  (port of R/working_fit.R, R/cell_shapes.R and
## R/curvature_adjustment.R of the paper's replication material)
## ===========================================================================

## ---- cell shapes from the fourth moment (single-frailty CNBG) -------------
## With a = 1/kappa: kappa_2(M) = a, kappa_3(M) = 2 a^2, kappa_4(M) = 6 a^3,
## A = (alpha + 1) Ybar / alpha, c_3 = Ybar^2 (alpha+1)(alpha+2)/alpha^2,
## c_4 = Ybar^3 (alpha+1)(alpha+2)(alpha+3)/alpha^3,
## E[(X - nu)^4] = 3 mu_z^2 + c_4 nu + (4 c_3 + 3 A^2) B nu^2 + 6 A kappa_3 nu^3 + kappa_4 nu^4,
## CV2(z) = E[(X - nu)^4] / mu_z^2 - 1, s_ij = 1 / CV2_ij (Proposition 5.2 of the paper).
.cnbg_fourth_central_moment <- function(nu, A, kappa, alpha) {
  Ybar <- A * alpha / (alpha + 1)
  c3 <- Ybar^2 * (alpha + 1) * (alpha + 2) / alpha^2
  c4 <- Ybar^3 * (alpha + 1) * (alpha + 2) * (alpha + 3) / alpha^3
  a  <- 1 / kappa
  B  <- a; k3 <- 2 * a^2; k4 <- 6 * a^3
  mu_z <- A * nu + B * nu^2
  3 * mu_z^2 + c4 * nu + (4 * c3 + 3 * A^2) * B * nu^2 + 6 * A * k3 * nu^3 + k4 * nu^4
}

.cnbg_cv2_z <- function(nu, A, kappa, alpha) {
  mu_z <- A * nu + nu^2 / kappa
  .cnbg_fourth_central_moment(nu, A, kappa, alpha) / mu_z^2 - 1
}

#' Working shapes of the squared residuals
#'
#' The Gamma shape of the working density of each squared residual,
#' `s_ij = 1 / CV^2(z_ij)`, from the fourth central moment of the
#' single-frailty CNBG cell at `(A, kappa_B)` and the severity-shape
#' convention `alpha` (Proposition 5.2 of the paper). `kappa = Inf` gives the
#' shapes of the absorbed working mean.
#'
#' @param nu Cell means.
#' @param A,kappa Variance-function parameters.
#' @param alpha Severity shape convention.
#' @return A numeric vector of shapes.
#' @examples
#' cnbg_cell_shapes(c(1e4, 1e5, 1e6), A = 6667, kappa = 5, alpha = 3)
#' @export
cnbg_cell_shapes <- function(nu, A, kappa, alpha = 3) 1 / .cnbg_cv2_z(nu, A, kappa, alpha)

## ---- curvature components in phi = (log kappa, log A) ---------------------
## Cell score g = s (z - mu_z) / mu_z^2 d mu_z / d phi, I_W = sum s / mu_z^2 D D',
## V_cell = sum s^2 CV2 / mu_z^2 D D', V_cross = empirical within-row cross terms,
## V = V_cell + V_cross, C = M^{-1} M_A with M'M = I_W, M_A' M_A = I_W V^{-1} I_W.
.cnbg_working_components <- function(z, nu, row, s_cell, kappa_hat, A_hat, alpha) {
  mu_z <- A_hat * nu + if (is.finite(kappa_hat)) nu^2 / kappa_hat else 0
  D    <- cbind(if (is.finite(kappa_hat)) -nu^2 / kappa_hat else 0 * nu, A_hat * nu)
  r    <- (z - mu_z) / mu_z^2
  g    <- (s_cell * r) * D
  I_W  <- crossprod(D * sqrt(s_cell) / mu_z)
  cv2  <- .cnbg_cv2_z(nu, A_hat, kappa_hat, alpha)
  V_cell <- crossprod(D * (s_cell * sqrt(cv2) / mu_z))
  S_row  <- rowsum(g, row)
  V_emp  <- crossprod(S_row)
  V_cross <- V_emp - crossprod(g)
  V <- V_cell + V_cross
  if (!is.finite(kappa_hat)) {
    r1 <- function(X) c(X[2, 2] / I_W[2, 2], NA_real_)
    return(list(I_W = I_W, V_cell = V_cell, V_cross = V_cross, V = V, V_emp = V_emp,
                eig_cell = r1(V_cell), eig_full = r1(V), eig_emp = r1(V_emp),
                trace_ratio = V[2, 2] / I_W[2, 2]))
  }
  list(I_W = I_W, V_cell = V_cell, V_cross = V_cross, V = V, V_emp = V_emp,
       eig_cell = eigen(solve(I_W, V_cell), only.values = TRUE)$values,
       eig_full = eigen(solve(I_W, V), only.values = TRUE)$values,
       eig_emp  = eigen(solve(I_W, V_emp), only.values = TRUE)$values,
       trace_ratio = sum(diag(solve(I_W, V))) / 2)
}

.cnbg_curvature_C <- function(I_W, V) {
  sym <- function(X) (X + t(X)) / 2
  V_use <- V
  ev <- eigen(sym(V), symmetric = TRUE)
  if (any(ev$values <= 0)) {
    ev$values <- pmax(ev$values, 1e-8 * max(ev$values))
    V_use <- ev$vectors %*% diag(ev$values) %*% t(ev$vectors)
  }
  M   <- chol(sym(I_W))
  M_A <- chol(sym(I_W %*% solve(V_use, I_W)))
  list(C = solve(M) %*% M_A, V_used = V_use)
}

## ---- the working log-likelihood and the log prior ---------------------------
## sd: the data list of the Stan model (nu, z, s_cell, kappa_a, kappa_b, A_a, A_b,
## kappa_lb, adjust, phi_hat, C, absorbed). phi = rbind(log kappa, log A), 2 x G.
.cnbg_working_loglik <- function(sd, phi) {
  if (sd$adjust == 1L) phi <- sd$phi_hat + sd$C %*% (phi - sd$phi_hat)
  A_v <- exp(phi[2, ]); k_v <- exp(phi[1, ])
  mu_z <- outer(sd$nu, A_v) + if (sd$absorbed == 1L) 0 else outer(sd$nu^2, 1 / k_v)
  ll <- dgamma(sd$z, shape = sd$s_cell, rate = sd$s_cell / mu_z, log = TRUE)
  colSums(matrix(ll, nrow = length(sd$z)))
}
.cnbg_working_logprior <- function(sd, phi) {
  dgamma(exp(phi[1, ]), sd$kappa_a, sd$kappa_b, log = TRUE) +
    dgamma(exp(phi[2, ]), sd$A_a, sd$A_b, log = TRUE)
}

## ---- engine "grid": deterministic quadrature of the same posterior -------------
## The posterior density in phi-space is prior x likelihood x kappa x A (Jacobian).
## Draws are grid cells sampled by mass, jittered uniformly within the cell.
.cnbg_grid_posterior <- function(sd, n_k = 200L, n_A = 240L, n_draw = 4000L, seed = 1L,
                                 kappa_ub = 2000, A_span = 30) {
  absorbed <- sd$absorbed == 1L
  A0 <- sd$A_a / sd$A_b
  lA <- log(A0) + seq(-log(A_span), log(A_span), length.out = n_A)
  lk <- if (absorbed) log(sd$kappa_a / sd$kappa_b) else seq(log(sd$kappa_lb), log(kappa_ub), length.out = n_k)
  nk <- length(lk)
  lp <- matrix(-Inf, nk, n_A)
  for (a in seq_len(n_A)) {
    phi <- rbind(lk, rep(lA[a], nk))
    ll <- .cnbg_working_loglik(sd, phi)
    lp[, a] <- ll + .cnbg_working_logprior(sd, phi) + if (absorbed) lA[a] else lk + lA[a]
  }
  lp[!is.finite(lp)] <- -Inf
  mx <- max(lp); if (!is.finite(mx)) return(NULL)
  P <- exp(lp - mx); P <- P / sum(P)
  pk <- rowSums(P); pA <- colSums(P)
  wq <- function(grid, p, probs) {
    cdf <- cumsum(p); if (length(grid) == 1L) return(rep(grid, length(probs)))
    sapply(probs, function(q) { i <- which(cdf >= q)[1]; if (is.na(i)) grid[length(grid)] else grid[i] })
  }
  set.seed(seed)
  cell <- sample.int(nk * n_A, n_draw, replace = TRUE, prob = as.numeric(P))
  ki <- ((cell - 1L) %% nk) + 1L; ai <- ((cell - 1L) %/% nk) + 1L
  hA <- (lA[2] - lA[1]) / 2; hk <- if (nk > 1L) (lk[2] - lk[1]) / 2 else 0
  A_draw <- exp(lA[ai] + runif(n_draw, -hA, hA))
  k_draw <- if (absorbed) rgamma(n_draw, sd$kappa_a, sd$kappa_b) else exp(lk[ki] + runif(n_draw, -hk, hk))
  q3 <- c(0.025, 0.5, 0.975)
  k_q <- if (absorbed) qgamma(q3, sd$kappa_a, sd$kappa_b) else exp(wq(lk, pk, q3))
  k_mean <- if (absorbed) sd$kappa_a / sd$kappa_b else sum(pk * exp(lk))
  list(kappa = unname(k_q), A = unname(exp(wq(lA, pA, q3))),
       kappa_mean = k_mean, A_mean = sum(pA * exp(lA)),
       diag = c(rhat_max = NA_real_, bulk_ess_min = NA_real_, tail_ess_min = NA_real_, n_div = 0, n_treedepth = 0),
       post = data.frame(kappa = k_draw, A = A_draw),
       grid_argmax = { ij <- which(lp == mx, arr.ind = TRUE)[1, ]; c(kappa = exp(lk[ij[1]]), A = exp(lA[ij[2]])) },
       edge_mass = if (absorbed) 0 else sum(pk[exp(lk) >= 0.95 * kappa_ub]))
}

## posterior mode in the constrained space (prior x likelihood, no Jacobian), the
## quantity rstan::optimizing() returns; refined from a start by L-BFGS-B
.cnbg_grid_mode <- function(sd, start) {
  absorbed <- sd$absorbed == 1L
  f <- function(p) {
    phi <- if (absorbed) rbind(log(start[["kappa"]]), p) else matrix(p, 2, 1)
    v <- .cnbg_working_loglik(sd, phi) + .cnbg_working_logprior(sd, phi)
    if (is.finite(v)) -v else 1e100
  }
  p0 <- if (absorbed) log(start[["A"]]) else log(c(start[["kappa"]], start[["A"]]))
  lower <- if (absorbed) -Inf else c(log(sd$kappa_lb), -Inf)
  o <- tryCatch(optim(p0, f, method = "L-BFGS-B", lower = lower), error = function(e) NULL)
  if (is.null(o)) return(NULL)
  if (absorbed) c(kappa = unname(start[["kappa"]]), A = unname(exp(o$par)))
  else c(kappa = unname(exp(o$par[1])), A = unname(exp(o$par[2])))
}

## ---- engine "nuts": the paper's sampler (rstan, Suggests) ----------------------
.cnbg_stan_env <- new.env(parent = emptyenv())

.cnbg_stan_model <- function(name) {
  key <- paste0("mod_", name)
  if (!is.null(.cnbg_stan_env[[key]])) return(.cnbg_stan_env[[key]])
  if (!requireNamespace("rstan", quietly = TRUE))
    stop("this engine requires the 'rstan' package; install.packages('rstan')")
  f <- system.file("stan", paste0(name, ".stan"), package = "hgr")
  if (!nzchar(f)) stop("Stan model ", name, " not found in the installed package")
  .cnbg_stan_env[[key]] <- rstan::stan_model(file = f, model_name = name)
  .cnbg_stan_env[[key]]
}

.cnbg_monitor_2 <- function(fit) {
  m0 <- rstan::monitor(as.array(fit), warmup = 0, print = FALSE)
  m  <- as.data.frame(unclass(m0))[match(c("kappa", "A"), rownames(m0)), ]
  sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  nd <- sum(sapply(sp, function(p) sum(p[, "divergent__"])))
  td <- fit@stan_args[[1]]$control$max_treedepth; if (is.null(td)) td <- 10
  nt <- sum(sapply(sp, function(p) sum(p[, "treedepth__"] >= td)))
  c(rhat_max = max(m$Rhat), bulk_ess_min = min(m$Bulk_ESS), tail_ess_min = min(m$Tail_ESS), n_div = nd, n_treedepth = nt)
}

.cnbg_sample_working <- function(sd, seed, chains = 2, iter = 1000, warmup = 500) {
  fit <- tryCatch(rstan::sampling(.cnbg_stan_model("cnbg_working_adj"), data = sd, chains = chains, iter = iter,
                                  warmup = warmup, refresh = 0, seed = seed, show_messages = FALSE),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  post <- tryCatch(as.data.frame(fit), error = function(e) NULL)
  if (is.null(post) || nrow(post) == 0) return(NULL)
  q <- function(v) unname(quantile(v, c(0.025, 0.5, 0.975)))
  list(kappa = q(post$kappa), A = q(post$A), kappa_mean = mean(post$kappa), A_mean = mean(post$A),
       diag = .cnbg_monitor_2(fit), post = post[, c("kappa", "A")], stanfit = fit)
}

.cnbg_nuts_mode <- function(sd, start, seed) {
  mode <- tryCatch(rstan::optimizing(.cnbg_stan_model("cnbg_working_adj"), data = sd, seed = seed, hessian = FALSE,
                                     init = list(kappa = start[["kappa"]], A = start[["A"]]), verbose = FALSE)$par,
                   error = function(e) NULL)
  if (is.null(mode)) NULL else c(kappa = unname(mode["kappa"]), A = unname(mode["A"]))
}

#' Working-likelihood posterior of (kappa_B, A)
#'
#' The composite-likelihood Bayes fit of the paper's Section 5: each squared
#' Stage-1 residual `z_ij = (x_ij - nu_hat_ij)^2 / (1 - h_ij)` (leverage
#' corrected) is given the working density `Gamma(s_ij, s_ij / mu_z_ij)` with
#' `mu_z = A nu_hat + nu_hat^2 / kappa_B` (`mean = "full"`) or `mu_z = A
#' nu_hat` (`mean = "absorbed"`, the free-level fit in which the row levels
#' have absorbed the frailty and `kappa_B` stays at its prior), priors
#' `kappa_B ~ Gamma(prior)` and `A ~ Gamma(2, 2 / A_init)` with `A_init =
#' mean(z) / mean(nu_hat)`. Three passes: (1) shapes `s_ij` from the fourth
#' moment at the posterior medians of a first fit with `s = 1/2`
#' (`shapes = "cell"`; `"half"` and `"legacy2"` fix `s = 1/2` or `s = 2`);
#' (2) the unadjusted posterior and its mode; (3) the posterior under the
#' curvature adjustment of Ribatet, Cooley and Davison (2012), evaluated on
#' `phi = (log kappa_B, log A)` as `phi* = phi_hat + C (phi - phi_hat)` with
#' `C` from the Godambe matrices at the mode (one-dimensional in `log A` for
#' the absorbed mean).
#'
#' `engine = "nuts"` is the sampler the paper ran (two chains of 1,000
#' iterations, 500 warm-up; requires `rstan`). `engine = "grid"` evaluates
#' the same posterior on a deterministic grid in `phi` and needs no Stan; it
#' agrees with the sampler up to Monte Carlo error and is the default here
#' because it is fast and reproducible.
#'
#' Cells with leverage `h >= 1 - 1e-8` and cells whose fitted mean is
#' numerically zero are dropped (`drop_zero_fit`), as in the paper.
#'
#' @param ms A `cnbg_stage1` object.
#' @param alpha Severity-shape convention for the fourth moment (`3` in the
#'   paper; amounts do not identify it).
#' @param mean `"absorbed"` or `"full"`.
#' @param prior Prior on `kappa_B`: `"ref"`, `"empb"` or `c(shape, rate)`.
#' @param engine `"grid"` or `"nuts"`.
#' @param shapes `"cell"`, `"half"` or `"legacy2"`.
#' @param adjust Apply the curvature adjustment (pass 3).
#' @param seed Seed for the sampler, or for the jitter of the grid draws.
#' @param chains,iter,warmup Sampler settings (`engine = "nuts"`).
#' @param kappa_lb Lower bound of `kappa_B` in the working posterior.
#' @param n_draw Number of posterior draws returned.
#' @param leverage_adjust,h_cap Divide the squared residuals by
#'   `1 - min(h, h_cap)`.
#' @param drop_zero_fit,zero_fit_tol Drop cells with `nu_hat < zero_fit_tol *
#'   max(nu_hat)`.
#' @param grid_n Grid sizes `c(n_kappa, n_A)` for `engine = "grid"`.
#' @return An object of class `cnbg_working` with the cell shapes
#'   (`s_cell`), the mode, the curvature eigenvalues (`eig_full`), the
#'   unadjusted (`unadj`) and adjusted (`adj`) posterior summaries (2.5, 50
#'   and 97.5 percent quantiles of `kappa` and `A`, means, sampler
#'   diagnostics), the adjusted posterior draws (`post`, a data frame with
#'   `kappa` and `A`), `C`, `n_cells`, `engine`, `mean`, `prior`, and
#'   `failed` (`NULL` on success).
#' @examples
#' set.seed(1)
#' ms <- cnbg_stage1(simulate_cnbg(kappa_r = 5)$triangle)
#' f <- fit_cnbg_working(ms, mean = "absorbed")
#' f$adj$A
#' @export
fit_cnbg_working <- function(ms, alpha = 3, mean = c("absorbed", "full"), prior = "ref",
                             engine = c("grid", "nuts"), shapes = c("cell", "half", "legacy2"),
                             adjust = TRUE, seed = 1L, chains = 2, iter = 1000, warmup = 500,
                             kappa_lb = 0.1, n_draw = 4000L, leverage_adjust = TRUE, h_cap = 0.95,
                             drop_zero_fit = TRUE, zero_fit_tol = 1e-6, grid_n = c(200L, 240L)) {
  mean <- match.arg(mean); engine <- match.arg(engine); shapes <- match.arg(shapes)
  absorbed <- mean == "absorbed"
  pr <- .cnbg_prior(prior)
  sq <- (ms$dat$x - ms$nu_hat)^2
  if (leverage_adjust) sq <- sq / (1 - pmin(ms$hat, h_cap))
  keep <- is.finite(ms$hat) & ms$hat <= 1 - 1e-8 & is.finite(sq)
  n_zero_fit <- 0L
  if (drop_zero_fit) {
    zf <- ms$nu_hat < zero_fit_tol * max(ms$nu_hat)
    n_zero_fit <- sum(keep & zf); keep <- keep & !zf
  }
  out <- list(n_cells = sum(keep), n_zero_fit_dropped = n_zero_fit, shapes = shapes, alpha = alpha,
              mean = mean, engine = engine, prior = pr, failed = NULL)
  if (sum(keep) < 10L) { out$failed <- "fewer than 10 usable cells"; return(structure(out, class = "cnbg_working")) }
  nu <- as.numeric(ms$nu_hat[keep]); z <- as.numeric(sq[keep]); row <- ms$dat$i[keep]
  A_init <- max(mean(z) / mean(nu), 1e-3)
  base <- list(n_obs = length(z), nu = nu, z = z, kappa_a = pr[["shape"]], kappa_b = pr[["rate"]],
               A_a = 2.0, A_b = 2.0 / A_init, kappa_lb = kappa_lb, adjust = 0L,
               phi_hat = c(0, 0), C = diag(2), absorbed = as.integer(absorbed))
  out$A_init <- A_init
  sampler <- function(sd, sd_seed) {
    if (engine == "nuts") .cnbg_sample_working(sd, sd_seed, chains, iter, warmup)
    else .cnbg_grid_posterior(sd, grid_n[1], grid_n[2], n_draw, sd_seed)
  }
  ## ---- pass 1: shapes ----
  if (shapes == "legacy2") s_cell <- rep(2, length(z))
  else if (shapes == "half") s_cell <- rep(0.5, length(z))
  else {
    p1 <- sampler(c(base, list(s_cell = rep(0.5, length(z)))), seed)
    if (is.null(p1)) { s_cell <- rep(0.5, length(z)); out$pass1_failed <- TRUE }
    else {
      s_cell <- cnbg_cell_shapes(nu, p1$A[2], if (absorbed) Inf else p1$kappa[2], alpha)
      out$pass1 <- p1[c("kappa", "A", "diag")]; out$pass1_failed <- FALSE
    }
    s_cell <- pmin(pmax(s_cell, 1e-4), 0.5)       # CV2 >= 2 by eq:cv_lower, so s <= 0.5
  }
  out$s_cell <- s_cell
  sd <- c(base, list(s_cell = s_cell))
  ## ---- pass 2: unadjusted ----
  u <- sampler(sd, seed + 1L)
  if (is.null(u)) { out$failed <- "unadjusted"; return(structure(out, class = "cnbg_working")) }
  out$unadj <- u[c("kappa", "A", "kappa_mean", "A_mean", "diag")]
  out$unadj_post <- u$post
  start <- c(kappa = u$kappa[2], A = u$A[2])
  mode <- if (engine == "nuts") .cnbg_nuts_mode(sd, start, seed)
          else .cnbg_grid_mode(sd, c(kappa = unname(u$grid_argmax["kappa"]), A = unname(u$grid_argmax["A"])))
  ## boundary guard: a mode at A ~ 0 or at the kappa lower bound is a degenerate
  ## optimiser result; fall back to the posterior median and flag it
  mode_ok <- !is.null(mode) && all(is.finite(mode[c("kappa", "A")])) &&
             mode["A"] > 1e-3 * u$A[2] && (absorbed || (mode["kappa"] > kappa_lb * 1.01 && mode["kappa"] < 100 * u$kappa[3]))
  kh <- if (mode_ok) unname(mode["kappa"]) else u$kappa[2]
  Ah <- if (mode_ok) unname(mode["A"])     else u$A[2]
  if (absorbed) kh <- Inf
  out$mode <- c(kappa = kh, A = Ah); out$mode_from_median <- !mode_ok
  ## ---- curvature components at the mode ----
  wc <- .cnbg_working_components(z, nu, row, s_cell, kh, Ah, alpha)
  out$eig_full <- wc$eig_full; out$eig_cell <- wc$eig_cell; out$eig_emp <- wc$eig_emp; out$trace_ratio <- wc$trace_ratio
  if (!adjust) {
    out$adj <- out$unadj; out$post <- u$post; out$C <- diag(2); out$adjusted <- FALSE
    return(structure(out, class = "cnbg_working"))
  }
  ## ---- pass 3: adjusted ----
  if (absorbed) {
    cA <- sqrt(wc$I_W[2, 2] / wc$V[2, 2]); cc <- list(C = diag(c(1, cA)))
    sd_adj <- sd; sd_adj$adjust <- 1L; sd_adj$phi_hat <- c(0, log(Ah)); sd_adj$C <- cc$C
  } else {
    cc <- .cnbg_curvature_C(wc$I_W, wc$V)
    sd_adj <- sd; sd_adj$adjust <- 1L; sd_adj$phi_hat <- c(log(kh), log(Ah)); sd_adj$C <- cc$C
  }
  a <- sampler(sd_adj, seed + 2L)
  if (is.null(a)) { out$failed <- "adjusted"; return(structure(out, class = "cnbg_working")) }
  out$adj <- a[c("kappa", "A", "kappa_mean", "A_mean", "diag")]
  out$post <- a$post; out$C <- cc$C; out$adjusted <- TRUE
  if (!is.null(a$edge_mass)) out$edge_mass <- a$edge_mass
  structure(out, class = "cnbg_working")
}

#' @export
print.cnbg_working <- function(x, ...) {
  cat(sprintf("CNBG working-likelihood posterior (%s mean, %s engine, prior Gamma(%g, %g))\n",
              x$mean, x$engine, x$prior[["shape"]], x$prior[["rate"]]))
  if (!is.null(x$failed)) { cat("  FAILED:", x$failed, "\n"); return(invisible(x)) }
  cat(sprintf("  cells %d (shapes median %.3f)   curvature eigenvalue %.2f   adjusted: %s\n",
              x$n_cells, median(x$s_cell), x$eig_full[1], if (isTRUE(x$adjusted)) "yes" else "no"))
  cat(sprintf("  A:       median %.4g   95%% [%.4g, %.4g]\n", x$adj$A[2], x$adj$A[1], x$adj$A[3]))
  if (x$mean == "full")
    cat(sprintf("  kappa_B: median %.2f   95%% [%.2f, %.2f]   mean %.2f\n",
                x$adj$kappa[2], x$adj$kappa[1], x$adj$kappa[3], x$adj$kappa_mean))
  else cat("  kappa_B: not in the likelihood (absorbed by the free levels); taken from the prior\n")
  if (x$engine == "nuts")
    cat(sprintf("  sampler: Rhat %.3f, bulk ESS %.0f, divergences %d\n",
                x$adj$diag["rhat_max"], x$adj$diag["bulk_ess_min"], x$adj$diag["n_div"]))
  invisible(x)
}

## usable posterior in the sense of the paper (sampler converged)
.cnbg_working_usable <- function(f) {
  if (is.null(f) || !is.null(f$failed)) return(FALSE)
  if (f$engine == "grid") return(TRUE)
  d <- f$adj$diag
  isTRUE(d["rhat_max"] < 1.05 && d["bulk_ess_min"] > 100)
}
