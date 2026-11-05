## ===========================================================================
## ODP BOOTSTRAP DECOMPOSITION
## A port of the library behind "The over-dispersed Poisson bootstrap: A
## decomposition of predictive variance" (Van Oirbeek, 2026, working paper):
## 24_correspondence_lib.R (chain-ladder fit, Algorithm 1 with its variants,
## the first-order level and pattern terms, the Gibbs sampler),
## 51_posterior_exact_lib.R (the closed-form quasi-posterior of Proposition
## 7.1) and the certificate of 67_process_convention.R of its replication
## material. With the same seed, odp_bootstrap() draws exactly what the
## paper's scripts drew.
## ===========================================================================

.odp_check <- function(triangle) {
  if (inherits(triangle, "odp_fit")) return(triangle$triangle)
  tri <- as.matrix(triangle); storage.mode(tri) <- "double"
  if (nrow(tri) < 2 || ncol(tri) < 2)
    stop("the triangle needs at least two accident years and two development years")
  obs <- !is.na(tri); k <- rowSums(obs); J <- ncol(tri)
  for (i in seq_len(nrow(tri)))
    if (k[i] < 1 || any(obs[i, ] != (seq_len(J) <= k[i])))
      stop("accident year ", i, " must be observed in a contiguous block of development years that starts at the first")
  if (any(diff(k) > 0)) stop("a later accident year cannot be observed further than an earlier one")
  tri
}

#' Chain-ladder fit of an incremental triangle as the ODP quasi-likelihood fit
#'
#' The (quasi-)maximum likelihood fit of the over-dispersed Poisson model
#' `X_ij = phi Y_ij`, `Y_ij ~ Poisson(mu_i pi_j / phi)`, which reproduces the
#' chain-ladder (Section 2.1 of the paper): development factors, cumulative
#' shares `F`, pattern `pi`, levels `mu_i = C_i / F_i`, fitted increments, the
#' Pearson residuals on the cells with a positive fitted increment, the
#' dispersion `phi = sum r^2 / (n - p)` and the adjusted residual pool
#' `r_ij {n / (n - p)}^{1/2}` of line 1 of Algorithm 1, together with the
#' residual-pool factor `sigma_r^2 / phi` of Remark 4.5.
#'
#' Columns whose chain-ladder fitted increment is exactly zero (an observed
#' column total of zero) are structural zeros: their parameter is removed,
#' their cells leave the residual pool and their pseudo-cells stay at zero,
#' as in Section 7 of the paper. `dof = "full"` keeps `p = I + J - 1` in the
#' degrees-of-freedom correction, which is what the paper's scripts do;
#' `dof = "reduced"` counts only the retained columns.
#'
#' @param triangle Incremental triangle (matrix), accident years in rows,
#'   development years in columns, unobserved cells `NA`.
#' @param dof `"full"` (the paper's scripts: `p = I + J - 1`) or
#'   `"reduced"` (`p = I + J' - 1` with `J'` the retained columns).
#' @return An object of class `odp_fit`: `f`, `F`, `pi`, `mu`, `C` (observed
#'   row totals), `lat` (last observed column of each row), `m` (fitted
#'   increments), `zero_columns`, `n`, `p`, `phi`, `residuals` (raw Pearson),
#'   `pool` (adjusted residuals), `pool_factor` (`sigma_r^2 / phi`) and
#'   `reserve` (the chain-ladder reserve `sum_i mu_i (1 - F_i)`).
#' @details The fit stops on a triangle with a development factor below one
#'   or a non-positive observed row total, where the paper's bootstrap and
#'   quasi-posterior are not defined (Sections 2.1 and 7).
#' @references Van Oirbeek, R. (2026). The over-dispersed Poisson bootstrap:
#'   A decomposition of predictive variance. Working paper.
#' @export
odp_fit <- function(triangle, dof = c("full", "reduced")) {
  dof <- match.arg(dof); tri <- .odp_check(triangle)
  I <- nrow(tri); J <- ncol(tri); obs <- !is.na(tri)
  cum <- t(apply(tri, 1, function(r) { s <- cumsum(ifelse(is.na(r), 0, r)); s[is.na(r)] <- NA; s }))
  f <- vapply(seq_len(J - 1), function(j) { rows <- which(obs[, j + 1])
    sum(cum[rows, j + 1]) / sum(cum[rows, j]) }, 0)
  if (any(!is.finite(f)))
    stop("a cumulative column total entering a development factor is zero: the chain-ladder fit is undefined")
  Fv <- c(rev(1 / cumprod(rev(f))), 1)
  lat <- vapply(seq_len(I), function(i) max(which(obs[i, ])), 1L)
  C <- cum[cbind(seq_len(I), lat)]
  mu <- C / Fv[lat]; pi <- diff(c(0, Fv))
  if (any(pi < 0))
    stop("a development factor below one gives a negative fitted increment; the ODP bootstrap of the paper is not defined on this triangle (Section 2.1)")
  if (any(C <= 0))
    stop("every observed row total must be positive (a zero latest diagonal gives a zero level)")
  m <- outer(mu, pi)
  zero_cols <- which(pi == 0)
  mo <- m[obs]; y <- tri[obs]; pos <- mo > 0
  r <- (y[pos] - mo[pos]) / sqrt(mo[pos])
  n <- sum(pos)
  p <- if (dof == "full") I + J - 1L else I + (J - length(zero_cols)) - 1L
  if (n <= p) stop("fewer cells with a positive fitted increment than parameters: the dispersion is not estimable")
  phi <- sum(r^2) / (n - p)
  radj <- r * sqrt(n / (n - p))
  structure(list(triangle = tri, I = I, J = J, obs = obs, f = f, F = Fv, pi = pi, mu = mu, C = C, lat = lat,
                 m = m, zero_columns = zero_cols, n = n, p = p, dof = dof, phi = phi, residuals = r, pool = radj,
                 pool_factor = mean((radj - mean(radj))^2) / phi, reserve = sum(mu * (1 - Fv[lat]))),
            class = "odp_fit")
}

#' @export
print.odp_fit <- function(x, ...) {
  cat(sprintf("ODP chain-ladder fit: %d accident years x %d development years, %d observed cells\n", x$I, x$J, sum(x$obs)))
  cat(sprintf("  chain-ladder reserve %.4g   phi (Pearson, n = %d, p = %d) %.4g   pool factor sigma_r^2 / phi %.6f\n",
              x$reserve, x$n, x$p, x$phi, x$pool_factor))
  if (length(x$zero_columns)) cat("  structural-zero columns:", paste(x$zero_columns, collapse = ", "), "\n")
  invisible(x)
}

#' First-order level, pattern and process terms of the ODP predictive variance
#'
#' The three terms of equation (5) of the paper at a chain-ladder fit:
#' `Pi = phi sum_i mu_i (1 - F_i)` (process), `L = phi sum_i (1 - F_i)^2
#' mu_i / F_i` (level, the inverse-information quadratic form of the levels)
#' and `P = h' S^{-1} h` (pattern, with `S` the Schur complement of the
#' pattern block of the expected Fisher information on the tangent space of
#' the simplex and `h` the projected gradient). To first order in the volume
#' the ODP bootstrap variance is `Pi + (sigma_r^2 / phi)(L + P)` and the
#' log-flat quasi-posterior predictive variance is `Pi + L + P`
#' (Theorem 4.1). Structural-zero columns carry no parameter.
#'
#' @param fit An `odp_fit` (or an incremental triangle, fitted first).
#' @param phi Dispersion at which the terms are evaluated (default: the
#'   fit's Pearson estimate).
#' @return A named vector `c(level, pattern, process)` with attributes
#'   `shares` (each term over their sum) and `phi`.
#' @export
odp_variance_terms <- function(fit, phi = NULL) {
  if (!inherits(fit, "odp_fit")) fit <- odp_fit(fit)
  if (is.null(phi)) phi <- fit$phi
  mu <- fit$mu; pi <- fit$pi; obs <- fit$obs
  keep <- pi > 0; pi <- pi[keep]; obs <- obs[, keep, drop = FALSE]
  I <- length(mu); J <- length(pi)
  Fi <- vapply(seq_len(I), function(i) sum(pi[obs[i, ]]), 0)
  d_mu <- 1 - Fi
  d_pi <- vapply(seq_len(J), function(j) sum(mu[!obs[, j]]), 0)
  Imp <- obs * (1 / phi)
  Ipp <- diag(vapply(seq_len(J), function(j) sum(mu[obs[, j]]) / pi[j], 0) / phi, J)
  Q <- rbind(diag(J - 1), -1)                       # basis of the tangent space sum(d pi) = 0
  Iee <- t(Q) %*% Ipp %*% Q; Ime <- Imp %*% Q; de <- t(Q) %*% d_pi
  Imm_inv <- diag(phi * mu / Fi, I)
  S <- Iee - t(Ime) %*% Imm_inv %*% Ime
  h <- de - t(Ime) %*% Imm_inv %*% d_mu
  L <- sum(d_mu^2 * phi * mu / Fi)
  P <- as.numeric(t(h) %*% solve(S, h))
  proc <- phi * sum(mu * (1 - Fi))
  out <- c(level = L, pattern = P, process = proc)
  attr(out, "shares") <- out / sum(out); attr(out, "phi") <- phi
  out
}

#' The ODP bootstrap of England and Verrall with its variants and its decomposition
#'
#' Algorithm 1 of the paper, vectorised over replications: Pearson residuals
#' and the adjusted pool of line 1, pseudo-cells by resampling the pool with
#' replacement (line 3), the chain-ladder refit on each pseudo-triangle
#' (line 4) and the process draw (line 5), one Gamma variable with shape
#' `g* / phi` and scale `phi` on the projected reserve mean `g*`, or a zero
#' reserve when `g*` is not positive or a pseudo column total is zero. The
#' three variants that isolate the operations are returned on the same
#' pseudo-triangles: the no-refit bootstrap (pseudo row totals projected
#' with the original factors), the diagonal-fixed bootstrap (observed row
#' totals projected with the refitted factors) and the process-only draw.
#'
#' The decomposition of Theorem 4.1 is reported as ratios: the variance of
#' the no-refit parameter part over `(sigma_r^2 / phi) L` (one exactly, up
#' to Monte Carlo error), the variance of the refit increment over
#' `(sigma_r^2 / phi) P`, the cross term as a share of the parameter
#' variance, the diagonal-fixed variance over `(sigma_r^2 / phi) P`, and
#' the bootstrap predictive variance over `Pi + (sigma_r^2 / phi)(L + P)`.
#'
#' `pool = "conventional"` is the England and Verrall pool, whose mean need
#' not vanish under heterogeneous levels (Remark 4.5); `pool = "centred"` is
#' the centred and rescaled pool the paper recommends, with mean zero and
#' variance `phi`, so that the factor `sigma_r^2 / phi` is one.
#'
#' @param triangle Incremental triangle or an `odp_fit`.
#' @param B Number of replications (the paper uses 4,000 and 20,000).
#' @param pool `"conventional"` or `"centred"`.
#' @param dof Degrees-of-freedom convention of [odp_fit()].
#' @param level Level of the central predictive intervals.
#' @param seed Optional seed; with the same seed the draws are those of the
#'   paper's scripts.
#' @return An object of class `odp_boot`: the `fit`, the pool used and its
#'   factor, the parameter parts `g_full`, `g_nore`, `g_cdiag` (B each), the
#'   predictive draws `R_full`, `R_nore`, `R_cdiag`, `R_cond`, the terms of
#'   [odp_variance_terms()], the `decomposition` ratios, the `intervals` of
#'   the four predictives, and the share of replications on which the zero
#'   fallback of line 5 fired.
#' @examples
#' tri <- rbind(c(1, 1, 1), c(10, 1, NA), c(1, NA, NA))
#' b <- odp_bootstrap(tri, B = 2000, seed = 1)
#' b$decomposition
#' @export
odp_bootstrap <- function(triangle, B = 4000L, pool = c("conventional", "centred"), dof = c("full", "reduced"),
                          level = 0.95, seed = NULL) {
  pool <- match.arg(pool); dof <- match.arg(dof)
  fit <- if (inherits(triangle, "odp_fit")) triangle else odp_fit(triangle, dof)
  if (!is.null(seed)) set.seed(seed)
  B <- as.integer(B); I <- fit$I; J <- fit$J; obs <- fit$obs; phi <- fit$phi
  m <- fit$m[obs]                                     # column-major order of the observed cells
  radj <- if (pool == "conventional") fit$pool else {
    r <- fit$residuals; (r - mean(r)) * sqrt(fit$n * phi / sum((r - mean(r))^2)) }
  factor <- mean((radj - mean(radj))^2) / phi
  ridx <- which(obs, arr.ind = TRUE)
  Ystar <- matrix(rep(m, each = B), B) +
    matrix(sample(radj, B * length(m), replace = TRUE), B) * matrix(rep(sqrt(pmax(m, 0)), each = B), B)
  cum <- array(0, c(B, I, J))
  for (k in seq_len(nrow(ridx))) { i <- ridx[k, 1]; j <- ridx[k, 2]
    cum[, i, j] <- Ystar[, k] + (if (j > 1) cum[, i, j - 1] else 0) }
  fstar <- matrix(NA_real_, B, J - 1)
  for (j in seq_len(J - 1)) { rows <- which(obs[, j + 1])
    fstar[, j] <- rowSums(cum[, rows, j + 1, drop = FALSE]) / rowSums(cum[, rows, j, drop = FALSE]) }
  Fs <- apply(fstar, 1, function(ff) rev(1 / cumprod(rev(ff))))
  if (!is.matrix(Fs)) Fs <- matrix(Fs, nrow = 1)
  Fstar <- cbind(t(Fs), 1)
  Cstar <- sapply(seq_len(I), function(i) cum[, i, fit$lat[i]]); if (!is.matrix(Cstar)) Cstar <- matrix(Cstar, nrow = B)
  open <- which(fit$lat < J)
  Fs_lat <- sapply(open, function(i) Fstar[, fit$lat[i]]); if (!is.matrix(Fs_lat)) Fs_lat <- matrix(Fs_lat, nrow = B)
  Fh_lat <- fit$F[fit$lat[open]]
  g_full  <- rowSums(Cstar[, open, drop = FALSE] * (1 / Fs_lat - 1))
  g_nore  <- as.vector(Cstar[, open, drop = FALSE] %*% (1 / Fh_lat - 1))
  g_cdiag <- as.vector((1 / Fs_lat - 1) %*% fit$C[open])
  g_hat   <- sum(fit$C[open] * (1 / Fh_lat - 1))
  procdraw <- function(g) { g[!is.finite(g)] <- 0
    ifelse(g > 0, stats::rgamma(B, shape = pmax(g, 1e-12) / phi, scale = phi), 0) }
  R_full <- procdraw(g_full); R_nore <- procdraw(g_nore); R_cdiag <- procdraw(g_cdiag); R_cond <- procdraw(rep(g_hat, B))
  terms <- odp_variance_terms(fit, phi)
  d <- g_full - g_nore
  decomposition <- c(level = stats::var(g_nore) / (factor * terms[["level"]]),
                     pattern = stats::var(d) / (factor * terms[["pattern"]]),
                     cross = 2 * stats::cov(g_nore, d) / stats::var(g_full),
                     diagonal_fixed = stats::var(g_cdiag) / (factor * terms[["pattern"]]),
                     predictive = stats::var(R_full) / (terms[["process"]] + factor * (terms[["level"]] + terms[["pattern"]])))
  a <- (1 - level) / 2
  q <- function(x) unname(stats::quantile(x, c(a, 1 - a)))
  intervals <- rbind(full = q(R_full), no_refit = q(R_nore), diagonal_fixed = q(R_cdiag), process_only = q(R_cond))
  colnames(intervals) <- c("lower", "upper")
  structure(list(fit = fit, B = B, pool_type = pool, pool = radj, pool_factor = factor, phi = phi, reserve = g_hat,
                 terms = terms, g_full = g_full, g_nore = g_nore, g_cdiag = g_cdiag,
                 R_full = R_full, R_nore = R_nore, R_cdiag = R_cdiag, R_cond = R_cond,
                 decomposition = decomposition, intervals = intervals, level = level,
                 fallback_share = mean(!(is.finite(g_full) & g_full > 0)),
                 mean = mean(R_full), sd = stats::sd(R_full)), class = "odp_boot")
}

#' @export
print.odp_boot <- function(x, ...) {
  cat(sprintf("ODP bootstrap (Algorithm 1), B = %d, %s pool, factor sigma_r^2 / phi = %.6f, phi = %.4g\n",
              x$B, x$pool_type, x$pool_factor, x$phi))
  cat(sprintf("  chain-ladder reserve %.4g   bootstrap mean %.4g   sd %.4g\n", x$reserve, x$mean, x$sd))
  cat(sprintf("  %.0f%% intervals:\n", 100 * x$level)); print(round(x$intervals, 2))
  sh <- attr(x$terms, "shares")
  cat(sprintf("  first-order shares: process %.2f  level %.2f  pattern %.2f\n", sh[["process"]], sh[["level"]], sh[["pattern"]]))
  cat("  decomposition ratios (Theorem 4.1): "); print(round(x$decomposition, 4))
  if (x$fallback_share > 0) cat(sprintf("  zero-reserve fallback fired on %.2f%% of the replications\n", 100 * x$fallback_share))
  invisible(x)
}

.odp_inv_moment <- function(s, k, m) {               # E[F_k^{-m}] under the Beta factorisation
  l <- k:(s$Jp - 1); if (any(s$e[l] <= m)) return(Inf)
  exp(sum(lbeta(s$e[l] - m, s$cj[l + 1]) - lbeta(s$e[l], s$cj[l + 1])))
}

.odp_exact_moments <- function(s) {
  op <- s$open; z <- c(pattern = 0, level = 0, process = 0)
  if (!length(op)) return(list(mean = 0, var = 0, var_g = 0, components = z, finite_mean = TRUE, finite_var = TRUE))
  k <- s$kk[op]; C <- s$C[op]; phi <- s$phi
  m1 <- vapply(k, function(kk) .odp_inv_moment(s, kk, 1), 0); m2 <- vapply(k, function(kk) .odp_inv_moment(s, kk, 2), 0)
  finite_mean <- all(is.finite(m1)); finite_var <- all(is.finite(m2))
  if (!finite_mean) return(list(mean = Inf, var = Inf, var_g = Inf, components = z + Inf, finite_mean = FALSE, finite_var = FALSE))
  mean_R <- sum(C * (m1 - 1))
  if (!finite_var) return(list(mean = mean_R, var = Inf, var_g = Inf, components = z + Inf, finite_mean = TRUE, finite_var = FALSE))
  cross <- function(ka, kb) { if (ka > kb) { t <- ka; ka <- kb; kb <- t }
    l1 <- if (ka < kb) ka:(kb - 1) else integer(0); l2 <- kb:(s$Jp - 1)
    exp(sum(lbeta(s$e[l1] - 1, s$cj[l1 + 1]) - lbeta(s$e[l1], s$cj[l1 + 1])) +
        sum(lbeta(s$e[l2] - 2, s$cj[l2 + 1]) - lbeta(s$e[l2], s$cj[l2 + 1]))) }
  EgF <- outer(seq_along(op), seq_along(op), Vectorize(function(a, b) cross(k[a], k[b])))
  E_g2 <- sum(outer(C, C) * EgF) - 2 * sum(C) * sum(C * m1) + sum(C)^2
  var_Eg <- E_g2 - mean_R^2                            # Var(E[g | pi]): the pattern component
  E_var_lev <- sum(phi * C * (m2 - 2 * m1 + 1))        # E[Var(g | pi)]: the level component
  E_var_proc <- sum(phi * C * (m1 - 1))                # phi E[R]: the process component
  list(mean = mean_R, var = var_Eg + E_var_lev + E_var_proc, var_g = var_Eg + E_var_lev,
       components = c(pattern = var_Eg, level = E_var_lev, process = E_var_proc), finite_mean = TRUE, finite_var = TRUE)
}

.odp_exact_draw <- function(s, n) {
  U <- sapply(seq_len(s$Jp - 1), function(k) stats::rbeta(n, s$e[k], s$cj[k + 1]))
  if (!is.matrix(U)) U <- matrix(U, n)
  Fm <- t(apply(U, 1, function(u) c(rev(cumprod(rev(u))), 1)))
  if (nrow(Fm) != n) Fm <- t(Fm)
  g <- numeric(n); R <- numeric(n)
  for (i in s$open) { Fi <- Fm[, s$kk[i]]; mu <- stats::rgamma(n, shape = s$C[i] / s$phi, rate = Fi / s$phi)
    g <- g + mu * (1 - Fi); R <- R + s$phi * stats::rpois(n, mu * (1 - Fi) / s$phi) }
  list(g = g, R = R)
}

#' The closed-form quasi-posterior of the ODP pattern and the exact reserve moments
#'
#' Proposition 7.1 of the paper. Under the quasi-posterior with priors flat
#' in the log-parameters and the dispersion plugged in, with `c_j` and `r_i`
#' the observed column and row totals divided by `phi` and the
#' structural-zero columns removed, the successive share ratios
#' `U_k = F_k / F_{k+1}` are independent `Beta(e_k, c_{k+1})` with
#' `e_k = sum_{j <= k} c_j - sum_{i: k_i <= k} r_i`, and given the pattern
#' the levels are independent `Gamma(C_i / phi, F_{k_i} / phi)`. The
#' predictive mean of the reserve is finite if and only if `min_{l >= k_i}
#' e_l > 1` and its variance if and only if `min_{l >= k_i} e_l > 2` on every
#' open row; where finite, both are computed exactly, the variance as the sum
#' of its pattern, level and process components. A `Dirichlet(alpha, ...,
#' alpha)` prior on the pattern keeps the factorisation with `c_j + alpha`
#' and `e_k + k alpha` (the paper's Table 3); with `alpha > 0` the
#' structural-zero columns may be retained.
#'
#' @param triangle Incremental triangle or an `odp_fit`.
#' @param phi Dispersion plugged in (default: the fit's Pearson estimate,
#'   as in the paper).
#' @param alpha Dirichlet prior mass per column (`0` is the log-flat prior).
#' @param zero_columns `"remove"` (the paper's target) or `"retain"`
#'   (requires `alpha > 0`).
#' @param n Number of independent draws of the parameter part and of the
#'   predictive to return (`0` for none); the draws need no Markov chain.
#' @param dof Degrees-of-freedom convention of [odp_fit()].
#' @param seed Optional seed for the draws.
#' @return An object of class `odp_posterior`: the exponents `e`, the scaled
#'   column totals `c`, the tail criterion `min_{l >= k_i} e_l` of each open
#'   row, the finiteness flags, the exact `mean`, `var` (predictive),
#'   `var_parameter` (`Var(g | D)`) and `components`, the chain-ladder
#'   reserve, and the `draws` (`g`, `R`) when `n > 0`.
#' @examples
#' tri <- rbind(c(1, 1, 1), c(10, 1, NA), c(1, NA, NA))
#' odp_posterior(tri)              # e = (5.06, 0.92): the mean is infinite
#' @export
odp_posterior <- function(triangle, phi = NULL, alpha = 0, zero_columns = c("remove", "retain"), n = 0L,
                          dof = c("full", "reduced"), seed = NULL) {
  zero_columns <- match.arg(zero_columns); dof <- match.arg(dof)
  fit <- if (inherits(triangle, "odp_fit")) triangle else odp_fit(triangle, dof)
  if (is.null(phi)) phi <- fit$phi
  if (alpha < 0) stop("alpha must be non-negative")
  tri <- fit$triangle; I <- fit$I; J <- fit$J; obs <- fit$obs; k_i <- rowSums(obs)
  ctot <- colSums(tri, na.rm = TRUE); rtot <- rowSums(tri, na.rm = TRUE)
  cols <- if (zero_columns == "remove") which(ctot > 0) else seq_len(J)
  cj <- ctot[cols] / phi + alpha; ri <- rtot / phi
  if (any(cj <= 0)) stop("a retained column has a non-positive total: the quasi-posterior is improper (retain zero columns only with alpha > 0)")
  Jp <- length(cols); if (Jp < 2) stop("fewer than two retained columns")
  kk <- vapply(seq_len(I), function(i) sum(cols <= k_i[i]), 0L)
  e <- vapply(seq_len(Jp - 1), function(k) sum(cj[seq_len(k)]) - sum(ri[kk <= k]), 0)
  s <- list(I = I, Jp = Jp, kk = kk, cj = cj, ri = ri, e = e, C = rtot, phi = phi, open = which(kk < Jp))
  mo <- .odp_exact_moments(s)
  tail <- vapply(s$open, function(i) min(e[kk[i]:(Jp - 1)]), 0); names(tail) <- s$open
  draws <- NULL
  if (n > 0) { if (!is.null(seed)) set.seed(seed); draws <- .odp_exact_draw(s, as.integer(n)) }
  structure(list(phi = phi, alpha = alpha, zero_columns = zero_columns, columns = cols, e = e, c = cj, r = ri,
                 open_rows = s$open, tail_criterion = tail, finite_mean = mo$finite_mean, finite_var = mo$finite_var,
                 mean = mo$mean, var = mo$var, var_parameter = mo$var_g, components = mo$components,
                 reserve = fit$reserve, draws = draws, setup = s), class = "odp_posterior")
}

#' @export
print.odp_posterior <- function(x, ...) {
  cat(sprintf("Log-flat quasi-posterior of the ODP pattern (Proposition 7.1), phi = %.4g, alpha = %g, zero columns %s\n",
              x$phi, x$alpha, x$zero_columns))
  cat("  exponents e_k:", paste(formatC(x$e, digits = 4, format = "g"), collapse = " "), "\n")
  cat(sprintf("  smallest tail exponent on an open row %.3f: mean %s, variance %s\n", min(x$tail_criterion),
              if (x$finite_mean) "finite" else "infinite", if (x$finite_var) "finite" else "infinite"))
  if (x$finite_mean) cat(sprintf("  predictive mean %.4g (chain-ladder reserve %.4g)\n", x$mean, x$reserve))
  if (x$finite_var) cat(sprintf("  predictive variance %.4g = pattern %.4g + level %.4g + process %.4g\n",
                                x$var, x$components[["pattern"]], x$components[["level"]], x$components[["process"]]))
  if (!is.null(x$draws)) cat(sprintf("  %d independent draws attached (draws$g, draws$R)\n", length(x$draws$R)))
  invisible(x)
}

#' Gibbs sampler of the log-flat quasi-posterior of the ODP model
#'
#' The sampler the paper's simulation cells used (Section S7.1 of its
#' Supplementary Material): scaled units `Y = X / phi`, Gamma full
#' conditionals for the levels and for the unnormalised pattern factors,
#' both sets of factors normalised after every sweep, and the predictive
#' `phi` times a Poisson variable with mean `sum_i mu_i (1 - F_i) / phi`.
#' [odp_posterior()] supersedes it with independent draws from the closed
#' form; the sampler is kept for the record.
#'
#' @inheritParams odp_posterior
#' @param n_iter,burn Sweeps and burn-in.
#' @return A list with the draws `g` (parameter part) and `R` (predictive).
#' @export
odp_gibbs <- function(triangle, phi = NULL, n_iter = 4000L, burn = 1000L, dof = c("full", "reduced"), seed = NULL) {
  dof <- match.arg(dof)
  fit <- if (inherits(triangle, "odp_fit")) triangle else odp_fit(triangle, dof)
  if (is.null(phi)) phi <- fit$phi
  if (!is.null(seed)) set.seed(seed)
  tri <- fit$triangle; I <- fit$I; J <- fit$J; obs <- fit$obs
  Y <- ifelse(obs, tri / phi, 0)
  rs <- rowSums(Y); cs <- colSums(Y); zc <- cs == 0
  if (any(rs <= 0) || any(cs < 0)) stop("the sampler needs positive row totals and non-negative column totals")
  mu <- fit$mu / phi; gam <- fit$pi
  keep <- n_iter - burn; g <- numeric(keep); R <- numeric(keep)
  Fi_open <- function(gm) vapply(seq_len(I), function(i) sum(gm[obs[i, ]]), 0)
  for (t in seq_len(n_iter)) {
    mu  <- stats::rgamma(I, shape = rs, rate = as.vector(obs %*% gam))
    gam <- numeric(J); gam[!zc] <- stats::rgamma(sum(!zc), shape = cs[!zc], rate = as.vector(t(obs[, !zc, drop = FALSE]) %*% mu))
    s <- sum(gam); gam <- gam / s; mu <- mu * s
    if (t > burn) { Fi <- Fi_open(gam); m_fut <- sum(mu * (1 - Fi))
      g[t - burn] <- phi * m_fut; R[t - burn] <- phi * stats::rpois(1, m_fut) }
  }
  list(g = g, R = R, phi = phi)
}

#' Deterministic certificates on the bootstrap refit from the observed triangle
#'
#' Every pseudo-cell lies in `[m + r_min m^{1/2}, m + r_max m^{1/2}]` with
#' `r_min` and `r_max` the extreme residuals of the pool, so every pseudo
#' cumulative column total has bounds that depend on the observed triangle
#' only (Section 7 of the paper). The weak certificate holds when the lower
#' bound of every pseudo cumulative column total entering a factor is
#' positive: the refit map then has no pole on the support of the resampling
#' and the bootstrap variance is finite. The strong certificate also requires
#' every pseudo column increment entering a factor and every pseudo row total
#' to be positive: every refitted factor exceeds one, every projected cell is
#' positive, the zero fallback cannot fire and the aggregate and cell-wise
#' process draws coincide. Neither certifies the Monte Carlo estimate of the
#' variance.
#'
#' @param boot An `odp_boot` (its fit and pool are used).
#' @return A list: `no_pole`, `factors_above_one`, the relative margins
#'   `min_den_rel`, `min_inc_rel`, `min_C_rel`, the largest attainable factor
#'   `f_max` and the largest attainable parameter part `g_max` (finite only
#'   under the weak certificate).
#' @export
odp_certificates <- function(boot) {
  if (!inherits(boot, "odp_boot")) stop("boot must be an odp_boot object")
  fit <- boot$fit; obs <- fit$obs; I <- fit$I; J <- fit$J; M <- fit$m
  rmin <- min(boot$pool); rmax <- max(boot$pool)
  lo <- M + rmin * sqrt(pmax(M, 0)); hi <- M + rmax * sqrt(pmax(M, 0)); lo[!obs] <- 0; hi[!obs] <- 0
  fmax <- den_min <- num_min <- numeric(J - 1)
  for (j in seq_len(J - 1)) { rows <- which(obs[, j + 1])
    den_min[j] <- sum(lo[rows, 1:j]); num_min[j] <- sum(lo[rows, 1:(j + 1)]); num_max <- sum(hi[rows, 1:(j + 1)])
    fmax[j] <- if (den_min[j] > 0) num_max / den_min[j] else Inf }
  no_pole <- all(den_min > 0) && all(num_min > 0)
  inc_min <- vapply(seq_len(J - 1), function(j) sum(lo[which(obs[, j + 1]), j + 1]), 0)
  Cmin <- vapply(seq_len(I), function(i) sum(lo[i, 1:fit$lat[i]]), 0)
  strong <- all(den_min > 0) && all(inc_min > 0) && all(Cmin > 0)
  Cmax <- vapply(seq_len(I), function(i) sum(hi[i, 1:fit$lat[i]]), 0); open <- which(fit$lat < J)
  g_max <- if (no_pole) sum(vapply(open, function(i) Cmax[i] * (prod(fmax[fit$lat[i]:(J - 1)]) - 1), 0)) else Inf
  list(no_pole = no_pole, factors_above_one = strong,
       min_den_rel = min(den_min / vapply(seq_len(J - 1), function(j) sum(M[which(obs[, j + 1]), 1:j]), 0)),
       min_inc_rel = min(inc_min / vapply(seq_len(J - 1), function(j) sum(M[which(obs[, j + 1]), j + 1]), 0)),
       min_C_rel = min(Cmin / fit$C), f_max = max(fmax), g_max = g_max)
}
