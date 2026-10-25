## ===========================================================================
## INFORMATION ON kappa_B FROM SQUARED RESIDUALS  (Section 6 of the paper;
## make_tables_split.R of the replication material, tables godambe*.tex)
## ===========================================================================

.cnbg_cov_z <- function(n1, n2, A, kappa)
  n1^2 * n2^2 * (6 / kappa^3 + 2 / kappa^2) + A^2 * n1 * n2 / kappa + 2 * A * n1 * n2 * (n1 + n2) / kappa^2

#' Information on kappa_B from the amount triangle
#'
#' At a calibration with the cell means and `A` known, computes for each
#' true `kappa_B` the working-independence (diagonal) Fisher information of
#' the squared-residual estimating equation, the within-row covariance
#' correction `Delta` from the shared frailty, the Godambe information
#' `G = I_diag^2 / (I_diag + Delta)` (Proposition 6.1 of the paper), the
#' implied coefficient of variation `1 / (kappa sqrt(G))`, the delta-method
#' coefficient of variation of the `nu^-2` weighted least-squares estimator
#' (the paper's Table godambe_wls), and the heterogeneity-limit per-cell
#' information `1 / (2 kappa^2 + 6 kappa)` summed over the observed cells.
#' Single-frailty fourth moment throughout.
#'
#' At the defaults (the paper's calibration) and averaged over `alpha` in
#' `c(1, 3, 10)`, `G` reproduces the `G_avg` column of the paper's Table
#' godambe. The paper states no van Trees bound from `G`, since that would
#' need the Fisher information of the full likelihood rather than the
#' Godambe information of one estimating equation; none is returned.
#'
#' @param kappa_grid True `kappa_B` values.
#' @param alpha Severity shape.
#' @param I,J,Y_bar,mu_base,mu_growth,w Calibration (see [simulate_cnbg()]).
#' @return A data frame: `kappa`, `I_diag`, `Delta`, `G`, `cv_G`, `cv_wls`,
#'   `I_lead`.
#' @examples
#' ## the paper's G_avg column:
#' Gs <- sapply(c(1, 3, 10), function(a) cnbg_information_kappa(alpha = a)$G)
#' round(rowMeans(Gs), 4)   # 0.1644 0.0610 0.0120 0.0017 0.0001
#' @export
cnbg_information_kappa <- function(kappa_grid = c(3, 5, 10, 20, 50), alpha = 3,
                                   I = 10L, J = 10L, Y_bar = 5000,
                                   mu_base = 30, mu_growth = 0.05, w = NULL) {
  if (is.null(w)) w <- .cnbg_paper_w[seq_len(J)]
  stopifnot(length(w) == J)
  mu_i <- mu_base * exp(mu_growth * (seq_len(I) - 1L))
  w    <- w / sum(w)
  nu   <- outer(mu_i, w) * Y_bar
  obs  <- .cnbg_obs_mask(I, J)
  A    <- (alpha + 1) * Y_bar / alpha
  rows <- lapply(kappa_grid, function(kappa) {
    mz <- A * nu + nu^2 / kappa
    cv <- .cnbg_cv2_z(nu, A, kappa, alpha)
    I_diag <- sum(((nu^2 / kappa^2)^2 / (cv * mz^2))[obs])
    Delta <- 0
    for (i in seq_len(I)) {
      js <- which(obs[i, ]); if (length(js) < 2L) next
      for (p in seq_along(js)[-length(js)]) for (q in js[(p + 1L):length(js)]) {
        j1 <- js[p]; n1 <- nu[i, j1]; n2 <- nu[i, q]
        Delta <- Delta + 2 * (n1^2 / kappa^2) * (n2^2 / kappa^2) * .cnbg_cov_z(n1, n2, A, kappa) /
                 (cv[i, j1] * mz[i, j1]^2 * cv[i, q] * mz[i, q]^2)
      }
    }
    G <- I_diag^2 / (I_diag + Delta)
    nuo <- nu[obs]; cvo <- cv[obs]; mzo <- mz[obs]
    var_B_wls <- sum(cvo * mzo^2) / sum(nuo^2)^2        # delta method on kappa = 1/B
    data.frame(kappa = kappa, I_diag = I_diag, Delta = Delta, G = G,
               cv_G = 1 / (kappa * sqrt(G)), cv_wls = kappa * sqrt(var_B_wls),
               I_lead = sum(obs) / (2 * kappa^2 + 6 * kappa))
  })
  do.call(rbind, rows)
}

## Internal regression test: the G_avg column of the paper's Table godambe
## (clean_run/prop33_constants.csv of the replication material).
.cnbg_selftest_information <- function(tol_G = 0.01) {
  target_G <- c(0.1644, 0.06096, 0.01199, 0.001653, 7.538e-5)
  Gs <- sapply(c(1, 3, 10), function(a) cnbg_information_kappa(alpha = a)$G)
  G_avg <- rowMeans(Gs)
  if (!all(abs(G_avg / target_G - 1) < tol_G))
    stop("cnbg information bound does not reproduce the paper's table\n  G_avg = ",
         paste(signif(G_avg, 4), collapse = ", "))
  invisible(TRUE)
}
