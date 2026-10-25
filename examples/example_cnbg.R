###############################################################################
# Example: The Compound Negative Binomial-Gamma (CNBG) Amount-Triangle Model
#
# Follows the pipeline of the paper's simulation and Schedule P studies:
#   1. a triangle drawn from the CNBG model at the paper's calibration,
#   2. the Stage-1 mean structure and the second-moment estimator, which sits
#      at its cap on most correctly specified triangles (its free-level target
#      is negative; see two_stage_target()),
#   3. the working-likelihood posterior of A from the absorbed fit,
#   4. the modular predictive CNBG-Mod against the point-estimate WLS-cond,
#   5. the kappa_B diagnostic from the full-mean fit,
#   6. the information on kappa_B at the calibration, and
#   7. a real triangle (Taylor-Ashe).
#
# Reference:
# Van Oirbeek, R. (2026). What Paid Triangles Identify about Reserve
# Uncertainty: Information Boundaries in Loss Reserving.
###############################################################################

library(hgr)
set.seed(2026)

# =============================================================================
# 1. SIMULATE A CNBG PAID-LOSS TRIANGLE
# =============================================================================
cat("================================================================\n")
cat("1. A CNBG paid-loss triangle (truth known)\n")
cat("================================================================\n\n")

# kappa_r = row-frailty dispersion (small = heterogeneous accident years);
# kappa_c = cell-shock dispersion (Inf = the single-frailty main design);
# alpha = severity shape. delta > 0 adds the calendar drift of Phase 16.
sim <- simulate_cnbg(I = 10, J = 10, kappa_r = 5, kappa_c = Inf, alpha = 3)
tri <- sim$triangle
cat(sprintf("True kappa_r: %g   true A: %.0f   realised reserve: %s\n",
            sim$kappa_r, sim$A, format(round(sim$R_true), big.mark = ",")))

# =============================================================================
# 2. STAGE 1 AND THE SECOND-MOMENT ESTIMATOR
# =============================================================================
cat("\n================================================================\n")
cat("2. Stage 1 and the two-stage estimator\n")
cat("================================================================\n\n")

ms <- cnbg_stage1(tri, levels = "free")
print(ms)
wls <- fit_cnbg_wls(ms)
cat(sprintf("second-moment estimate: kappa_B = %.1f%s, A = %.0f\n",
            wls$kappa, if (wls$at_cap) " (at the cap)" else "", wls$A))
cat("\nIts population target under free levels (Proposition 4.2 / Remark 4.5):\n")
print(two_stage_target(kappa_r = c(3, 5, 10, 20, 50), kappa_c = c(5, 20, Inf))[, c(1, 2, 3, 5, 6)],
      digits = 3, row.names = FALSE)
cat("B_target <= 0 means the constrained estimator converges to the cap.\n")

# =============================================================================
# 3. THE WORKING-LIKELIHOOD POSTERIOR OF A (ABSORBED FIT)
# =============================================================================
cat("\n================================================================\n")
cat("3. Working-likelihood posterior (absorbed mean, grid engine)\n")
cat("================================================================\n\n")

# engine = "grid" needs no Stan; engine = "nuts" is the sampler the paper ran.
f_abs <- fit_cnbg_working(ms, alpha = 3, mean = "absorbed", prior = "ref", engine = "grid")
print(f_abs)

# =============================================================================
# 4. CNBG-Mod AGAINST WLS-cond
# =============================================================================
cat("\n================================================================\n")
cat("4. Predictive reserve: CNBG-Mod and WLS-cond\n")
cat("================================================================\n\n")

fit <- fit_cnbg(tri, prior = "ref", alpha = 3, engine = "grid")
print(fit)
res_mod <- reserve_cnbg(fit, B = 4000)
print(res_mod)
res_wls <- bootstrap_cnbg_cond(tri, B = 4000)
print(res_wls)
cat(sprintf("\nRealised reserve inside the CNBG-Mod 95%% interval: %s; inside WLS-cond: %s\n",
            ifelse(sim$R_true >= res_mod$pi[1] && sim$R_true <= res_mod$pi[2], "yes", "no"),
            ifelse(sim$R_true >= res_wls$pi[1] && sim$R_true <= res_wls$pi[2], "yes", "no")))
cat("WLS-cond carries neither parameter nor mean-structure uncertainty and is\n")
cat("the lower edge of the paper's design space.\n")

# =============================================================================
# 5. THE kappa_B DIAGNOSTIC
# =============================================================================
cat("\n================================================================\n")
cat("5. kappa_B posterior diagnostic (full-mean fit)\n")
cat("================================================================\n\n")

print(diagnose_kappa(fit))
cat("\nUnder free levels kappa_r is absorbed by the levels and kappa_c survives\n")
cat("only below a threshold, so 'not_informative' is the expected verdict on\n")
cat("most paid triangles; the predictive takes kappa_B from the prior.\n")

cat("\nAcross true kappa_r the posterior median hardly moves:\n")
cat(sprintf("%-12s %20s %16s\n", "true kappa_r", "post. median kappa_B", "regime"))
for (k_true in c(3, 10, 50)) {
  tri_k <- simulate_cnbg(kappa_r = k_true, alpha = 3)$triangle
  fit_k <- fit_cnbg(tri_k, engine = "grid")
  dg_k  <- diagnose_kappa(fit_k)
  cat(sprintf("%-12g %20.2f %16s\n", k_true,
              if (is.null(fit_k$full) || !is.null(fit_k$full$failed)) NA else fit_k$full$adj$kappa[2], dg_k$regime))
}

# =============================================================================
# 6. INFORMATION ON kappa_B AT THE CALIBRATION
# =============================================================================
cat("\n================================================================\n")
cat("6. Fisher / Godambe information on kappa_B\n")
cat("================================================================\n\n")

info <- cnbg_information_kappa(kappa_grid = c(3, 5, 10, 20, 50), alpha = 3)
print(info, digits = 3, row.names = FALSE)
cat("\nG falls by three orders of magnitude from kappa_B = 3 to 50: the second\n")
cat("moment of paid amounts carries less and less on the dispersion.\n")

# =============================================================================
# 7. A REAL TRIANGLE: TAYLOR-ASHE
# =============================================================================
cat("\n================================================================\n")
cat("7. Taylor-Ashe (real data)\n")
cat("================================================================\n\n")

taylor_ashe <- matrix(c(
  357848, 766940, 610542, 482940, 527326, 574398, 146342, 139950, 227229, 67948,
  352118, 884021, 933894, 1183289, 445745, 320996, 527804, 266172, 425046, NA,
  290507, 1001799, 926219, 1016654, 750816, 146923, 495992, 280405, NA, NA,
  310608, 1108250, 776189, 1562400, 272482, 352053, 206286, NA, NA, NA,
  443160, 693190, 991983, 769488, 504851, 470639, NA, NA, NA, NA,
  396132, 937085, 847498, 805037, 705960, NA, NA, NA, NA, NA,
  440832, 847631, 1131398, 1063269, NA, NA, NA, NA, NA, NA,
  359480, 1061648, 1443370, NA, NA, NA, NA, NA, NA, NA,
  376686, 986608, NA, NA, NA, NA, NA, NA, NA, NA,
  344014, NA, NA, NA, NA, NA, NA, NA, NA, NA
), nrow = 10, ncol = 10, byrow = TRUE)

ta <- prepare_triangle(taylor_ashe)            # the paper's data rules (negatives to NA)
fit_ta <- fit_cnbg(ta$triangle, engine = "grid")
print(fit_ta)
print(diagnose_kappa(fit_ta))
print(reserve_cnbg(fit_ta, B = 4000))

# =============================================================================
# 8. (OPTIONAL) THE EXACT FRAILTY-AUGMENTED POSTERIOR -- requires rstan
# =============================================================================
# The paper's per-portfolio reference, at roughly 200 times the cost of the
# working fit. levels = "anchored" holds the cell means at the Stage-1 fit;
# levels = "free" makes the accident-year levels parameters (Phase 17).
#
# if (requireNamespace("rstan", quietly = TRUE)) {
#   fx <- fit_cnbg_exact(tri, prior = "ref", alpha = 3, levels = "anchored",
#                        chains = 2, iter = 1000, warmup = 500, seed = 2026)
#   print(fx); print(diagnose_kappa(fx)); print(reserve_cnbg(fx))
# }
