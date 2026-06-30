###############################################################################
# Example: The Compound Negative Binomial-Gamma (CNBG) Amount-Triangle Model
#
# Demonstrates what a PAID-LOSS triangle can and cannot identify about reserve
# uncertainty. The headline of the method is not a better point reserve but a
# DIAGNOSTIC: whether amount-only inference is trustworthy on a given portfolio
# at all. This example shows
#   - the kappa posterior diagnostic separating a supported portfolio from an
#     unsupported one (the analogue of c_hat in the multinomial example),
#   - kappa failing to track the truth even on data drawn from the model (the
#     information ceiling), and
#   - the conditioning-respecting posterior predictive reserve.
#
# Unlike the other examples, this one leads with SIMULATION, because the
# diagnostic only means something when the truth is known. A real triangle
# (Taylor-Ashe) is included at the end -- it is the canonical no-Bayes case,
# which is itself the instructive outcome.
#
# Reference:
# Van Oirbeek, R. (2026). What Paid Triangles Can and Cannot Identify:
# An Information Ceiling for Loss-Reserve Uncertainty.
###############################################################################

library(hgr)

set.seed(2026)

# =============================================================================
# 1. SIMULATE A CNBG PAID-LOSS TRIANGLE
# =============================================================================

cat("================================================================\n")
cat("1. A CNBG paid-loss triangle (truth known)\n")
cat("================================================================\n\n")

# kappa = frailty dispersion (small = heterogeneous accident years);
# alpha = severity shape. This is the single-frailty regime of the paper:
# all heterogeneity is the reserve-dominant row frailty kappa_r.
sim <- simulate_cnbg(I = 10, J = 10, kappa = 5, alpha = 3)
tri <- sim$triangle

cat(sprintf("True kappa:        %d\n", 5))
cat(sprintf("Realised reserve:  %s\n",
            format(round(sim$R_true), big.mark = ",")))


# =============================================================================
# 2. FIT THE WORKING-LIKELIHOOD POSTERIOR (recommended default)
# =============================================================================

cat("\n================================================================\n")
cat("2. Fit CNBG-Bayes (deterministic grid; no Stan dependency)\n")
cat("================================================================\n\n")

fit <- fit_cnbg(tri, method = "bayes", prior = "ref")
print(fit)

# The WLS point estimate (lower edge of the design space) for comparison.
fit_wls <- fit_cnbg(tri, method = "wls")
print(fit_wls)


# =============================================================================
# 3. THE TRUST DIAGNOSTIC
# =============================================================================
# Does THIS triangle inform kappa, or is the posterior just returning the
# prior? Regimes: informative | prior_driven | flight | no_bayes.
# Read it as a routing decision, not a goodness-of-fit test.

cat("\n================================================================\n")
cat("3. Kappa posterior diagnostic\n")
cat("================================================================\n\n")

dg <- diagnose_kappa(fit)
print(dg)

cat(sprintf("\n  contraction c = %.3f  (c < 1 => posterior tighter than prior)\n",
            dg$c))
cat(sprintf("  mean shift  s = %.2f  (prior-SD units)\n", dg$s))
cat(sprintf("  KL          = %.2f\n", dg$KL))
cat(sprintf("  verdict: amount-only inference %s on this portfolio\n",
            ifelse(isTRUE(dg$informative), "is SUPPORTED",
                   "is NOT supported -> escalate to counts (NB-CL) or joint model")))


# =============================================================================
# 4. CONDITIONING-RESPECTING POSTERIOR PREDICTIVE RESERVE
# =============================================================================
# The observed triangle is held fixed; the per-row frailty is restored by a
# conjugate update; kappa_r is integrated over the prior (it is absorbed by the
# free row levels and not identified from amounts alone). Contrast with the
# WLS-cond point bootstrap, which omits parameter and frailty uncertainty.

cat("\n================================================================\n")
cat("4. Posterior predictive reserve\n")
cat("================================================================\n\n")

res_bayes <- reserve_cnbg(fit, B = 4000)
print(res_bayes)

res_wls <- bootstrap_cnbg_cond(tri, B = 4000)
print(res_wls)

cat(sprintf("\nRealised reserve falls in the CNBG-Bayes 95%% PI: %s\n",
            ifelse(sim$R_true >= res_bayes$pi[1] &&
                   sim$R_true <= res_bayes$pi[2], "yes", "no")))
cat("WLS-cond intervals are narrower -- they drop parameter, mean-structure,\n")
cat("and shared-frailty uncertainty; that gap is the value of the full draw.\n")


# =============================================================================
# 5. THE INFORMATION CEILING: kappa does not track the truth
# =============================================================================
# Even on data drawn exactly from the model, the amount triangle barely
# identifies kappa. Across a 10x wide range of true kappa, the posterior
# median hardly moves -- the signature of the ceiling, not an estimator bug.

cat("\n================================================================\n")
cat("5. kappa fails to track truth (the information ceiling)\n")
cat("================================================================\n\n")

cat(sprintf("%-12s %18s %14s\n", "true kappa", "post. median kappa", "regime"))
cat(paste(rep("-", 46), collapse = ""), "\n")

for (k_true in c(3, 5, 10, 20, 50)) {
  tri_k <- simulate_cnbg(I = 10, J = 10, kappa = k_true, alpha = 3)$triangle
  fit_k <- fit_cnbg(tri_k, method = "bayes", prior = "ref")
  dg_k  <- diagnose_kappa(fit_k)
  cat(sprintf("%-12d %18.2f %14s\n",
              k_true, fit_k$kappa_median, dg_k$regime))
}

cat("\nThe posterior median clusters near the prior region regardless of truth:\n")
cat("amounts identify the variance SCALE (A), not the frailty LOCATION (kappa).\n")


# =============================================================================
# 6. THE INFORMATION BOUND (Fisher / Godambe / van Trees)
# =============================================================================
# The theoretical counterpart to Section 5: the data Fisher information for
# kappa collapses and the prior fraction approaches 1 as kappa grows.

cat("\n================================================================\n")
cat("6. Fisher / Godambe / van Trees bound on kappa\n")
cat("================================================================\n\n")

info <- cnbg_information_kappa(kappa_grid = c(3, 5, 10, 20, 50), alpha = 3)
print(info, digits = 3, row.names = FALSE)

cat("\nprior_frac -> 1 and vantrees_sd flattens as kappa grows: at large kappa\n")
cat("the prior carries essentially all the information, at any triangle size.\n")


# =============================================================================
# 7. A REAL TRIANGLE: TAYLOR-ASHE (the canonical no-Bayes case)
# =============================================================================
# Taylor-Ashe is the field's most-cited 10x10 paid benchmark. The paper
# documents that the working-likelihood fit does not yield a usable kappa
# posterior on it -- the no-Bayes regime. That is the instructive outcome:
# the diagnostic flags a portfolio whose amount triangle cannot support a
# stable kappa inference, and the routing decision is to obtain counts.

cat("\n\n================================================================\n")
cat("7. Taylor-Ashe (real data) -- the diagnostic in the field\n")
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

fit_ta <- fit_cnbg(taylor_ashe, method = "bayes", prior = "ref")
print(fit_ta)

dg_ta <- diagnose_kappa(fit_ta)
print(dg_ta)

cat("\nWhatever the regime here, the lesson is the same: the diagnostic decides\n")
cat("WHETHER to trust amount-only inference, before any predictive procedure is\n")
cat("selected. When it does not engage, the count triangle (NB-CL) or the joint\n")
cat("frequency-severity model is the principled next step, not another\n")
cat("amount-only method.\n")


# =============================================================================
# 8. (OPTIONAL) THE EXACT FRAILTY-AUGMENTED POSTERIOR -- requires rstan
# =============================================================================
# A theoretical reference, fitted by HMC. Fragile on 10x10 triangles (the
# kappa-lambda funnel) and not for population-scale deployment; the working
# likelihood above is the recommended default. Uncomment to run with rstan.

# if (requireNamespace("rstan", quietly = TRUE)) {
#   fit_x <- fit_cnbg_exact(tri, prior = "ref", alpha = 3,
#                           chains = 4, iter = 2000, seed = 2026)
#   print(fit_x)               # reports Rhat / divergences as findings
#   print(diagnose_kappa(fit_x))
#   print(reserve_cnbg(fit_x))
# }
