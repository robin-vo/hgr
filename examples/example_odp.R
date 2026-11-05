###############################################################################
# Example: The ODP bootstrap and the decomposition of its predictive variance
#
# Follows Sections 3, 4 and 7 of the paper on the Taylor-Ashe triangle:
#   1. the chain-ladder fit as the ODP quasi-likelihood fit, with the Pearson
#      dispersion and the residual-pool factor,
#   2. the first-order level, pattern and process terms,
#   3. Algorithm 1 with its variants and the decomposition ratios,
#   4. the closed-form quasi-posterior, its moment conditions and exact
#      moments, with independent draws,
#   5. the two certificates on the refit, and
#   6. the centred and rescaled pool.
#
# Reference:
# Van Oirbeek, R. (2026). The over-dispersed Poisson bootstrap: A
# decomposition of predictive variance. Working paper.
###############################################################################

library(hgr)
library(ChainLadder)

# Incremental Taylor-Ashe triangle, future cells NA
inc <- cum2incr(GenIns)
inc[!(row(inc) + col(inc) <= 11)] <- NA

# =============================================================================
# 1. THE FIT
# =============================================================================
fit <- odp_fit(inc)
fit

# =============================================================================
# 2. FIRST-ORDER TERMS (equation (5) of the paper)
# =============================================================================
odp_variance_terms(fit)

# =============================================================================
# 3. ALGORITHM 1 WITH ITS VARIANTS (same seed, same draws as the paper's scripts)
# =============================================================================
b <- odp_bootstrap(fit, B = 4000, seed = 2026)
b
b$decomposition          # level ratio one up to Monte Carlo error; pattern and predictive ratios near one
b$intervals              # full, no-refit, diagonal-fixed and process-only predictives

# =============================================================================
# 4. THE CLOSED-FORM QUASI-POSTERIOR (Proposition 7.1)
# =============================================================================
po <- odp_posterior(fit, n = 20000, seed = 1)
po
po$tail_criterion        # min_{l >= k_i} e_l per open row: > 1 for a finite mean, > 2 for a finite variance
var(b$R_full) / po$var   # bootstrap variance over the exact quasi-posterior variance

# The 3 x 3 example of Section 7: an infinite posterior mean under the log-flat prior
ex3 <- rbind(c(1, 1, 1), c(10, 1, NA), c(1, NA, NA))
odp_posterior(ex3)
odp_posterior(ex3, alpha = 1)       # Dirichlet(1) prior on the pattern: e_k + k, finite mean and variance

# =============================================================================
# 5. CERTIFICATES FROM THE OBSERVED TRIANGLE
# =============================================================================
odp_certificates(b)

# =============================================================================
# 6. THE CENTRED AND RESCALED POOL (Remark 4.5)
# =============================================================================
bc <- odp_bootstrap(fit, B = 4000, pool = "centred", seed = 2026)
bc$pool_factor           # one by construction
