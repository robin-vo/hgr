## Reference values: the paper's replication library (24_correspondence_lib.R, 51_posterior_exact_lib.R),
## run on the 3 x 3 example of Section 7 and on the Taylor-Ashe triangle with the seeds below.

ex3 <- rbind(c(1, 1, 1), c(10, 1, NA), c(1, NA, NA))

test_that("the chain-ladder fit and the first-order terms reproduce the paper's 3 x 3 example", {
  fit <- odp_fit(ex3)
  expect_equal(fit$f, c(13 / 11, 1.5))
  expect_equal(fit$F, c(22 / 39, 2 / 3, 1))
  expect_equal(fit$n, 6L); expect_equal(fit$p, 5L)
  expect_equal(fit$phi, 2.1756198347, tolerance = 1e-9)
  tt <- odp_variance_terms(fit)
  expect_equal(unname(tt[["level"]]), 7.2820333311, tolerance = 1e-9)
  expect_equal(unname(tt[["pattern"]]), 121.1661694122, tolerance = 1e-9)
  expect_equal(unname(tt[["process"]]), 13.6470698723, tolerance = 1e-9)
  expect_equal(sum(attr(tt, "shares")), 1)
})

test_that("L + P equals the full inverse-information quadratic form of the reserve gradient", {
  fit <- odp_fit(ex3); phi <- fit$phi; mu <- fit$mu; pi <- fit$pi; obs <- fit$obs
  I <- length(mu); J <- length(pi)
  Fi <- vapply(1:I, function(i) sum(pi[obs[i, ]]), 0)
  Q <- rbind(diag(J - 1), -1)
  Imm <- diag(Fi / (phi * mu), I); Imp <- obs / phi
  Ipp <- diag(vapply(1:J, function(j) sum(mu[obs[, j]]) / pi[j], 0) / phi, J)
  Ifull <- rbind(cbind(Imm, Imp %*% Q), cbind(t(Q) %*% t(Imp), t(Q) %*% Ipp %*% Q))
  d <- c(1 - Fi, as.vector(t(Q) %*% vapply(1:J, function(j) sum(mu[!obs[, j]]), 0)))
  tt <- odp_variance_terms(fit)
  expect_equal(unname(tt[["level"]] + tt[["pattern"]]), as.numeric(t(d) %*% solve(Ifull, d)), tolerance = 1e-10)
})

test_that("Theorem 4.1(i) holds exactly: the exhaustive no-refit variance equals (sigma_r^2 / phi) L", {
  fit <- odp_fit(ex3); radj <- fit$pool; obs <- fit$obs; m <- fit$m[obs]
  idx <- as.matrix(expand.grid(rep(list(1:6), 6)))                     # all 6^6 equally likely resamples
  Ystar <- matrix(rep(m, each = nrow(idx)), nrow(idx)) + matrix(radj[idx], nrow(idx)) * matrix(rep(sqrt(m), each = nrow(idx)), nrow(idx))
  ridx <- which(obs, arr.ind = TRUE)
  Cstar <- sapply(1:3, function(i) rowSums(Ystar[, ridx[, 1] == i, drop = FALSE]))
  open <- which(fit$lat < 3)
  g_nore <- as.vector(Cstar[, open] %*% (1 / fit$F[fit$lat[open]] - 1))
  v <- mean((g_nore - mean(g_nore))^2)
  expect_equal(v, 7.1879674106, tolerance = 1e-9)
  expect_equal(v, fit$pool_factor * unname(odp_variance_terms(fit)[["level"]]), tolerance = 1e-12)
})

test_that("the closed-form quasi-posterior gives the paper's exponents and moment conditions on the 3 x 3 example", {
  po <- odp_posterior(ex3)
  expect_equal(po$e, c(5.0560303894, 0.9192782526), tolerance = 1e-9)
  expect_equal(unname(po$c), c(5.5156695157, 0.9192782526, 0.4596391263), tolerance = 1e-9)
  expect_false(po$finite_mean); expect_false(po$finite_var)
  expect_equal(po$mean, Inf)
  po1 <- odp_posterior(ex3, alpha = 1)                                   # Dirichlet(1): c_j + 1, e_k + k
  expect_equal(po1$e, po$e + c(1, 2), tolerance = 1e-12)
  expect_true(po1$finite_mean)
})

test_that("the Beta factorisation is a density identity in the share ratios (log-flat and Dirichlet(alpha) priors)", {
  set.seed(3)
  I <- 6; J <- 6; mu <- 400 * exp(0.1 * (0:5)); pi <- c(.35, .25, .18, .12, .06, .04)
  tri <- outer(mu, pi); tri[] <- rpois(length(tri), tri); tri[row(tri) + col(tri) > I + 1] <- NA
  for (alpha in c(0, 1, 2.5)) {
    po <- odp_posterior(tri, alpha = alpha); s <- po$setup; Jp <- s$Jp
    vals <- replicate(300, {
      u <- runif(Jp - 1)
      Fk <- c(rev(cumprod(rev(u))), 1); pik <- diff(c(0, Fk))
      target <- sum((s$cj - 1) * log(pik)) - sum(s$ri * log(Fk[s$kk]))   # prod pi_j^{c_j - 1} prod F_i^{-r_i}
      jac <- sum((2:(Jp - 1) - 1) * log(u[2:(Jp - 1)]))                    # |d pi / d U| = prod_{l >= 2} U_l^{l - 1}
      target + jac - sum(dbeta(u, s$e, s$cj[-1], log = TRUE))
    })
    expect_lt(sd(vals), 1e-8)
  }
})

test_that("exact moments agree with independent draws from the closed form", {
  set.seed(4)
  I <- 8; J <- 5; mu <- rep(2000, I); pi <- c(.45, .25, .15, .10, .05)
  tri <- outer(mu, pi); tri[] <- rpois(length(tri), tri); tri[row(tri) + col(tri) > I + 1] <- NA
  po <- odp_posterior(tri, n = 200000, seed = 5)
  expect_true(po$finite_var)
  expect_equal(mean(po$draws$R), po$mean, tolerance = 0.005)
  expect_equal(var(po$draws$R), po$var, tolerance = 0.02)
  expect_equal(var(po$draws$g), po$var_parameter, tolerance = 0.02)
  expect_equal(unname(sum(po$components)), po$var)
})

test_that("the bootstrap, the closed form and the Gibbs sampler reproduce the paper's library on Taylor-Ashe", {
  skip_if_not_installed("ChainLadder")
  inc <- ChainLadder::cum2incr(ChainLadder::GenIns); inc[!(row(inc) + col(inc) <= 11)] <- NA
  b <- odp_bootstrap(inc, B = 500, seed = 2026)
  expect_equal(b$phi, 52601.36151147, tolerance = 1e-10)
  expect_equal(b$pool_factor * b$phi, 52600.80727413, tolerance = 1e-10)
  expect_equal(b$reserve, 18680855.611924, tolerance = 1e-10)
  expect_equal(sum(b$g_full), 9481133054.665970, tolerance = 1e-10)
  expect_equal(sum(b$g_nore), 9399645798.782644, tolerance = 1e-10)
  expect_equal(sum(b$g_cdiag), 9416192355.748899, tolerance = 1e-10)
  expect_equal(sum(b$R_full), 9491912889.359859, tolerance = 1e-10)
  expect_equal(sum(b$R_nore), 9412843562.946671, tolerance = 1e-10)
  expect_equal(sum(b$R_cdiag), 9422778913.529692, tolerance = 1e-10)
  expect_equal(sum(b$R_cond), 9317224057.439194, tolerance = 1e-10)
  expect_equal(as.vector(b$terms), c(4382271506832.196, 3311921772229.811, 982638439386.3625), tolerance = 1e-10)
  po <- odp_posterior(inc, phi = b$phi)
  expect_equal(po$e, c(63.256366, 194.885621, 286.073280, 350.709382, 341.497986, 303.318328, 242.258235, 161.979172, 72.878627), tolerance = 1e-7)
  expect_equal(po$mean, 18846372.832312, tolerance = 1e-10)
  expect_equal(po$var, 9007874177925.3105, tolerance = 1e-10)
  expect_equal(po$var_parameter, 8016529307392.9648, tolerance = 1e-10)
  d <- odp_posterior(inc, phi = b$phi, n = 2000, seed = 7)$draws
  expect_equal(sum(d$g), 37740668104.020859, tolerance = 1e-10)
  expect_equal(sum(d$R), 37706286573.626953, tolerance = 1e-10)
  gb <- odp_gibbs(inc, phi = b$phi, n_iter = 600, burn = 100, seed = 11)
  expect_equal(sum(gb$g), 9526807947.281364, tolerance = 1e-10)
  expect_equal(sum(gb$R), 9539467315.550734, tolerance = 1e-10)
  expect_equal(b$decomposition[["level"]], 1, tolerance = 0.15)            # exactly one up to Monte Carlo error at B = 500
})

test_that("the centred and rescaled pool has mean zero and variance phi, so the factor is one", {
  skip_if_not_installed("ChainLadder")
  inc <- ChainLadder::cum2incr(ChainLadder::GenIns); inc[!(row(inc) + col(inc) <= 11)] <- NA
  b <- odp_bootstrap(inc, B = 50, pool = "centred", seed = 1)
  expect_equal(mean(b$pool), 0, tolerance = 1e-12)
  expect_equal(mean(b$pool^2), b$phi, tolerance = 1e-12)
  expect_equal(b$pool_factor, 1, tolerance = 1e-12)
})

test_that("the certificates reproduce the supplement's 3 x 3 example and the fit refuses inadmissible triangles", {
  tri <- rbind(c(100, 1, 1), c(100, 5, NA), c(100, NA, NA))               # Supplement S5: every cumulative bound positive, a factor below one is attainable
  b <- odp_bootstrap(tri, B = 200, seed = 1)
  ce <- odp_certificates(b)
  expect_true(ce$no_pole); expect_false(ce$factors_above_one)
  expect_true(is.finite(ce$g_max))
  bad <- rbind(c(1, 1, 1), c(-2, 3, NA), c(2, NA, NA))                    # positive margins, first factor -3 (Section 2.1)
  expect_error(odp_fit(bad), "below one")
  expect_error(odp_fit(rbind(c(1, NA), c(1, 1))), "contiguous|further")
})

test_that("structural-zero columns are removed from the pool and from the posterior, with the paper's degrees-of-freedom convention", {
  set.seed(6)
  tri <- outer(rep(500, 6), c(.5, .3, .15, .05, 0, 0)); tri[] <- rpois(length(tri), tri); tri[row(tri) + col(tri) > 7] <- NA
  tri[, 5:6] <- ifelse(is.na(tri[, 5:6]), NA, 0)
  fit <- odp_fit(tri)
  expect_equal(fit$zero_columns, 5:6)
  expect_equal(fit$n, sum(!is.na(tri)) - 3L)                             # the three observed cells of the zero columns leave the pool
  expect_equal(fit$p, 6L + 6L - 1L)
  expect_equal(odp_fit(tri, dof = "reduced")$p, 6L + 4L - 1L)
  po <- odp_posterior(tri)
  expect_equal(po$columns, 1:4)
  expect_error(odp_posterior(tri, zero_columns = "retain"), "improper")
  expect_true(is.list(odp_posterior(tri, zero_columns = "retain", alpha = 1)))
})
