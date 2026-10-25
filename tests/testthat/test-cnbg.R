test_that("the data-generating process reproduces the paper's calibration", {
  set.seed(1)
  s <- simulate_cnbg(kappa_r = 5, alpha = 3)
  expect_equal(dim(s$triangle), c(10L, 10L))
  expect_equal(sum(is.na(s$triangle)), 45L)
  expect_equal(s$A, (1 + 1 / 3) * 5000)
  expect_equal(round(s$mu[10], 2), round(30 * exp(0.05 * 9), 2))
  expect_true(all(s$eps == 1))
  set.seed(1)
  s2 <- simulate_cnbg(kappa_r = 5, kappa_c = 5, delta = 0.1)
  expect_true(any(s2$eps != 1))
  expect_true(s2$R_true >= 0)
})

test_that("the second-moment estimator is the nu^-2 NNLS with the cap and df factor", {
  set.seed(2)
  ms <- cnbg_stage1(simulate_cnbg(kappa_r = 5)$triangle)
  w <- fit_cnbg_wls(ms)
  nu <- as.numeric(ms$nu_hat); z <- (ms$dat$x - nu)^2; sw <- 1 / nu
  nn <- nnls::nnls(cbind(sw * nu, sw * nu^2), sw * z)
  B <- max(nn$x[2], 1 / 1000)
  expect_equal(w$kappa, (1 / B) * (ms$n_obs - ms$q) / ms$n_obs)
  expect_true(w$kappa <= 1000)
})

test_that("the two-stage target reproduces the paper's table", {
  tt <- two_stage_target(kappa_r = c(3, 5, 10, 20, 50), kappa_c = 5)
  expect_equal(round(tt$kappa_target, 1), c(15.2, 17.9, 22.2, 26.3, 30.4))
  expect_true(all(tt$B_star < 0))
  expect_equal(round(tt$c_q, 3), c(0.461, 0.440, 0.417, 0.402, 0.392))        # results/two_stage_cellshock_report.txt
  expect_equal(round(tt$kappa_c_threshold, 1), c(17.5, 13.7, 10.9, 9.5, 8.6))
  single <- two_stage_target(kappa_r = 5, kappa_c = Inf)
  expect_equal(single$kappa_target, Inf)
})

test_that("the information bound reproduces the paper's Godambe column", {
  expect_true(hgr:::.cnbg_selftest_information())
  info <- cnbg_information_kappa()
  expect_true(all(c("I_diag", "Delta", "G", "cv_G", "cv_wls", "I_lead") %in% names(info)))
  expect_false("vantrees_sd" %in% names(info))
})

test_that("cell shapes respect the lower bound CV2 >= 2", {
  s <- cnbg_cell_shapes(c(1e3, 1e4, 1e5, 1e6), A = 6667, kappa = 5, alpha = 3)
  expect_true(all(s > 0 & s <= 0.5))
})

test_that("the grid working fit and the CNBG-Mod predictive run on a simulated triangle", {
  set.seed(3)
  sim <- simulate_cnbg(kappa_r = 5, alpha = 3)
  fit <- fit_cnbg(sim$triangle, engine = "grid")
  expect_true(fit$ok)
  expect_equal(fit$abs$mean, "absorbed")
  expect_true(is.finite(fit$abs$adj$A[2]))
  expect_true(fit$abs$adj$A[1] < sim$A * 3 && fit$abs$adj$A[3] > sim$A / 3)
  r <- reserve_cnbg(fit, B = 300)
  expect_s3_class(r, "cnbg_reserve")
  expect_equal(r$method, "CNBG-Mod-Ref")
  expect_true(r$pi[1] < r$pi[2])
  d <- diagnose_kappa(fit)
  expect_true(d$regime %in% c("informative", "not_informative"))
  expect_true(d$tier %in% c("prior-driven", "weak", "moderate", "strong"))
  r2 <- reserve_cnbg(fit, B = 200, prior = "empb")
  expect_equal(r2$method, "CNBG-Mod-EmpB")
  b <- bootstrap_cnbg_cond(sim$triangle, B = 300)
  expect_equal(b$method, "WLS-cond")
})

test_that("prepare_triangle applies the paper's negative-increment rule", {
  x <- simulate_cnbg()$full
  x[2, 3] <- -5
  p <- prepare_triangle(x)
  expect_equal(p$n_negative, 1L)
  expect_true(is.na(p$triangle[2, 3]))
  expect_equal(sum(is.na(p$triangle)), 46L)
  q <- prepare_triangle(x, negatives = "floor")
  expect_equal(q$triangle[2, 3], 1)
})
