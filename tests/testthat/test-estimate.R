test_that("estimators are unbiased under MCAR and EC is robust to non-stationarity", {
  set.seed(11)
  mc <- pc_montecarlo(R = 60, n_old = 3000, n_new = 800, k = 4, regime = "MCAR")
  expect_true(all(abs(mc$bias) < 0.03))
  mc2 <- pc_montecarlo(R = 60, n_old = 3000, n_new = 800, k = 4, regime = "MNAR_nonstat")
  expect_lt(abs(mc2$bias[mc2$estimator == "ec"]), 0.03)
  expect_gt(abs(mc2$bias[mc2$estimator == "sm"]), 0.04)
})

test_that("SSM repairs state-dependent attrition and T2 detects it", {
  set.seed(12)
  mc <- pc_montecarlo(R = 60, n_old = 3000, n_new = 800, k = 4, regime = "MNAR_state")
  expect_lt(abs(mc$bias[mc$estimator == "ssm"]), 0.03)
  expect_gt(attr(mc, "tests")[["reject_T2"]], 0.8)
})

test_that("decomposition identity and bounds hold", {
  set.seed(13)
  s <- pc_simulate(n_old = 2000, n_new = 500, k = 4, regime = "MNAR_trait")
  f <- pc_estimate(s$old, s$new, k = 4)
  d <- f$decomposition
  expect_equal(unname(d["naive"]), unname(d["conditioning_sm"] + d["attrition"]))
  b <- pc_bounds(f, y_range = range(c(s$old$y, s$new$y), na.rm = TRUE))
  expect_lt(b$bounds[["lower"]], b$bounds[["upper"]])
  expect_true(is.finite(b$breakdown[["ec_rho_star"]]))
})

test_that("pc_from_wide reproduces pc_simulate inputs", {
  set.seed(14)
  T <- 10; n <- 400
  R <- matrix(rbinom(n * T, 1, 0.9), n, T); R[1:200, 1:4] <- 0L   # persons 1:200 enter at wave 5
  wide <- data.frame(cohort = c(rep(5L, 200), rep(1L, 200)))
  for (w in 1:T) wide[[paste0("r_w", w)]] <- R[, w]
  wide$y1 <- rnorm(n); wide$y5 <- rnorm(n)
  inp <- pc_from_wide(wide, cohort = "cohort", resp = paste0("r_w", 1:T), y_cols = c("1" = "y1", "5" = "y5"), t = 5, c_old = 1)
  expect_equal(inp$k, 4); expect_equal(inp$m, 4)
  expect_equal(nrow(inp$old), 200); expect_equal(nrow(inp$new), 200)
  f <- pc_estimate(inp$old, inp$new, k = inp$k, m = inp$m)
  expect_s3_class(f, "pc_estimate")
})
