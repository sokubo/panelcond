test_that("estimators are unbiased under MCAR and EC is robust to non-stationarity", {
  set.seed(11)
  mc <- pc_montecarlo(R = 60, n_old = 3000, n_new = 800, k = 4, regime = "MCAR")
  expect_true(all(abs(mc$bias) < 0.03))
  expect_true(all(abs(mc$mean_se / mc$mc_sd - 1)[mc$estimator != "ipw"] < 0.3))
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
  b <- pc_bounds(f, y_range = range(c(s$old$y, s$new$y), na.rm = TRUE), sd_new = sd(s$new$y))
  expect_lt(b$bounds[["lower"]], b$bounds[["upper"]])
  expect_true(is.finite(b$breakdown[["ec_rho_star"]]))
  expect_equal(b$breakdown[["ec_eta_star"]], f$estimates$estimate[f$estimates$estimator == "ec"])
  expect_equal(pc_bounds(f, y_range = c(4, -4))$bounds, pc_bounds(f, y_range = c(-4, 4))$bounds)
  expect_error(pc_bounds(f, y_range = c(-4, 4), sd_new = 0), "positive")
  expect_error(pc_bounds(f, y_range = c(0, 1, 2)), "two")
})

test_that("the EC standard error is the first-order variance over disjoint groups", {
  set.seed(16)
  s <- pc_simulate(n_old = 1500, n_new = 400, k = 4, regime = "MNAR_trait")
  f <- pc_estimate(s$old, s$new, k = 4)
  old <- s$old; io <- old$s_prior == 1 & !is.na(old$y)
  yo <- old$y[io]; es <- old$y_entry[io]; ens <- old$y_entry[!io]
  p <- sum(io) / nrow(old)
  v <- (var(yo) - 2 * (1 - p) * cov(yo, es) + (1 - p)^2 * var(es)) / sum(io) + (1 - p)^2 * var(ens) / sum(!io) + var(s$new$y) / nrow(s$new)
  expect_equal(f$estimates$se[f$estimates$estimator == "ec"], sqrt(v))
  ## a perfectly persistent outcome and a constant fresh outcome: the variance is that of the cohort's entry mean
  old2 <- old; old2$y[io] <- old2$y_entry[io]; new2 <- s$new; new2$y <- 1
  f2 <- pc_estimate(old2, new2, k = 4)
  expect_equal(f2$estimates$se[f2$estimates$estimator == "ec"]^2, (p * var(es) + (1 - p) * var(ens)) / nrow(old))
  expect_equal(f2$estimates$estimate[f2$estimates$estimator == "ec"], mean(old$y_entry) - 1)
})

test_that("item non-response is reported, T1 stays finite and indicators are checked", {
  set.seed(17)
  s <- pc_simulate(n_old = 1000, n_new = 300, k = 3)
  new <- s$new; new$y[1:20] <- NA
  old <- s$old; old$y_entry[which(old$s_prior == 1)[1:10]] <- NA
  expect_message(f <- pc_estimate(old, new, k = 3), "item non-response")
  expect_true(all(is.finite(f$tests$statistic)))
  expect_equal(unname(f$info$dropped[c("new_y", "y_entry")]), c(20, 10))
  expect_equal(f$estimates$n_new[f$estimates$estimator == "naive"], 280)
  bad <- s$old; bad$s_prior[1:5] <- NA
  expect_error(pc_estimate(bad, s$new, k = 3), "contains NA")
  bad <- s$old; bad$y <- factor(bad$y > 0)
  expect_error(pc_estimate(bad, s$new, k = 3), "numeric")
  expect_error(pc_estimate(s$old, s$new, k = 3, m = 5), "between 0 and k")
  expect_error(pc_estimate(s$old, s$new[0, ], k = 3), "at least one row")
})

test_that("the joint person bootstrap and EC-adj run and agree with the analytic results", {
  set.seed(18)
  s <- pc_simulate(n_old = 3000, n_new = 800, k = 4, regime = "MAR_X")
  f <- pc_estimate(s$old, s$new, k = 4, adjust = "x", nboot = 100)
  e <- f$estimates
  expect_equal(e$estimator, c("naive", "sm", "ssm", "ec", "ec_adj", "ipw"))
  expect_true(all(is.finite(e$se_boot)))
  ratio <- e$se_boot / e$se
  expect_true(all(abs(ratio[!is.na(ratio)] - 1) < 0.3))
  expect_lt(abs(e$estimate[e$estimator == "ec_adj"] - e$estimate[e$estimator == "ec"]), 0.06)
  expect_true(all(is.finite(f$tests$statistic_boot)))
  expect_error(pc_estimate(s$old, s$new, k = 4, adjust = "u"), "present in both")
  f1 <- pc_estimate(s$old, s$new, k = 4, nboot = 1)
  expect_true(all(is.na(f1$estimates$se_boot)))
})

test_that("pc_from_wide reproduces the inputs and checks its arguments", {
  set.seed(14)
  nw <- 10; n <- 400
  resp <- matrix(rbinom(n * nw, 1, 0.9), n, nw); resp[1:200, 1:4] <- 0L   # persons 1:200 enter at wave 5
  wide <- data.frame(id = seq_len(n), cohort = c(rep(5L, 200), rep(1L, 200)))
  for (w in 1:nw) wide[[paste0("r_w", w)]] <- resp[, w]
  wide$y1 <- rnorm(n); wide$y5 <- rnorm(n)
  expect_message(inp <- pc_from_wide(wide, cohort = "cohort", resp = paste0("r_w", 1:nw), y_cols = c("1" = "y1", "5" = "y5"), t = 5, c_old = 1, id = "id"),
                 "not first-wave respondents")
  expect_equal(inp$k, 4); expect_equal(inp$m, 4)
  expect_equal(nrow(inp$old), 200); expect_equal(nrow(inp$new), sum(resp[1:200, 5]))
  expect_false(anyNA(inp$new$y))
  expect_false(anyNA(inp$new$s_m1))              # t + k + 1 = 10 is the last wave, so SSM is feasible
  expect_false(anyNA(inp$old$s_next))
  f <- pc_estimate(inp$old, inp$new, k = inp$k, m = inp$m)
  expect_s3_class(f, "pc_estimate")
  ## columns for unfielded waves are all zero: m and s_m1 follow the last observed wave
  w2 <- wide; w2[paste0("r_w", 7:nw)] <- 0L
  inp2 <- suppressMessages(pc_from_wide(w2, cohort = "cohort", resp = paste0("r_w", 1:nw), y_cols = c("1" = "y1", "5" = "y5"), t = 5, c_old = 1))
  expect_equal(inp2$m, 1); expect_true(all(is.na(inp2$new$s_m1)))
  expect_error(pc_from_wide(rbind(wide, wide[1:3, ]), cohort = "cohort", resp = paste0("r_w", 1:nw), y_cols = c("1" = "y1", "5" = "y5"), t = 5, c_old = 1, id = "id"), "duplicated")
  expect_error(pc_from_wide(wide, cohort = "cohort", resp = paste0("r_w", 1:nw), y_cols = c("1" = "y1"), t = 5, c_old = 1), "comparison wave")
  expect_error(pc_from_wide(wide, cohort = "cohort", resp = paste0("r_w", 1:nw), y_cols = c("y1", "y5"), t = 5, c_old = 1), "named")
  w3 <- wide; w3[paste0("r_w", 5:nw)] <- 0L
  expect_error(pc_from_wide(w3, cohort = "cohort", resp = paste0("r_w", 1:nw), y_cols = c("1" = "y1", "5" = "y5"), t = 5, c_old = 1), "no response")
  w4 <- wide; w4$cohort[1:10] <- 9L
  expect_message(expect_message(pc_from_wide(w4, cohort = "cohort", resp = paste0("r_w", 1:nw), y_cols = c("1" = "y1", "5" = "y5"), t = 5, c_old = 1),
                                "not first-wave respondents"), "neither cohort")
})

test_that("pc_increment recovers the dose increment and its SE has no fresh-arm term", {
  set.seed(15)
  a <- pc_simulate(n_old = 6000, n_new = 10, k = 12, regime = "MNAR_trait")
  b <- pc_simulate(n_old = 6000, n_new = 10, k = 8,  regime = "MNAR_trait")
  A <- data.frame(y = a$old$y, y_entry = a$old$y_entry, s = a$old$s_prior)
  B <- data.frame(y = b$old$y, y_entry = b$old$y_entry, s = b$old$s_prior)
  inc <- pc_increment(A, B, kA = 12, kB = 8, nboot = 50)
  truth <- a$truth[["pce_surv"]] - b$truth[["pce_surv"]]
  expect_lt(abs(inc$ec_diff[["estimate"]] - truth), 3 * inc$ec_diff[["se"]] + 0.02)
  vc <- function(D) { i <- D$s == 1 & !is.na(D$y); p <- mean(i)
    (var(D$y[i]) - 2 * (1 - p) * cov(D$y[i], D$y_entry[i]) + (1 - p)^2 * var(D$y_entry[i])) / sum(i) + (1 - p)^2 * var(D$y_entry[!i]) / sum(!i) }
  expect_equal(inc$ec_diff[["se"]], sqrt(vc(A) + vc(B)))
  expect_lt(abs(inc$ec_diff[["se_boot"]] / inc$ec_diff[["se"]] - 1), 0.3)
  expect_warning(pc_increment(A, B, kA = 8, kB = 12), "should exceed")
  A2 <- A; A2$s[1:3] <- NA
  expect_error(pc_increment(A2, B, kA = 12, kB = 8), "contains NA")
})

test_that("pc_simulate marks the symmetric design infeasible when m < k and checks regimes", {
  set.seed(19)
  s <- pc_simulate(n_old = 300, n_new = 100, k = 4, m = 2)
  expect_true(all(is.na(s$new$s_m1)))
  f <- pc_estimate(s$old, s$new, k = 4, m = 2)
  expect_true(is.na(f$estimates$estimate[f$estimates$estimator == "ssm"]))
  expect_true(f$info$sm_partial)
  expect_error(pc_simulate(n_old = 50, n_new = 20, k = 2, regime = list(aX = 0)), "needs elements")
  expect_error(pc_simulate(n_old = 50, n_new = 20, k = 2, m = 3), "between 0 and k")
  expect_error(pc_montecarlo(R = 0, n_old = 50, n_new = 20, k = 2), "positive integer")
})
