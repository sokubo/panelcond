test_that("the JLPS-type schedule has a four-dimensional identified set", {
  d <- pc_design_rank(c(1, 5, 13), 19)
  s <- d$summary
  expect_equal(c(s$cells, s$parameters, s$rank), c(41, 39, 35))
  expect_equal(s$stride_d, 4)
  expect_equal(s$kernel_dim, 4)
  expect_equal(s$identified_dim, 4)
  expect_true(s$identified_equals_d)
  expect_equal(s$components, 1)
  expect_true(s$trapezoid)
  expect_true(s$followup_condition)
  expect_equal(ncol(d$tau_null_basis), s$identified_dim)
  expect_equal(crossprod(d$tau_null_basis), diag(s$identified_dim), tolerance = 1e-10)
})

test_that("identified functionals are the ones the theory names", {
  d <- pc_design_rank(c(1, 5, 13), 19)
  lam <- function(i, w) { v <- numeric(19); v[i] <- w; v }
  expect_true(pc_is_identified(d, lam(c(1, 5, 9), c(1, -2, 1))))    # centred lag-4 second difference at s = 5
  expect_true(pc_is_identified(d, lam(c(5, 9, 13), c(1, -2, 1))))   # and at s = 9
  expect_false(pc_is_identified(d, lam(c(2, 3, 4), c(1, -2, 1))))   # ordinary second difference
  expect_false(pc_is_identified(d, lam(5, 1)))                      # a level
  expect_false(pc_is_identified(d, lam(c(1, 5), c(-1, 1))))         # a lag-4 increment
  # the identified space is larger than the span of the lag-4 second differences
  expect_true(pc_is_identified(d, lam(c(2, 3, 4), c(1, -2, 1)) - lam(c(6, 7, 8), c(1, -2, 1))))
  expect_error(pc_is_identified(d, numeric(3)), "one weight per observed tenure")
})

test_that("stride one with too short a follow-up leaves more than the affine line", {
  d <- pc_design_rank(c(1, 4, 8), 8)                 # w = 0 < (e_2 - e_1) - d = 2
  expect_equal(d$summary$stride_d, 1)
  expect_equal(d$summary$identified_dim, 3)
  expect_false(d$summary$identified_equals_d)
  expect_false(d$summary$followup_condition)
  expect_true(d$summary$trapezoid)
})

test_that("the three-cohort formula holds on trapezoids", {
  for (E in list(c(1, 5, 13), c(1, 4, 8), c(1, 7, 10), c(1, 3, 4), c(1, 6, 11))) {
    for (T in seq(max(E), max(E) + 8)) {
      d <- pc_design_rank(E, T); s <- d$summary
      expect_equal(s$identified_dim, max(s$stride_d, (E[2] - E[1]) - s$followup_w))
      expect_equal(s$three_cohort_dimension, s$identified_dim)
      expect_equal(s$components, 1)
      expect_equal(s$kernel_dim, s$identified_dim)
    }
  }
})

test_that("a disconnected support has a kernel larger than its tenure projection", {
  # cohorts 1, 4, 6, each observed at entry and once more: no two share a period
  six <- pc_design_rank(cells = data.frame(e = rep(c(1, 4, 6), each = 2), t = c(1, 2, 4, 5, 6, 7)))
  s <- six$summary
  expect_equal(c(s$cells, s$parameters, s$rank), c(6, 9, 6))
  expect_equal(s$kernel_dim, 3)
  expect_equal(s$identified_dim, 1)
  expect_equal(s$components, 3)
  expect_false(s$trapezoid)
  expect_true(is.na(s$followup_condition))
  expect_equal(s$kernel_dim, s$identified_dim + s$components - 1)
  expect_false(pc_is_identified(six, c(0, 1)))       # tau(2) is not identified

  # two components, but the cohorts 1 and 2 share period 2
  five <- pc_design_rank(cells = data.frame(e = rep(c(1, 2, 5), each = 2), t = c(1, 2, 2, 3, 5, 6)))
  expect_equal(five$summary$components, 2)
  expect_equal(five$summary$kernel_dim, 2)
  expect_equal(five$summary$identified_dim, 1)
  expect_equal(five$summary$kernel_dim, five$summary$identified_dim + five$summary$components - 1)
})

test_that("kernel_dim = identified_dim + (components - 1) on every support in a sweep", {
  mk <- function(E, T) do.call(rbind, lapply(E, function(e) data.frame(e = e, t = seq.int(e, T))))
  set.seed(11); bad <- 0L; disc <- 0L
  for (rep in 1:200) {
    E <- sort(sample(1:6, sample(2:3, 1))); T <- sample(max(E):8, 1)
    cl <- mk(E, T)
    drop <- sample(nrow(cl), sample(0:2, 1))
    if (length(drop)) cl <- cl[-drop, , drop = FALSE]
    if (!all(E %in% cl$e) || !any(cl$t - cl$e + 1 == 1)) next
    s <- pc_design_rank(cells = cl)$summary
    disc <- disc + (s$components > 1L)
    if (s$kernel_dim != s$identified_dim + s$components - 1L) bad <- bad + 1L
  }
  expect_equal(bad, 0L)
  expect_gt(disc, 0L)
})

test_that("a single cohort identifies nothing about the path", {
  d <- pc_design_rank(3, 8)
  expect_equal(d$summary$identified_dim, length(d$tenures) - 1L)
  expect_true(is.na(d$summary$stride_d))
  expect_false(pc_is_identified(d, c(0, 1, rep(0, 4))))
})

test_that("input checks fire", {
  expect_error(pc_design_rank(c(5, 1), 9), "strictly increasing")
  expect_error(pc_design_rank(c(1, 5), 3), "at least the last entry wave")
  expect_error(pc_design_rank(c(1, 5)), "supply")
  expect_error(pc_design_rank(cells = data.frame(e = 1, t = 0)), "not admissible")
  expect_error(pc_is_identified(list(), 1), "must come from")
})
