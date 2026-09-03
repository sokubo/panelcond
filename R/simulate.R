#' Simulate a continuing and a fresh cohort
#'
#' Data-generating process of the paper's Monte Carlo study: a fixed observed
#' covariate `x`, a persistent unobserved trait `u`, independent per-wave
#' transients, a saturating dose-response `tau(k) = tau_inf * (1 - exp(-k/kappa))`
#' with heterogeneity `tau_i = tau(k) * (1 + gam * u)`, and per-wave response
#' decisions with dropout logit `a0 + aX*x + (aU + shift)*u + aE*eps_w`.
#'
#' @param n_old,n_new cohort sizes.
#' @param k dose of the continuing cohort at the comparison wave.
#' @param m matched follow-up waves for the fresh cohort (default `k`).
#' @param regime one of `"MCAR"`, `"MAR_X"`, `"MNAR_trait"`, `"MNAR_nonstat"`,
#'   `"MNAR_state"`, or a named list with elements `aX`, `aU`, `aE`, `shift_old`.
#' @param tau_inf,kappa,gam,a0 dose-response and attrition parameters.
#' @return A list with `old`, `new` (ready for [pc_estimate()]), `truth` (the
#'   population PCE, the survivor PCE among the continuing cohort's respondents,
#'   and the survivor PCE among those who also respond at `t + 1`), `k`, `m`.
#' @export
pc_simulate <- function(n_old = 4800, n_new = 960, k = 4, m = k, regime = "MNAR_trait",
                        tau_inf = 0.30, kappa = 4, gam = 0.30, a0 = -2.2) {
  regimes <- list(MCAR = list(aX = 0, aU = 0, aE = 0, shift_old = 0),
                  MAR_X = list(aX = 0.8, aU = 0, aE = 0, shift_old = 0),
                  MNAR_trait = list(aX = 0.4, aU = 0.8, aE = 0, shift_old = 0),
                  MNAR_nonstat = list(aX = 0.4, aU = 0.8, aE = 0, shift_old = 0.6),
                  MNAR_state = list(aX = 0.4, aU = 0.4, aE = 0.8, shift_old = 0))
  reg <- if (is.character(regime)) regimes[[match.arg(regime, names(regimes))]] else regime
  tau_k <- function(kk) tau_inf * (1 - exp(-kk / kappa))
  sim_cohort <- function(n, k_prior, shift, k_future) {
    x <- stats::rbinom(n, 1, 0.5); u <- stats::rnorm(n)
    W <- k_prior + 1 + k_future
    eps <- matrix(stats::rnorm(n * W), n, W)
    Y0 <- 0.5 * x + 0.5 * u + eps
    t <- k_prior + 1
    tau_i <- tau_k(k_prior) * (1 + gam * u)
    lp <- a0 + reg$aX * x + reg$aE * eps
    if (k_prior > 0) lp[, seq_len(k_prior)] <- lp[, seq_len(k_prior), drop = FALSE] + (reg$aU + shift) * u
    lp[, t:W] <- lp[, t:W, drop = FALSE] + reg$aU * u
    keep <- matrix(stats::rbinom(n * W, 1, 1 - stats::plogis(lp)), n, W)
    s_prior <- if (k_prior > 0) as.integer(rowSums(keep[, seq_len(k_prior), drop = FALSE]) == k_prior) else rep(1L, n)
    s_next <- keep[, t]
    s_fut <- function(j) if (j > 0) as.integer(rowSums(keep[, t:(t + j - 1), drop = FALSE]) == j) else rep(1L, n)
    list(x = x, u = u, y_entry = Y0[, 1], y = Y0[, t] + tau_i, tau_i = tau_i, s_prior = s_prior, s_next = s_next,
         s_m = s_fut(min(k_future, max(k_future - 1, 0))), s_m1 = s_fut(k_future))
  }
  o <- sim_cohort(n_old, k, reg$shift_old, 1)
  nw <- sim_cohort(n_new, 0, 0, m + 1)
  old <- data.frame(y = ifelse(o$s_prior == 1, o$y, NA), y_entry = o$y_entry, s_prior = o$s_prior, s_next = o$s_next, x = o$x)
  new <- data.frame(y = nw$y, s_m = nw$s_m, s_m1 = nw$s_m1, x = nw$x)
  io <- o$s_prior == 1
  list(old = old, new = new, k = k, m = m,
       truth = c(pce_pop = mean(o$tau_i), pce_surv = mean(o$tau_i[io]), pce_surv_next = mean(o$tau_i[io & o$s_next == 1]),
                 tau_k = tau_k(k)))
}

#' Monte Carlo evaluation of the estimators
#'
#' @param R number of replications.
#' @param ... arguments passed to [pc_simulate()].
#' @return A data frame with one row per estimator: bias and RMSE relative to the
#'   estimator's own estimand (population PCE for IPW, survivor PCE otherwise),
#'   coverage of nominal 95% intervals, and the rejection rates of the two
#'   diagnostic tests.
#' @export
pc_montecarlo <- function(R = 200, ...) {
  args <- list(...)
  rows <- vector("list", R)
  for (r in seq_len(R)) {
    s <- do.call(pc_simulate, args)
    f <- pc_estimate(s$old, s$new, k = s$k, m = s$m)
    e <- f$estimates
    rows[[r]] <- c(naive = e$estimate[1], sm = e$estimate[2], ssm = e$estimate[3], ec = e$estimate[4], ipw = e$estimate[5],
                   se_naive = e$se[1], se_sm = e$se[2], se_ssm = e$se[3], se_ec = e$se[4],
                   T1 = f$tests$statistic[1], T2 = f$tests$statistic[2],
                   pop = s$truth[["pce_pop"]], surv = s$truth[["pce_surv"]], surv_next = s$truth[["pce_surv_next"]])
  }
  res <- do.call(rbind, rows)
  tg <- list(naive = "surv", sm = "surv", ssm = "surv_next", ec = "surv", ipw = "pop")
  out <- do.call(rbind, lapply(names(tg), function(nm) {
    err <- res[, nm] - res[, tg[[nm]]]
    se <- if (nm == "ipw") NA else res[, paste0("se_", nm)]
    data.frame(estimator = nm, estimand = tg[[nm]], bias = mean(err, na.rm = TRUE), rmse = sqrt(mean(err^2, na.rm = TRUE)),
               coverage95 = if (nm == "ipw") NA else mean(abs(err) <= 1.96 * se, na.rm = TRUE))
  }))
  attr(out, "tests") <- c(reject_T1 = mean(abs(res[, "T1"]) > 1.96, na.rm = TRUE), reject_T2 = mean(abs(res[, "T2"]) > 1.96, na.rm = TRUE))
  out
}
