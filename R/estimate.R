#' Estimate the panel conditioning effect with a refreshment sample
#'
#' Given the continuing cohort's respondents at the comparison wave and the fresh
#' cohort's first-wave respondents, computes five estimators of the panel
#' conditioning effect (PCE) at dose \eqn{k}: the naive contrast, survival matching
#' (SM), symmetric survival matching (SSM), the entry-wave correction (EC), and,
#' when `old$x` is supplied, attrition inverse-probability weighting on `x` (IPW).
#' Also returns the decomposition of the naive contrast, the two diagnostic
#' statistics \eqn{T_1} (SM vs EC; detects non-stationary attrition) and
#' \eqn{T_2} (SSM vs SM; detects state-dependent attrition), and the entry-wave
#' selection differential.
#'
#' @param old data frame for the continuing cohort with columns `y` (outcome at the
#'   comparison wave `t`; `NA` for non-respondents), `y_entry` (outcome at the
#'   cohort's entry wave, observed for every entrant), `s_prior` (1 if the unit
#'   responded at every wave up to and including `t`), `s_next` (1 if the unit also
#'   responds at `t + 1`; may be `NA` if unavailable) and, optionally, `x` (a
#'   discrete covariate for IPW).
#' @param new data frame for the fresh cohort with columns `y` (first-wave outcome),
#'   `s_m` (1 if the unit responds at `t + 1, ..., t + m` with `m = k`) and,
#'   optionally, `s_m1` (1 if the unit responds through `t + k + 1`).
#' @param k number of prior interviews of the continuing cohort at the comparison
#'   wave (the dose).
#' @param m number of matched follow-up waves actually available for SM; defaults
#'   to `k`. When `m < k` the SM estimate is flagged as partial.
#' @return An object of class `pc_estimate`: a list with `estimates` (data frame:
#'   estimator, estimate, se, n_old, n_new), `decomposition` (naive, conditioning =
#'   SM, attrition = naive - SM), `delta_entry` (entry-wave selection differential),
#'   `tests` (T1, T2 with two-sided p-values), and `info`.
#' @details The estimators are described in Okubo (2026), "Panel conditioning as a
#'   dose-response causal effect". All variances are analytic: independent-sample
#'   formulas for the contrasts, the paired formula for EC, and the share-weighted
#'   formulas of Appendix A.7 for the diagnostic statistics.
#' @examples
#' set.seed(1)
#' sim <- pc_simulate(n_old = 2000, n_new = 500, k = 4, regime = "MNAR_trait")
#' fit <- pc_estimate(sim$old, sim$new, k = 4)
#' fit
#' @export
pc_estimate <- function(old, new, k, m = k) {
  stopifnot(is.data.frame(old), is.data.frame(new), all(c("y", "y_entry", "s_prior") %in% names(old)),
            all(c("y", "s_m") %in% names(new)))
  io <- old$s_prior == 1 & !is.na(old$y)
  yo <- old$y[io]
  yn <- new$y[!is.na(new$y)]
  nm <- new$s_m == 1 & !is.na(new$y)
  yn_m <- new$y[nm]
  vm <- function(v) stats::var(v) / length(v)
  naive <- mean(yo) - mean(yn); se_naive <- sqrt(vm(yo) + vm(yn))
  sm <- mean(yo) - mean(yn_m); se_sm <- sqrt(vm(yo) + vm(yn_m))
  ## entry-wave correction
  ent_all <- old$y_entry[!is.na(old$y_entry)]
  ent_s <- old$y_entry[io]
  delta <- mean(ent_s, na.rm = TRUE) - mean(ent_all)
  dpair <- old$y[io] - old$y_entry[io]
  n_s <- sum(!is.na(dpair)); n_all <- length(ent_all)
  cov_term <- 2 * stats::var(ent_s, na.rm = TRUE) * (n_s / n_all) / n_all
  se_ec <- sqrt(max(stats::var(dpair, na.rm = TRUE) / n_s + stats::var(ent_all) / n_all - cov_term, 0) + vm(yn))
  ec <- (mean(yo) - delta) - mean(yn)
  ## symmetric survival matching
  ssm <- NA_real_; se_ssm <- NA_real_; T2 <- NA_real_
  if ("s_next" %in% names(old) && "s_m1" %in% names(new) && any(!is.na(old$s_next)) && any(!is.na(new$s_m1))) {
    is <- io & old$s_next == 1 & !is.na(old$s_next)
    yo_s <- old$y[is]; m1 <- nm & new$s_m1 == 1 & !is.na(new$s_m1); yn_m1 <- new$y[m1]
    ssm <- mean(yo_s) - mean(yn_m1); se_ssm <- sqrt(vm(yo_s) + vm(yn_m1))
    ## T2: SSM vs SM, Appendix A.7
    r_s <- mean(old$s_next[io] == 1, na.rm = TRUE); yo_ns <- old$y[io & old$s_next == 0 & !is.na(old$s_next)]
    v_old2 <- (1 - r_s)^2 * (vm(yo_s) + stats::var(yo_ns) / max(length(yo_ns), 2))
    drop1 <- nm & !m1; q1 <- mean(m1[nm])
    v_fr2 <- (1 - q1)^2 * (stats::var(new$y[drop1]) / max(sum(drop1), 2) + vm(yn_m1))
    T2 <- (ssm - sm) / sqrt(v_old2 + v_fr2)
  }
  ## T1: SM vs EC
  ent_ns <- old$y_entry[!io & !is.na(old$y_entry)]; p_s <- mean(io)
  v_delta <- (1 - p_s)^2 * (vm(ent_s[!is.na(ent_s)]) + stats::var(ent_ns) / max(length(ent_ns), 2))
  q_m <- mean(nm); v_fresh <- (1 - q_m)^2 * (stats::var(new$y[!nm]) / max(sum(!nm), 2) + vm(yn_m))
  T1 <- (sm - ec) / sqrt(v_delta + v_fresh)
  ## IPW on x
  ipw <- NA_real_
  if ("x" %in% names(old)) {
    pS <- tapply(old$s_prior, old$x, mean)[as.character(old$x)]
    w <- 1 / pmax(as.numeric(pS), 0.02)
    ipw <- sum(w[io] * old$y[io]) / sum(w[io]) - mean(yn)
  }
  est <- data.frame(estimator = c("naive", "sm", "ssm", "ec", "ipw"),
                    estimate = c(naive, sm, ssm, ec, ipw),
                    se = c(se_naive, se_sm, se_ssm, se_ec, NA),
                    n_old = c(length(yo), length(yo), if (is.na(ssm)) NA else sum(io & old$s_next == 1, na.rm = TRUE), length(yo), length(yo)),
                    n_new = c(length(yn), length(yn_m), if (is.na(ssm)) NA else sum(nm & new$s_m1 == 1, na.rm = TRUE), length(yn), length(yn)),
                    stringsAsFactors = FALSE)
  out <- list(estimates = est,
              decomposition = c(naive = naive, conditioning_sm = sm, attrition = naive - sm),
              delta_entry = delta,
              tests = data.frame(test = c("T1_sm_vs_ec", "T2_ssm_vs_sm"), statistic = c(T1, T2),
                                 p = 2 * stats::pnorm(-abs(c(T1, T2))), detects = c("non-stationary attrition", "state-dependent attrition")),
              info = list(k = k, m = m, sm_partial = m < k, p_survive = mean(old$s_prior == 1, na.rm = TRUE),
                          mean_old_surv = mean(yo), mean_new = mean(yn)))
  class(out) <- "pc_estimate"
  out
}

#' @export
print.pc_estimate <- function(x, ...) {
  cat("panelcond: refreshment-sample estimates of the conditioning effect (dose k =", x$info$k, ")\n")
  e <- x$estimates; e$estimate <- round(e$estimate, 4); e$se <- round(e$se, 4)
  print(e, row.names = FALSE)
  if (isTRUE(x$info$sm_partial)) cat("  note: SM uses", x$info$m, "<", x$info$k, "matched waves (partial matching)\n")
  cat(sprintf("\ndecomposition: naive %.4f = conditioning (SM) %.4f + attrition %.4f\n",
              x$decomposition["naive"], x$decomposition["conditioning_sm"], x$decomposition["attrition"]))
  cat(sprintf("entry-wave selection differential: %.4f\n", x$delta_entry))
  t <- x$tests; t$statistic <- round(t$statistic, 3); t$p <- signif(t$p, 3)
  cat("\ndiagnostic tests:\n"); print(t, row.names = FALSE)
  invisible(x)
}

#' Worst-case bounds and breakdown values
#'
#' Worst-case (Manski-type) bounds on the population conditioning effect that
#' make no assumption about attrition, and breakdown values for the survival-
#' matched and entry-wave estimators.
#'
#' @param fit a `pc_estimate` object.
#' @param y_range numeric length-2 support of the outcome (e.g. `c(0, 1)` for a
#'   proportion, `c(0, 100)` for a feeling thermometer).
#' @param sd_new standard deviation of the fresh cohort's outcome, used to express
#'   the SM breakdown value in standard-deviation units; optional.
#' @return A list with `bounds` (lower, upper bound on the population PCE), and
#'   `breakdown`: `sm_std` (the SM estimate in SD units, which equals the violation
#'   of attrition exchangeability needed to explain it away) and `ec_rho_star`
#'   (the ratio of the wave-`t` to the entry-wave selection differential at which
#'   the EC estimate vanishes; values far from 1 are robust).
#' @export
pc_bounds <- function(fit, y_range, sd_new = NULL) {
  stopifnot(inherits(fit, "pc_estimate"), length(y_range) == 2)
  p <- fit$info$p_survive; mo <- fit$info$mean_old_surv; mn <- fit$info$mean_new
  e <- fit$estimates
  sm <- e$estimate[e$estimator == "sm"]; naive <- e$estimate[e$estimator == "naive"]
  list(bounds = c(lower = p * mo + (1 - p) * y_range[1] - mn, upper = p * mo + (1 - p) * y_range[2] - mn),
       breakdown = c(sm_std = if (is.null(sd_new)) NA_real_ else sm / sd_new,
                     ec_rho_star = if (abs(fit$delta_entry) > 1e-12) naive / fit$delta_entry else NA_real_))
}
