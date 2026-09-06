#' Estimate the panel conditioning effect with a refreshment sample
#'
#' Given the continuing cohort's respondents at the comparison wave and the fresh
#' cohort's first-wave respondents, computes five estimators of the panel
#' conditioning effect (PCE) at dose \eqn{k}: the naive contrast, survival matching
#' (SM), symmetric survival matching (SSM), the entry-wave correction (EC), and,
#' when `old$x` is supplied, attrition inverse-probability weighting on `x` (IPW).
#' With `adjust` the regression-standardised entry-wave correction (EC-adj) is
#' added. Also returns the decomposition of the naive contrast, the two diagnostic
#' statistics \eqn{T_1} (SM vs EC; detects non-stationary attrition) and
#' \eqn{T_2} (SSM vs SM; detects state-dependent attrition), and the entry-wave
#' selection differential.
#'
#' @param old data frame for the continuing cohort with columns `y` (outcome at the
#'   comparison wave `t`; `NA` for non-respondents), `y_entry` (outcome at the
#'   cohort's entry wave, observed for every entrant; `NA` for item non-response),
#'   `s_prior` (1 if the unit responded at every wave up to and including `t`),
#'   `s_next` (1 if the unit also responds at `t + 1`; may be `NA` if unavailable)
#'   and, optionally, `x` (a discrete covariate for IPW).
#' @param new data frame for the fresh cohort with columns `y` (first-wave outcome;
#'   `NA` for item non-response), `s_m` (1 if the unit responds at
#'   `t + 1, ..., t + m` with `m = k`) and, optionally, `s_m1` (1 if the unit
#'   responds through `t + k + 1`; may be `NA` if unavailable).
#' @param k number of prior interviews of the continuing cohort at the comparison
#'   wave (the dose); a positive integer.
#' @param m number of matched follow-up waves actually available for SM; defaults
#'   to `k`. When `m < k` the SM estimate is flagged as partial.
#' @param adjust optional character vector naming numeric covariate columns present
#'   in both `old` and `new` (for example sex, birth year, education at entry). The
#'   entry-wave differential and the fresh-cohort mean are then standardised by
#'   linear regression to the survivors' covariate distribution (EC-adj). Rows with
#'   a missing covariate are not allowed.
#' @param nboot number of replications of the joint person bootstrap (persons
#'   resampled with replacement within each cohort; every estimator and both
#'   diagnostic differences recomputed in each replicate). `0` (the default) gives
#'   analytic standard errors only. Call `set.seed()` beforehand for reproducibility.
#' @return An object of class `pc_estimate`: a list with `estimates` (data frame:
#'   estimator, estimate, se, n_old, n_new, and `se_boot` when `nboot > 0`),
#'   `decomposition` (naive, conditioning = SM, attrition = naive - SM),
#'   `delta_entry` (entry-wave selection differential), `tests` (T1, T2 with
#'   two-sided p-values from the analytic denominators and, when `nboot > 0`, from
#'   the bootstrap denominators), and `info` (k, m, whether SM is partial, the share
#'   of the continuing cohort with an observed wave-`t` answer, the arm means, the
#'   number of entry-wave pairs and the counts of rows dropped for item
#'   non-response). IPW and EC-adj have no analytic standard error; use `nboot`.
#' @details The estimators are those of Okubo (2026), "Panel conditioning as a
#'   dose-response causal effect". Analytic variances are first-order:
#'   independent-sample formulas for the contrasts; for EC the variance of
#'   \eqn{\bar Y_t^S - (1-p)\bar Y_c^S + (1-p)\bar Y_c^{NS} - \bar Y_t^0}, with
#'   \eqn{p} the survivors' share, which keeps the within-person covariance of the
#'   survivors' two answers; and share-weighted formulas over disjoint groups for
#'   the diagnostic statistics. IPW uses the survival share within each level of
#'   `x`, floored at 0.02. Rows of `old` with `s_prior = 1` but `y` missing are
#'   treated as non-survivors; rows with `y_entry` missing enter neither the
#'   differential nor its variance. A message reports these counts.
#' @examples
#' set.seed(1)
#' sim <- pc_simulate(n_old = 2000, n_new = 500, k = 4, regime = "MNAR_trait")
#' fit <- pc_estimate(sim$old, sim$new, k = 4)
#' fit
#' ## joint person bootstrap and covariate adjustment on x
#' pc_estimate(sim$old, sim$new, k = 4, adjust = "x", nboot = 50)
#' @export
pc_estimate <- function(old, new, k, m = k, adjust = NULL, nboot = 0) {
  stopifnot("`old` must be a data frame" = is.data.frame(old),
            "`new` must be a data frame" = is.data.frame(new),
            "`old` needs columns y, y_entry, s_prior" = all(c("y", "y_entry", "s_prior") %in% names(old)),
            "`new` needs columns y, s_m" = all(c("y", "s_m") %in% names(new)),
            "`k` must be a single positive integer" = is.numeric(k) && length(k) == 1 && k >= 1 && k == round(k),
            "`m` must be a single integer between 0 and k" = is.numeric(m) && length(m) == 1 && m >= 0 && m <= k && m == round(m),
            "`nboot` must be a single non-negative integer" = is.numeric(nboot) && length(nboot) == 1 && nboot >= 0 && nboot == round(nboot))
  if (nrow(old) == 0 || nrow(new) == 0) stop("`old` and `new` must each have at least one row", call. = FALSE)
  old$y <- check_outcome(old$y, "old$y"); old$y_entry <- check_outcome(old$y_entry, "old$y_entry"); new$y <- check_outcome(new$y, "new$y")
  old$s_prior <- check_indicator(old$s_prior, "old$s_prior"); new$s_m <- check_indicator(new$s_m, "new$s_m")
  if ("s_next" %in% names(old)) old$s_next <- check_indicator(old$s_next, "old$s_next", allow_na = TRUE)
  if ("s_m1" %in% names(new)) new$s_m1 <- check_indicator(new$s_m1, "new$s_m1", allow_na = TRUE)
  if (!is.null(adjust)) {
    stopifnot("`adjust` must name columns present in both `old` and `new`" = is.character(adjust) && all(adjust %in% names(old)) && all(adjust %in% names(new)))
    for (v in adjust) {
      if (!(is.numeric(old[[v]]) || is.logical(old[[v]])) || !(is.numeric(new[[v]]) || is.logical(new[[v]]))) stop(sprintf("covariate `%s` must be numeric or logical in both cohorts", v), call. = FALSE)
      if (anyNA(old[[v]]) || anyNA(new[[v]])) stop(sprintf("covariate `%s` has missing values; impute or drop those rows first", v), call. = FALSE)
    }
    for (D in list(old, new)) if (qr(cbind(1, as.matrix(D[adjust]) * 1))$rank < length(adjust) + 1)
      warning("`adjust` covariates are constant or collinear within a cohort; the redundant terms are dropped from that cohort's regression", call. = FALSE)
  }
  io <- old$s_prior == 1L & !is.na(old$y)
  if (!any(io)) stop("no continuing-cohort respondent with an observed wave-t answer (s_prior = 1 and y not NA)", call. = FALSE)
  if (!any(!is.na(new$y))) stop("no fresh-cohort respondent with an observed answer", call. = FALSE)
  if (!any(new$s_m == 1L & !is.na(new$y))) warning("no fresh-cohort respondent with s_m = 1: the survival-matched estimate is NA", call. = FALSE)
  n_item_old <- sum(old$s_prior == 1L & is.na(old$y)); n_item_new <- sum(is.na(new$y)); n_item_entry <- sum(is.na(old$y_entry))
  if (n_item_old + n_item_new + n_item_entry > 0)
    message(sprintf("item non-response dropped: %d continuing survivors without y, %d fresh respondents without y, %d entrants without y_entry",
                    n_item_old, n_item_new, n_item_entry))
  if ("x" %in% names(old)) {
    if (anyNA(old$x)) message(sprintf("%d rows of `old` have x = NA and are excluded from IPW", sum(is.na(old$x))))
    pS <- tapply(old$s_prior[!is.na(old$x)], old$x[!is.na(old$x)], mean)
    if (any(pS < 0.02)) warning(sprintf("survival share below 0.02 in %d level(s) of x; IPW weights are floored at 1/0.02", sum(pS < 0.02)), call. = FALSE)
  }
  r <- pc_point(old, new, adjust)
  est_names <- c("naive", "sm", "ssm", "ec", if (!is.null(adjust)) "ec_adj", "ipw")
  est <- data.frame(estimator = est_names, estimate = r$est[est_names], se = r$se[est_names],
                    n_old = r$n_old[est_names], n_new = r$n_new[est_names], stringsAsFactors = FALSE, row.names = NULL)
  tests <- data.frame(test = c("T1_sm_vs_ec", "T2_ssm_vs_sm"), statistic = c(r$T1, r$T2),
                      p = 2 * stats::pnorm(-abs(c(r$T1, r$T2))),
                      detects = c("non-stationary attrition", "state-dependent attrition"), stringsAsFactors = FALSE)
  if (nboot > 0) {
    bm <- matrix(NA_real_, nboot, length(est_names) + 2, dimnames = list(NULL, c(est_names, "d1", "d2")))
    for (b in seq_len(nboot)) {
      rb <- pc_point(old[sample.int(nrow(old), replace = TRUE), , drop = FALSE], new[sample.int(nrow(new), replace = TRUE), , drop = FALSE], adjust)
      bm[b, ] <- c(rb$est[est_names], rb$est["sm"] - rb$est["ec"], rb$est["ssm"] - rb$est["sm"])
    }
    sdb <- function(v) if (sum(!is.na(v)) >= 2) stats::sd(v, na.rm = TRUE) else NA_real_
    est$se_boot <- apply(bm[, est_names, drop = FALSE], 2, sdb)
    tests$se_boot <- c(sdb(bm[, "d1"]), sdb(bm[, "d2"]))
    tests$statistic_boot <- c(r$est["sm"] - r$est["ec"], r$est["ssm"] - r$est["sm"]) / tests$se_boot
    tests$p_boot <- 2 * stats::pnorm(-abs(tests$statistic_boot))
  }
  out <- list(estimates = est,
              decomposition = c(naive = unname(r$est["naive"]), conditioning_sm = unname(r$est["sm"]), attrition = unname(r$est["naive"] - r$est["sm"])),
              delta_entry = r$delta,
              tests = tests,
              info = list(k = k, m = m, sm_partial = m < k, p_survive = mean(io), mean_old_surv = r$mean_old, mean_new = r$mean_new,
                          n_pairs = r$n_pairs, dropped = c(old_y = n_item_old, new_y = n_item_new, y_entry = n_item_entry), nboot = nboot))
  class(out) <- "pc_estimate"
  out
}

## point estimates, analytic standard errors and diagnostic statistics; no checks, no messages
pc_point <- function(old, new, adjust = NULL) {
  io <- old$s_prior == 1L & !is.na(old$y); yo <- old$y[io]
  yobs <- !is.na(new$y); yn <- new$y[yobs]
  nm <- yobs & new$s_m == 1L; yn_m <- new$y[nm]
  mean0 <- function(v) if (length(v)) mean(v) else NA_real_
  var0 <- function(v) if (length(v) >= 2) stats::var(v) else NA_real_
  vm <- function(v) var0(v) / length(v)
  naive <- mean0(yo) - mean0(yn); se_naive <- sqrt(vm(yo) + vm(yn))
  sm <- mean0(yo) - mean0(yn_m); se_sm <- sqrt(vm(yo) + vm(yn_m))
  ## entry-wave correction; p_s is the survivors' share of the entrants with an observed entry answer
  eobs <- !is.na(old$y_entry)
  ent_s <- old$y_entry[io & eobs]; ent_ns <- old$y_entry[!io & eobs]
  n_s <- length(ent_s); n_ns <- length(ent_ns); p_s <- n_s / (n_s + n_ns)
  delta <- mean0(ent_s) - mean0(c(ent_s, ent_ns))
  ec <- (mean0(yo) - delta) - mean0(yn)
  se_ec <- sqrt(ec_var_cont(old$y[io & eobs], ent_s, ent_ns) + vm(yn))
  ## symmetric survival matching and T2 (SSM vs SM)
  ssm <- se_ssm <- T2 <- NA_real_; n_ssm_old <- n_ssm_new <- NA_integer_
  if ("s_next" %in% names(old) && "s_m1" %in% names(new) && any(!is.na(old$s_next)) && any(!is.na(new$s_m1))) {
    io1 <- io & !is.na(old$s_next) & old$s_next == 1L; yo_s <- old$y[io1]
    m1 <- nm & !is.na(new$s_m1) & new$s_m1 == 1L; yn_m1 <- new$y[m1]
    ssm <- mean0(yo_s) - mean0(yn_m1); se_ssm <- sqrt(vm(yo_s) + vm(yn_m1))
    n_ssm_old <- length(yo_s); n_ssm_new <- length(yn_m1)
    yo_ns <- old$y[io & !is.na(old$s_next) & old$s_next == 0L]
    r_s <- length(yo_s) / (length(yo_s) + length(yo_ns))
    v_old2 <- (1 - r_s)^2 * (vm(yo_s) + var0(yo_ns) / max(length(yo_ns), 2))
    y_d1 <- new$y[nm & !m1]; q1 <- sum(m1) / sum(nm)
    v_fr2 <- (1 - q1)^2 * (var0(y_d1) / max(length(y_d1), 2) + vm(yn_m1))
    T2 <- (ssm - sm) / sqrt(v_old2 + v_fr2)
  }
  ## T1 (SM vs EC); the fresh shares are taken among respondents with an observed answer
  v_delta <- (1 - p_s)^2 * (var0(ent_s) / max(n_s, 2) + var0(ent_ns) / max(n_ns, 2))
  y_nm <- new$y[yobs & !nm]; q_m <- sum(nm) / sum(yobs)
  v_fresh <- (1 - q_m)^2 * (var0(y_nm) / max(length(y_nm), 2) + vm(yn_m))
  T1 <- (sm - ec) / sqrt(v_delta + v_fresh)
  ## IPW on x: survivors reweighted by the inverse survival share of their x level
  ipw <- NA_real_
  if ("x" %in% names(old)) {
    xok <- !is.na(old$x)
    pS <- tapply(old$s_prior[xok], old$x[xok], mean)
    w <- 1 / pmax(as.numeric(pS[as.character(old$x)]), 0.02)
    iw <- io & xok
    ipw <- sum(w[iw] * old$y[iw]) / sum(w[iw]) - mean0(yn)
  }
  ## EC-adj: regression standardisation of the entry differential and the fresh mean to the survivors' covariates
  ec_adj <- NA_real_
  if (!is.null(adjust)) {
    Xo <- cbind(1, as.matrix(old[adjust]) * 1); Xn <- cbind(1, as.matrix(new[adjust]) * 1)
    ip <- io & eobs
    m_c <- ols_pred(old$y_entry, Xo, Xo[ip, , drop = FALSE])
    m_0 <- ols_pred(new$y, Xn, Xo[io, , drop = FALSE])
    ec_adj <- mean0(yo) - (mean0(ent_s) - m_c) - m_0
  }
  list(est = c(naive = naive, sm = sm, ssm = ssm, ec = ec, ec_adj = ec_adj, ipw = ipw),
       se = c(naive = se_naive, sm = se_sm, ssm = se_ssm, ec = se_ec, ec_adj = NA_real_, ipw = NA_real_),
       n_old = c(naive = length(yo), sm = length(yo), ssm = n_ssm_old, ec = length(yo), ec_adj = length(yo), ipw = length(yo)),
       n_new = c(naive = length(yn), sm = length(yn_m), ssm = n_ssm_new, ec = length(yn), ec_adj = length(yn), ipw = length(yn)),
       T1 = unname(T1), T2 = unname(T2), delta = delta, mean_old = mean0(yo), mean_new = mean0(yn), n_pairs = n_s)
}

## first-order variance of the continuing-cohort part of the entry-wave estimator,
## Ybar_t^S - Ybar_c^S + Ybar_c = Ybar_t^S - (1 - p) Ybar_c^S + (1 - p) Ybar_c^NS,
## from the survivors' pairs (y_s, e_s) and the non-survivors' entry answers e_ns
ec_var_cont <- function(y_s, e_s, e_ns) {
  ok <- !is.na(y_s) & !is.na(e_s); ys <- y_s[ok]; es <- e_s[ok]; ens <- e_ns[!is.na(e_ns)]
  n_s <- length(ys); n_ns <- length(ens); p <- n_s / (n_s + n_ns)
  if (n_s < 2) return(NA_real_)
  v_s <- (stats::var(ys) - 2 * (1 - p) * stats::cov(ys, es) + (1 - p)^2 * stats::var(es)) / n_s
  v_ns <- if (n_ns >= 2) (1 - p)^2 * stats::var(ens) / n_ns else 0
  v_s + v_ns
}

## least-squares fit of y on Xfit (complete cases), averaged prediction over the rows of Xpred
ols_pred <- function(y, Xfit, Xpred) {
  ok <- !is.na(y)
  if (sum(ok) <= ncol(Xfit) || nrow(Xpred) == 0) return(NA_real_)
  b <- stats::lm.fit(Xfit[ok, , drop = FALSE], y[ok])$coefficients; b[is.na(b)] <- 0
  mean(Xpred %*% b)
}

check_outcome <- function(v, what) {
  if (is.logical(v)) v <- as.numeric(v)
  if (!is.numeric(v)) stop(sprintf("`%s` must be numeric or logical; recode factors and character vectors first", what), call. = FALSE)
  v
}

check_indicator <- function(v, what, allow_na = FALSE) {
  if (is.logical(v)) v <- as.integer(v)
  if (!is.numeric(v) || !all(v[!is.na(v)] %in% c(0, 1))) stop(sprintf("`%s` must be a 0/1 indicator", what), call. = FALSE)
  if (!allow_na && anyNA(v)) stop(sprintf("`%s` contains NA; code a unit with no response as 0", what), call. = FALSE)
  as.integer(v)
}

#' @export
print.pc_estimate <- function(x, ...) {
  cat(sprintf("panelcond: refreshment-sample estimates of the conditioning effect (dose k = %s)\n", format(x$info$k)))
  e <- x$estimates; e$estimate <- round(e$estimate, 4); e$se <- round(e$se, 4)
  if (!is.null(e$se_boot)) e$se_boot <- round(e$se_boot, 4)
  print(e, row.names = FALSE)
  if (isTRUE(x$info$sm_partial)) cat("  note: SM uses", x$info$m, "<", x$info$k, "matched waves (partial matching)\n")
  if (x$info$nboot > 0) cat(sprintf("  se_boot: joint person bootstrap, %d replications\n", x$info$nboot))
  cat(sprintf("\ndecomposition: naive %.4f = conditioning (SM) %.4f + attrition %.4f\n",
              x$decomposition["naive"], x$decomposition["conditioning_sm"], x$decomposition["attrition"]))
  cat(sprintf("entry-wave selection differential: %.4f (%d pairs)\n", x$delta_entry, x$info$n_pairs))
  tt <- x$tests; tt$statistic <- round(tt$statistic, 3); tt$p <- signif(tt$p, 3)
  if (!is.null(tt$se_boot)) { tt$se_boot <- round(tt$se_boot, 4); tt$statistic_boot <- round(tt$statistic_boot, 3); tt$p_boot <- signif(tt$p_boot, 3) }
  cat("\ndiagnostic tests:\n"); print(tt, row.names = FALSE)
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
#'   the SM breakdown value in standard-deviation units; optional, positive.
#' @return A list with `bounds` (lower, upper bound on the population PCE, obtained
#'   by assigning the continuing-cohort members without an observed wave-`t`
#'   answer the endpoints of `y_range`) and `breakdown`: `sm_std` (the SM estimate
#'   in SD units, which equals the violation of attrition exchangeability, in SD
#'   units, needed to explain it away), `ec_eta_star` (the additive drift in the
#'   selection differential between entry and `t` at which the EC estimate
#'   vanishes; equal to the estimate itself), `ec_eta_ci95` (the drift at which the
#'   95\% interval reaches zero; uses the bootstrap SE when available) and
#'   `ec_rho_star` (the ratio of the wave-`t` to the entry-wave selection
#'   differential at which the EC estimate vanishes; the further from 1, the
#'   larger the proportional change needed).
#' @examples
#' set.seed(1)
#' sim <- pc_simulate(n_old = 2000, n_new = 500, k = 4)
#' fit <- pc_estimate(sim$old, sim$new, k = 4)
#' pc_bounds(fit, y_range = c(-4, 4), sd_new = sd(sim$new$y))
#' @export
pc_bounds <- function(fit, y_range, sd_new = NULL) {
  stopifnot("`fit` must be a pc_estimate object" = inherits(fit, "pc_estimate"),
            "`y_range` must be two finite numbers" = is.numeric(y_range) && length(y_range) == 2 && all(is.finite(y_range)),
            "`sd_new` must be a single positive number" = is.null(sd_new) || (is.numeric(sd_new) && length(sd_new) == 1 && is.finite(sd_new) && sd_new > 0))
  y_range <- sort(y_range)
  p <- fit$info$p_survive; mo <- fit$info$mean_old_surv; mn <- fit$info$mean_new
  if (mo < y_range[1] || mo > y_range[2] || mn < y_range[1] || mn > y_range[2]) warning("an arm mean lies outside `y_range`", call. = FALSE)
  e <- fit$estimates
  sm <- e$estimate[e$estimator == "sm"]; naive <- e$estimate[e$estimator == "naive"]; ec <- e$estimate[e$estimator == "ec"]
  se_ec <- if (!is.null(e$se_boot) && !is.na(e$se_boot[e$estimator == "ec"])) e$se_boot[e$estimator == "ec"] else e$se[e$estimator == "ec"]
  list(bounds = c(lower = p * mo + (1 - p) * y_range[1] - mn, upper = p * mo + (1 - p) * y_range[2] - mn),
       breakdown = c(sm_std = if (is.null(sd_new)) NA_real_ else sm / sd_new,
                     ec_eta_star = ec,
                     ec_eta_ci95 = ec - sign(ec) * 1.96 * se_ec,
                     ec_rho_star = if (abs(fit$delta_entry) > 1e-12) naive / fit$delta_entry else NA_real_))
}
