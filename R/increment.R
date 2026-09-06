#' Incremental dose contrast between two continuing cohorts
#'
#' When two cohorts were drawn from the same target population at different
#' entry waves (for example an original sample and a cohort-matched refreshment
#' sample), their survivors at a later wave carry different doses \eqn{k_A > k_B}.
#' Applying the entry-wave correction to each cohort with its own entry wave
#' identifies the increment \eqn{\tau^{S_A}(k_A) - \tau^{S_B}(k_B)} without any
#' fresh cohort at the comparison wave. This is the design used when a later
#' refreshment sample targets a different population and cannot serve as a
#' dose-0 control.
#'
#' @param A,B data frames with columns `y` (outcome at the comparison wave;
#'   `NA` for non-respondents), `y_entry` (outcome at the cohort's own entry wave,
#'   observed for all entrants; `NA` for item non-response) and `s` (1 if the unit
#'   responded at every wave up to and including the comparison wave).
#' @param kA,kB doses of the two cohorts at the comparison wave, `kA > kB`; used
#'   for labelling only.
#' @param nboot number of replications of the joint person bootstrap (persons
#'   resampled within each cohort); `0` gives analytic standard errors only.
#' @return An object of class `pc_increment`: a list with `naive_diff` (difference
#'   of survivors' means, with SE and p), `ec_diff` (entry-corrected increment,
#'   with SE and p), the two selection differentials `delta`, and `info`. With
#'   `nboot > 0` each contrast also carries `se_boot` and `p_boot`.
#' @details The analytic SE of the increment adds the two cohorts' first-order
#'   variances of \eqn{\bar Y_t^S - \bar Y_c^S + \bar Y_c} (see [pc_estimate()]);
#'   there is no fresh-arm term.
#' @examples
#' set.seed(2)
#' a <- pc_simulate(n_old = 3000, n_new = 10, k = 12, regime = "MNAR_trait")$old
#' b <- pc_simulate(n_old = 1000, n_new = 10, k = 8,  regime = "MNAR_trait")$old
#' pc_increment(data.frame(y = a$y, y_entry = a$y_entry, s = a$s_prior),
#'              data.frame(y = b$y, y_entry = b$y_entry, s = b$s_prior), kA = 12, kB = 8)
#' @export
pc_increment <- function(A, B, kA, kB, nboot = 0) {
  stopifnot("`A` and `B` must be data frames" = is.data.frame(A) && is.data.frame(B),
            "`A` needs columns y, y_entry, s" = all(c("y", "y_entry", "s") %in% names(A)),
            "`B` needs columns y, y_entry, s" = all(c("y", "y_entry", "s") %in% names(B)),
            "`kA` and `kB` must be single non-negative integers" = all(vapply(list(kA, kB), function(v) is.numeric(v) && length(v) == 1 && v >= 0 && v == round(v), logical(1))),
            "`nboot` must be a single non-negative integer" = is.numeric(nboot) && length(nboot) == 1 && nboot >= 0 && nboot == round(nboot))
  if (kA <= kB) warning("`kA` should exceed `kB`; the contrast is reported as A minus B", call. = FALSE)
  prep <- function(D, what) {
    D$y <- check_outcome(D$y, paste0(what, "$y")); D$y_entry <- check_outcome(D$y_entry, paste0(what, "$y_entry"))
    D$s <- check_indicator(D$s, paste0(what, "$s"))
    if (!any(D$s == 1L & !is.na(D$y))) stop(sprintf("cohort %s has no survivor with an observed answer", what), call. = FALSE)
    n_drop <- c(sum(D$s == 1L & is.na(D$y)), sum(is.na(D$y_entry)))
    if (sum(n_drop) > 0) message(sprintf("cohort %s: item non-response dropped: %d survivors without y, %d entrants without y_entry", what, n_drop[1], n_drop[2]))
    D
  }
  A <- prep(A, "A"); B <- prep(B, "B")
  one <- function(D) {
    i <- D$s == 1L & !is.na(D$y); eobs <- !is.na(D$y_entry)
    y <- D$y[i]; e_s <- D$y_entry[i & eobs]; e_ns <- D$y_entry[!i & eobs]
    delta <- mean(e_s) - mean(c(e_s, e_ns))
    list(m = mean(y), v_m = if (length(y) >= 2) stats::var(y) / length(y) else NA_real_, delta = delta,
         v_ec = ec_var_cont(D$y[i & eobs], e_s, e_ns), n = length(y), n_pairs = length(e_s))
  }
  point <- function(A, B) { a <- one(A); b <- one(B); c(naive = a$m - b$m, ec = (a$m - a$delta) - (b$m - b$delta)) }
  a <- one(A); b <- one(B)
  naive <- a$m - b$m; se_naive <- sqrt(a$v_m + b$v_m)
  ec <- (a$m - a$delta) - (b$m - b$delta); se_ec <- sqrt(a$v_ec + b$v_ec)
  out <- list(naive_diff = c(estimate = naive, se = se_naive, p = 2 * stats::pnorm(-abs(naive / se_naive))),
              ec_diff = c(estimate = ec, se = se_ec, p = 2 * stats::pnorm(-abs(ec / se_ec))),
              delta = c(A = a$delta, B = b$delta),
              info = list(kA = kA, kB = kB, nA = a$n, nB = b$n, pairsA = a$n_pairs, pairsB = b$n_pairs, nboot = nboot))
  if (nboot > 0) {
    bm <- t(replicate(nboot, point(A[sample.int(nrow(A), replace = TRUE), , drop = FALSE], B[sample.int(nrow(B), replace = TRUE), , drop = FALSE])))
    sb <- apply(bm, 2, stats::sd, na.rm = TRUE)
    out$naive_diff <- c(out$naive_diff, se_boot = unname(sb["naive"]), p_boot = 2 * stats::pnorm(-abs(naive / unname(sb["naive"]))))
    out$ec_diff <- c(out$ec_diff, se_boot = unname(sb["ec"]), p_boot = 2 * stats::pnorm(-abs(ec / unname(sb["ec"]))))
  }
  class(out) <- "pc_increment"
  out
}

#' @export
print.pc_increment <- function(x, ...) {
  cat(sprintf("panelcond: incremental dose contrast, k = %s vs %s (n = %d, %d)\n", format(x$info$kA), format(x$info$kB), x$info$nA, x$info$nB))
  line <- function(lab, v) {
    cat(sprintf("  %s %.4f (se %.4f, p = %.3g)", lab, v["estimate"], v["se"], v["p"]))
    if (!is.na(v["se_boot"])) cat(sprintf("; bootstrap se %.4f, p = %.3g", v["se_boot"], v["p_boot"]))
    cat("\n")
  }
  line("naive difference of survivors:", x$naive_diff)
  line("entry-corrected increment:     ", x$ec_diff)
  cat(sprintf("  selection differentials: A %.4f (%d pairs), B %.4f (%d pairs)\n", x$delta["A"], x$info$pairsA, x$delta["B"], x$info$pairsB))
  if (x$info$nboot > 0) cat(sprintf("  bootstrap: persons resampled within cohort, %d replications\n", x$info$nboot))
  invisible(x)
}
