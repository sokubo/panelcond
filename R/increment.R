#' Incremental dose contrast between two continuing cohorts
#'
#' When two cohorts were drawn from the same target population at different
#' entry waves (for example an original sample and a cohort-matched refreshment
#' sample), their survivors at a later wave carry different doses \eqn{k_A > k_B}.
#' Applying the entry-wave correction to each cohort with its own entry wave
#' identifies the increment \eqn{\tau^S(k_A) - \tau^S(k_B)} without any fresh
#' cohort at the comparison wave. This is the design used when a later
#' refreshment sample targets a different population and cannot serve as a
#' dose-0 control.
#'
#' @param A,B data frames with columns `y` (outcome at the comparison wave;
#'   `NA` for non-respondents), `y_entry` (outcome at the cohort's own entry wave,
#'   observed for all entrants) and `s` (1 if the unit responded at every wave up
#'   to and including the comparison wave).
#' @param kA,kB doses of the two cohorts at the comparison wave.
#' @return A list with `naive_diff` (difference of survivors' means, with SE),
#'   `ec_diff` (entry-corrected increment, with SE), the two selection
#'   differentials, and `info`.
#' @examples
#' set.seed(2)
#' a <- pc_simulate(n_old = 3000, n_new = 10, k = 12, regime = "MNAR_trait")$old
#' b <- pc_simulate(n_old = 1000, n_new = 10, k = 8,  regime = "MNAR_trait")$old
#' pc_increment(data.frame(y = a$y, y_entry = a$y_entry, s = a$s_prior),
#'              data.frame(y = b$y, y_entry = b$y_entry, s = b$s_prior), kA = 12, kB = 8)
#' @export
pc_increment <- function(A, B, kA, kB) {
  stopifnot(all(c("y", "y_entry", "s") %in% names(A)), all(c("y", "y_entry", "s") %in% names(B)))
  one <- function(D) {
    i <- D$s == 1 & !is.na(D$y)
    y <- D$y[i]; e_s <- D$y_entry[i]; e_all <- D$y_entry[!is.na(D$y_entry)]
    delta <- mean(e_s, na.rm = TRUE) - mean(e_all)
    n_s <- sum(!is.na(y - e_s))
    v <- stats::var(y - e_s, na.rm = TRUE) / n_s + stats::var(e_all) / length(e_all)
    list(m = mean(y), v_m = stats::var(y) / length(y), delta = delta, v_ec = v, n = length(y))
  }
  a <- one(A); b <- one(B)
  naive <- a$m - b$m; se_naive <- sqrt(a$v_m + b$v_m)
  ec <- (a$m - a$delta) - (b$m - b$delta); se_ec <- sqrt(a$v_ec + b$v_ec)
  out <- list(naive_diff = c(estimate = naive, se = se_naive, p = 2 * stats::pnorm(-abs(naive / se_naive))),
              ec_diff = c(estimate = ec, se = se_ec, p = 2 * stats::pnorm(-abs(ec / se_ec))),
              delta = c(A = a$delta, B = b$delta), info = list(kA = kA, kB = kB, nA = a$n, nB = b$n))
  class(out) <- "pc_increment"; out
}

#' @export
print.pc_increment <- function(x, ...) {
  cat(sprintf("panelcond: incremental dose contrast, k = %d vs %d (n = %d, %d)\n", x$info$kA, x$info$kB, x$info$nA, x$info$nB))
  cat(sprintf("  naive difference of survivors: %.4f (se %.4f, p = %.3g)\n", x$naive_diff["estimate"], x$naive_diff["se"], x$naive_diff["p"]))
  cat(sprintf("  entry-corrected increment:      %.4f (se %.4f, p = %.3g)\n", x$ec_diff["estimate"], x$ec_diff["se"], x$ec_diff["p"]))
  cat(sprintf("  selection differentials: A %.4f, B %.4f\n", x$delta["A"], x$delta["B"]))
  invisible(x)
}
