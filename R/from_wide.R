#' Build estimation inputs from a wide panel
#'
#' Converts a wide data frame (one row per person) with per-wave response
#' indicators into the `old` and `new` inputs of [pc_estimate()].
#'
#' @param data wide data frame.
#' @param cohort name of the column giving the entry wave of each person (integer
#'   wave number).
#' @param resp character vector of column names of response indicators (0/1),
#'   one per wave in wave order, e.g. `paste0("r_w", 1:19)`.
#' @param y_cols named character vector of outcome columns by wave, e.g.
#'   `c("5" = "dq10_1", "1" = "zq13_1")`; must contain the comparison wave and the
#'   continuing cohort's entry wave.
#' @param t comparison wave.
#' @param c_old entry wave of the continuing cohort.
#' @param c_new entry wave of the fresh cohort (normally `t`).
#' @param x optional name of a discrete covariate column for IPW.
#' @return A list with `old`, `new`, `k = t - c_old` and `m` (the number of
#'   follow-up waves available for the fresh cohort, capped at `k`).
#' @export
pc_from_wide <- function(data, cohort, resp, y_cols, t, c_old, c_new = t, x = NULL) {
  stopifnot(cohort %in% names(data), all(resp %in% names(data)))
  R <- as.matrix(data[, resp]); R[is.na(R)] <- 0L
  surv <- function(from, to) as.integer(rowSums(R[, from:to, drop = FALSE]) == (to - from + 1))
  k <- t - c_old
  last <- max(which(colSums(R) > 0))
  m <- min(k, last - t)
  yt <- data[[y_cols[[as.character(t)]]]]; ye <- data[[y_cols[[as.character(c_old)]]]]
  old_i <- data[[cohort]] == c_old & !is.na(data[[cohort]])
  new_i <- data[[cohort]] == c_new & !is.na(data[[cohort]])
  old <- data.frame(y = ifelse(surv(c_old, t)[old_i] == 1, yt[old_i], NA), y_entry = ye[old_i],
                    s_prior = surv(c_old, t)[old_i],
                    s_next = if (t + 1 <= ncol(R)) R[old_i, t + 1] else NA_integer_)
  if (!is.null(x)) old$x <- data[[x]][old_i]
  new <- data.frame(y = ifelse(R[new_i, t] == 1, yt[new_i], NA),
                    s_m = if (m >= 1) surv(t, t + m)[new_i] else R[new_i, t],
                    s_m1 = if (t + k + 1 <= ncol(R)) surv(t, t + k + 1)[new_i] else NA_integer_)
  list(old = old, new = new, k = k, m = m)
}
