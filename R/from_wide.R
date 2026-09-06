#' Build estimation inputs from a wide panel
#'
#' Converts a wide data frame (one row per person) with per-wave response
#' indicators into the `old` and `new` inputs of [pc_estimate()].
#'
#' @param data wide data frame, one row per person.
#' @param cohort name of the column giving the entry wave of each person (integer
#'   wave number). Rows whose entry wave is neither `c_old` nor `c_new` are
#'   ignored, as are fresh-cohort rows without a response at `t`.
#' @param resp character vector of column names of response indicators (0/1;
#'   `NA` is read as 0), one per wave in wave order, e.g. `paste0("r_w", 1:19)`.
#' @param y_cols named character vector of outcome columns by wave, e.g.
#'   `c("1" = "y_w1", "5" = "y_w5")`; must contain the comparison wave and the
#'   continuing cohort's entry wave.
#' @param t comparison wave.
#' @param c_old entry wave of the continuing cohort, `c_old < t`.
#' @param c_new entry wave of the fresh cohort (normally `t`).
#' @param x optional name of a discrete covariate column for IPW.
#' @param id optional name of a person identifier column; duplicated identifiers
#'   are an error.
#' @return A list with `old`, `new`, `k = t - c_old` and `m` (the number of
#'   follow-up waves available for the fresh cohort, capped at `k`). Waves after
#'   the last wave with any response are treated as not yet fielded: `s_next` and
#'   `s_m1` are `NA` when the required waves are unavailable.
#' @examples
#' set.seed(3)
#' n <- 400; R <- matrix(rbinom(n * 8, 1, 0.9), n, 8); R[1:200, 1:3] <- 0L
#' wide <- data.frame(id = 1:n, entry = c(rep(4L, 200), rep(1L, 200)))
#' for (w in 1:8) wide[[paste0("r_w", w)]] <- R[, w]
#' wide$y_w1 <- rnorm(n); wide$y_w4 <- rnorm(n)
#' inp <- pc_from_wide(wide, cohort = "entry", resp = paste0("r_w", 1:8),
#'                     y_cols = c("1" = "y_w1", "4" = "y_w4"), t = 4, c_old = 1, id = "id")
#' pc_estimate(inp$old, inp$new, k = inp$k, m = inp$m)
#' @export
pc_from_wide <- function(data, cohort, resp, y_cols, t, c_old, c_new = t, x = NULL, id = NULL) {
  stopifnot("`data` must be a data frame" = is.data.frame(data),
            "`cohort` column not found in `data`" = is.character(cohort) && length(cohort) == 1 && cohort %in% names(data),
            "some `resp` columns are not in `data`" = is.character(resp) && all(resp %in% names(data)),
            "`y_cols` must be a named character vector" = is.character(y_cols) && !is.null(names(y_cols)) && all(nzchar(names(y_cols))),
            "some `y_cols` columns are not in `data`" = all(y_cols %in% names(data)),
            "`t`, `c_old`, `c_new` must be single wave numbers" = all(vapply(list(t, c_old, c_new), function(v) is.numeric(v) && length(v) == 1, logical(1))),
            "`c_old` must be an earlier wave than `t`" = c_old < t,
            "`t` exceeds the number of waves in `resp`" = t <= length(resp),
            "`c_old` must be at least wave 1" = c_old >= 1)
  if (!as.character(t) %in% names(y_cols)) stop(sprintf("`y_cols` has no entry for the comparison wave %s", format(t)), call. = FALSE)
  if (!as.character(c_old) %in% names(y_cols)) stop(sprintf("`y_cols` has no entry for the continuing cohort's entry wave %s", format(c_old)), call. = FALSE)
  if (!is.null(x) && !(is.character(x) && length(x) == 1 && x %in% names(data))) stop("`x` must name a column of `data`", call. = FALSE)
  if (!is.null(id)) {
    if (!(is.character(id) && length(id) == 1 && id %in% names(data))) stop("`id` must name a column of `data`", call. = FALSE)
    if (anyDuplicated(data[[id]])) stop(sprintf("%d duplicated values in `%s`; one row per person is required", sum(duplicated(data[[id]])), id), call. = FALSE)
  }
  if (!all(vapply(data[resp], function(v) is.numeric(v) || is.logical(v), logical(1)))) stop("response indicator columns must be numeric or logical", call. = FALSE)
  if (c_new < 1 || c_new > length(resp)) stop("`c_new` must be a wave covered by `resp`", call. = FALSE)
  R <- as.matrix(data[, resp, drop = FALSE]) * 1; R[is.na(R)] <- 0
  if (!all(R %in% c(0, 1))) stop("response indicators must be coded 0/1 (NA is read as 0)", call. = FALSE)
  storage.mode(R) <- "integer"
  surv <- function(from, to) as.integer(rowSums(R[, from:to, drop = FALSE]) == (to - from + 1))
  k <- t - c_old
  fielded <- colSums(R) > 0
  if (!fielded[t]) stop(sprintf("no response recorded at the comparison wave %s", format(t)), call. = FALSE)
  last <- max(which(fielded))
  if (any(!fielded[seq_len(last)])) message(sprintf("no response recorded at wave(s) %s; survival through these waves is impossible", paste(which(!fielded[seq_len(last)]), collapse = ", ")))
  m <- min(k, last - t)
  yt <- data[[y_cols[[as.character(t)]]]]; ye <- data[[y_cols[[as.character(c_old)]]]]
  ent <- data[[cohort]]
  old_i <- !is.na(ent) & ent == c_old
  new_i <- !is.na(ent) & ent == c_new
  if (!any(old_i)) stop(sprintf("no row with `%s` == %s (continuing cohort)", cohort, format(c_old)), call. = FALSE)
  if (!any(new_i)) stop(sprintf("no row with `%s` == %s (fresh cohort)", cohort, format(c_new)), call. = FALSE)
  n_other <- sum(!(old_i | new_i))
  if (n_other > 0) message(sprintf("%d rows belong to neither cohort and are ignored", n_other))
  n_nonresp <- sum(new_i & R[, t] == 0L)
  if (n_nonresp > 0) message(sprintf("%d fresh-cohort rows without a response at wave %s are not first-wave respondents and are dropped", n_nonresp, format(t)))
  new_i <- new_i & R[, t] == 1L
  if (!any(new_i)) stop(sprintf("no fresh-cohort respondent at wave %s", format(t)), call. = FALSE)
  s_old <- surv(c_old, t)
  old <- data.frame(y = ifelse(s_old[old_i] == 1L, yt[old_i], NA), y_entry = ye[old_i], s_prior = s_old[old_i],
                    s_next = if (t + 1 <= last) R[old_i, t + 1] else NA_integer_)
  if (!is.null(x)) old$x <- data[[x]][old_i]
  new <- data.frame(y = yt[new_i],
                    s_m = if (m >= 1) surv(t, t + m)[new_i] else R[new_i, t],
                    s_m1 = if (t + k + 1 <= last) surv(t, t + k + 1)[new_i] else NA_integer_)
  if (!is.null(id)) { old$id <- data[[id]][old_i]; new$id <- data[[id]][new_i] }
  list(old = old, new = new, k = k, m = m)
}
