#' What a staggered panel design identifies about the conditioning path
#'
#' Builds the cell-mean design of a staggered panel, \eqn{\mu(e,t) = \alpha(t) +
#' g(e) + \tau(s)} with tenure \eqn{s = t - e + 1} and the normalisations
#' \eqn{\tau(1) = 0}, \eqn{g(e_1) = 0}, on the support \eqn{\{(e,t): e \in E,\ e
#' \le t \le T\}} (or on an arbitrary set of fielded cells), and reports how much
#' of the conditioning path the design can identify from cell means alone.
#'
#' The identified set of \eqn{\tau} is an affine family whose direction space is
#' the *tenure projection* of the design's kernel; `identified_dim` is its
#' dimension. Two kinds of direction are in that projection for every staggered
#' design: the affine function \eqn{m(s-1)}, and, when the entry waves share a
#' common stride \eqn{d} (the greatest common divisor of their spacings), the
#' \eqn{d-1} functions of \eqn{s \bmod d} that vanish at \eqn{s = 1}. On a
#' common-end trapezoid whose last cohort is followed long enough these exhaust
#' it and `identified_dim` equals \eqn{d} (see Details); on other supports it can
#' be larger, and the rank computation is then the only check.
#'
#' A linear functional \eqn{\lambda'\tau} is identified exactly when it
#' annihilates that projection ([pc_is_identified()]). The centred lag-\eqn{d}
#' second differences \eqn{\Delta_d^2\tau(s) = \tau(s+d) - 2\tau(s) + \tau(s-d)}
#' are identified, and are the natural examples rather than the whole identified
#' space, which has dimension (number of free tenure coordinates) minus
#' `identified_dim` and contains contrasts across residue classes as well.
#' Ordinary second differences are identified only when \eqn{d = 1}, and no level
#' of \eqn{\tau} and no lag-\eqn{d} increment ever is.
#'
#' @param entries integer vector of entry waves \eqn{e_1 < \dots < e_K}.
#' @param t_max last observed wave (all cohorts are observed from entry to
#'   `t_max`). Ignored when `cells` is supplied.
#' @param cells optional data frame with columns `e` and `t` giving the fielded
#'   cells of an arbitrary support.
#' @param tol tolerance for the rank computations.
#' @return An object of class `pc_design`: a list with `summary` (one row: cells,
#'   free parameters, rank, `kernel_dim` = dimension of the whole kernel,
#'   `identified_dim` = dimension of its tenure projection and hence of the
#'   identified set of \eqn{\tau}, the stride `d`, whether `identified_dim`
#'   equals `d`, `components` = number of connected components of the
#'   cohort--period incidence graph, whether the support is a common-end
#'   trapezoid, the follow-up `w = t_max - e_K` of the last cohort, whether the
#'   follow-up condition `w >= (e_2 - e_1) - d` holds, and, for a trapezoid with
#'   three cohorts, the exact dimension `max(d, (e_2 - e_1) - w)`), `tenures`
#'   (the observed tenure values), and `tau_null_basis` (an orthonormal matrix
#'   whose `identified_dim` columns span the directions of the identified set in
#'   the \eqn{\tau} coordinates, one row per observed tenure).
#' @details `kernel_dim` and `identified_dim` differ when the cohort--period
#'   incidence graph is disconnected:
#'   `kernel_dim = identified_dim + (components - 1)`, the extra directions
#'   moving \eqn{\alpha} against \eqn{g} on the components that do not contain
#'   the reference cohort, and \eqn{\tau} not at all. Only `identified_dim`
#'   describes what the design leaves unknown about the conditioning path; on a
#'   disconnected support, knowledge of the cohort effects \eqn{g} therefore need
#'   not pin the path down. Every staggered trapezoid is connected.
#'
#'   The follow-up condition and the three-cohort formula are Lemma 1 of the
#'   companion paper on identification in staggered panels (Okubo 2026,
#'   arXiv:2609.28871, version 1)
#'   and hold for common-end trapezoids: there the
#'   identified set has dimension exactly `d` whenever the last cohort is
#'   observed for at least `(e_2 - e_1) - d` waves after its entry, and with three
#'   cohorts the dimension is `max(d, (e_2 - e_1) - w)` for every `t_max`. With
#'   more than three cohorts the follow-up condition is sufficient but not
#'   necessary (earlier cohorts can supply what the last one lacks), and
#'   `max(d, (e_2 - e_1) - w)` is an upper bound. Both
#'   are reported as `NA` for supports that are not trapezoids. Use
#'   [pc_is_identified()] to test a specific functional.
#' @examples
#' pc_design_rank(c(1, 5, 13), 19)          # JLPS-type schedule: stride 4, identified_dim 4
#' pc_design_rank(c(1, 5, 13, 20), 23)      # a fourth entry at an even wave restores stride 1
#' pc_design_rank(c(1, 7, 10), 12)          # three cohorts, follow-up too short: dimension 4 > d = 3
#'
#' # a support on which no two cohorts share a period: the kernel is larger than
#' # its tenure projection, and knowing g would not pin tau(2)
#' cells <- data.frame(e = rep(c(1, 4, 6), each = 2), t = c(1, 2, 4, 5, 6, 7))
#' pc_design_rank(cells = cells)
#' @export
pc_design_rank <- function(entries, t_max = NULL, cells = NULL, tol = 1e-9) {
  if (is.null(cells)) {
    if (is.null(t_max)) stop("supply `t_max` or `cells`")
    entries <- as.integer(entries)
    if (any(diff(entries) <= 0)) stop("`entries` must be strictly increasing")
    if (t_max < max(entries)) stop("`t_max` must be at least the last entry wave")
    cells <- do.call(rbind, lapply(entries, function(e) data.frame(e = e, t = seq.int(e, t_max))))
  } else {
    if (!all(c("e", "t") %in% names(cells))) stop("`cells` needs columns `e` and `t`")
    cells <- as.data.frame(cells)[, c("e", "t")]
  }
  cells <- unique(cells[order(cells$e, cells$t), ]); cells$s <- cells$t - cells$e + 1L
  if (any(cells$s < 1)) stop("cells with t < e are not admissible")
  E <- sort(unique(cells$e)); Tt <- sort(unique(cells$t)); S <- sort(unique(cells$s)); e0 <- E[1]
  cols <- c(paste0("a", Tt), if (length(E) > 1) paste0("g", E[-1]), paste0("u", S[S != 1]))
  X <- matrix(0, nrow(cells), length(cols), dimnames = list(NULL, cols))
  X[cbind(seq_len(nrow(cells)), match(paste0("a", cells$t), cols))] <- 1
  ge <- cells$e != e0; if (any(ge)) X[cbind(which(ge), match(paste0("g", cells$e[ge]), cols))] <- 1
  us <- cells$s != 1; if (any(us)) X[cbind(which(us), match(paste0("u", cells$s[us]), cols))] <- 1
  rk <- qr(X, tol = tol)$rank; kernel_dim <- ncol(X) - rk
  gcd2 <- function(a, b) { while (b) { r <- a %% b; a <- b; b <- r }; a }
  d <- if (length(E) > 1) Reduce(gcd2, abs(diff(E))) else NA_integer_
  sv <- svd(X, nu = 0, nv = ncol(X))
  Nb <- if (kernel_dim > 0) sv$v[, seq.int(rk + 1, ncol(X)), drop = FALSE] else matrix(0, ncol(X), 0)
  tau_basis <- matrix(0, length(S), ncol(Nb), dimnames = list(paste0("tau", S), NULL))
  if (ncol(Nb) > 0) tau_basis[S != 1, ] <- Nb[grep("^u", cols), , drop = FALSE]
  ## the identified set of tau is the COLUMN SPACE of tau_basis; its dimension can be
  ## smaller than the kernel's when the incidence graph is disconnected, so reduce to
  ## an orthonormal basis and report its rank rather than the kernel's nullity.
  if (ncol(tau_basis) > 0) {
    s2 <- svd(tau_basis)
    keep <- if (length(s2$d)) s2$d > tol * max(1, s2$d[1]) else logical(0)
    tau_basis <- s2$u[, keep, drop = FALSE]
    dimnames(tau_basis) <- list(paste0("tau", S), NULL)
  }
  identified_dim <- ncol(tau_basis)
  ## connected components of the cohort-period incidence graph
  nodes <- unique(c(paste0("e", cells$e), paste0("t", cells$t)))
  uf <- seq_along(nodes)
  find <- function(i) { while (uf[i] != i) { uf[i] <<- uf[uf[i]]; i <- uf[i] }; i }
  for (k in seq_len(nrow(cells))) {
    a <- find(match(paste0("e", cells$e[k]), nodes)); b <- find(match(paste0("t", cells$t[k]), nodes))
    if (a != b) uf[a] <- b
  }
  components <- length(unique(vapply(seq_along(nodes), find, integer(1))))
  w <- max(Tt) - max(E); D2 <- if (length(E) > 1) E[2] - E[1] else NA_integer_
  is_trap <- nrow(cells) == sum(max(Tt) - E + 1) &&
    all(paste(cells$e, cells$t) %in% paste(rep(E, max(Tt) - E + 1), unlist(lapply(E, function(e) seq.int(e, max(Tt))))))
  summary <- data.frame(cells = nrow(cells), parameters = ncol(X), rank = rk,
                        kernel_dim = kernel_dim, identified_dim = identified_dim, stride_d = d,
                        identified_equals_d = if (is.na(d)) NA else identified_dim == d,
                        components = components, trapezoid = is_trap, followup_w = w,
                        followup_condition = if (is.na(D2) || !is_trap) NA else w >= D2 - d,
                        three_cohort_dimension = if (is_trap && length(E) == 3) max(d, D2 - w) else NA_integer_)
  structure(list(summary = summary, tenures = S, tau_null_basis = tau_basis, cells = cells, X = X), class = "pc_design")
}

#' @export
print.pc_design <- function(x, ...) {
  s <- x$summary
  cat(sprintf("Staggered panel design: %d cells, %d free parameters, rank %d, kernel dimension %d\n",
              s$cells, s$parameters, s$rank, s$kernel_dim))
  free <- length(x$tenures) - 1L
  cat(sprintf("Identified set of tau has dimension %d (of %d free tenure coordinate%s)%s.\n",
              s$identified_dim, free, if (free == 1L) "" else "s",
              if (!is.na(s$stride_d)) sprintf("; stride d = %d%s", s$stride_d,
                if (isTRUE(s$identified_equals_d)) ", which it equals" else ", which it exceeds") else ""))
  if (s$identified_dim >= free) {
    cat("No nonzero functional of the conditioning path is identified on this support.\n")
  } else if (!is.na(s$stride_d) && isTRUE(s$identified_equals_d) &&
             any(x$tenures - s$stride_d %in% x$tenures & x$tenures + s$stride_d %in% x$tenures)) {
    cat(sprintf("Centred lag-%d second differences tau(s+%d) - 2 tau(s) + tau(s-%d) are identified%s;\n  they are examples, not the whole identified space (use pc_is_identified() to test a functional).\n",
                s$stride_d, s$stride_d, s$stride_d,
                if (s$stride_d == 1) "" else "; ordinary second differences are not"))
  }
  if (s$components > 1L)
    cat(sprintf("The cohort-period incidence graph has %d components, so the kernel is %d dimensions larger than\n  its tenure projection: those directions move alpha against g and tau not at all. Knowing g would not\n  pin the path down on this support.\n",
                s$components, s$components - 1L))
  if (!s$trapezoid) {
    cat("The support is not a common-end trapezoid; the follow-up rule of Lemma 1 does not apply and the rank\n  computation above is the check.\n")
  } else if (!is.na(s$followup_condition)) {
    cat(sprintf("Follow-up of last cohort w = %d; condition w >= (e_2 - e_1) - d %s%s.\n", s$followup_w,
                if (isTRUE(s$followup_condition)) "holds" else "fails",
                if (isTRUE(s$followup_condition)) ""
                else if (isTRUE(s$identified_equals_d)) ", but the identified set already has dimension d (the condition is sufficient, not necessary, with more than three cohorts)"
                else ": following the last cohort longer is one way to shrink the identified set"))
  }
  invisible(x)
}

#' Is a linear functional of the conditioning path identified by a design?
#'
#' @param design an object returned by [pc_design_rank()].
#' @param lambda numeric vector of weights, one per observed tenure (in the order
#'   of `design$tenures`), defining the functional \eqn{\sum_s \lambda_s \tau(s)}.
#' @param tol tolerance.
#' @return `TRUE` if the functional annihilates every direction of the identified
#'   set (hence is point identified from cell means), `FALSE` otherwise.
#' @examples
#' d <- pc_design_rank(c(1, 5, 13), 19)
#' lam <- numeric(19); lam[c(2, 3, 4)] <- c(1, -2, 1); pc_is_identified(d, lam)   # D^2 tau(3): FALSE
#' lam <- numeric(19); lam[c(1, 5, 9)] <- c(1, -2, 1); pc_is_identified(d, lam)   # D_4^2 tau(5): TRUE
#' @export
pc_is_identified <- function(design, lambda, tol = 1e-8) {
  if (!inherits(design, "pc_design")) stop("`design` must come from pc_design_rank()")
  if (length(lambda) != length(design$tenures)) stop("`lambda` must have one weight per observed tenure")
  if (ncol(design$tau_null_basis) == 0) return(TRUE)
  max(abs(crossprod(design$tau_null_basis, lambda))) < tol
}
