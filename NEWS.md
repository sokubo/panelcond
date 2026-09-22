# panelcond 0.1.3

* New `pc_design_rank()` and `pc_is_identified()`: build the cell-mean design of
  a staggered panel (or of an arbitrary set of fielded cells), report its rank
  and nullity, the *stride* of the refreshment schedule (the greatest common
  divisor of the spacings between entry waves), whether the follow-up of the
  last cohort satisfies the condition of Lemma 1 of the companion identification
  paper, and a basis of the directions of the identified set of the conditioning
  path; test whether a given linear functional of the path is identified. On a
  schedule with stride `d` the centred lag-`d` second differences
  `tau(s+d) - 2 tau(s) + tau(s-d)` are identified and ordinary second differences
  are not, but they are examples rather than the whole identified space: use
  `pc_is_identified()` to test a functional.
* `pc_design_rank()` distinguishes the dimension of the whole kernel
  (`kernel_dim`) from the dimension of its tenure projection (`identified_dim`),
  which is the dimension of the identified set of the conditioning path. The two
  differ when the cohort-period incidence graph is disconnected, by one dimension
  per extra component; an earlier draft of this function reported the kernel's
  dimension as the identified set's and, on such a support, both overstated it and
  blamed the follow-up of the last cohort. `components` and `trapezoid` are new
  columns of the summary, and the follow-up rule of Lemma 1 is reported only for
  common-end trapezoids, where it applies. `tau_null_basis` is now an orthonormal
  basis of the identified set's direction space, with `identified_dim` columns.
* `inst/CITATION` no longer points at a working paper that is not going to appear
  under that title; it cites the package, and records that the two companion
  manuscripts are in preparation.
* `tests/testthat/test-design_rank.R` covers the schedules quoted in the
  documentation, the three-cohort formula on a grid of trapezoids, the
  disconnected-support case, and the identity
  `kernel_dim = identified_dim + (components - 1)` on a sweep of supports.
* Documentation now credits the entry-wave correction to Das, Toepoel and van
  Soest (2011, Sociological Methods & Research 40: 32-56, Assumption 3): the
  observable combination is theirs; `pc_estimate()` reinterprets it, under the
  restriction that the survivors' selection differential on the *latent*
  outcome is the same at entry and at the comparison wave, as an estimator of
  the conditioning effect *among survivors* rather than of a population effect.
  The two coincide under a homogeneous effect and differ otherwise.
* The diagnostics `T1` and `T2` are documented as compatibility checks among the
  maintained restrictions (including survivor-effect homogeneity), not as
  classifiers of which restriction fails; when non-stationary and
  state-dependent attrition are both present, `T2` rejects and `T1` may not,
  and the least-biased estimator can be the one whose restriction is most
  violated.

# panelcond 0.1.2

* The analytic standard error of the entry-wave correction in `pc_estimate()`
  and of the increment in `pc_increment()` now uses the first-order variance of
  the survivors' change plus the cohort's entry mean, which keeps the covariance
  between the survivors' two answers. The earlier expression omitted it and was
  conservative (variance up to 10% too large at dose 4, about 18% for the
  increment).
* `T1` is no longer `NA` when the fresh cohort has item non-response.
* `pc_estimate()` and `pc_increment()` gain `nboot`, a joint person bootstrap
  (persons resampled within each cohort) that supplies standard errors for every
  estimator, including IPW and EC-adj, and bootstrap denominators for `T1`, `T2`.
* `pc_estimate()` gains `adjust`, the regression-standardised entry-wave
  correction (EC-adj).
* `pc_bounds()` returns the additive drift at which the EC estimate vanishes
  (`ec_eta_star`) and at which its 95% interval reaches zero (`ec_eta_ci95`);
  `y_range` is sorted and `sd_new` must be positive.
* `pc_from_wide()` gains `id` (duplicated identifiers are an error), treats waves
  after the last observed response as not fielded, and checks its arguments.
* `pc_simulate()` returns `s_m1 = NA` when `m < k`; `pc_montecarlo()` reports
  the Monte Carlo standard deviation, the mean analytic standard error and the
  Monte Carlo standard error of the bias.
* Input checks with informative messages throughout; rows dropped for item
  non-response are reported.
* Vignette and `inst/CITATION` added.

# panelcond 0.1.1

* New `pc_increment()`: incremental dose contrast between two continuing cohorts
  drawn from the same population (entry-wave correction on both sides), for panels
  whose later refreshment sample targets a different population.
* `pc_bounds()` gains `sd_new`; documentation clarifies the survivor estimand.

# panelcond 0.1.0

* Initial development version: `pc_estimate()`, `pc_from_wide()`, `pc_simulate()`,
  `pc_montecarlo()`, `pc_bounds()`.
