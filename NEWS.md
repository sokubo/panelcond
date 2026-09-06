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
