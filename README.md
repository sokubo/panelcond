# panelcond: refreshment-sample designs for panel conditioning

`panelcond` estimates the panel conditioning effect, the causal effect of prior
interviews on a wave-specific survey answer, from a panel that has added a
refreshment sample. It implements the naive fresh-versus-continuing contrast,
survival matching, symmetric survival matching, the entry-wave correction of Das, Toepoel and van Soest (2011), reinterpreted as a survivor-effect estimator (with
a regression-standardised version), attrition IPW and the incremental dose
contrast between two continuing cohorts; the decomposition of the naive contrast
into conditioning and attrition components; two diagnostic tests that
distinguish state-dependent from non-stationary attrition; a joint person
bootstrap; worst-case bounds and breakdown values; and a simulation engine for
design planning. There are no dependencies beyond base R.

```r
remotes::install_github("sokubo/panelcond")
library(panelcond)
set.seed(1)
sim <- pc_simulate(n_old = 4800, n_new = 960, k = 4, regime = "MNAR_trait")
fit <- pc_estimate(sim$old, sim$new, k = 4, adjust = "x", nboot = 200)
fit
pc_bounds(fit, y_range = c(-4, 4), sd_new = sd(sim$new$y))
mc <- pc_montecarlo(R = 200, k = 4, regime = "MNAR_state")
mc
attr(mc, "tests")
```

From a wide panel with one row per person, per-wave response indicators
`r_w1, ..., r_w8`, an entry-wave column and outcome columns for waves 1 and 4:

```r
set.seed(2)
n <- 600; R <- matrix(rbinom(n * 8, 1, 0.9), n, 8); R[1:300, 1:3] <- 0L
wide <- data.frame(id = seq_len(n), entry_wave = c(rep(4L, 300), rep(1L, 300)))
for (w in 1:8) wide[[paste0("r_w", w)]] <- R[, w]
wide$y_w1 <- rnorm(n); wide$y_w4 <- rnorm(n)
inp <- pc_from_wide(wide, cohort = "entry_wave", resp = paste0("r_w", 1:8),
                    y_cols = c("1" = "y_w1", "4" = "y_w4"), t = 4, c_old = 1, id = "id")
pc_estimate(inp$old, inp$new, k = inp$k, m = inp$m)
```

Two continuing cohorts from the same population (an original sample and a
cohort-matched refreshment) with doses 12 and 8 at the same wave, each corrected
by its own entry wave:

```r
set.seed(3)
a <- pc_simulate(n_old = 4800, n_new = 10, k = 12, regime = "MNAR_trait")$old
b <- pc_simulate(n_old = 960,  n_new = 10, k = 8,  regime = "MNAR_trait")$old
pc_increment(data.frame(y = a$y, y_entry = a$y_entry, s = a$s_prior),
             data.frame(y = b$y, y_entry = b$y_entry, s = b$s_prior), kA = 12, kB = 8, nboot = 200)
```

The vignette (`vignette("panelcond")`) walks through the designs on simulated
data. `pc_design_rank()` / `pc_is_identified()` implement the identification
results of Okubo (2026), *Panel Conditioning in Fixed-Effects Models: Identification and Bias Propagation*,
[arXiv:2609.28871](https://arxiv.org/abs/2609.28871). The estimators come from a
companion manuscript on refreshment-sample designs, which is in preparation.
`citation("panelcond")` gives the package and the posted paper, and will add the
second once it is posted.

## Design planning

`pc_design_rank(entries, t_max)` builds the cell-mean design of a staggered
panel and reports what the schedule identifies about the conditioning path
before any data are collected: the stride `d` (greatest common divisor of the
spacings between entry waves), the dimension of the identified set, whether the
last cohort is followed long enough, and a basis of the unidentified directions;
`pc_is_identified()` tests a specific functional. A schedule entering cohorts at
waves 1, 5 and 13 has stride 4, so its identified set is four-dimensional: the
centred lag-4 second differences `tau(s+4) - 2 tau(s) + tau(s-4)` are identified
and ordinary second differences are not. Those are examples rather than the whole
identified space, which is larger; `pc_is_identified()` is the test. Adding a
fourth entry at an even wave (20 or 22) restores stride 1.

`pc_design_rank()` reports the dimension of the whole kernel and the dimension of
its tenure projection separately. They differ when two cohorts never share a
period: the extra directions move the period and cohort effects against each
other and the conditioning path not at all, so on such a support even exact
knowledge of the cohort effects does not pin the path down.

```r
d <- pc_design_rank(c(1, 5, 13), 19); d
lam <- numeric(19); lam[c(1, 5, 9)] <- c(1, -2, 1); pc_is_identified(d, lam)   # TRUE
lam <- numeric(19); lam[2:4] <- c(1, -2, 1);        pc_is_identified(d, lam)   # FALSE
```

