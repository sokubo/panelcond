# panelcond: refreshment-sample designs for panel conditioning

`panelcond` estimates the panel conditioning effect, the causal effect of prior
interviews on a wave-specific survey answer, from a panel that has added a
refreshment sample. It implements the naive fresh-versus-continuing contrast,
survival matching, symmetric survival matching, the entry-wave correction (with
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
data. Methods are described in Okubo, S. (2026), *Panel conditioning as a
dose-response causal effect: identification with refreshment samples*, working
paper; see `citation("panelcond")`.
