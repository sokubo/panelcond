# panelcond — refreshment-sample designs for panel conditioning

`panelcond` estimates the **panel conditioning effect** — the causal effect of
prior interviews on a wave-specific survey answer — from a panel that has added a
refreshment sample. It implements the naive fresh-versus-continuing contrast,
survival matching, **symmetric survival matching**, the **entry-wave correction**,
and attrition IPW; the decomposition of the naive contrast into conditioning and
attrition components; two diagnostic tests that distinguish state-dependent from
non-stationary attrition; worst-case bounds and breakdown values; and a simulation
engine for design planning. The core has no dependencies.

```r
remotes::install_github("sokubo/panelcond")
library(panelcond)
sim <- pc_simulate(n_old = 4800, n_new = 960, k = 4, regime = "MNAR_trait")
fit <- pc_estimate(sim$old, sim$new, k = 4)
fit
pc_bounds(fit, y_range = c(-4, 4))
pc_montecarlo(R = 200, k = 4, regime = "MNAR_state")
```

From a wide panel with per-wave response indicators:

```r
inp <- pc_from_wide(wide, cohort = "entry_wave", resp = paste0("r_w", 1:19),
                    y_cols = c("1" = "zq13_1", "5" = "dq10_1"), t = 5, c_old = 1)
pc_estimate(inp$old, inp$new, k = inp$k, m = inp$m)
```

Methods are described in Okubo, S. (2026), *Panel conditioning as a dose–response
causal effect: identification with refreshment samples*, working paper.
