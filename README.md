# hvtiRpropensity

<!-- badges: start -->
[![R-CMD-check](https://github.com/ehrlinger/hvtiRpropensity/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ehrlinger/hvtiRpropensity/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/ehrlinger/hvtiRpropensity/graph/badge.svg)](https://app.codecov.io/gh/ehrlinger/hvtiRpropensity)
[![active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/badges/latest/active.svg)
[![pkgdown](https://github.com/ehrlinger/hvtiRpropensity/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/ehrlinger/hvtiRpropensity/actions/workflows/pkgdown.yaml)

[![R package version](https://img.shields.io/github/r-package/v/ehrlinger/hvtiRpropensity)](https://github.com/ehrlinger/hvtiRpropensity)

[![lint](https://github.com/ehrlinger/hvtiRpropensity/actions/workflows/lint.yaml/badge.svg)](https://github.com/ehrlinger/hvtiRpropensity/actions/workflows/lint.yaml)
<!-- badges: end -->

Propensity score methods for cardiac surgery comparative-effectiveness
research.  Ports the balancing logic from SAS programs into a tidy R API
whose outputs are compatible with
[hvtiPlotR](https://ehrlinger.github.io/hvtiPlotR/) and
[hvtiRutilities](https://github.com/ehrlinger/hvtiRutilities).

**Status: the methods are implemented, the SAS acceptance check is not.**
Score estimation (`ps_logistic()`, `ps_ordinal()`, `ps_nominal()`,
`bs_continuous()`, `bs_count()`), both balancing methods (`ps_match()`,
`ps_weight()`) and four sensitivity analyses (`sa_rosenbaum()`, `sa_evalue()`,
`sa_overlap()`, `sa_trim_sweep()`) are exported and tested.

What is not yet done is parity. Every test in this package runs against the
synthetic generators; no test reproduces a real run of the SAS programs whose
balancing logic this ports. Until one does, treat agreement with SAS as
unverified rather than assumed.

**Full documentation:** <https://ehrlinger.github.io/hvtiRpropensity/>

## Installation

```r
# Install the development version from GitHub:
# install.packages("pak")
pak::pak("ehrlinger/hvtiRpropensity")
```

## Quick start

```r
library(hvtiRpropensity)

# Reproducible synthetic dataset
dta <- sample_ps_data(n = 500, seed = 42)
dta$prob_t <- NULL  # Preserve the generator's true score only when it is an input.

# --- Propensity score estimation ---
obj <- ps_logistic(tavr ~ age + female + ef + diabetes + hypertension, data = dta)

# --- Nearest-neighbour 1:1 matching ---
m <- ps_match(obj$data)
print(m)
summary(m)
matched <- m$data[m$data$match == 1L, ]

# --- IPTW weighting (ATE) ---
w <- ps_weight(dta, estimand = "ATE")
print(w)
summary(w)
```

## Logistic model bundles

`fit_logistic()` fits binary, proportional-odds ordinal, and generalized-logit
nominal outcome models. Clinically meaningful levels are always explicit;
factor storage order never decides the event, order, or reference category.

```r
binary_model <- fit_logistic(
  tavr ~ age + female + ef + diabetes + hypertension,
  data = dta,
  family = "binary",
  outcome_levels = c(0, 1),
  event_level = 1
)
```

Every model-producing result has four slots:

- `$data`: patient-level data with predictions or scores;
- `$meta`: formula, declared levels, model family, row accounting, package
  versions, and `bundle_version = 1L`;
- `$tables`: `estimates`, `covariance`, `by_imputation`, and `fit_status`, plus
  workflow-specific diagnostics;
- `$models`: every fitted model, including every imputation-specific fit.

The `estimates` table has columns `term`, `estimate`, `std.error`, `statistic`,
`df`, `p.value`, `conf.low`, `conf.high`, `odds_ratio`, and `pooled`.
`by_imputation` has `imputation`, `term`, `estimate`, and `std.error`;
`fit_status` has `imputation`, `converged`, `n_input`, `n_analyzed`, and
`n_excluded`.

For stacked imputations, patient predictions are averaged across fits while
coefficients and covariance use Rubin's rules. Missing patients, changed
categorical predictor levels, aliased terms, non-convergence, and unusable
covariance matrices stop the analysis instead of producing a partial pool.

`validate_logistic()` applies a saved version-1 binary bundle to a declared
validation cohort without refitting. It reports calibration, observed versus
expected events, AUC, and Brier score. Ordinal and nominal validation are not
part of this first interface.

## Functions

**Score estimation**

| Function | Description |
|---|---|
| `fit_logistic()` | Binary, ordinal, or nominal outcome-model bundle, with optional MI |
| `validate_logistic()` | No-refit validation of a saved binary model bundle |
| `ps_logistic()` | Binary propensity score via logistic regression (with optional MI) |
| `ps_ordinal()` | Ordered treatment propensity via cumulative logit (MASS::polr) |
| `ps_nominal()` | Nominal treatment propensity via generalised logit (nnet::multinom) |
| `bs_continuous()` | Balancing score for continuous exposures via linear regression |
| `bs_count()` | Balancing score for count exposures via negative-binomial / Poisson |

**Balancing methods**

| Function | Description |
|---|---|
| `ps_match()` | Greedy nearest-neighbour 1:1 matching with optional caliper |
| `ps_weight()` | IPTW weighting — ATE, ATT, or ATC estimand |

**Sensitivity analysis**

| Function | Description |
|---|---|
| `sa_rosenbaum()` | Rosenbaum gamma bounds for matched analyses |
| `sa_evalue()` | E-values (VanderWeele & Ding 2017) for unmeasured confounding |
| `sa_overlap()` | PS overlap / positivity diagnostics |
| `sa_trim_sweep()` | IPTW trim-threshold sensitivity sweep |

**Utilities**

| Function | Description |
|---|---|
| `sample_ps_data()` | Synthetic cardiac-surgery dataset (binary treatment) |
| `sample_ps_data_ordinal()` | Synthetic dataset with ordered treatment |
| `sample_ps_data_nominal()` | Synthetic dataset with nominal treatment |
| `sample_ps_data_count()` | Synthetic dataset with count exposure |
| `is_ps_data()` | Predicate for `ps_data` objects |

## Related packages

- **hvtiPlotR** — mirror histograms, balance plots, survival and hazard curves
- **hvtiRutilities** — data wrangling utilities for cardiac surgery datasets
