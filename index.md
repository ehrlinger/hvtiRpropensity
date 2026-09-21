# hvtiRpropensity

[![R package
version](https://img.shields.io/github/r-package/v/ehrlinger/hvtiRpropensity)](https://github.com/ehrlinger/hvtiRpropensity)

[![lint](https://github.com/ehrlinger/hvtiRpropensity/actions/workflows/lint.yaml/badge.svg)](https://github.com/ehrlinger/hvtiRpropensity/actions/workflows/lint.yaml)

Propensity score methods for cardiac surgery comparative-effectiveness
research. Ports the balancing logic from SAS programs into a tidy R API
whose outputs are compatible with
[hvtiPlotR](https://ehrlinger.github.io/hvtiPlotR/) and
[hvtiRutilities](https://github.com/ehrlinger/hvtiRutilities).

**Status: the methods are implemented, the SAS acceptance check is
not.** Score estimation
([`ps_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_logistic.md),
[`ps_ordinal()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_ordinal.md),
[`ps_nominal()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_nominal.md),
[`bs_continuous()`](https://ehrlinger.github.io/hvtiRpropensity/reference/bs_continuous.md),
[`bs_count()`](https://ehrlinger.github.io/hvtiRpropensity/reference/bs_count.md)),
both balancing methods
([`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md),
[`ps_weight()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_weight.md))
and four sensitivity analyses
([`sa_rosenbaum()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_rosenbaum.md),
[`sa_evalue()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_evalue.md),
[`sa_overlap()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_overlap.md),
[`sa_trim_sweep()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_trim_sweep.md))
are exported and tested.

What is not yet done is parity. Every test in this package runs against
the synthetic generators; no test reproduces a real run of the SAS
programs whose balancing logic this ports. Until one does, treat
agreement with SAS as unverified rather than assumed.

**Full documentation:** <https://ehrlinger.github.io/hvtiRpropensity/>

## Installation

``` r

# Install the development version from GitHub:
# install.packages("pak")
pak::pak("ehrlinger/hvtiRpropensity")
```

## Quick start

``` r

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

[`fit_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/fit_logistic.md)
fits binary, proportional-odds ordinal, and generalized-logit nominal
outcome models. Clinically meaningful levels are always explicit; factor
storage order never decides the event, order, or reference category.

``` r

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
- `$meta`: formula, declared levels, model family, row accounting,
  package versions, and `bundle_version = 1L`;
- `$tables`: `estimates`, `covariance`, `by_imputation`, and
  `fit_status`, plus workflow-specific diagnostics;
- `$models`: every fitted model, including every imputation-specific
  fit.

The `estimates` table has columns `term`, `estimate`, `std.error`,
`statistic`, `df`, `p.value`, `conf.low`, `conf.high`, `odds_ratio`, and
`pooled`. `by_imputation` has `imputation`, `term`, `estimate`, and
`std.error`; `fit_status` has `imputation`, `converged`, `n_input`,
`n_analyzed`, and `n_excluded`.

For stacked imputations, patient predictions are averaged across fits
while coefficients and covariance use Rubin’s rules. Missing patients,
changed categorical predictor levels, aliased terms, non-convergence,
and unusable covariance matrices stop the analysis instead of producing
a partial pool.

[`validate_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/validate_logistic.md)
applies a saved version-1 binary bundle to a declared validation cohort
without refitting. It reports calibration, observed versus expected
events, AUC, and Brier score. Ordinal and nominal validation are not
part of this first interface.

## Functions

**Score estimation**

| Function | Description |
|----|----|
| [`fit_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/fit_logistic.md) | Binary, ordinal, or nominal outcome-model bundle, with optional MI |
| [`validate_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/validate_logistic.md) | No-refit validation of a saved binary model bundle |
| [`ps_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_logistic.md) | Binary propensity score via logistic regression (with optional MI) |
| [`ps_ordinal()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_ordinal.md) | Ordered treatment propensity via cumulative logit (MASS::polr) |
| [`ps_nominal()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_nominal.md) | Nominal treatment propensity via generalised logit (nnet::multinom) |
| [`bs_continuous()`](https://ehrlinger.github.io/hvtiRpropensity/reference/bs_continuous.md) | Balancing score for continuous exposures via linear regression |
| [`bs_count()`](https://ehrlinger.github.io/hvtiRpropensity/reference/bs_count.md) | Balancing score for count exposures via negative-binomial / Poisson |

**Balancing methods**

| Function | Description |
|----|----|
| [`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md) | Greedy nearest-neighbour 1:1 matching with optional caliper |
| [`ps_weight()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_weight.md) | IPTW weighting — ATE, ATT, or ATC estimand |

**Sensitivity analysis**

| Function | Description |
|----|----|
| [`sa_rosenbaum()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_rosenbaum.md) | Rosenbaum gamma bounds for matched analyses |
| [`sa_evalue()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_evalue.md) | E-values (VanderWeele & Ding 2017) for unmeasured confounding |
| [`sa_overlap()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_overlap.md) | PS overlap / positivity diagnostics |
| [`sa_trim_sweep()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_trim_sweep.md) | IPTW trim-threshold sensitivity sweep |

**Utilities**

| Function | Description |
|----|----|
| [`sample_ps_data()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sample_ps_data.md) | Synthetic cardiac-surgery dataset (binary treatment) |
| [`sample_ps_data_ordinal()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sample_ps_data_ordinal.md) | Synthetic dataset with ordered treatment |
| [`sample_ps_data_nominal()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sample_ps_data_nominal.md) | Synthetic dataset with nominal treatment |
| [`sample_ps_data_count()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sample_ps_data_count.md) | Synthetic dataset with count exposure |
| [`is_ps_data()`](https://ehrlinger.github.io/hvtiRpropensity/reference/is_ps_data.md) | Predicate for `ps_data` objects |

## Related packages

- **hvtiPlotR** — mirror histograms, balance plots, survival and hazard
  curves
- **hvtiRutilities** — data wrangling utilities for cardiac surgery
  datasets
