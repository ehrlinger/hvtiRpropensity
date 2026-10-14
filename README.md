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

**Status: the methods are implemented, and SAS parity is checked for the
model fits and the CCF macros.** Score estimation (`ps_logistic()`,
`ps_forest()`, `ps_ordinal()`, `ps_nominal()`, `bs_continuous()`,
`bs_count()`), the outcome-model bundles (`fit_logistic()`, `fit_linear()`,
`validate_logistic()`), matching and weighting (`ps_match()`, `ps_weight()`),
balance statistics (`ps_stddiff()`, `ps_stddiff_perm()`, `ps_mw_var()`), the
support and survival summaries (`ps_support()`, `ps_rmst()`) and four
sensitivity analyses (`sa_rosenbaum()`, `sa_evalue()`, `sa_overlap()`,
`sa_trim_sweep()`) are exported and tested.

Since 0.1.10, tests compare the package against SAS 9.4 runs on synthetic
fixtures (`tests/testthat/fixtures/`): the binary, ordinal and nominal
logistic fits, the Poisson balancing score, `fit_linear()`, `bs_continuous()`,
Rubin pooling against `PROC MIANALYZE`, and the CCF macros `%stddiff`,
`%stddiffci` and `%mw_var`. Matching, weighting, `ps_rmst()`, `ps_forest()`
and the `sa_*()` functions have no SAS oracle yet, so treat their agreement
with SAS as unverified.

**Full documentation:** <https://ehrlinger.github.io/hvtiRpropensity/>

## Installation

hvtiRpropensity is part of the HVTI R family. Install it with the rest of
the family through [hvtiR](https://github.com/ehrlinger/hvtiR):

```r
# install.packages("pak")
pak::pak("ehrlinger/hvtiR")
hvtiR::install()
```

`hvtiR::install()` installs every member in one call, so the versions the
members require of one another resolve together. Update later with
`hvtiR::update()`, from a fresh R session: it will not overwrite a package
that is already loaded.

To install hvtiRpropensity on its own:

```r
pak::pak("ehrlinger/hvtiRpropensity")
```

## Quick start

```r
library(hvtiRpropensity)

# Reproducible synthetic dataset. Drop the generator's true score so the
# example estimates its own.
dta <- sample_ps_data(n = 500, seed = 42)
dta$prob_t <- NULL

# --- Propensity score estimation ---
obj <- ps_logistic(tavr ~ age + female + ef + diabetes + hypertension, data = dta)

# --- 1:1 nearest-neighbour matching on the logit, caliper 0.2 SD ---
m <- ps_match(obj$data, score_col = "logit_t", scale = "logit",
              caliper_sd = 0.2, seed = 42)
print(m)
summary(m)
matched <- m$data[m$data$match == 1L, ]

# --- IPTW weighting (ATE) from the estimated score ---
w <- ps_weight(obj$data, estimand = "ATE")
print(w)
summary(w)
```

The caliper of 0.2 standard deviations of the logit of the propensity score
is the Austin (2011) recommendation. Give `seed` on every `ps_match()` call.
Matching shuffles the order of the treated patients, and the seed makes the
pairs reproducible.

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
`df`, `p.value`, `conf.low`, `conf.high`, `odds_ratio`, and `pooled`. Count models
(`bs_count()`) name the exponentiated column `rate_ratio`, since a log-link
coefficient exponentiates to a rate ratio.
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
| `ps_logistic()` | Binary propensity score via logistic regression (with optional MI) |
| `ps_forest()` | Binary propensity score from a random forest, out-of-bag (`randomForestSRC`) |
| `ps_ordinal()` | Ordered treatment propensity via cumulative logit (`MASS::polr`) |
| `ps_nominal()` | Nominal treatment propensity via generalised logit (`nnet::multinom`) |
| `bs_continuous()` | Balancing score for continuous exposures via linear regression |
| `bs_count()` | Balancing score for count exposures via negative-binomial or Poisson |

**Outcome models**

| Function | Description |
|---|---|
| `fit_logistic()` | Binary, ordinal, or nominal outcome-model bundle, with optional MI |
| `fit_linear()` | Linear outcome-model bundle, with optional MI and Rubin pooling |
| `validate_logistic()` | No-refit validation of a saved binary model bundle |

**Balancing methods**

| Function | Description |
|---|---|
| `ps_match()` | Greedy 1:1 nearest-neighbour matching, on the probability or logit scale, with an optional caliper (absolute or `caliper_sd`) |
| `ps_weight()` | IPTW weighting for the ATE, ATT, or ATC estimand, optionally trimmed |

**Balance and treatment effects**

| Function | Description |
|---|---|
| `ps_stddiff()` | Standardized differences by variable type, optionally weighted (CCF `%stddiff`) |
| `ps_stddiff_perm()` | Permutation reference for those differences (CCF `%stddiffci`) |
| `ps_mw_var()` | Matching-weight difference in means with a bootstrap SD (CCF `%mw_var`) |
| `ps_rmst()` | Restricted mean survival contrast, weighted or on matched pairs |
| `ps_support()` | Compares common-support and trimming rules and how well they agree |

**Sensitivity analysis**

| Function | Description |
|---|---|
| `sa_rosenbaum()` | Rosenbaum gamma bounds for matched analyses (unmeasured confounding) |
| `sa_evalue()` | E-values (VanderWeele & Ding 2017) for unmeasured confounding |
| `sa_overlap()` | Overlap and positivity diagnostics on the score distribution |
| `sa_trim_sweep()` | How effective sample size and weights move across IPTW trim thresholds |

**Utilities**

| Function | Description |
|---|---|
| `sample_ps_data()` | Synthetic cardiac-surgery dataset (binary treatment) |
| `sample_ps_data_ordinal()` | Synthetic dataset with ordered treatment |
| `sample_ps_data_nominal()` | Synthetic dataset with nominal treatment |
| `sample_ps_data_count()` | Synthetic dataset with count exposure |
| `is_ps_data()` | Predicate for `ps_data` objects |

The four generators return plain data frames and `is_ps_data()` returns
`TRUE` or `FALSE`. Every other export returns a `ps_data` object with
`$data`, `$meta` and `$tables`.

## Documentation

The [Get started](https://ehrlinger.github.io/hvtiRpropensity/articles/hvtiRpropensity.html)
vignette walks the workflow from score to sensitivity analysis, including a
worked example that reshapes the balance tables for `hvtiPlotR::hv_balance()`.

The SAS programs this package ports ship with it:
`system.file("sas-templates", package = "hvtiRpropensity")`.

## Related packages

- **hvtiPlotR**: mirror histograms, balance plots, survival and hazard curves.
  `hv_mirror_hist()` reads the `$data` of a `ps_match()` or `ps_weight()`
  object directly. `hv_balance()` needs the SMD tables reshaped to long form
  in percent; the vignette shows how.
- **hvtiRutilities**: data wrangling utilities for cardiac surgery datasets.
