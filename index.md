# hvtiRpropensity

[![R package
version](https://img.shields.io/github/r-package/v/ehrlinger/hvtiRpropensity)](https://github.com/ehrlinger/hvtiRpropensity)

[![lint](https://github.com/ehrlinger/hvtiRpropensity/actions/workflows/lint.yaml/badge.svg)](https://github.com/ehrlinger/hvtiRpropensity/actions/workflows/lint.yaml)

Propensity score methods for cardiac surgery comparative-effectiveness
research. Ports the balancing logic from SAS programs into a tidy R API
whose outputs are compatible with
[hvtiPlotR](https://ehrlinger.github.io/hvtiPlotR/) and
[hvtiRutilities](https://github.com/ehrlinger/hvtiRutilities).

**Status: the methods are implemented, and SAS parity is checked for the
model fits and the CCF macros.** Score estimation
([`ps_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_logistic.md),
[`ps_forest()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_forest.md),
[`ps_ordinal()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_ordinal.md),
[`ps_nominal()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_nominal.md),
[`bs_continuous()`](https://ehrlinger.github.io/hvtiRpropensity/reference/bs_continuous.md),
[`bs_count()`](https://ehrlinger.github.io/hvtiRpropensity/reference/bs_count.md)),
the outcome-model bundles
([`fit_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/fit_logistic.md),
[`fit_linear()`](https://ehrlinger.github.io/hvtiRpropensity/reference/fit_linear.md),
[`validate_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/validate_logistic.md)),
matching and weighting
([`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md),
[`ps_weight()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_weight.md)),
balance statistics
([`ps_stddiff()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_stddiff.md),
[`ps_stddiff_perm()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_stddiff_perm.md),
[`ps_mw_var()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_mw_var.md)),
the support and survival summaries
([`ps_support()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_support.md),
[`ps_rmst()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_rmst.md))
and four sensitivity analyses
([`sa_rosenbaum()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_rosenbaum.md),
[`sa_evalue()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_evalue.md),
[`sa_overlap()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_overlap.md),
[`sa_trim_sweep()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_trim_sweep.md))
are exported and tested.

Since 0.1.10, tests compare the package against SAS 9.4 runs on
synthetic fixtures (`tests/testthat/fixtures/`): the binary, ordinal and
nominal logistic fits, the Poisson balancing score,
[`fit_linear()`](https://ehrlinger.github.io/hvtiRpropensity/reference/fit_linear.md),
[`bs_continuous()`](https://ehrlinger.github.io/hvtiRpropensity/reference/bs_continuous.md),
Rubin pooling against `PROC MIANALYZE`, and the CCF macros `%stddiff`,
`%stddiffci` and `%mw_var`. Matching, weighting,
[`ps_rmst()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_rmst.md),
[`ps_forest()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_forest.md)
and the `sa_*()` functions have no SAS oracle yet, so treat their
agreement with SAS as unverified.

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

The caliper of 0.2 standard deviations of the logit of the propensity
score is the Austin (2011) recommendation. Give `seed` on every
[`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md)
call. Matching shuffles the order of the treated patients, and the seed
makes the pairs reproducible.

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
`pooled`. Count models
([`bs_count()`](https://ehrlinger.github.io/hvtiRpropensity/reference/bs_count.md))
name the exponentiated column `rate_ratio`, since a log-link coefficient
exponentiates to a rate ratio. `by_imputation` has `imputation`, `term`,
`estimate`, and `std.error`; `fit_status` has `imputation`, `converged`,
`n_input`, `n_analyzed`, and `n_excluded`.

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
| [`ps_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_logistic.md) | Binary propensity score via logistic regression (with optional MI) |
| [`ps_forest()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_forest.md) | Binary propensity score from a random forest, out-of-bag (`randomForestSRC`) |
| [`ps_ordinal()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_ordinal.md) | Ordered treatment propensity via cumulative logit ([`MASS::polr`](https://rdrr.io/pkg/MASS/man/polr.html)) |
| [`ps_nominal()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_nominal.md) | Nominal treatment propensity via generalised logit ([`nnet::multinom`](https://rdrr.io/pkg/nnet/man/multinom.html)) |
| [`bs_continuous()`](https://ehrlinger.github.io/hvtiRpropensity/reference/bs_continuous.md) | Balancing score for continuous exposures via linear regression |
| [`bs_count()`](https://ehrlinger.github.io/hvtiRpropensity/reference/bs_count.md) | Balancing score for count exposures via negative-binomial or Poisson |

**Outcome models**

| Function | Description |
|----|----|
| [`fit_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/fit_logistic.md) | Binary, ordinal, or nominal outcome-model bundle, with optional MI |
| [`fit_linear()`](https://ehrlinger.github.io/hvtiRpropensity/reference/fit_linear.md) | Linear outcome-model bundle, with optional MI and Rubin pooling |
| [`validate_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/validate_logistic.md) | No-refit validation of a saved binary model bundle |

**Balancing methods**

| Function | Description |
|----|----|
| [`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md) | Greedy 1:1 nearest-neighbour matching, on the probability or logit scale, with an optional caliper (absolute or `caliper_sd`) |
| [`ps_weight()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_weight.md) | IPTW weighting for the ATE, ATT, or ATC estimand, optionally trimmed |

**Balance and treatment effects**

| Function | Description |
|----|----|
| [`ps_stddiff()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_stddiff.md) | Standardized differences by variable type, optionally weighted (CCF `%stddiff`) |
| [`ps_stddiff_perm()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_stddiff_perm.md) | Permutation reference for those differences (CCF `%stddiffci`) |
| [`ps_mw_var()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_mw_var.md) | Matching-weight difference in means with a bootstrap SD (CCF `%mw_var`) |
| [`ps_rmst()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_rmst.md) | Restricted mean survival contrast, weighted or on matched pairs |
| [`ps_support()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_support.md) | Compares common-support and trimming rules and how well they agree |

**Sensitivity analysis**

| Function | Description |
|----|----|
| [`sa_rosenbaum()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_rosenbaum.md) | Rosenbaum gamma bounds for matched analyses (unmeasured confounding) |
| [`sa_evalue()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_evalue.md) | E-values (VanderWeele & Ding 2017) for unmeasured confounding |
| [`sa_overlap()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_overlap.md) | Overlap and positivity diagnostics on the score distribution |
| [`sa_trim_sweep()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_trim_sweep.md) | How effective sample size and weights move across IPTW trim thresholds |

**Utilities**

| Function | Description |
|----|----|
| [`sample_ps_data()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sample_ps_data.md) | Synthetic cardiac-surgery dataset (binary treatment) |
| [`sample_ps_data_ordinal()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sample_ps_data_ordinal.md) | Synthetic dataset with ordered treatment |
| [`sample_ps_data_nominal()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sample_ps_data_nominal.md) | Synthetic dataset with nominal treatment |
| [`sample_ps_data_count()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sample_ps_data_count.md) | Synthetic dataset with count exposure |
| [`is_ps_data()`](https://ehrlinger.github.io/hvtiRpropensity/reference/is_ps_data.md) | Predicate for `ps_data` objects |

The four generators return plain data frames and
[`is_ps_data()`](https://ehrlinger.github.io/hvtiRpropensity/reference/is_ps_data.md)
returns `TRUE` or `FALSE`. Every other export returns a `ps_data` object
with `$data`, `$meta` and `$tables`.

## Documentation

The [Get
started](https://ehrlinger.github.io/hvtiRpropensity/articles/hvtiRpropensity.html)
vignette walks the workflow from score to sensitivity analysis,
including a worked example that reshapes the balance tables for
[`hvtiPlotR::hv_balance()`](https://ehrlinger.github.io/hvtiPlotR/reference/hv_balance.html).

The SAS programs this package ports ship with it:
`system.file("sas-templates", package = "hvtiRpropensity")`.

## Related packages

- **hvtiPlotR**: mirror histograms, balance plots, survival and hazard
  curves. `hv_mirror_hist()` reads the `$data` of a
  [`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md)
  or
  [`ps_weight()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_weight.md)
  object directly. `hv_balance()` needs the SMD tables reshaped to long
  form in percent; the vignette shows how.
- **hvtiRutilities**: data wrangling utilities for cardiac surgery
  datasets.
