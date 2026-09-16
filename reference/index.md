# Package index

## Package Overview

Package-level documentation and the base S3 class shared by all
hvtiRpropensity data objects.

- [`hvtiRpropensity-package`](https://ehrlinger.github.io/hvtiRpropensity/reference/hvtiRpropensity-package.md)
  [`hvtiRpropensity`](https://ehrlinger.github.io/hvtiRpropensity/reference/hvtiRpropensity-package.md)
  : hvtiRpropensity: Propensity Score Methods for Cardiac Surgery
  Studies
- [`is_ps_data()`](https://ehrlinger.github.io/hvtiRpropensity/reference/is_ps_data.md)
  : Test whether an object is a ps_data instance

## Data

Synthetic datasets for examples and unit tests.

- [`sample_ps_data()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sample_ps_data.md)
  : Generate sample data for propensity score examples
- [`sample_ps_data_ordinal()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sample_ps_data_ordinal.md)
  : Generate sample data with an ordinal treatment for ps_ordinal()
  examples
- [`sample_ps_data_nominal()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sample_ps_data_nominal.md)
  : Generate sample data with a nominal treatment for ps_nominal()
  examples
- [`sample_ps_data_count()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sample_ps_data_count.md)
  : Generate sample data with a count outcome for bs_count() examples

## Propensity Score Estimation

Functions for estimating propensity scores and balancing scores from raw
patient data. Supports binary, ordinal, nominal, continuous, and count
outcomes, with optional multiple-imputation averaging.

- [`ps_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_logistic.md)
  : Estimate propensity scores via binary logistic regression
- [`ps_ordinal()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_ordinal.md)
  : Estimate propensity scores via ordinal (cumulative logit) logistic
  regression
- [`ps_nominal()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_nominal.md)
  : Estimate propensity scores via nominal (generalised logit) logistic
  regression
- [`bs_continuous()`](https://ehrlinger.github.io/hvtiRpropensity/reference/bs_continuous.md)
  : Estimate a balancing score via linear regression
- [`bs_count()`](https://ehrlinger.github.io/hvtiRpropensity/reference/bs_count.md)
  : Estimate a balancing score via negative-binomial or Poisson
  regression

## Propensity Score Methods

Functions for balancing treated and control groups via matching or
inverse-probability-of-treatment weighting.

- [`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md)
  : Propensity Score Nearest-Neighbour Matching
- [`ps_weight()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_weight.md)
  : Inverse-Probability-of-Treatment Weighting (IPTW)

## Balance Diagnostics

Standardized differences between treated and control groups, by variable
type and optionally weighted.

- [`ps_stddiff()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_stddiff.md)
  : Standardized differences between two groups

## Sensitivity Analysis

Functions for assessing robustness of propensity score analyses to
unmeasured confounding and model assumptions.

- [`sa_rosenbaum()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_rosenbaum.md)
  : Rosenbaum Sensitivity Bounds for Matched Analyses
- [`sa_evalue()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_evalue.md)
  : E-Value Sensitivity Analysis
- [`sa_overlap()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_overlap.md)
  : Propensity Score Overlap and Positivity Diagnostics
- [`sa_trim_sweep()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_trim_sweep.md)
  : IPTW Weight-Trimming Sensitivity Sweep

## Internal

S3 methods for printing and summarising ps_data objects.

- [`print(`*`<ps_data>`*`)`](https://ehrlinger.github.io/hvtiRpropensity/reference/print.ps_data.md)
  : Print a ps_data object
- [`summary(`*`<ps_data>`*`)`](https://ehrlinger.github.io/hvtiRpropensity/reference/summary.ps_data.md)
  : Summarise a ps_data object
