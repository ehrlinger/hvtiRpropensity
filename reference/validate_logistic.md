# Validate a saved binary logistic model

Scores a declared validation cohort with every fitted model retained in
a binary
[`fit_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/fit_logistic.md)
bundle. It does not refit or modify the source bundle.

## Usage

``` r
validate_logistic(
  model,
  data,
  outcome_col = NULL,
  prediction_col = "predicted",
  groups = 10L
)
```

## Arguments

- model:

  A binary `lm_fit` bundle with bundle version 1.

- data:

  Validation data containing the stored predictors and outcome.

- outcome_col:

  Outcome column. Defaults to the bundle's outcome column.

- prediction_col:

  Output column for averaged predicted probabilities.

- groups:

  Number of equal-frequency calibration groups.

## Value

An object of class `c("lm_validation", "ps_data")`.
