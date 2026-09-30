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

## Examples

``` r
dta <- sample_ps_data(n = 150, seed = 42)
odd <- seq_len(nrow(dta)) %% 2 == 1
# Fit on one half, then score the other half without refitting.
fit <- fit_logistic(
  tavr ~ age + female + ef,
  data           = dta[odd, ],
  family         = "binary",
  outcome_levels = c(0, 1),
  event_level    = 1
)
val <- validate_logistic(fit, dta[!odd, ])
val$tables$performance
#>     n observed expected  oe_ratio       auc    brier
#> 1 150       75 77.24495 0.9709373 0.7607111 0.199376
val$tables$calibration
#>    group  n observed  expected
#> 1      1 15        2  2.426137
#> 2      2 15        4  4.038097
#> 3      3 15        3  5.137954
#> 4      4 15        7  6.241516
#> 5      5 15        6  7.155727
#> 6      6 15        9  8.128755
#> 7      7 15       10  9.049067
#> 8      8 15       11 10.100560
#> 9      9 15        9 11.685894
#> 10    10 15       14 13.281238
```
