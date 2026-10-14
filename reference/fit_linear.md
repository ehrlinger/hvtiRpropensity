# Fit an ordinary least-squares model bundle

Fits [`stats::lm()`](https://rdrr.io/r/stats/lm.html) to one dataset, or
once per imputation of stacked long-form data, and combines imputations
with Rubin's rules. A single fit is inferred on the t distribution with
its residual degrees of freedom, as `PROC REG` reports it.

## Usage

``` r
fit_linear(
  formula,
  data,
  outcome_col = NULL,
  id_col = "id",
  imputation_col = NULL,
  prediction_col = "fitted"
)
```

## Arguments

- formula:

  Model formula. The response may be transformed, e.g. `log(los) ~ age`.
  A `.` on the right-hand side expands to the data's other columns,
  excluding the response, `id_col` and `imputation_col`.

- data:

  A data frame, or stacked imputations when `imputation_col` is set.

- outcome_col:

  Outcome column. By default the one variable in the formula's response;
  if supplied, it must name that variable. The response must read
  exactly one column.

- id_col:

  Patient identifier column. Required only when `imputation_col` is set.

- imputation_col:

  Imputation-index column, or `NULL` for one dataset.

- prediction_col:

  Name of the fitted-value column added to `$data`, averaged across
  imputations.

## Value

An object of class `c("rm_fit", "ps_data")`. `$tables$fit_quality` gives
each fit's residual degrees of freedom, residual standard error and
R-squared; R-squared is `NA` for imputed data, where Rubin's rules
define none.

## See also

[`fit_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/fit_logistic.md),
[`bs_continuous()`](https://ehrlinger.github.io/hvtiRpropensity/reference/bs_continuous.md)

## Examples

``` r
dta <- sample_ps_data(n = 150, seed = 42)
fit <- fit_linear(ef ~ age + female, data = dta)
fit$tables$estimates[, c("term", "estimate", "std.error", "p.value")]
#>          term   estimate  std.error      p.value
#> 1 (Intercept) 68.7280984 5.12425021 2.156296e-32
#> 2         age -0.1867446 0.06779147 6.237331e-03
#> 3      female  0.1259720 1.20495525 9.168075e-01
```
