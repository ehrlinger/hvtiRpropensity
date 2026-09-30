# Fit a general logistic model bundle

Fits binary, proportional-odds ordinal, or generalized-logit nominal
models under an explicit outcome-level contract. Stacked long-form
imputations are fitted separately and combined with Rubin's rules.

## Usage

``` r
fit_logistic(
  formula,
  data,
  family = c("binary", "ordinal", "nominal"),
  outcome_col = NULL,
  id_col = "id",
  imputation_col = NULL,
  outcome_levels,
  event_level = NULL,
  reference_level = NULL,
  prediction_prefix = "prob",
  trace = FALSE
)
```

## Arguments

- formula:

  Model formula.

- data:

  A data frame containing the outcome and predictors.

- family:

  Model family.

- outcome_col:

  Outcome column. By default it is read from `formula`.

- id_col:

  Patient identifier column.

- imputation_col:

  Imputation-index column, or `NULL` for one dataset.

- outcome_levels:

  Complete outcome levels in their declared order.

- event_level:

  Event level for a binary model.

- reference_level:

  Reference level for a nominal model.

- prediction_prefix:

  Name for a binary probability column, or prefix for ordinal and
  nominal level-specific probability columns.

- trace:

  Whether model engines should print fitting progress.

## Value

An object of class `c("lm_fit", "ps_data")`.

## Examples

``` r
dta <- sample_ps_data(n = 150, seed = 42)
fit <- fit_logistic(
  tavr ~ age + female + ef,
  data           = dta,
  family         = "binary",
  outcome_levels = c(0, 1),
  event_level    = 1
)
fit$tables$estimates[, c("term", "estimate", "odds_ratio")]
#>          term    estimate odds_ratio
#> 1 (Intercept)  3.55713646 35.0626500
#> 2         age -0.09113276  0.9128965
#> 3      female  0.18166589  1.1992135
#> 4          ef  0.05839068  1.0601291
head(fit$data[, c("id", "tavr", "prob")])
#>   id tavr      prob
#> 1  1    0 0.1657534
#> 2  2    0 0.5493138
#> 3  3    0 0.2856728
#> 4  4    0 0.2041096
#> 5  5    0 0.5254897
#> 6  6    0 0.3978214
```
