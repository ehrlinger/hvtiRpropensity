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
