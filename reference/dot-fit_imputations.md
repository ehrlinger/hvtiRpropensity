# Fit one model per imputation with strict patient alignment

Fit one model per imputation with strict patient alignment

## Usage

``` r
.fit_imputations(
  data,
  response_col,
  id_col,
  imputation_col,
  fit_fn,
  predict_fn,
  predictor_cols = NULL,
  require_stable_response = TRUE
)
```

## Arguments

- data:

  A complete data frame or stacked long-form imputations.

- response_col:

  Response column whose patient values must remain stable.

- id_col:

  Patient identifier column.

- imputation_col:

  Imputation-index column, or `NULL` for one dataset.

- fit_fn:

  Function accepting one data-frame partition.

- predict_fn:

  Function accepting a fit and its aligned partition.

- predictor_cols:

  Predictor columns whose categorical levels must agree across
  imputations.

- require_stable_response:

  Whether each patient's response must agree across imputations.

## Value

A list containing base data, averaged predictions, fitted models, fit
status and imputation labels.
