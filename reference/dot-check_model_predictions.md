# Validate predictions against their input partition

Validate predictions against their input partition

## Usage

``` r
.check_model_predictions(predictions, n, imputation, columns = NULL)
```

## Arguments

- predictions:

  Numeric vector or matrix.

- n:

  Expected number of rows.

- imputation:

  Imputation label used in errors.

- columns:

  Expected matrix column names, or `NULL` for the first fit.

## Value

Matrix column names, or `NULL` for vector predictions.
