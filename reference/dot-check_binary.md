# Assert that a column contains only 0/1 (or logical) values

Emits a warning (not an error) if NA values are present, because NA
treatment assignments silently drop patients from downstream matching
and weighting steps.

## Usage

``` r
.check_binary(data, col, call_env = rlang::caller_env(), warn_na = TRUE)
```

## Arguments

- data:

  A data frame.

- col:

  Name of the column to check.

- call_env:

  Environment for the error call.

- warn_na:

  Whether to warn about missing values. A caller that counts missing
  treatment and missing score together, and states the total in one
  warning, sets this to `FALSE`.

## Value

Invisible NULL.
