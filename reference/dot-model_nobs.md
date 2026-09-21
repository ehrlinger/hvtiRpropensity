# Count rows used by a fitted model

Count rows used by a fitted model

## Usage

``` r
.model_nobs(fit, fallback)
```

## Arguments

- fit:

  A fitted model object.

- fallback:

  Count to return when the model has no
  [`nobs()`](https://rdrr.io/r/stats/nobs.html) method.

## Value

A single integer.
