# Stop when a fitted linear model has an aliased term

[`stats::lm()`](https://rdrr.io/r/stats/lm.html) returns `NA` for a term
it cannot estimate, which would otherwise surface later as an anonymous
non-finite coefficient.

## Usage

``` r
.check_aliased(models)
```

## Arguments

- models:

  Named list of fitted models.

## Value

Invisible `NULL`.
