# Count rows used by a fitted model

Every model this package fits uses `na.action = na.exclude`, so the rows
dropped for a missing value are recorded in the fit whatever the family.
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) has no
[`nobs()`](https://rdrr.io/r/stats/nobs.html) method, so counting from
[`nobs()`](https://rdrr.io/r/stats/nobs.html) would report every input
row as analysed.

## Usage

``` r
.model_nobs(fit, n_input)
```

## Arguments

- fit:

  A fitted model object.

- n_input:

  Number of rows the model was given.

## Value

A single integer.
