# Build a tidy SMD table across covariates, through ps_stddiff()

Every covariate is treated as Gaussian, so the denominator is
`sqrt((var1 + var0) / 2)`, and with `weight_col` the variances are the
`PROC MEANS` weighted ones, divided by `n - 1`. This is the one formula
the package uses; see
[`ps_stddiff()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_stddiff.md).
A subset with a group absent, such as the matched rows of a match that
found no pairs, gets `NA` for every covariate rather than an error.

## Usage

``` r
.smd_table(data, treatment, covariates = NULL, weight_col = NULL)
```

## Arguments

- data:

  A data frame.

- treatment:

  Name of the binary treatment column.

- covariates:

  Character vector of covariate column names. If `NULL`, all numeric
  columns other than `treatment` are used.

- weight_col:

  Optional name of a weight column.

## Value

A data frame with columns `variable` and `smd`, rounded to 4 places.
