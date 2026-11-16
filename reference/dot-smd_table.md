# Build a tidy SMD table across covariates, through ps_stddiff()

Each covariate is routed to the
[`ps_stddiff()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_stddiff.md)
formula for its type. A logical column, or a numeric column whose
non-missing values are all 0 or 1, takes the binary formula. A factor or
character column takes the categorical formula of Yang and Dalton
(2012), one row per column. With two levels left that distance is the
absolute value of the binary formula, and it is reported that way,
unsigned, so a factor whose third level is absent from the matched rows
keeps the same scale and sign convention as before the match (#122).
Every other numeric column is Gaussian, with the denominator
`sqrt((var1 + var0) / 2)`. With `weight_col`, proportions and variances
are weighted as
[`ps_stddiff()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_stddiff.md)
describes. A subset with a group absent, such as the matched rows of a
match that found no pairs, gets `NA` for every covariate rather than an
error.

## Usage

``` r
.smd_table(
  data,
  treatment,
  covariates = NULL,
  weight_col = NULL,
  reserved = character(),
  keep = character(),
  types = NULL
)
```

## Arguments

- data:

  A data frame.

- treatment:

  Name of the binary treatment column.

- covariates:

  Character vector of covariate column names, used as given. If `NULL`,
  the default set from
  [`.default_covariates()`](https://ehrlinger.github.io/hvtiRpropensity/reference/dot-default_covariates.md)
  is used, excluding `treatment`, `weight_col` and `reserved`.

- weight_col:

  Optional name of a weight column.

- reserved:

  Further columns never used as default covariates.

- keep:

  Model predictors protected in the default set; see
  [`.default_covariates()`](https://ehrlinger.github.io/hvtiRpropensity/reference/dot-default_covariates.md).

- types:

  Named character vector from
  [`.covariate_types()`](https://ehrlinger.github.io/hvtiRpropensity/reference/dot-covariate_types.md).
  A caller that builds two tables from different rows, such as before
  and after a match, types the covariates once on the larger set and
  passes the result to both, so a column is never Gaussian in one table
  and binary in the other. `NULL` types them on `data`.

## Value

A data frame with columns `variable` and `smd`, rounded to 4 places, one
row per covariate.
