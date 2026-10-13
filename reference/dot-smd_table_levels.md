# Standardized differences between pairs of treatment levels

Each pair is compared on the patients in those two levels only, with
`level` as the treated group and `versus` as the comparison group.

## Usage

``` r
.smd_table_levels(
  data,
  treatment,
  levels,
  scheme,
  covariates = NULL,
  reserved = character(),
  keep = character()
)
```

## Arguments

- data:

  A data frame.

- treatment:

  Name of the treatment column.

- levels:

  Treatment levels, reference first.

- scheme:

  Any of `"reference"`, `"adjacent"`, `"all"`.

- covariates:

  Covariate columns. `NULL` uses
  [`.default_covariates()`](https://ehrlinger.github.io/hvtiRpropensity/reference/dot-default_covariates.md),
  excluding `treatment` and `reserved`.

- reserved:

  Columns never treated as covariates.

- keep:

  Model predictors protected in the default set; see
  [`.default_covariates()`](https://ehrlinger.github.io/hvtiRpropensity/reference/dot-default_covariates.md).

## Value

A data frame with columns `variable`, `level`, `versus`, `smd`.
