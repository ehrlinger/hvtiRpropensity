# Default covariates for a balance table

Every numeric column except those in `reserved`, the fixed output
columns this package appends (`prob_t`, `logit_t`, `mt_wt`, `iptw`,
`match`, `pair_id`, `quintile`, `decile`, `bs`, `cluster`), and constant
columns. A column is constant when its non-missing values take at most
one value. A function that generates columns under other names, such as
the `prob_<level>` scores of
[`ps_nominal()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_nominal.md),
passes those exact names in `reserved`.

## Usage

``` r
.default_covariates(data, reserved = character(), keep = character())
```

## Arguments

- data:

  A data frame.

- reserved:

  Columns never used as covariates, such as the treatment, score and id
  columns.

- keep:

  Columns protected from the fixed output names and the constant rule,
  typically [`all.vars()`](https://rdrr.io/r/base/allnames.html) of a
  formula's right-hand side.

## Value

A character vector of column names.

## Details

Columns in `keep`, the predictors a model formula declares, are never
dropped by the fixed output names or the constant rule; only `reserved`
can drop them.
