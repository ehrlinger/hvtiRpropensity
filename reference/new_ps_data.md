# Construct a validated ps_data object

Internal entry point used by every `ps_*()` function to create its
return value. Enforces the four-slot contract (`$data`, `$meta`,
`$tables`, `$models`) and attaches the two-level S3 class vector.

## Usage

``` r
new_ps_data(data, meta, tables = list(), models = list(), subclass)
```

## Arguments

- data:

  A data frame — the original data with propensity scores or weights
  appended.

- meta:

  A named list of metadata (column names, formula, method parameters,
  computed statistics, etc.).

- tables:

  A named list of diagnostic objects (SMD tables, group counts,
  effective N, etc.). May be
  [`list()`](https://rdrr.io/r/base/list.html).

- models:

  A named list of fitted model objects. May be
  [`list()`](https://rdrr.io/r/base/list.html).

- subclass:

  A single string naming the specific subclass (e.g. `"ps_match"`).

## Value

A named list of class `c(subclass, "ps_data")`.
