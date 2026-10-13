# Expand a dot formula over the model columns only

[`stats::terms()`](https://rdrr.io/r/stats/terms.html) expands `.` to
every other column, which would fit the patient identifier and
imputation index as covariates. The expansion is done once, against
`data` without those columns and without the response variables, and the
explicit formula is returned.

## Usage

``` r
.expand_dot_formula(formula, data, exclude = NULL)
```

## Arguments

- formula:

  Two-sided model formula.

- data:

  Model data.

- exclude:

  Column names that `.` must not expand to, such as the identifier and
  imputation columns. `NULL` entries are ignored.

## Value

`formula` unchanged when its right-hand side has no `.`, otherwise the
expanded formula, in the original formula's environment.
