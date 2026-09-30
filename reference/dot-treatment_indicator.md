# Recode a two-valued treatment column to a 0/1 indicator

With `treated_level = NULL` the column must be 0/1 or logical, and 1 or
`TRUE` is the treated value. With a `treated_level`, the column may hold
any two values and `treated_level` names the treated one. Both values
must be present.

## Usage

``` r
.treatment_indicator(
  data,
  col,
  treated_level = NULL,
  call_env = rlang::caller_env()
)
```

## Arguments

- data:

  A data frame.

- col:

  Name of the treatment column.

- treated_level:

  The treated value, or `NULL`.

- call_env:

  Environment for the error call.

## Value

A list: `trt`, an integer 0/1 vector (`NA` where the column is missing),
and `labels`, the reference then the treated value as character.
