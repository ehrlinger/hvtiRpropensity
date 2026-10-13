# Warn once about rows excluded for a missing treatment or score

Warn once about rows excluded for a missing treatment or score

## Usage

``` r
.warn_excluded(excluded, treatment_col, score_col, what)
```

## Arguments

- excluded:

  Logical vector, `TRUE` for a row left out of the computation.

- treatment_col, score_col:

  Column names, for the message.

- what:

  What the excluded rows are left out of, for the message.

## Value

`excluded`, invisibly.
