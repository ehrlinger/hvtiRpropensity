# Stop when `outcome_col` is not the formula's response

A mismatched `outcome_col` would fit one column but run the
stable-response check, and record the metadata, against another. A
transformed response such as `log(ef)` is matched by the column it
reads, `ef`.

## Usage

``` r
.check_outcome_response(formula, outcome_col)
```

## Arguments

- formula:

  Two-sided model formula.

- outcome_col:

  Supplied outcome column name.

## Value

Invisible `NULL`.
