# Resolve `outcome_col` from the formula's response

The response must read exactly one column, whether or not `outcome_col`
is supplied: `I(count + other)` or `cbind(count, other)` would fit one
quantity while the stable-response check, and the metadata, used
another. A transformed response such as `log(ef)` reads one column,
`ef`. A supplied `outcome_col` must name that column.

## Usage

``` r
.check_outcome_response(formula, outcome_col = NULL)
```

## Arguments

- formula:

  Two-sided model formula.

- outcome_col:

  Supplied outcome column name, or `NULL`.

## Value

The outcome column name.
