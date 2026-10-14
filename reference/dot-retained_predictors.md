# Variables the model's retained terms read

[`all.vars()`](https://rdrr.io/r/base/allnames.html) on the right-hand
side also returns variables a formula subtracts, such as `dropme` in
`y ~ . - dropme`. These are the variables of the terms the model keeps,
plus any [`offset()`](https://rdrr.io/r/stats/offset.html).

## Usage

``` r
.retained_predictors(formula)
```

## Arguments

- formula:

  Model formula.

## Value

Character vector of column names, possibly empty.
