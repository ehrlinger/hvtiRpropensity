# Test whether an object is a ps_data instance

Test whether an object is a ps_data instance

## Usage

``` r
is_ps_data(x)
```

## Arguments

- x:

  Any R object.

## Value

`TRUE` if `x` inherits from `"ps_data"`, `FALSE` otherwise.

## Examples

``` r
dta <- sample_ps_data()
obj <- ps_match(dta, treatment_col = "tavr")
is_ps_data(obj)   # TRUE
#> [1] TRUE
is_ps_data(list()) # FALSE
#> [1] FALSE
```
