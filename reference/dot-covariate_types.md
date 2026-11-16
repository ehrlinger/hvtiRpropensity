# Standardized-difference type of each covariate

A logical column, or a numeric column whose non-missing values are all 0
or 1, is `"binary"`; any other numeric column is `"gaussian"`; a factor
or character column is `"categorical"`. Any other type aborts.

## Usage

``` r
.covariate_types(data, covariates)
```

## Arguments

- data:

  A data frame.

- covariates:

  Covariate column names.

## Value

A character vector named by `covariates`.
