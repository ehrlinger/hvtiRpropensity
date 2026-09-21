# Combine coefficient vectors and covariance matrices with Rubin's rules

Combine coefficient vectors and covariance matrices with Rubin's rules

## Usage

``` r
.rubin_pool(estimates, covariances, conf_level = 0.95)
```

## Arguments

- estimates:

  List of identically named coefficient vectors.

- covariances:

  List of identically named covariance matrices.

- conf_level:

  Confidence level.

## Value

A list with final estimates, covariance and per-imputation values.
