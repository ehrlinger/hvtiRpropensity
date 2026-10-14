# Combine coefficient vectors and covariance matrices with Rubin's rules

The degrees of freedom are Barnard and Rubin's (1999) small-sample
adjustment, computed as `mice::pool()` computes them. With
`df_com = Inf` they reduce to Rubin's (1987) large-sample degrees of
freedom.

## Usage

``` r
.rubin_pool(estimates, covariances, conf_level = 0.95, df_com = Inf)
```

## Arguments

- estimates:

  List of identically named coefficient vectors.

- covariances:

  List of identically named covariance matrices.

- conf_level:

  Confidence level.

- df_com:

  Complete-data degrees of freedom, from
  [`.complete_df()`](https://ehrlinger.github.io/hvtiRpropensity/reference/dot-complete_df.md).
  Must be positive; `Inf` gives Rubin's (1987) degrees of freedom.

## Value

A list with final estimates, covariance and per-imputation values.
