# Covariance of a proportional-odds fit from the expected information

[`MASS::polr()`](https://rdrr.io/pkg/MASS/man/polr.html) inverts the
observed information (the numerical Hessian of the log likelihood). SAS
`PROC LOGISTIC` reports the inverse of the expected (Fisher)
information, and for the cumulative logit, which is not a canonical
link, the two differ at the estimate. This computes the expected
information analytically, so ordinal standard errors match SAS: \\I =
\sum_i w_i \sum_j g\_{ij} g\_{ij}^T / p\_{ij}\\, where \\p\_{ij}\\ is
the fitted probability of category \\j\\ and \\g\_{ij}\\ its gradient in
the coefficients and thresholds.

## Usage

``` r
.polr_expected_vcov(fit)
```

## Arguments

- fit:

  A [`MASS::polr()`](https://rdrr.io/pkg/MASS/man/polr.html) fit with
  the logistic link and its model frame.

## Value

The covariance matrix, named and ordered as `vcov(fit)`.

## Details

The category probabilities are computed without cancellation, so extreme
linear predictors give finite contributions, and an
[`offset()`](https://rdrr.io/r/stats/offset.html) enters through the
fitted linear predictor.
