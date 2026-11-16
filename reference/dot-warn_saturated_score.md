# Warn when a propensity score is numerically 0 or 1

A finite, converged logistic fit can still give a patient with an
extreme covariate value a probability within `10 * .Machine$double.eps`
of 0 or 1, the threshold at which
[`glm.fit()`](https://rdrr.io/r/stats/glm.html) warns.
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) caps the linear
predictor there, so the logit is near +/-36 rather than infinite, but
the value reflects the cap, not the model, and an inverse-probability
weight on the wrong side of it is of order `1e15`. Name those patients
when the score is made, rather than leave it to a later step.

## Usage

``` r
.warn_saturated_score(p, score_col, logit_col, id = NULL)
```

## Arguments

- p:

  Propensity scores.

- score_col, logit_col:

  Column names, for the message.

- id:

  Patient identifiers, or `NULL` for row numbers.

## Value

Invisible `NULL`.
