# Stop when a binary logistic fit predicts the outcome perfectly

Under complete or quasi-complete separation the maximum likelihood
estimate does not exist: the likelihood keeps rising as some
coefficients grow without bound, and
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) stops only because
its deviance criterion is met, with a warning at most. The coefficients
and their standard errors are then meaningless. The test is whether the
fit has reached a maximum: 25 further Newton steps are taken from the
reported coefficients, and if any coefficient moves by more than `1e-4`
of its size (or of 1, if smaller) the fit is diverging. At a finite
maximum the steps leave the coefficients where they are, to rounding.

## Usage

``` r
.check_separation(models, outcome_col)
```

## Arguments

- models:

  Named list of fitted models.

- outcome_col:

  The outcome column, for the message.

## Value

Invisible `NULL`.

## Details

A fitted probability that rounds to exactly 0 or 1 is not, alone,
evidence of separation: one extreme but valid covariate value does it
with finite, well-estimated coefficients (PR \#122 review). That case is
left to the callers:
[`ps_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_logistic.md)
warns and names those patients. Only
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) fits are checked.
