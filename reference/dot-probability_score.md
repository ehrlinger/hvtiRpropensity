# The propensity score of a ps_data object on the probability scale

A
[`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md)
object built with `scale = "logit"` names a logit column in
`$meta$score_col`. Readers that need a probability call this helper,
which applies [`stats::plogis()`](https://rdrr.io/r/stats/Logistic.html)
when `$meta$scale` is `"logit"` and otherwise returns the column
unchanged, then checks that every value lies in \[0, 1\].

## Usage

``` r
.probability_score(
  x,
  score_col = x$meta$score_col,
  call_env = rlang::caller_env()
)
```

## Arguments

- x:

  A `ps_data` object.

- score_col:

  Name of the score column. Default `x$meta$score_col`.

- call_env:

  Environment for the error call.

## Value

A numeric vector of probabilities, one per row of `x$data`.
