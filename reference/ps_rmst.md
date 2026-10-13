# Observed-outcome restricted mean survival contrast, weighted

Design-based comparison of restricted mean survival time (RMST) between
treated and control patients using only the observed follow-up: weighted
Kaplan-Meier curves (weights from the propensity score) integrated to
`tau`, differenced, with a stratified bootstrap interval.

## Usage

``` r
ps_rmst(
  x,
  time_col,
  event_col,
  tau,
  weights = c("unweighted", "ato"),
  subsets = list(),
  n_boot = 200L,
  refit = FALSE,
  seed = NULL,
  time_unit_days = 365.2425,
  clip = 0.001
)
```

## Arguments

- x:

  A `ps_data` object with a score and treatment column. A
  [`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md)
  object is analysed on its matched pairs only, and a
  [`ps_weight()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_weight.md)
  object contributes its stored weight column; see Details. A score on
  the logit scale (`scale = "logit"`, as a
  [`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md)
  object can carry) is back-transformed with
  [`stats::plogis()`](https://rdrr.io/r/stats/Logistic.html) before any
  weights are built.

- time_col, event_col:

  Follow-up time and event indicator (0/1) columns of `x$data`.

- tau:

  Horizon, in the units of `time_col`. It must not exceed the largest
  observed time among the rows that enter either arm's curve, for any
  analysed subset and weighting.

- weights:

  Weighting schemes to run: any of `"unweighted"`, `"ato"` (overlap),
  `"att"`, `"ate"`. For a
  [`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md)
  input the default is `"unweighted"` alone, and no other scheme is
  accepted.

- subsets:

  Optional named list of logical vectors (one value per row); each is
  run as a restricted analysis (for example a common-support subset).
  The full data is always included as `"all"`.

- n_boot:

  Bootstrap draws. Resampling is by matched pair for a
  [`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md)
  input and by row, stratified by arm, otherwise.

- refit:

  If `FALSE` (default) the score is held fixed in the bootstrap, so the
  intervals ignore the uncertainty in estimating the propensity score
  and understate the total uncertainty. If `TRUE` the score model is
  refit on each draw (supported for `ps_logistic` and `ps_forest`
  inputs).

- seed:

  Bootstrap seed: `NULL` (default) or one whole number. With a seed, the
  random number stream is set to `abs(seed)` immediately before the
  draws, the result is reproducible whatever ran earlier in the session,
  the caller's stream is restored afterwards, and the seed is kept in
  `$meta$seed`. With `NULL` the draws come from the caller's stream and
  advance it, and `$meta$seed` is `NA`.

- time_unit_days:

  Days per unit of `time_col`, for the `diff_days` column. Default
  `365.2425` (time in years).

- clip:

  Scores are clipped to `[clip, 1 - clip]` before weighting.

## Value

A `ps_rmst` / `ps_data` object. `$tables$estimates` has one row per
subset-and-weighting estimator: `estimator`, `subset`, `weighting`, `n`
(rows analysed; for a matched input, twice the matched pairs),
`ess_treated`, `ess_control`, `rmst_treated`, `rmst_control`, `diff`
(treated minus control, positive favours treated), `diff_days`,
`lo_days`, `hi_days` (2.5 and 97.5 percent bootstrap) and `n_failed`
(bootstrap draws that gave no estimate). `$tables$curves` has the
weighted Kaplan-Meier steps (`estimator`, `arm`, `time`, `surv`) for
plotting; `$data` is `x$data` plus one weight column `w_<weighting>` per
scheme. `$meta$seed` is the seed used (`NA` when none was given).
`$meta$scale` is the scale of `score_col`, carried from `x`
(`"probability"` when `x` records none). `$meta$resampling` is `"pair"`
for a matched input and `"row"` (stratified by arm) otherwise, and
`$meta$stored_weighting` names the scheme that used a
[`ps_weight()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_weight.md)
object's stored weights (`NA` when none did).

## Details

**Weights.** `"unweighted"` gives every row weight 1. `"ato"` weights
treated rows by `1 - p` and controls by `p`, targeting the overlap
population; `"att"` weights controls by `p / (1 - p)`, targeting the
treated; `"ate"` weights by `1 / p` and `1 / (1 - p)`, targeting the
whole cohort. Here `p` is the score clipped to `[clip, 1 - clip]`.

**Matched input.** For a
[`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md)
object only the rows with `match_col == 1` are analysed, so `n` is twice
the number of matched pairs, and the estimand is the one the match
targets. Weighting a matched set by the score adjusts twice for the same
confounders, so only `"unweighted"` is accepted. The bootstrap resamples
whole pairs: it draws `pair_id` values with replacement, as many as
there are pairs, and takes both members of each drawn pair, so a pair
drawn twice enters twice. This keeps the within-pair correlation that a
row bootstrap ignores. A subset that splits a pair keeps the member it
contains, which is then resampled alone.

**Weighted input.** For a
[`ps_weight()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_weight.md)
object, the scheme that matches its estimand (`"ate"` for
`estimand = "ATE"`, `"att"` for `"ATT"`) uses the stored column
`x$meta$weight_col` instead of recomputing weights from the score, so
any trimming is kept. Stabilisation multiplies each arm by a constant,
which leaves each arm's Kaplan-Meier curve and effective sample size
unchanged. The other schemes are computed from the score, and an
`estimand = "ATC"` object has no matching scheme. In the bootstrap each
row keeps its stored weight.

**Horizon.** The weighted Kaplan-Meier curve is not defined past the
last observed time in an arm, counting only the rows that enter the
curve (time and event present, weight above zero), so a `tau` beyond it
stops with an error rather than carrying the last survival value
forward. A bootstrap draw can still end short of `tau`, and its curve is
then carried flat to `tau`.

**Failures.** A bootstrap draw that errors, or that leaves an arm empty,
gives `NA`; the interval is taken over the other draws and the number of
`NA` draws is reported as `n_failed`.

## Examples

``` r
if (requireNamespace("survival", quietly = TRUE)) {
  dta <- sample_ps_data(n = 200, seed = 5)[, c("id", "tavr", "age", "ef")]
  dta$t <- stats::rexp(nrow(dta), 0.2)
  dta$e <- as.integer(dta$t < 5)
  dta$t <- pmin(dta$t, 5)
  res <- ps_rmst(ps_logistic(tavr ~ age + ef, dta), "t", "e", tau = 4, n_boot = 20,
                  seed = 1024)
  res$tables$estimates
}
#>        estimator subset  weighting   n ess_treated ess_control rmst_treated
#> 1 all/unweighted    all unweighted 400    200.0000    200.0000     2.728810
#> 2        all/ato    all        ato 400    157.8078    155.2728     2.737144
#>   rmst_control        diff diff_days   lo_days  hi_days n_failed
#> 1     2.828556 -0.09974563 -36.43134 -158.0537 51.11612        0
#> 2     2.791022 -0.05387833 -19.67865 -109.9290 74.75615        0
```
