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
  seed = 1024L,
  time_unit_days = 365.2425,
  clip = 0.001
)
```

## Arguments

- x:

  A `ps_data` object with a score and treatment column.

- time_col, event_col:

  Follow-up time and event indicator (0/1) columns of `x$data`.

- tau:

  Horizon, in the units of `time_col`.

- weights:

  Weighting schemes to run: any of `"unweighted"`, `"ato"` (overlap),
  `"att"`, `"ate"`.

- subsets:

  Optional named list of logical vectors (one value per row); each is
  run as a restricted analysis (for example a common-support subset).
  The full data is always included as `"all"`.

- n_boot:

  Bootstrap draws. Resampling is stratified by arm.

- refit:

  If `FALSE` (default) the score is held fixed in the bootstrap, which
  understates uncertainty. If `TRUE` the score model is refit on each
  draw (supported for `ps_logistic` and `ps_forest` inputs).

- seed:

  Bootstrap seed.

- time_unit_days:

  Days per unit of `time_col`, for the `diff_days` column. Default
  `365.2425` (time in years).

- clip:

  Scores are clipped to `[clip, 1 - clip]` before weighting.

## Value

A `ps_rmst` / `ps_data` object. `$tables$estimates` has one row per
subset-and-weighting estimator: `estimator`, `subset`, `weighting`, `n`,
`ess_treated`, `ess_control`, `rmst_treated`, `rmst_control`, `diff`
(treated minus control, positive favours treated), `diff_days`,
`lo_days`, `hi_days` (2.5 and 97.5 percent bootstrap) and `n_failed`.
`$tables$curves` has the weighted Kaplan-Meier steps (`estimator`,
`arm`, `time`, `surv`) for plotting; `$data` is `x$data` plus one weight
column `w_<weighting>` per scheme.

## Examples

``` r
if (requireNamespace("survival", quietly = TRUE)) {
  dta <- sample_ps_data(n = 200, seed = 5)[, c("id", "tavr", "age", "ef")]
  dta$t <- stats::rexp(nrow(dta), 0.2)
  dta$e <- as.integer(dta$t < 5)
  dta$t <- pmin(dta$t, 5)
  res <- ps_rmst(ps_logistic(tavr ~ age + ef, dta), "t", "e", tau = 4, n_boot = 20)
  res$tables$estimates
}
#>        estimator subset  weighting   n ess_treated ess_control rmst_treated
#> 1 all/unweighted    all unweighted 400    200.0000    200.0000     2.730575
#> 2        all/ato    all        ato 400    157.8078    155.2728     2.727508
#>   rmst_control        diff diff_days    lo_days  hi_days n_failed
#> 1     2.701130  0.02944547  10.75474  -63.96384 89.13301        0
#> 2     2.820534 -0.09302584 -33.97699 -151.78152 69.04643        0
```
