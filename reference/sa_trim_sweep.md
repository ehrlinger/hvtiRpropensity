# IPTW Weight-Trimming Sensitivity Sweep

Re-computes IPTW weights across a range of winsorisation thresholds and
reports how the weight distribution and effective sample size (ESS)
change as trimming increases. This helps identify a trimming threshold
that adequately controls extreme weights without an unacceptable loss of
effective sample size.

## Usage

``` r
sa_trim_sweep(
  x,
  trim_seq = seq(0, 0.1, by = 0.01),
  estimand = c("ATE", "ATT", "ATC"),
  stabilise = TRUE,
  score_col = "prob_t",
  treatment_col = "tavr",
  treated_level = NULL
)
```

## Arguments

- x:

  A `ps_weight` object (returned by
  [`ps_weight()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_weight.md))
  **or** a `ps_match` object **or** a plain data frame that contains a
  propensity score column and a binary treatment column.

- trim_seq:

  Numeric vector of trimming values to evaluate. Each value is the lower
  tail quantile; weights outside the (`trim`, `1 - trim`) range are
  winsorised. `0` means no trimming. Default `seq(0, 0.10, by = 0.01)`.

- estimand:

  One of `"ATE"`, `"ATT"`, or `"ATC"` – the causal estimand used to
  compute IPTW weights. Taken from `x$meta$estimand` when `x` is a
  `ps_weight` object; defaults to `"ATE"` otherwise.

- stabilise:

  Logical. Whether to stabilise weights. Taken from `x$meta$stabilised`
  when available; defaults to `TRUE`. As in
  [`ps_weight()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_weight.md),
  stabilisation applies to the ATE only.

- score_col:

  Name of the propensity score column. Taken from `x$meta$score_col`
  when `x` is a `ps_data` object. Default `"prob_t"`. The score of a
  [`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md)
  object built with `scale = "logit"` is a logit, and is
  back-transformed with
  [`stats::plogis()`](https://rdrr.io/r/stats/Logistic.html) before use.

- treatment_col:

  Name of the binary treatment column. Taken from `x$meta$treatment_col`
  when `x` is a `ps_data` object. Default `"tavr"`.

- treated_level:

  The value of `treatment_col` that is treated. Read from the metadata
  of a `ps_data` object when it records one. `NULL` (default) requires a
  0/1 or logical column and takes 1 or `TRUE`.

## Value

An object of class `c("sa_trim_sweep", "ps_data")`:

- `$data`:

  The input data frame (`x$data` for a `ps_data` object), unchanged.

- `$meta`:

  `score_col`, `treatment_col`, `treated_level`, `scale`, `estimand`,
  `stabilised` (whether stabilisation was applied, so `FALSE` for the
  ATT and ATC), `trim_seq` (sorted), `n_total` and `n_used`, the
  patients with a non-missing treatment and score.

- `$tables$sweep`:

  Data frame with one row per trim value and columns `trim` (the
  winsorisation threshold), `ess_control` and `ess_treated` (effective
  sample sizes), `max_weight` (maximum weight after trimming),
  `sd_weight` (standard deviation of the weights after trimming) and
  `pct_trimmed` (percentage of patients whose weight was winsorised).

- `$models`:

  An empty list.

## Details

For a `ps_data` input, a supplied `estimand`, `stabilise`, `score_col`,
`treatment_col` or `treated_level` must agree with the object's
metadata; a conflict aborts. Arguments not supplied are read from the
metadata, then take their defaults.

## See also

[`ps_weight()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_weight.md),
[`sa_rosenbaum()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_rosenbaum.md),
[`sa_evalue()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_evalue.md),
[`sa_overlap()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_overlap.md)

## Examples

``` r
dta    <- sample_ps_data(n = 400, seed = 3)
w_obj  <- ps_weight(dta, estimand = "ATE")
sweep  <- sa_trim_sweep(w_obj)
head(sweep$tables$sweep)
#>   trim ess_control ess_treated max_weight sd_weight pct_trimmed
#> 1 0.00       303.1       305.0     4.4440    0.5578           0
#> 2 0.01       312.6       312.7     3.1193    0.5209           2
#> 3 0.02       316.2       315.5     2.8898    0.5068           4
#> 4 0.03       320.7       318.9     2.7022    0.4895           6
#> 5 0.04       328.8       325.9     2.4054    0.4558           8
#> 6 0.05       339.1       334.9     2.0931    0.4126          10

# Directly from a data frame
sweep2 <- sa_trim_sweep(dta, score_col = "prob_t", treatment_col = "tavr",
                        estimand = "ATE")
```
