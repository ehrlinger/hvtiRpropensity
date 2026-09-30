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
  when available; defaults to `TRUE`.

- score_col:

  Name of the propensity score column. Taken from `x$meta$score_col`
  when `x` is a `ps_data` object. Default `"prob_t"`.

- treatment_col:

  Name of the binary treatment column. Taken from `x$meta$treatment_col`
  when `x` is a `ps_data` object`. Default `"tavr"\`.

- treated_level:

  The value of `treatment_col` that is treated. Read from the metadata
  of a `ps_data` object when it records one. `NULL` (default) requires a
  0/1 or logical column and takes 1 or `TRUE`.

## Value

A data frame with one row per trim value and columns:

- `trim`:

  The winsorisation threshold.

- `ess_control`:

  Effective sample size for the control group.

- `ess_treated`:

  Effective sample size for the treated group.

- `max_weight`:

  Maximum weight across all patients after trimming.

- `sd_weight`:

  Standard deviation of weights after trimming.

- `pct_trimmed`:

  Percentage of patients whose weight was winsorised.

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
head(sweep)
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
