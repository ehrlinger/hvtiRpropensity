# Rosenbaum Sensitivity Bounds for Matched Analyses

Assesses how sensitive a matched-pair analysis is to unmeasured
confounding using the Rosenbaum bounds framework. For each value of the
sensitivity parameter Gamma – the maximum odds ratio for differential
treatment assignment within a matched pair due to hidden bias – the
function reports the upper and lower bounds on the p-value of the
Wilcoxon signed-rank test.

## Usage

``` r
sa_rosenbaum(
  x,
  outcome_col,
  gamma_max = 3,
  gamma_inc = 0.25,
  alpha = 0.05,
  treatment_col = NULL,
  pair_id_col = "pair_id",
  treated_level = NULL,
  alternative = c("auto", "greater", "less")
)
```

## Arguments

- x:

  A `ps_match` object (returned by
  [`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md))
  **or** a plain data frame whose rows are matched patients.

- outcome_col:

  Name of the numeric outcome column. Required.

- gamma_max:

  Maximum Gamma to test. Default `3`.

- gamma_inc:

  Increment between Gamma values. Default `0.25`.

- alpha:

  Significance threshold. Default `0.05`.

- treatment_col:

  Name of the binary treatment column (0/1). Taken from
  `x$meta$treatment_col` when `x` is a `ps_match` object. Required when
  `x` is a plain data frame.

- pair_id_col:

  Name of the pair-ID column created by
  [`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md).
  Taken from `x$meta$pair_id_col` when recorded. Default `"pair_id"`.

- treated_level:

  The value of `treatment_col` that is treated. Read from the metadata
  of a `ps_data` object when it records one. `NULL` (default) requires a
  0/1 or logical column and takes 1 or `TRUE`.

- alternative:

  Direction of the test: `"auto"` (default; chosen from the sign of the
  median difference, see Details), `"greater"` (treated outcome above
  control) or `"less"`.

## Value

An object of class `c("sa_rosenbaum", "ps_data")`:

- `$data`:

  The input data frame (`x$data` for a `ps_match` object).

- `$meta`:

  Column names, `treated_level`, the settings (`gamma`, the grid tested,
  `gamma_max`, `gamma_inc`, `alpha`, `alternative`, the direction
  tested, and `direction_from_data`), and the results:
  `sensitivity_value`, the largest Gamma for which the upper p-value is
  below `alpha` (`NA` if the analysis is not significant at Gamma = 1),
  `n_pairs`, the complete matched pairs used, `n_dropped`, `n_nonzero`
  and the statistic `t_plus`.

- `$tables$bounds`:

  Data frame with columns `gamma`, `p_upper`, `p_lower`, `reject_upper`
  (logical – is the upper bound \< `alpha`?).

- `$tables$pairs`:

  One row per complete pair: `pair_id`, `treated` and `control`
  outcomes, and `difference`.

- `$models`:

  An empty list.

## Details

A Gamma of 1 corresponds to no hidden bias (equivalent to perfect
randomisation within pairs). The *sensitivity value* is the largest
Gamma at which the upper p-value bound remains below `alpha`; a larger
sensitivity value indicates a more robust result.

The statistic is the Wilcoxon signed-rank sum \\T^+\\ of the ranks of
\\\|d_i\|\\ over pairs with \\d_i \> 0\\, where \\d_i\\ is the treated
minus the control outcome and zero differences are dropped. Tied
\\\|d_i\|\\ take their average rank. Under hidden bias Gamma each pair
is positive with probability at most \\p^+ = \Gamma / (1 + \Gamma)\\, so
the worst-case null distribution has mean \\p^+ \sum r_i\\ and variance
\\p^+ (1 - p^+) \sum r_i^2\\, and the best case uses \\p^- = 1 / (1 +
\Gamma)\\. Using \\\sum r_i^2\\ rather than \\n(n+1)(2n+1)/6\\ corrects
the variance for ties; without ties the two agree. The p-values use the
normal approximation with a continuity correction of 0.5, so at Gamma =
1 they equal the one-sided
[`stats::wilcox.test()`](https://rdrr.io/r/stats/wilcox.test.html)
p-value with `exact = FALSE` and `correct = TRUE`.

Rosenbaum's method assumes the direction of the test is fixed before the
outcome is seen. Set `alternative` to `"greater"` (treated above
control) or `"less"` to state it. The default `"auto"` keeps the earlier
behaviour and picks the direction from the data: it tests the lower
tail, with a warning, when the median nonzero difference is negative. A
direction chosen from the data makes the bounds slightly optimistic, so
prefer an explicit `alternative` in a pre-specified analysis.

Column names and the treated level are read from a
[`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md)
object's metadata. A supplied `treatment_col`, `pair_id_col` or
`treated_level` must agree with that metadata; a conflict aborts.

## References

Rosenbaum, P. R. (2002). *Observational Studies* (2nd ed.). Springer.

## See also

[`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md),
[`sa_evalue()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_evalue.md),
[`sa_overlap()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_overlap.md),
[`sa_trim_sweep()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_trim_sweep.md)

## Examples

``` r
dta <- sample_ps_data(n = 300, seed = 1)
m   <- ps_match(dta)
res <- sa_rosenbaum(m, outcome_col = "ef", gamma_max = 2, alternative = "greater")
res
#> <sa_rosenbaum>
#>   Outcome     : ef (greater, 300 pairs)
#>   Sensitivity : Gamma = 1.75 at alpha = 0.05
#>  gamma      p_upper      p_lower reject_upper
#>   1.00 5.990619e-11 5.990619e-11         TRUE
#>   1.25 7.884045e-07 2.220446e-16         TRUE
#>   1.50 2.257868e-04 0.000000e+00         TRUE
#>   1.75 7.391183e-03 0.000000e+00         TRUE
#>   2.00 6.394389e-02 0.000000e+00        FALSE
res$tables$bounds
#>   gamma      p_upper      p_lower reject_upper
#> 1  1.00 5.990619e-11 5.990619e-11         TRUE
#> 2  1.25 7.884045e-07 2.220446e-16         TRUE
#> 3  1.50 2.257868e-04 0.000000e+00         TRUE
#> 4  1.75 7.391183e-03 0.000000e+00         TRUE
#> 5  2.00 6.394389e-02 0.000000e+00        FALSE
res$meta$sensitivity_value
#> [1] 1.75
```
