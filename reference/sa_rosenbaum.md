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
  treated_level = NULL
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
  Default `"pair_id"`.

- treated_level:

  The value of `treatment_col` that is treated. Read from the metadata
  of a `ps_data` object when it records one. `NULL` (default) requires a
  0/1 or logical column and takes 1 or `TRUE`.

## Value

A list with three elements:

- `$bounds`:

  Data frame with columns `gamma`, `p_upper`, `p_lower`, `reject_upper`
  (logical – is the upper bound \< `alpha`?).

- `$sensitivity_value`:

  Largest Gamma for which the analysis remains significant (upper
  p-value \< `alpha`), or `NA` if the initial analysis is not
  significant.

- `$n_pairs`:

  Number of complete matched pairs used.

## Details

A Gamma of 1 corresponds to no hidden bias (equivalent to perfect
randomisation within pairs). The *sensitivity value* is the largest
Gamma at which the upper p-value bound remains below `alpha`; a larger
sensitivity value indicates a more robust result.

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
res <- sa_rosenbaum(m, outcome_col = "ef", gamma_max = 2)
res$bounds
#>   gamma      p_upper      p_lower reject_upper
#> 1  1.00 6.990297e-12 6.990297e-12         TRUE
#> 2  1.25 1.513417e-07 0.000000e+00         TRUE
#> 3  1.50 6.328064e-05 0.000000e+00         TRUE
#> 4  1.75 2.813541e-03 0.000000e+00         TRUE
#> 5  2.00 3.141844e-02 0.000000e+00         TRUE
res$sensitivity_value
#> [1] 2
```
