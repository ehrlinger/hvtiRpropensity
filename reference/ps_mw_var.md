# Bootstrap variance of a matching-weight treatment effect

Estimates the weighted difference in means, group 1 minus group 0, for
each outcome, with a bootstrap standard deviation, percentiles, z and p.
A port of the CCF `%mw_var` SAS macro (Rajeswaran 2014), for matching
weights as described by Li and Greene (2013).

## Usage

``` r
ps_mw_var(
  data,
  treatment_col = "tavr",
  outcomes,
  weight_col,
  n_rep = 1000L,
  seed = NULL
)
```

## Arguments

- data:

  A data frame.

- treatment_col:

  Name of the group column, coded 0/1 or logical, with both groups
  present. Default `"tavr"`.

- outcomes:

  Character vector of numeric outcome columns.

- weight_col:

  Name of the column of matching weights, which must be positive and
  finite.

- n_rep:

  Number of bootstrap replicates. Default `1000`, the macro's
  `RESAMPL=`.

- seed:

  Optional integer seed.

## Value

An object of class `c("ps_mw_var", "ps_data")` with:

- `$data`:

  The input data frame, unchanged.

- `$meta`:

  Named list: `treatment_col`, `weight_col`, `outcomes`, `method`,
  `n_rep`, `seed`, `n_total`, `n_dropped`.

- `$tables`:

  Named list with `mw_var`, one row per outcome: `outcome`, `label`,
  `sumwgt_0`, `mean_0`, `sd_0`, `sumwgt_1`, `mean_1`, `sd_1`,
  `estimate`, `n_rep`, `sd`, `p2_5`, `p16`, `p50`, `p84`, `p97_5`, `z`
  and `p`.

## Details

**Estimate.** For each group, `PROC MEANS` weighted summaries: the mean
`sum(w y) / sum(w)`, the SD with an `n - 1` denominator, and the sum of
weights. The macro labels the sum of weights `N`; it is reported here as
`sumwgt_0` and `sumwgt_1` because it is not a count. For a 0/1 outcome
the means are weighted proportions.

**Bootstrap.** Each replicate resamples rows with replacement within
each group, to that group's size, and recomputes the weighted means. The
treatment group is the only stratum, as in the macro. `sd` is the
standard deviation of the replicate differences, and the percentiles use
`quantile(type = 4)`, which is SAS `PCTLDEF=1` as the macro sets it.
Replicates where an outcome has no non-missing value in a group are left
out for that outcome, and `n_rep` counts those that remain.

**The weights are held fixed.** Each resampled row keeps its original
weight; the propensity model and the weights are not re-estimated. This
reproduces the macro, and it ignores the uncertainty in estimating the
weights, so `sd` can be too small.

**Missing values.** Rows with a missing group are dropped. A missing
outcome is dropped for that outcome only.

**Random numbers.** With `seed`, the result is reproducible and the
caller's random number stream is restored afterwards.

## References

Li L, Greene T. A weighting analogue to pair matching in propensity
score analysis. The International Journal of Biostatistics.
2013;9(2):215-234.

## See also

[`ps_weight()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_weight.md),
[`ps_stddiff()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_stddiff.md)

## Examples

``` r
dta <- sample_ps_data(n = 200, seed = 42)
dta$mt_wt <- pmin(dta$prob_t, 1 - dta$prob_t) /
  ifelse(dta$tavr == 1, dta$prob_t, 1 - dta$prob_t)
obj <- ps_mw_var(dta, treatment_col = "tavr", outcomes = c("age", "ef"),
                 weight_col = "mt_wt", n_rep = 200, seed = 1)
obj$tables$mw_var
#>   outcome label sumwgt_0   mean_0     sd_0 sumwgt_1   mean_1     sd_1
#> 1     age  <NA> 116.0528 74.50040 5.446871 115.3951 74.37045 5.632652
#> 2      ef  <NA> 116.0528 54.52772 6.903189 115.3951 54.55226 7.246387
#>      estimate n_rep        sd      p2_5        p16         p50       p84
#> 1 -0.12994823   200 0.8665785 -1.939877 -0.9633344 -0.14053587 0.6971630
#> 2  0.02453471   200 1.0395746 -2.040412 -1.1955832 -0.05329846 0.9746622
#>      p97_5           z         p
#> 1 1.442813 -0.14995552 0.8807997
#> 2 1.780236  0.02360072 0.9811711
```
