# Standardized differences between two groups

Computes the standardized difference, group 1 minus group 0, for each
named variable, choosing the formula by the variable's type. A port of
the CCF `%stddiff` SAS macro (Artis 2019, from Yang 2012), and intended
to reproduce its `est.stddiff` output.

## Usage

``` r
ps_stddiff(
  data,
  treatment_col = "tavr",
  gaussian = NULL,
  nong_ord = NULL,
  binary = NULL,
  categorical = NULL,
  weight_col = NULL
)
```

## Arguments

- data:

  A data frame.

- treatment_col:

  Name of the group column, coded 0/1 or logical, with both groups
  present. Default `"tavr"`.

- gaussian, nong_ord, binary, categorical:

  Character vectors of column names, by type. See Details. At least one
  variable is required, and a variable may appear under one type only. A
  `binary` column must be coded 0/1 or logical; a `categorical` column
  may have any coding.

- weight_col:

  Optional name of a column of positive weights, for example matching
  weights from
  [`ps_weight()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_weight.md).
  Weights must be finite.

## Value

An object of class `c("ps_stddiff", "ps_data")` with:

- `$data`:

  The input data frame, unchanged.

- `$meta`:

  Named list: `treatment_col`, `weight_col`, `variables` (a list by
  type), `method`, `n_total`, `n_dropped`.

- `$tables`:

  Named list with `stddiff`, a data frame with one row per variable in
  the order given: `variable`, `label` (the column's `"label"`
  attribute, or `NA`), `type` (`"gaussian"`, `"nong_ord"`, `"binary"`,
  `"categorical"` or `"onelevel"`) and `stddiff`, which is `NA` when the
  pooled variance is zero.

## Details

**Variable types.** The type decides the formula, so it is given, not
inferred: an ordinal NYHA class stored as an integer would otherwise be
treated as Gaussian.

- `gaussian`: difference in means over `sqrt((var1 + var0) / 2)`. The
  denominator averages the two variances; it does not pool them by
  sample size, so it differs from a Cohen's d when the groups are
  unequal.

- `nong_ord`: the same formula on ranks of the pooled sample, ties
  averaged, as `PROC RANK` does by default. Ranks are unweighted; only
  their means and variances use `weight_col`.

- `binary` and `categorical` are routed by the number of levels
  observed, as the macro routes its `BINARY=` and `CATG=` lists: one
  level gives 0, two use the binary formula on the proportion at the
  last level in sort order (`1` for a 0/1 variable, the last level for a
  factor), and more than two use the Mahalanobis form of Yang and Dalton
  (2012). For that form the last level in sort order is the one dropped,
  and `S` has diagonal `(t_i (1 - t_i) + c_i (1 - c_i)) / 2` and
  off-diagonal `-(t_i t_j + c_i c_j) / 2`. A categorical difference is
  never negative. When `S` is singular, as when the two groups share no
  levels, the difference is `NA` with a warning naming the variable.

**Weights.** Means are `sum(w x) / sum(w)`; variances are
`sum(w (x - m)^2) / (n - 1)`, which is what `PROC MEANS` reports with a
`WEIGHT` statement. Weights must be positive.

**Missing values.** Rows with a missing group are dropped. Otherwise a
missing value is dropped for that variable only, and a variable with no
non-missing values gets `NA`.

**Units.** The result is a proportion of a standard deviation, not a
percent. The SAS output labels the column "(%)" but does not multiply by
100.

## References

Yang D, Dalton JE. A unified approach to measuring the effect size
between two groups using SAS. SAS Global Forum 2012, paper 335-2012.

## See also

[`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md),
[`ps_weight()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_weight.md)

## Examples

``` r
dta <- sample_ps_data(n = 200, seed = 42)
obj <- ps_stddiff(dta, treatment_col = "tavr",
                  gaussian = c("age", "ef"), binary = c("female", "diabetes"))
obj$tables$stddiff
#>   variable label     type     stddiff
#> 1      age  <NA> gaussian -0.79935800
#> 2       ef  <NA> gaussian  0.61385077
#> 3   female  <NA>   binary  0.39993164
#> 4 diabetes  <NA>   binary -0.05555041
```
