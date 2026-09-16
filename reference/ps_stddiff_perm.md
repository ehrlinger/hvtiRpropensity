# Permutation reference for standardized differences

Computes
[`ps_stddiff()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_stddiff.md)
on the observed groups and on `n_perm` random permutations of the group
labels, and reports the observed value beside percentiles of its
permutation distribution. A port of the CCF `%stddiffci` SAS macro
(Artis 2020).

## Usage

``` r
ps_stddiff_perm(
  data,
  treatment_col = "tavr",
  gaussian = NULL,
  nong_ord = NULL,
  binary = NULL,
  categorical = NULL,
  weight_col = NULL,
  reweight = NULL,
  n_perm = 1000L,
  seed = NULL
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

  Optional name of the column
  [`ps_stddiff()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_stddiff.md)
  weights by. Requires `reweight`, which supplies its values.

- reweight:

  A function of a data frame returning one positive, finite weight per
  row, for example one that refits a propensity model and returns
  matching weights. Required when `weight_col` is given, and rejected
  without it.

- n_perm:

  Number of permutations. Default `1000`, the macro's.

- seed:

  Optional integer seed.

## Value

An object of class `c("ps_stddiff_perm", "ps_data")` with:

- `$data`:

  The input data frame, unchanged.

- `$meta`:

  Named list: `treatment_col`, `weight_col`, `variables`, `method`,
  `n_perm`, `seed`, `n_total`.

- `$tables`:

  Named list with `stddiff_perm`, one row per variable: `variable`,
  `label`, `type`, `observed` (as
  [`ps_stddiff()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_stddiff.md)
  reports it) and the permutation percentiles `p2_5`, `p16`, `p50`,
  `p84`, `p97_5`.

## Details

**What the percentiles are.** They describe the standardized difference
you would see if group membership carried no information about the
variable. They are not a confidence interval for the observed value,
whatever the SAS macro's name suggests: an observed difference outside
the 2.5 to 97.5 range is larger than label shuffling alone tends to
produce. The 16 and 84 percentiles bound the central 68%, the range the
macro's 2021 correction fixed.

**Percentile definition.** `quantile(type = 2)`, which is SAS
`PCTLDEF=5`, the `PROC UNIVARIATE` default the macro relies on.
Permutations where a variable's difference is `NA` are left out of that
variable's percentiles.

**Weights.** A matching or IPTW weight depends on the group, so a weight
carried over unchanged onto permuted groups describes neither the
observed design nor the null. `reweight` must therefore be given
whenever `weight_col` is: a function taking a data frame and returning
one weight per row. It is called once on the unpermuted data, for the
observed row, and once on each permuted data set, and its result is
written to `weight_col` before
[`ps_stddiff()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_stddiff.md)
runs. The macro instead expects the permuted weights to exist as columns
before it is called.

**Missing groups.** Rows with a missing group keep it and are not
permuted into either group.

**Random numbers.** With `seed`, the permutations are reproducible and
the caller's random number stream is restored afterwards. Without it,
the current stream is used and advanced.

## See also

[`ps_stddiff()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_stddiff.md)

## Examples

``` r
dta <- sample_ps_data(n = 200, seed = 42)
obj <- ps_stddiff_perm(dta, treatment_col = "tavr", gaussian = c("age", "ef"),
                       binary = "female", n_perm = 200, seed = 1)
obj$tables$stddiff_perm
#>   variable label     type   observed       p2_5         p16          p50
#> 1      age  <NA> gaussian -0.7993580 -0.1920948 -0.10274954 -0.004489759
#> 2       ef  <NA> gaussian  0.6138508 -0.1746507 -0.08978417  0.003163732
#> 3   female  <NA>   binary  0.3999316 -0.2225974 -0.09059306 -0.010055712
#>          p84     p97_5
#> 1 0.09089048 0.1726795
#> 2 0.11432990 0.1990921
#> 3 0.10068703 0.1919339
```
