# Generate sample data with a count outcome for bs_count() examples

Produces a synthetic data frame with a Poisson-distributed count outcome
(`rbc_tot`, total intra-operative RBC units transfused), along with
baseline covariates and an optional multiply-imputed stacked format.
Mirrors the dataset used in `tp.pm.count.balncing_score.sas`.

## Usage

``` r
sample_ps_data_count(n = 500L, seed = 42L, m = 1L, n_imputations)
```

## Arguments

- n:

  Number of patients. Default `500L`.

- seed:

  Random seed. Default `42L`. The fixed default gives the same demo
  dataset on every call. `NULL` draws from the caller's random-number
  stream, as
  [`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md)
  does: the data differ from call to call, and
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before the call
  reproduces them. A negative seed gives the same data as its absolute
  value.

- m:

  Number of imputations. If `> 1`, returns `m` datasets stacked in the
  long format, with an `_IMPUTATION_` index column. Each is an
  independent draw of `n` patients, covariates and outcome alike, that
  reuses the `id` values `1` to `n`. They exercise the stacked format;
  they do not resemble real imputations, in which one patient's observed
  values and outcome agree across datasets. Default `1L` (single
  complete dataset).

- n_imputations:

  Deprecated: the old name of `m`. It still works, and warns once per
  session, for one release.

## Value

A data frame with columns: `id`, `rbc_tot` (count), `age`, `female`,
`diabetes`, `hypertension`, `hct_pr` (pre-op haematocrit), `cpb_time`
(cardiopulmonary bypass time, minutes), and optionally `_IMPUTATION_`
when `m > 1`.

## See also

[`bs_count()`](https://ehrlinger.github.io/hvtiRpropensity/reference/bs_count.md),
[`bs_continuous()`](https://ehrlinger.github.io/hvtiRpropensity/reference/bs_continuous.md)

## Examples

``` r
dta <- sample_ps_data_count(n = 200, seed = 7)
hist(dta$rbc_tot, main = "RBC units distribution")


# Stacked multiply-imputed version (3 imputations)
dta_mi <- sample_ps_data_count(n = 200, seed = 7, m = 3)
table(dta_mi[["_IMPUTATION_"]])
#> 
#>   1   2   3 
#> 200 200 200 
```
