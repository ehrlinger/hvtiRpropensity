# Generate sample data with a count outcome for bs_count() examples

Produces a synthetic data frame with a Poisson-distributed count outcome
(`rbc_tot`, total intra-operative RBC units transfused), along with
baseline covariates and an optional multiply-imputed stacked format.
Mirrors the dataset used in `tp.pm.count.balncing_score.sas`.

## Usage

``` r
sample_ps_data_count(n = 500L, seed = 42L, n_imputations = 1L)
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

- n_imputations:

  If `> 1`, returns a stacked multiply-imputed data frame with an
  `_IMPUTATION_` index column containing `n_imputations` replicates
  (each with slightly jittered covariate values to mimic imputed
  datasets). Default `1L` (single complete dataset).

## Value

A data frame with columns: `id`, `rbc_tot` (count), `age`, `female`,
`diabetes`, `hypertension`, `hct_pr` (pre-op haematocrit), `cpb_time`
(cardiopulmonary bypass time, minutes), and optionally `_IMPUTATION_`
when `n_imputations > 1`.

## See also

[`bs_count()`](https://ehrlinger.github.io/hvtiRpropensity/reference/bs_count.md),
[`bs_continuous()`](https://ehrlinger.github.io/hvtiRpropensity/reference/bs_continuous.md)

## Examples

``` r
dta <- sample_ps_data_count(n = 200, seed = 7)
hist(dta$rbc_tot, main = "RBC units distribution")


# Stacked multiply-imputed version (3 imputations)
dta_mi <- sample_ps_data_count(n = 200, seed = 7, n_imputations = 3)
table(dta_mi[["_IMPUTATION_"]])
#> 
#>   1   2   3 
#> 200 200 200 
```
