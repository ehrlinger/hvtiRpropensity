# Generate sample data with a nominal treatment for ps_nominal() examples

Produces a synthetic data frame with a four-level nominal treatment
variable (`rtyp`: tricuspid repair types COS, PER, DEV, CE), mimicking
the dataset used in `tp.lm.logistic_propensity.polytomous.sas`.

## Usage

``` r
sample_ps_data_nominal(n = 200L, seed = 42L)
```

## Arguments

- n:

  Number of patients per treatment level. Total rows = `4 * n`. Default
  `200L`.

- seed:

  Random seed. Default `42L`. The fixed default gives the same demo
  dataset on every call. `NULL` draws from the caller's random-number
  stream, as
  [`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md)
  does: the data differ from call to call, and
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before the call
  reproduces them.

## Value

A data frame with `4 * n` rows and columns: `id`, `rtyp` (unordered
factor), `age`, `female`, `ef`, `diabetes`, `hypertension`.

## See also

[`ps_nominal()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_nominal.md)

## Examples

``` r
dta <- sample_ps_data_nominal(n = 100, seed = 7)
table(dta$rtyp)
#> 
#> COS PER DEV  CE 
#> 100 100 100 100 
str(dta)
#> 'data.frame':    400 obs. of  7 variables:
#>  $ id          : int  1 2 3 4 5 6 7 8 9 10 ...
#>  $ rtyp        : Factor w/ 4 levels "COS","PER","DEV",..: 1 1 1 1 1 1 1 1 1 1 ...
#>  $ age         : num  80.9 46 51.1 53.9 48.3 ...
#>  $ female      : int  0 0 0 1 1 0 1 1 0 1 ...
#>  $ ef          : num  52.8 65.9 47.8 54.1 62.1 ...
#>  $ diabetes    : int  0 0 0 0 0 0 0 0 0 0 ...
#>  $ hypertension: int  1 0 0 1 1 1 0 1 0 1 ...
```
