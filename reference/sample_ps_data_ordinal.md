# Generate sample data with an ordinal treatment for ps_ordinal() examples

Produces a synthetic data frame with a three-level ordered treatment
variable (`nyha_grp`: NYHA functional classes I, II, III), mimicking a
typical cardiac-surgery dataset used with
`tp.lm.logistic_propensity_ordinal.sas`.

## Usage

``` r
sample_ps_data_ordinal(n = 300L, seed = 42L)
```

## Arguments

- n:

  Number of patients per treatment level. Total rows = `3 * n`. Default
  `300L`.

- seed:

  Random seed. Default `42L`. The fixed default gives the same demo
  dataset on every call. `NULL` draws from the caller's random-number
  stream, as
  [`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md)
  does: the data differ from call to call, and
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before the call
  reproduces them.

## Value

A data frame with `3 * n` rows and columns: `id`, `nyha_grp` (ordered
factor: I \< II \< III), `age`, `female`, `ef`, `diabetes`,
`hypertension`.

## See also

[`ps_ordinal()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_ordinal.md)

## Examples

``` r
dta <- sample_ps_data_ordinal(n = 100, seed = 7)
table(dta$nyha_grp)
#> 
#>   I  II III 
#> 100 100 100 
str(dta)
#> 'data.frame':    300 obs. of  7 variables:
#>  $ id          : int  1 2 3 4 5 6 7 8 9 10 ...
#>  $ nyha_grp    : Ord.factor w/ 3 levels "I"<"II"<"III": 1 1 1 1 1 1 1 1 1 1 ...
#>  $ age         : num  82.6 51.2 55.8 58.3 53.3 ...
#>  $ female      : int  1 1 1 0 0 1 0 0 1 0 ...
#>  $ ef          : num  60.2 71.8 55.7 61.3 68.4 ...
#>  $ diabetes    : int  0 0 0 0 0 0 0 0 0 0 ...
#>  $ hypertension: int  1 0 0 1 1 1 0 1 0 1 ...
```
