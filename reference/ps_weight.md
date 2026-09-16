# Inverse-Probability-of-Treatment Weighting (IPTW)

Computes IPTW weights from a pre-computed propensity score column and
appends them to the dataset. Supports ATE, ATT, and ATC estimands.
Returns a `ps_weight` object whose `$tables` slot contains SMD tables
(unweighted and weighted) and effective-N diagnostics compatible with
`hvtiPlotR::hv_mirror_hist()` (weighted mode) and
`hvtiPlotR::hv_balance()`.

## Usage

``` r
ps_weight(
  data,
  treatment_col = "tavr",
  score_col = "prob_t",
  estimand = c("ATE", "ATT", "ATC"),
  stabilise = TRUE,
  trim = NULL,
  covariates = NULL,
  weight_col = "iptw"
)
```

## Arguments

- data:

  A data frame. Must contain `treatment_col` and `score_col`.

- treatment_col:

  Name of the binary treatment column (0/1 or logical). Default
  `"tavr"`.

- score_col:

  Name of the numeric propensity score column (values in \[0, 1\]).
  Default `"prob_t"`.

- estimand:

  Causal estimand. One of:

  - `"ATE"` (default) — average treatment effect. Treated weights =
    `1/ps`; control weights = `1/(1-ps)`.

  - `"ATT"` — average treatment effect on the treated. Treated weights =
    `1`; control weights = `ps/(1-ps)`.

  - `"ATC"` — average treatment effect on the controls. Treated weights
    = `(1-ps)/ps`; control weights = `1`.

- stabilise:

  Logical. If `TRUE` (default), weights are stabilised by multiplying by
  the marginal treatment probability.

- trim:

  Optional numeric in (0, 0.5). If supplied, weights outside the
  (`trim`, `1 - trim`) quantile range are winsorised. Default `NULL` (no
  trimming).

- covariates:

  Character vector of covariate columns for SMD diagnostics. If `NULL`
  (default), all numeric columns other than `treatment_col` and
  `score_col` are used.

- weight_col:

  Name of the output weight column appended to `$data`. Default
  `"iptw"`. If the column already exists it is overwritten.

## Value

An object of class `c("ps_weight", "ps_data")` with:

- `$data`:

  The input data frame with `weight_col` appended.

- `$meta`:

  Named list: `treatment_col`, `score_col`, `weight_col`, `estimand`,
  `stabilised`, `trim`, `method`, `n_total`.

- `$tables`:

  Named list: `smd_unweighted`, `smd_weighted`, `group_counts`,
  `effective_n`.

## See also

[`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md),
[`sample_ps_data()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sample_ps_data.md),
`hvtiPlotR::hv_mirror_hist()`, `hvtiPlotR::hv_balance()`

## Examples

``` r
dta <- sample_ps_data(n = 200, seed = 42)

# 1. ATE weights (default)
obj <- ps_weight(dta)
print(obj)
#> <ps_weight>
#>   N total     : 400
#>   Treatment   : tavr
#>   PS column   : prob_t
#>   Method      : IPTW-ATE
#>   Tables      : smd_unweighted, smd_weighted, group_counts, effective_n 
summary(obj)
#> Summary of <ps_weight>
#> 
#> Smd unweighted:
#>                  variable     smd
#> id                     id  3.4555
#> age                   age -0.7994
#> female             female  0.3989
#> ef                     ef  0.6139
#> diabetes         diabetes -0.0554
#> hypertension hypertension -0.0310
#> match               match      NA
#> 
#> Smd weighted:
#>                  variable     smd
#> id                     id  3.4514
#> age                   age  0.1130
#> female             female -0.0943
#> ef                     ef -0.0295
#> diabetes         diabetes -0.0278
#> hypertension hypertension -0.0973
#> match               match      NA
#> 
#> Group counts:
#>     group   n
#> 1 control 200
#> 2 treated 200
#> 
#> Effective n:
#>     group n_effective
#> 1 control       127.9
#> 2 treated       107.6
#> 

# 2. ATT weights, stabilised
obj_att <- ps_weight(dta, estimand = "ATT", stabilise = TRUE)
summary(obj_att)
#> Summary of <ps_weight>
#> 
#> Smd unweighted:
#>                  variable     smd
#> id                     id  3.4555
#> age                   age -0.7994
#> female             female  0.3989
#> ef                     ef  0.6139
#> diabetes         diabetes -0.0554
#> hypertension hypertension -0.0310
#> match               match      NA
#> 
#> Smd weighted:
#>                  variable     smd
#> id                     id  4.7721
#> age                   age  0.1031
#> female             female -0.1383
#> ef                     ef  0.0429
#> diabetes         diabetes -0.1745
#> hypertension hypertension -0.3172
#> match               match      NA
#> 
#> Group counts:
#>     group   n
#> 1 control 200
#> 2 treated 200
#> 
#> Effective n:
#>     group n_effective
#> 1 control        62.5
#> 2 treated       200.0
#> 

# 3. Extract data with weights appended
head(obj$data[, c("id", "tavr", "prob_t", "iptw")])
#>   id tavr    prob_t      iptw
#> 1  1    0 0.1112484 0.5625869
#> 2  2    0 0.5270896 1.0572826
#> 3  3    0 0.3951929 0.8267098
#> 4  4    0 0.2836606 0.6979932
#> 5  5    0 0.3751620 0.8002074
#> 6  6    0 0.3148051 0.7297193
```
