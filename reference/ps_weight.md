# Inverse-Probability-of-Treatment Weighting (IPTW)

Computes IPTW weights from a pre-computed propensity score column and
appends them to the dataset. `estimand` picks the target: the average
treatment effect in the whole sample (ATE, the default), in the treated
(ATT) or in the controls (ATC). For the matched population (ATM), use
the `mt_wt` column from
[`ps_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_logistic.md)
or
[`ps_forest()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_forest.md)
instead. Returns a `ps_weight` object whose `$tables` slot contains SMD
tables (unweighted and weighted) and effective-N diagnostics. Its
`$data` feeds
[`hvtiPlotR::hv_mirror_hist()`](https://ehrlinger.github.io/hvtiPlotR/reference/hv_mirror_hist.html)
in weighted mode, and
[`hvtiPlotR::hv_balance()`](https://ehrlinger.github.io/hvtiPlotR/reference/hv_balance.html)
reads the object itself (hvtiPlotR 2.8.22 or later).

## Usage

``` r
ps_weight(
  data,
  treatment_col,
  score_col = "prob_t",
  estimand = c("ATE", "ATT", "ATC"),
  stabilise = TRUE,
  trim = NULL,
  covariate_cols = NULL,
  weight_col = "iptw",
  treated_level = NULL,
  id_col = "id",
  covariates
)
```

## Arguments

- data:

  A data frame. Must contain `treatment_col` and `score_col`.

- treatment_col:

  Name of the binary treatment column (0/1 or logical). Required: it has
  no default, because a study's column name never is one. Omitting it is
  deprecated, warns once per session and still means `"tavr"`, for one
  release.

- score_col:

  Name of the numeric propensity score column (values in \[0, 1\]).
  Default `"prob_t"`.

- estimand:

  Causal estimand. One of:

  - `"ATE"` (default): average treatment effect. Treated weights =
    `1/ps`; control weights = `1/(1-ps)`.

  - `"ATT"`: average treatment effect on the treated. Treated weights =
    `1`; control weights = `ps/(1-ps)`.

  - `"ATC"`: average treatment effect on the controls. Treated weights =
    `(1-ps)/ps`; control weights = `1`.

- stabilise:

  Logical. If `TRUE` (default) and `estimand = "ATE"`, each weight is
  multiplied by the marginal probability of the patient's own arm
  (`P(treated)` or `P(control)`). Stabilisation is defined for the ATE
  only: for `"ATT"` and `"ATC"` it has no effect, so the reference arm
  keeps weights of exactly 1. `$meta$stabilised` records whether it was
  applied.

- trim:

  Optional numeric in (0, 0.5). If supplied, weights outside the
  (`trim`, `1 - trim`) quantile range are winsorised, within the arm
  that is reweighted. For `"ATT"` only the control weights are trimmed,
  at quantiles of the control weights, and every treated weight stays
  exactly 1; `"ATC"` is the mirror image. For `"ATE"` both arms are
  reweighted, so the quantiles are taken over all weights pooled across
  both arms (Lee, Lessler and Stuart, 2011). Default `NULL` (no
  trimming).

- covariate_cols:

  Character vector of covariate columns for SMD diagnostics, used as
  given. A 0/1 or logical column takes the binary standardized
  difference, a factor or character column the categorical one of
  [`ps_stddiff()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_stddiff.md)
  (one row per column), any other numeric column the Gaussian one. If
  `NULL` (default), every numeric, logical, factor or character column
  that varies is used, except identifier-like text columns (every value
  different), `treatment_col`, `score_col`, `weight_col`, `id_col` and
  the fixed columns this package appends (`prob_t`, `logit_t`, `mt_wt`,
  `iptw`, `match`, `pair_id`, `quintile`, `decile`, `bs`, `cluster`).

- weight_col:

  Name of the output weight column appended to `$data`. Default
  `"iptw"`. If the column already exists it is overwritten.

- treated_level:

  The value of `treatment_col` that is treated. `NULL` (default)
  requires a 0/1 or logical column and takes 1 or `TRUE`. Give a value,
  for example `"transcatheter"`, to use a column holding any two values.

- id_col:

  The patient identifier column, left out of the default covariates.
  Default `"id"`; it need not exist in `data`.

- covariates:

  Deprecated: the old name of `covariate_cols`. It still works, and
  warns once per session, for one release.

## Value

An object of class `c("ps_weight", "ps_data")` with:

- `$data`:

  The input data frame with `weight_col` appended.

- `$meta`:

  Named list: `treatment_col`, `treated_level` (the value of
  `treatment_col` taken as treated), `score_col`, `weight_col`,
  `estimand`, `stabilised`, `trim`, `method`, `n_total`, `n_excluded`
  (rows with a missing treatment or score).

- `$tables`:

  Named list: `smd_unweighted`, `smd_weighted` (`variable`, `smd`),
  `group_counts` (`group`, `treated`, `n`) and `effective_n` (`group`,
  `treated`, `n_effective`, the Kish effective sample size). Group
  tables name each group by its value in `treatment_col` (`0`/`1`,
  `FALSE`/`TRUE`, or the two values when `treated_level` is given) and
  carry a logical `treated` column.

## Details

A patient with a missing treatment or a missing score cannot be
weighted. Those rows stay in `$data` with a missing `weight_col`, are
left out of the stabilising treatment proportion, the SMD tables, the
group counts and the effective sample sizes, and one warning states how
many there are. `$meta$n_excluded` keeps the count.

## References

Lee BK, Lessler J, Stuart EA (2011). Weight trimming and propensity
score weighting. *PLoS ONE*, 6(3), e18174.
[doi:10.1371/journal.pone.0018174](https://doi.org/10.1371/journal.pone.0018174)

## See also

[`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md),
[`sample_ps_data()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sample_ps_data.md),
[`hvtiPlotR::hv_mirror_hist()`](https://ehrlinger.github.io/hvtiPlotR/reference/hv_mirror_hist.html),
[`hvtiPlotR::hv_balance()`](https://ehrlinger.github.io/hvtiPlotR/reference/hv_balance.html)

## Examples

``` r
dta <- sample_ps_data(n = 200, seed = 42)

# 1. ATE weights (default)
obj <- ps_weight(dta, treatment_col = "tavr")
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
#>       variable     smd
#> 1          age -0.7994
#> 2       female  0.3999
#> 3           ef  0.6139
#> 4     diabetes -0.0556
#> 5 hypertension -0.0310
#> 
#> Smd weighted:
#>       variable     smd
#> 1          age  0.1143
#> 2       female -0.0956
#> 3           ef -0.0298
#> 4     diabetes -0.0281
#> 5 hypertension -0.0987
#> 
#> Group counts:
#>   group treated   n
#> 1     0   FALSE 200
#> 2     1    TRUE 200
#> 
#> Effective n:
#>   group treated n_effective
#> 1     0   FALSE       127.9
#> 2     1    TRUE       107.6
#> 

# 2. ATT weights: treated patients keep a weight of 1
obj_att <- ps_weight(dta, treatment_col = "tavr", estimand = "ATT")
summary(obj_att)
#> Summary of <ps_weight>
#> 
#> Smd unweighted:
#>       variable     smd
#> 1          age -0.7994
#> 2       female  0.3999
#> 3           ef  0.6139
#> 4     diabetes -0.0556
#> 5 hypertension -0.0310
#> 
#> Smd weighted:
#>       variable     smd
#> 1          age  0.0734
#> 2       female -0.0987
#> 3           ef  0.0305
#> 4     diabetes -0.1246
#> 5 hypertension -0.2262
#> 
#> Group counts:
#>   group treated   n
#> 1     0   FALSE 200
#> 2     1    TRUE 200
#> 
#> Effective n:
#>   group treated n_effective
#> 1     0   FALSE        62.5
#> 2     1    TRUE       200.0
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
