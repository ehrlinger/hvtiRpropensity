# Propensity Score Nearest-Neighbour Matching

Performs 1:1 nearest-neighbour matching without replacement on a
pre-computed propensity score column. Returns a `ps_match` object whose
`$data` slot is the **full** original dataset with a `match` column set
to `1` for matched pairs and `0` for unmatched patients.

## Usage

``` r
ps_match(
  data,
  treatment_col = "tavr",
  score_col = "prob_t",
  id_col = "id",
  caliper = NULL,
  covariates = NULL,
  match_col = "match",
  ratio = 1L,
  seed = NULL,
  treated_level = NULL,
  scale = c("probability", "logit"),
  caliper_sd = NULL
)
```

## Arguments

- data:

  A data frame. Must contain `treatment_col`, `score_col`, and `id_col`.

- treatment_col:

  Name of the binary treatment column (0/1 or logical). Default
  `"tavr"`.

- score_col:

  Name of the numeric propensity score column, on the scale named by
  `scale`: probabilities in \[0, 1\] for `scale = "probability"`, or
  logits (any finite value, such as the `logit_t` column
  [`ps_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_logistic.md)
  appends) for `scale = "logit"`. Default `"prob_t"`.

- id_col:

  Name of the patient identifier column. Default `"id"`.

- caliper:

  Maximum allowable difference in `score_col`, on its own scale, for a
  match to be accepted. `NULL` (default) imposes no caliper unless
  `caliper_sd` is given. Give at most one of `caliper` and `caliper_sd`.

- covariates:

  Character vector of covariate columns for SMD diagnostics, used as
  given. If `NULL` (default), every numeric column that varies is used,
  except `treatment_col`, `score_col`, `id_col`, `match_col` and the
  fixed columns this package appends (`prob_t`, `logit_t`, `mt_wt`,
  `iptw`, `pair_id`, `quintile`, `decile`, `bs`, `cluster`).

- match_col:

  Name of the output match-indicator column appended to `$data`. Default
  `"match"`. If the column already exists it is overwritten.

- ratio:

  Matching ratio (controls : treated). Currently only the number `1`
  (1:1) is supported. Default `1`.

- seed:

  Seed for the random shuffling of treated-patient order that reduces
  order-dependence. `NULL` (default) or one whole number. With a seed,
  the random number stream is set to `abs(seed)` immediately before the
  draws, the result is reproducible whatever ran earlier in the session,
  the caller's stream is restored afterwards, and the seed is kept in
  `$meta$seed`. With `NULL` the draws come from the caller's stream and
  advance it, and `$meta$seed` is `NA`.

- treated_level:

  The value of `treatment_col` that is treated. `NULL` (default)
  requires a 0/1 or logical column and takes 1 or `TRUE`. Give a value,
  for example `"transcatheter"`, to use a column holding any two values.

- scale:

  The scale `score_col` is on, and so the scale distances and the
  caliper are measured on: `"probability"` (default) or `"logit"`. The
  recommended caliper, 0.2 standard deviations of the logit of the
  propensity score (Austin 2011), needs `scale = "logit"`.

- caliper_sd:

  Caliper in pooled standard-deviation units of `score_col`: the caliper
  applied is `caliper_sd * sqrt((var_t + var_c) / 2)`, where `var_t` and
  `var_c` are the variances of `score_col` among treated and control
  patients (missing values dropped). `NULL` (default) uses `caliper`
  instead. `scale = "logit", caliper_sd = 0.2` is the Austin (2011)
  recommendation.

## Value

An object of class `c("ps_match", "ps_data")` with:

- `$data`:

  The full input data frame with `match_col` and `pair_id` appended.
  `match_col` is `1` for matched patients and `0` otherwise. `pair_id`
  is an integer identifying which matched pair each patient belongs to
  (`NA` for unmatched patients). Filter to `match == 1` to obtain the
  matched subset.

- `$meta`:

  Named list: `treatment_col`, `score_col`, `id_col`, `match_col`,
  `pair_id_col`, `scale`, `caliper` (the absolute caliper applied, on
  `scale`), `caliper_sd`, `method`, `seed`, `n_total`, `n_matched`,
  `n_unmatched`.

- `$tables`:

  Named list: `smd_before`, `smd_after`, `group_counts_before`,
  `group_counts_after`. Group tables name each group by its value in
  `treatment_col` (`0`/`1` or `FALSE`/`TRUE`) and carry a logical
  `treated` column.

## Details

The `$tables` slot contains standardised mean difference (SMD) tables
before and after matching, and group counts. These are the diagnostics
that `hvtiPlotR::hv_mirror_hist()` and `hvtiPlotR::hv_balance()` expect.

## References

Austin PC (2011). Optimal caliper widths for propensity-score matching
when estimating differences in means and differences in proportions in
observational studies. *Pharmaceutical Statistics*, 10(2), 150–161.
[doi:10.1002/pst.433](https://doi.org/10.1002/pst.433)

## See also

[`ps_weight()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_weight.md),
[`sample_ps_data()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sample_ps_data.md),
`hvtiPlotR::hv_mirror_hist()`, `hvtiPlotR::hv_balance()`

## Examples

``` r
dta <- sample_ps_data(n = 200, seed = 42)

# 1. Build the matched object
obj <- ps_match(dta)

# 2. Inspect diagnostics
print(obj)
#> <ps_match>
#>   N total     : 400
#>   Treatment   : tavr
#>   PS column   : prob_t
#>   Method      : nearest-neighbour 1:1
#>   Tables      : smd_before, smd_after, group_counts_before, group_counts_after 
summary(obj)
#> Summary of <ps_match>
#> 
#> Smd before:
#>                  variable     smd
#> age                   age -0.7994
#> female             female  0.3989
#> ef                     ef  0.6139
#> diabetes         diabetes -0.0554
#> hypertension hypertension -0.0310
#> 
#> Smd after:
#>                  variable     smd
#> age                   age -0.7994
#> female             female  0.3989
#> ef                     ef  0.6139
#> diabetes         diabetes -0.0554
#> hypertension hypertension -0.0310
#> 
#> Group counts before:
#>   group treated   n
#> 1     0   FALSE 200
#> 2     1    TRUE 200
#> 
#> Group counts after:
#>   group treated   n
#> 1     0   FALSE 200
#> 2     1    TRUE 200
#> 

# 3. Extract the full data (with match indicator)
head(obj$data)
#>   id tavr      age female       ef diabetes hypertension    prob_t match
#> 1  1    0 88.96767      0 51.95379        1            1 0.1112484     1
#> 2  2    0 73.48241      0 59.60242        0            1 0.5270896     1
#> 3  3    0 80.90503      1 52.38991        1            0 0.3951929     1
#> 4  4    0 83.06290      0 59.35072        0            1 0.2836606     1
#> 5  5    0 81.23415      1 50.53527        0            0 0.3751620     1
#> 6  6    0 77.15100      0 51.42113        0            0 0.3148051     1
#>   pair_id
#> 1     186
#> 2      45
#> 3     104
#> 4     129
#> 5       3
#> 6      27
table(obj$data$match)
#> 
#>   1 
#> 400 

# 4. Matched subset only
matched <- obj$data[obj$data$match == 1L, ]
nrow(matched)
#> [1] 400

# 5. Match on the logit with a 0.2 SD caliper (Austin 2011)
dta$logit_t <- stats::qlogis(dta$prob_t)
obj_logit <- ps_match(dta, score_col = "logit_t", scale = "logit", caliper_sd = 0.2)
obj_logit$meta$caliper
#> [1] 0.2170102
```
