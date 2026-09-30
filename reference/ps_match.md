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
  seed = NULL
)
```

## Arguments

- data:

  A data frame. Must contain `treatment_col`, `score_col`, and `id_col`.

- treatment_col:

  Name of the binary treatment column (0/1 or logical). Default
  `"tavr"`.

- score_col:

  Name of the numeric propensity score column (values in \[0, 1\]).
  Default `"prob_t"`.

- id_col:

  Name of the patient identifier column. Default `"id"`.

- caliper:

  Maximum allowable difference in propensity score for a match to be
  accepted. `NULL` (default) imposes no caliper. Typical values:
  `0.1 * sd(score)` (20% of SD rule) or `0.2`.

- covariates:

  Character vector of covariate columns for SMD diagnostics. If `NULL`
  (default), all numeric columns other than `treatment_col` and
  `score_col` are used.

- match_col:

  Name of the output match-indicator column appended to `$data`. Default
  `"match"`. If the column already exists it is overwritten.

- ratio:

  Matching ratio (controls : treated). Currently only `1` (1:1) is
  supported. Default `1`.

- seed:

  Integer seed for the random shuffling of treated- patient order that
  reduces order-dependence. Pass an integer for a reproducible matched
  set; `NULL` (default) uses the session's current RNG state (not
  reproducible across sessions).

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
  `pair_id_col`, `caliper`, `method`, `n_total`, `n_matched`,
  `n_unmatched`.

- `$tables`:

  Named list: `smd_before`, `smd_after`, `group_counts_before`,
  `group_counts_after`.

## Details

The `$tables` slot contains standardised mean difference (SMD) tables
before and after matching, and group counts. These are the diagnostics
that `hvtiPlotR::hv_mirror_hist()` and `hvtiPlotR::hv_balance()` expect.

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
#>     group   n
#> 1 control 200
#> 2 treated 200
#> 
#> Group counts after:
#>     group   n
#> 1 control 200
#> 2 treated 200
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
#> 1     185
#> 2      44
#> 3     100
#> 4     127
#> 5       2
#> 6     118
table(obj$data$match)
#> 
#>   1 
#> 400 

# 4. Matched subset only
matched <- obj$data[obj$data$match == 1L, ]
nrow(matched)
#> [1] 400
```
