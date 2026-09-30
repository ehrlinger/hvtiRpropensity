# Propensity Score Overlap and Positivity Diagnostics

Examines the propensity score (PS) distribution for potential positivity
violations – regions of covariate space where the probability of
receiving one treatment is near zero – and quantifies the common support
region.

## Usage

``` r
sa_overlap(
  x,
  score_col = "prob_t",
  treatment_col = "tavr",
  trim_threshold = 0.05,
  treated_level = NULL
)
```

## Arguments

- x:

  A `ps_data` object **or** a plain data frame.

- score_col:

  Name of the propensity score column. Taken from `x$meta$score_col`
  when `x` is a `ps_data` object. Default `"prob_t"`.

- treatment_col:

  Name of the binary treatment column (0/1). Taken from
  `x$meta$treatment_col` when `x` is a `ps_data`
  object`. Default `"tavr"\`.

- trim_threshold:

  Numeric in (0, 0.5). PS values below this threshold or above
  `1 - trim_threshold` are flagged as near-positivity violations.
  Default `0.05`.

- treated_level:

  The value of `treatment_col` that is treated. Read from the metadata
  of a `ps_data` object when it records one. `NULL` (default) requires a
  0/1 or logical column and takes 1 or `TRUE`.

## Value

A list with:

- `$summary`:

  Data frame: one row per treatment group with `group` (the value in
  `treatment_col`), a logical `treated`, `n`, `ps_min`, `ps_q25`,
  `ps_median`, `ps_mean`, `ps_q75`, `ps_max`.

- `$overlap_region`:

  Named numeric vector: `lower` and `upper` bounds of the empirical
  common support region.

- `$outside_overlap`:

  Data frame: `group`, `treated`, `n_outside`, `pct_outside` – patients
  whose PS is outside the overlap region.

- `$positivity_flags`:

  Data frame: `group`, `treated`, `n_near_zero` (PS \<
  `trim_threshold`), `n_near_one` (PS \> `1 - trim_threshold`),
  `pct_near_zero`, `pct_near_one`.

## Details

The function reports:

- PS summary statistics by treatment group.

- The empirical overlap region: \[max(min_treated, min_control),
  min(max_treated, max_control)\].

- The proportion of each group whose PS lies outside the overlap region.

- The number and proportion of patients with PS below `trim_threshold`
  or above `1 - trim_threshold` (near-positivity-violation flags).

## See also

[`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md),
[`ps_weight()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_weight.md),
[`sa_rosenbaum()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_rosenbaum.md),
[`sa_evalue()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_evalue.md),
[`sa_trim_sweep()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_trim_sweep.md)

## Examples

``` r
dta <- sample_ps_data(n = 400, seed = 7)
res <- sa_overlap(dta, score_col = "prob_t", treatment_col = "tavr")
res$overlap_region
#>  lower  upper 
#> 0.0856 0.9082 
res$positivity_flags
#>   group treated n_near_zero n_near_one pct_near_zero pct_near_one
#> 1     0   FALSE           1          0           0.2          0.0
#> 2     1    TRUE           0          2           0.0          0.5

# From a ps_match object
m   <- ps_match(dta)
res <- sa_overlap(m)
res$outside_overlap
#>   group treated n_outside pct_outside
#> 1     0   FALSE        11         2.8
#> 2     1    TRUE        19         4.8
```
