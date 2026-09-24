# Random-forest propensity score (out-of-bag)

Estimates P(treated) with a random-forest classifier
([`randomForestSRC::rfsrc()`](https://www.randomforestsrc.org//reference/rfsrc.html))
and returns the **out-of-bag** class probability as the score, so no
patient is scored by trees that saw them. The result has the same
contract as
[`ps_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_logistic.md),
so it feeds
[`ps_weight()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_weight.md),
[`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md),
[`sa_overlap()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_overlap.md)
and `hvtiPlotR::hv_balance()` / `hvtiPlotR::hv_mirror_hist()` unchanged.

## Usage

``` r
ps_forest(
  formula,
  data,
  treatment_col = NULL,
  id_col = "id",
  score_col = "prob_t",
  logit_col = "logit_t",
  weight_col = "mt_wt",
  treated_level = NULL,
  ntree = 500L,
  seed = NULL,
  clip = 0.001,
  ...
)
```

## Arguments

- formula:

  Two-sided formula, `treatment ~ x1 + x2 + ...`.

- data:

  Data frame with no missing values in the formula variables.

- treatment_col:

  Treatment column name; defaults to the formula response. Binary, 0/1
  or logical.

- id_col:

  Optional patient identifier column, carried through.

- score_col, logit_col, weight_col:

  Names of the appended score, logit and overlap-weight columns.

- treated_level:

  Value of `treatment_col` that is "treated" (the score is P(treated)).
  Defaults to `1` / `TRUE`.

- ntree:

  Number of trees.

- seed:

  Optional integer seed passed to
  [`randomForestSRC::rfsrc()`](https://www.randomforestsrc.org//reference/rfsrc.html).

- clip:

  Scores are clipped to `[clip, 1 - clip]` before the logit and weight
  are computed. Forest probabilities can be exactly 0 or 1.

- ...:

  Further arguments for
  [`randomForestSRC::rfsrc()`](https://www.randomforestsrc.org//reference/rfsrc.html).

## Value

A `ps_forest` / `ps_data` object. `$data` is `data` plus `score_col`,
`logit_col`, `weight_col` (overlap weights), `quintile` and `decile`;
`$meta` follows
[`ps_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_logistic.md)
with `method = "forest-oob"`; `$tables` holds `smd` and `group_counts`;
`$models$forest` is the fitted forest.

## Examples

``` r
# \donttest{
if (requireNamespace("randomForestSRC", quietly = TRUE)) {
  dta <- sample_ps_data(n = 150, seed = 1)[, c("id", "tavr", "age", "ef")]
  obj <- ps_forest(tavr ~ age + ef, dta, ntree = 100, seed = 1)
  summary(obj$data$prob_t)
}
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.0010  0.2548  0.4615  0.5019  0.7744  0.9990 
# }
```
