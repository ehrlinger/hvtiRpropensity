# Random-forest propensity score (out-of-bag)

Estimates P(treated) with a random-forest classifier
([`randomForestSRC::rfsrc()`](https://www.randomforestsrc.org//reference/rfsrc.html))
and returns the **out-of-bag** class probability as the score, so no
patient is scored by trees that saw them. The result appends the same
score columns as
[`ps_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_logistic.md),
so it feeds
[`ps_weight()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_weight.md),
[`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md)
and
[`sa_overlap()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_overlap.md)
unchanged, and its `$data` feeds
[`hvtiPlotR::hv_mirror_hist()`](https://ehrlinger.github.io/hvtiPlotR/reference/hv_mirror_hist.html).

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

  Two-sided formula, `treatment ~ x1 + x2 + ...`. A `.` on the
  right-hand side expands to the data's other columns, excluding the
  treatment and `id_col`; terms the formula subtracts, as in `. - x3`,
  are left out of the forest.

- data:

  Data frame with no missing values in the formula variables.

- treatment_col:

  Treatment column name; defaults to the formula response. Binary, 0/1
  or logical.

- id_col:

  Optional patient identifier column, carried through.

- score_col, logit_col, weight_col:

  Names of the appended score, logit and matching-weight columns.

- treated_level:

  Value of `treatment_col` that is "treated" (the score is P(treated)).
  Defaults to `1` / `TRUE`.

- ntree:

  Number of trees.

- seed:

  `NULL` (default) or one whole number. With a seed, R's random number
  stream is set to `abs(seed)` immediately before the fit and
  [`randomForestSRC::rfsrc()`](https://www.randomforestsrc.org//reference/rfsrc.html)
  is given `seed = -abs(seed)`, the sign it requires; the caller's
  stream is restored afterwards and the seed is kept in `$meta$seed`.
  With `NULL`, `rfsrc()` takes its seed from the caller's stream, so the
  scores depend on whatever drew random numbers earlier in the session.

- clip:

  Scores are clipped to `[clip, 1 - clip]` before the logit and weight
  are computed. Forest probabilities can be exactly 0 or 1.

- ...:

  Further arguments for
  [`randomForestSRC::rfsrc()`](https://www.randomforestsrc.org//reference/rfsrc.html).

## Value

A `ps_forest` / `ps_data` object. `$data` is `data` plus `score_col`,
`logit_col`, `weight_col`, `quintile` and `decile`. `weight_col` holds
the Li-Greene matching weight computed as in
[`ps_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_logistic.md),
`min(p, 1 - p) / p` for treated patients and `min(p, 1 - p) / (1 - p)`
for controls, which targets the ATM estimand (average treatment effect
among the matched; Li & Greene, 2013,
[doi:10.1515/ijb-2012-0030](https://doi.org/10.1515/ijb-2012-0030) );
`$meta` holds `formula`, `treatment_col`, `id_col`, `score_col`,
`logit_col`, `weight_col`, `treated_level`, `treatment_levels`, `method`
(`"forest-oob"`), `ntree`, `seed` (the seed used, `NA` when none was
given) and `n_total`; `$tables` holds `smd` and `group_counts` (`group`,
the value in `treatment_col`; a logical `treated`; `n`);
`$models$forest` is the fitted forest.

## Threads and reproducibility

A seed makes the forest reproducible: seeded this way, fits have matched
to the last digit on a 192-thread OpenMP Linux server and a
single-threaded Mac. `randomForestSRC` uses every core unless
`options(rf.cores = )` says otherwise, and some `randomForestSRC`
procedures have not been stable from run to run under OpenMP. When two
runs with the same seed disagree, set `options(rf.cores = 1L)` before
calling `ps_forest()`; it costs speed and buys a single-threaded,
repeatable fit. This function leaves the option alone.

## Examples

``` r
# \donttest{
if (requireNamespace("randomForestSRC", quietly = TRUE)) {
  op <- options(rf.cores = 2L)
  dta <- sample_ps_data(n = 150, seed = 1)[, c("id", "tavr", "age", "ef")]
  obj <- ps_forest(tavr ~ age + ef, dta, ntree = 100, seed = 1)
  summary(obj$data$prob_t)
  options(op)
}
# }
```
