# hvtiRpropensity: Propensity Score Methods for Cardiac Surgery Studies

Provides functions for propensity score and balancing score analysis
used in cardiac surgery comparative-effectiveness research. Translates
the Cleveland Clinic HVTI SAS template library into a consistent R
workflow.

**Step 1: estimate scores**:

- [`ps_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_logistic.md):
  binary treatment propensity score (logistic regression; supports
  multiply-imputed data via a stacked data frame).

- [`ps_forest()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_forest.md):
  binary treatment propensity score from a random forest
  (`randomForestSRC`), scored out-of-bag.

- [`ps_ordinal()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_ordinal.md):
  ordered treatment propensity score (proportional-odds model via
  [`MASS::polr`](https://rdrr.io/pkg/MASS/man/polr.html); equivalent to
  SAS `PROC LOGISTIC` default for ordinal responses).

- [`ps_nominal()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_nominal.md):
  nominal treatment propensity score (multinomial logistic via
  [`nnet::multinom`](https://rdrr.io/pkg/nnet/man/multinom.html);
  equivalent to SAS `PROC LOGISTIC link=glogit`).

- [`bs_continuous()`](https://ehrlinger.github.io/hvtiRpropensity/reference/bs_continuous.md):
  balancing score for a continuous exposure (linear regression;
  equivalent to SAS `PROC REG` + `%grp` macro).

- [`bs_count()`](https://ehrlinger.github.io/hvtiRpropensity/reference/bs_count.md):
  balancing score for a count exposure (negative-binomial or Poisson
  regression; equivalent to SAS `PROC GENMOD dist=nb`).

**Outcome models**:

- [`fit_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/fit_logistic.md):
  binary, ordinal or nominal outcome-model bundle, with Rubin pooling
  over stacked imputations and Barnard-Rubin degrees of freedom,
  matching `mice::pool()`.

- [`fit_linear()`](https://ehrlinger.github.io/hvtiRpropensity/reference/fit_linear.md):
  linear outcome-model bundle, with Rubin pooling.

- [`validate_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/validate_logistic.md):
  applies a saved binary bundle to a validation cohort without
  refitting.

**Step 2: balance and use scores**:

- [`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md):
  1:1 nearest-neighbour matching on a pre-computed score, on the
  probability or the logit scale. `scale = "logit", caliper_sd = 0.2` is
  the Austin (2011) caliper. Returns a `pair_id` column in `$data` for
  use with
  [`sa_rosenbaum()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_rosenbaum.md)
  and
  [`ps_rmst()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_rmst.md).

- [`ps_weight()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_weight.md):
  IPTW weights (ATE / ATT / ATC) from a pre-computed score.

- [`ps_support()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_support.md):
  compares common-support and trimming rules.

- [`ps_rmst()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_rmst.md):
  restricted mean survival contrast, weighted, or on the matched pairs
  of a
  [`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md)
  object.

**Balance statistics**:

- [`ps_stddiff()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_stddiff.md):
  standardized differences by variable type, optionally weighted (CCF
  `%stddiff`).

- [`ps_stddiff_perm()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_stddiff_perm.md):
  permutation reference for those differences (CCF `%stddiffci`).

- [`ps_mw_var()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_mw_var.md):
  matching-weight difference in means with a bootstrap SD (CCF
  `%mw_var`).

**Step 3: sensitivity analysis** (optional):

- [`sa_rosenbaum()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_rosenbaum.md):
  Rosenbaum gamma bounds for matched analyses; how much hidden bias
  would overturn the result.

- [`sa_evalue()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_evalue.md):
  E-values (VanderWeele & Ding 2017) for unmeasured confounding. Odds
  and hazard ratios need `rare`.

- [`sa_overlap()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_overlap.md):
  overlap and positivity diagnostics on the score distribution.

- [`sa_trim_sweep()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_trim_sweep.md):
  how the effective sample size and the weights change across IPTW trim
  thresholds.

**Data**:
[`sample_ps_data()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sample_ps_data.md),
[`sample_ps_data_ordinal()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sample_ps_data_ordinal.md),
[`sample_ps_data_nominal()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sample_ps_data_nominal.md)
and
[`sample_ps_data_count()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sample_ps_data_count.md)
generate synthetic datasets;
[`is_ps_data()`](https://ehrlinger.github.io/hvtiRpropensity/reference/is_ps_data.md)
tests for a `ps_data` object.

Every function except the four data generators (which return data
frames) and
[`is_ps_data()`](https://ehrlinger.github.io/hvtiRpropensity/reference/is_ps_data.md)
returns a `ps_data` S3 object. Its `$data` slot holds the input data
with any score, weight or flag columns appended, `$meta` records the
columns and parameters used, and `$tables` holds the diagnostics: SMD
tables, group counts, effective sample sizes, or a sensitivity
function's bounds and sweeps.

[`hvtiPlotR::hv_mirror_hist()`](https://ehrlinger.github.io/hvtiPlotR/reference/hv_mirror_hist.html)
reads the `$data` of a
[`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md)
or
[`ps_weight()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_weight.md)
object directly.
[`hvtiPlotR::hv_balance()`](https://ehrlinger.github.io/hvtiPlotR/reference/hv_balance.html)
does not read `$tables` as they are: it needs one long data frame with
`variable`, `group` and `std_diff` columns, the difference in percent,
so the SMD tables are stacked and multiplied by 100 first. The vignette
shows the reshape.

The SAS programs this package ports ship with it, in
`system.file("sas-templates", package = "hvtiRpropensity")`.

## Typical two-step workflow (binary treatment)

    library(hvtiRpropensity)

    # Step 1: estimate propensity score
    dta   <- sample_ps_data(n = 500, seed = 42)
    dta$prob_t <- NULL
    score <- ps_logistic(
      tavr ~ age + female + ef + diabetes + hypertension,
      data = dta
    )
    print(score)

    # Step 2a: match on the logit of the score, caliper 0.2 SD
    m <- ps_match(score$data, score_col = "logit_t", scale = "logit",
                  caliper_sd = 0.2, seed = 42)
    matched <- m$data[m$data$match == 1L, ]

    # Step 2b: or weight
    w <- ps_weight(score$data, score_col = score$meta$score_col)

## References

Austin PC (2011). Optimal caliper widths for propensity-score matching
when estimating differences in means and differences in proportions in
observational studies. *Pharmaceutical Statistics*, 10(2), 150–161.
[doi:10.1002/pst.433](https://doi.org/10.1002/pst.433)

## See also

- [`ps_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_logistic.md),
  [`ps_forest()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_forest.md),
  [`ps_ordinal()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_ordinal.md),
  [`ps_nominal()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_nominal.md):
  propensity score estimation

- [`bs_continuous()`](https://ehrlinger.github.io/hvtiRpropensity/reference/bs_continuous.md),
  [`bs_count()`](https://ehrlinger.github.io/hvtiRpropensity/reference/bs_count.md):
  balancing scores for continuous and count exposures

- [`fit_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/fit_logistic.md),
  [`fit_linear()`](https://ehrlinger.github.io/hvtiRpropensity/reference/fit_linear.md),
  [`validate_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/validate_logistic.md):
  outcome models

- [`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md),
  [`ps_weight()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_weight.md):
  matching and IPTW weighting

- [`ps_stddiff()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_stddiff.md),
  [`ps_stddiff_perm()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_stddiff_perm.md),
  [`ps_mw_var()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_mw_var.md):
  balance statistics

- [`ps_support()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_support.md),
  [`ps_rmst()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_rmst.md):
  support rules and survival contrasts

- [`sa_rosenbaum()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_rosenbaum.md),
  [`sa_evalue()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_evalue.md),
  [`sa_overlap()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_overlap.md),
  [`sa_trim_sweep()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_trim_sweep.md):
  sensitivity analysis

- [`sample_ps_data()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sample_ps_data.md):
  synthetic dataset for examples / tests

- [`is_ps_data()`](https://ehrlinger.github.io/hvtiRpropensity/reference/is_ps_data.md):
  predicate for `ps_data` objects

## Author

**Maintainer**: John Ehrlinger <john.ehrlinger@gmail.com>

Authors:

- John Ehrlinger <john.ehrlinger@gmail.com>
