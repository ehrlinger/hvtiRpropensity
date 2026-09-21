# hvtiRpropensity: Propensity Score Methods for Cardiac Surgery Studies

Provides functions for propensity score and balancing score analysis
used in cardiac surgery comparative-effectiveness research. Translates
the Cleveland Clinic HVTI SAS template library into a consistent R
workflow.

**Step 1 — estimate scores** (new in this version):

- [`ps_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_logistic.md)
  — binary treatment propensity score (logistic regression; supports
  multiply-imputed data via a stacked data frame).

- [`ps_ordinal()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_ordinal.md)
  — ordered treatment propensity score (proportional-odds model via
  [`MASS::polr`](https://rdrr.io/pkg/MASS/man/polr.html); equivalent to
  SAS `PROC LOGISTIC` default for ordinal responses).

- [`ps_nominal()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_nominal.md)
  — nominal treatment propensity score (multinomial logistic via
  [`nnet::multinom`](https://rdrr.io/pkg/nnet/man/multinom.html);
  equivalent to SAS `PROC LOGISTIC link=glogit`).

- [`bs_continuous()`](https://ehrlinger.github.io/hvtiRpropensity/reference/bs_continuous.md)
  — balancing score for a continuous exposure (linear regression;
  equivalent to SAS `PROC REG` + `%grp` macro).

- [`bs_count()`](https://ehrlinger.github.io/hvtiRpropensity/reference/bs_count.md)
  — balancing score for a count exposure (negative-binomial or Poisson
  regression; equivalent to SAS `PROC GENMOD dist=nb`).

**Step 2 — balance and use scores**:

- [`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md)
  — 1:1 nearest-neighbour matching on a pre-computed score. Returns a
  `pair_id` column in `$data` for use with
  [`sa_rosenbaum()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_rosenbaum.md).

- [`ps_weight()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_weight.md)
  — IPTW weights (ATE / ATT / ATC) from a pre-computed score.

**Step 3 — sensitivity analysis** (optional):

- [`sa_rosenbaum()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_rosenbaum.md)
  — Rosenbaum gamma bounds for matched analyses (requires `rbounds`).

- [`sa_evalue()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_evalue.md)
  — E-values (VanderWeele & Ding 2017); no extra dependency.

- [`sa_overlap()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_overlap.md)
  — Overlap and positivity diagnostics on the PS distribution.

- [`sa_trim_sweep()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_trim_sweep.md)
  — IPTW trim-threshold sensitivity sweep.

All functions return a `ps_data` S3 object whose `$data` slot holds the
original dataset with score / weight columns appended, and whose
`$tables` slot holds balance diagnostics compatible with
`hvtiPlotR::hv_mirror_hist()` and `hvtiPlotR::hv_balance()`.

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

    # Step 2a: match on the estimated score
    m <- ps_match(score$data, score_col = score$meta$score_col)
    matched <- m$data[m$data$match == 1L, ]

    # Step 2b: or weight
    w <- ps_weight(score$data, score_col = score$meta$score_col)

## See also

- [`ps_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_logistic.md)
  — binary propensity score estimation

- [`ps_ordinal()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_ordinal.md)
  — ordinal propensity score estimation

- [`ps_nominal()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_nominal.md)
  — nominal propensity score estimation

- [`bs_continuous()`](https://ehrlinger.github.io/hvtiRpropensity/reference/bs_continuous.md)
  — continuous-outcome balancing score

- [`bs_count()`](https://ehrlinger.github.io/hvtiRpropensity/reference/bs_count.md)
  — count-outcome balancing score

- [`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md)
  — nearest-neighbour 1:1 matching

- [`ps_weight()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_weight.md)
  — IPTW weighting (ATE / ATT / ATC)

- [`sa_rosenbaum()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_rosenbaum.md)
  — Rosenbaum sensitivity bounds (matched analyses)

- [`sa_evalue()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_evalue.md)
  — E-value sensitivity analysis

- [`sa_overlap()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_overlap.md)
  — overlap and positivity diagnostics

- [`sa_trim_sweep()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_trim_sweep.md)
  — IPTW trim-threshold sensitivity sweep

- [`sample_ps_data()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sample_ps_data.md)
  — synthetic dataset for examples / tests

- [`is_ps_data()`](https://ehrlinger.github.io/hvtiRpropensity/reference/is_ps_data.md)
  — predicate for `ps_data` objects

## Author

**Maintainer**: John Ehrlinger <john.ehrlinger@gmail.com>

Authors:

- John Ehrlinger <john.ehrlinger@gmail.com>
