# Fit a general logistic model bundle

Fits binary, proportional-odds ordinal, or generalized-logit nominal
models under an explicit outcome-level contract. Stacked long-form
imputations are fitted separately, one model per imputation. The
coefficients are combined with Rubin's rules; the predicted
probabilities are averaged per patient.

## Usage

``` r
fit_logistic(
  formula,
  data,
  family = c("binary", "ordinal", "nominal"),
  outcome_col = NULL,
  id_col = "id",
  imputation_col = NULL,
  outcome_levels,
  event_level = NULL,
  reference_level = NULL,
  prediction_prefix = "prob",
  trace = FALSE
)
```

## Arguments

- formula:

  Model formula.

- data:

  A data frame containing the outcome and predictors.

- family:

  Model family: `"binary"` (default,
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html) with a logit link),
  `"ordinal"` (proportional-odds cumulative logit,
  [`MASS::polr()`](https://rdrr.io/pkg/MASS/man/polr.html); a two-level
  outcome is fitted with
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html)) or `"nominal"`
  (generalized logit,
  [`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html)).

- outcome_col:

  Outcome column. By default it is read from `formula`.

- id_col:

  Patient identifier column. Default `"id"`. Required, and checked, only
  when `imputation_col` is set, where it aligns each patient's rows
  across imputations.

- imputation_col:

  Imputation-index column, or `NULL` for one dataset.

- outcome_levels:

  Complete outcome levels in their declared order.

- event_level:

  Event level for a binary model.

- reference_level:

  Reference level for a nominal model.

- prediction_prefix:

  Name for a binary probability column, or prefix for ordinal and
  nominal level-specific probability columns.

- trace:

  Whether model engines should print fitting progress.

## Value

An object of class `c("lm_fit", "ps_data")` with:

- `$data`:

  `data` (the lowest-numbered imputation when `imputation_col` is set,
  without that column) with the predicted probabilities appended: one
  column named `prediction_prefix` for a binary model (`prob` by
  default), the probability of `event_level`; or one
  `<prediction_prefix>_<level>` column per outcome level for an ordinal
  or nominal model, whose outcome column is returned as a factor
  (ordered for an ordinal model). A row missing a predictor has missing
  probabilities.

- `$meta`:

  Named list: `bundle_version` (1), `formula`, `outcome_col`,
  `model_family`, `outcome_levels`, `event_level` (binary only, else
  `NULL`), `reference_level` (nominal only, else `NULL`),
  `cumulative_direction` (ordinal only, else `NULL`), `predictors`,
  `id_col`, `imputation_col`, `n_imputations`, `n_input`, `n_analyzed`,
  `n_excluded` (summed over imputations), `prediction_prefix` and
  `package_versions` (R, this package and the model engine).

- `$tables`:

  Named list: `estimates` (one row per term: `term`, `estimate`,
  `std.error`, `statistic`, `df`, `p.value`, `conf.low`, `conf.high`,
  `odds_ratio` (`NA` for an ordinal threshold), and `pooled`, `TRUE`
  when Rubin's rules combined imputations; a pooled `df` is Barnard and
  Rubin's, see "Pooled degrees of freedom"), `covariance` (the
  coefficient covariance matrix, pooled when imputed), `by_imputation`
  (`imputation`, `term`, `estimate`, `std.error`) and `fit_status` (one
  row per imputation: `imputation`, `converged`, `n_input`,
  `n_analyzed`, `n_excluded`).

- `$models`:

  The fitted models, one per imputation, named by imputation index
  (`"1"` for a single dataset). For binary bundles,
  [`validate_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/validate_logistic.md)
  scores new data with them; it rejects ordinal and nominal bundles.

## Details

**Multiple imputation.** With `imputation_col` set, one model is fitted
to each imputed dataset. The two outputs are then combined differently.
`$tables$estimates` and `$tables$covariance` pool the coefficients and
their variances with Rubin's rules, with Barnard and Rubin's degrees of
freedom (see "Pooled degrees of freedom" below). The probability columns
in `$data` are the per-patient mean, on the probability scale, of the
predictions from each imputation's model. Rubin's rules play no part
there: averaging the predictions is what Leyrat et al. (2019) call the
"across" approach. For a propensity score they favour the "within"
approach instead; see
[`ps_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_logistic.md)
for the recipe.

`$data` is built from the lowest-numbered imputation: its rows, ordered
as they appear there, its covariate values, and no `imputation_col`.

## Pooled degrees of freedom

With `imputation_col` set, each pooled coefficient is tested against a t
distribution, and the `df` column of `$tables$estimates` gives its
degrees of freedom. We use Barnard and Rubin's (1999) small-sample
degrees of freedom, computed the way `mice::pool()` computes them, so a
pooled table here and one from `mice` agree.

Rubin's (1987) original degrees of freedom are \\\nu\_{old} = (m - 1) /
\lambda^2\\, where \\m\\ is the number of imputations and \\\lambda =
(1 + 1/m) B / T\\ is the proportion of a coefficient's total variance
\\T\\ due to missing data, through the between-imputation variance
\\B\\. (This is `lambda` in `mice`, not its finite-sample fraction of
missing information, `fmi`.) When little is missing, \\\lambda\\ is
small and \\\nu\_{old}\\ can run into the millions, larger than the
sample could ever support, and the intervals come out too narrow.
Barnard and Rubin cap it with the complete-data degrees of freedom
\\\nu\_{com}\\, the value the analysis would have had with no missing
data: \$\$\nu\_{obs} = \frac{\nu\_{com} + 1}{\nu\_{com} + 3} \nu\_{com}
(1 - \lambda), \qquad \nu = \left(\frac{1}{\nu\_{old}} +
\frac{1}{\nu\_{obs}}\right)^{-1}.\$\$ So \\\nu\\ never exceeds
\\\nu\_{com}\\, and it approaches \\\nu\_{old}\\ as \\\nu\_{com}\\
grows.

\\\nu\_{com}\\ is \\n - k\\, the residual degrees of freedom of the
first imputation's model: \\n\\ analysed patients minus \\k\\ estimated
parameters. That holds for every family, linear, logistic, ordinal,
nominal and count, as it does in `mice`.
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) reports
no residual degrees of freedom, so for a nominal model \\k\\ is its
effective degrees of freedom, again as `mice` does. Pooling stops with
an error when \\\nu\_{com}\\ is zero or negative, as it is for a
saturated model, because no small-sample degrees of freedom exist then.

Two consequences are worth knowing. First, a single logistic, ordinal,
nominal or count fit is tested against the normal distribution (`df` is
`Inf`), but the same model pooled over imputations reports a finite `df`
and a t test. A single linear fit uses its residual degrees of freedom.
Second, SAS `PROC MIANALYZE` reports \\\nu\_{old}\\ unless you give it
the complete-data degrees of freedom with `EDF=`. To reproduce these
`df`, and the p-values and intervals that depend on them, run
`PROC MIANALYZE` with `EDF=` set to \\n - k\\.

Barnard J, Rubin DB (1999). Small-sample degrees of freedom with
multiple imputation. *Biometrika*, 86(4), 948–955.
[doi:10.1093/biomet/86.4.948](https://doi.org/10.1093/biomet/86.4.948)

Rubin DB (1987). *Multiple Imputation for Nonresponse in Surveys*.
Wiley.
[doi:10.1002/9780470316696](https://doi.org/10.1002/9780470316696)

## References

Leyrat C, Seaman SR, White IR, Douglas I, Smeeth L, Kim J, Resche-Rigon
M, Carpenter JR, Williamson EJ (2019). Propensity score analysis with
partially observed covariates: how should multiple imputation be used?
*Statistical Methods in Medical Research*, 28(1), 3–19.
[doi:10.1177/0962280217713032](https://doi.org/10.1177/0962280217713032)

## Examples

``` r
dta <- sample_ps_data(n = 150, seed = 42)
fit <- fit_logistic(
  tavr ~ age + female + ef,
  data           = dta,
  family         = "binary",
  outcome_levels = c(0, 1),
  event_level    = 1
)
fit$tables$estimates[, c("term", "estimate", "odds_ratio")]
#>          term    estimate odds_ratio
#> 1 (Intercept)  3.55713646 35.0626500
#> 2         age -0.09113276  0.9128965
#> 3      female  0.18166589  1.1992135
#> 4          ef  0.05839068  1.0601291
head(fit$data[, c("id", "tavr", "prob")])
#>   id tavr      prob
#> 1  1    0 0.1657534
#> 2  2    0 0.5493138
#> 3  3    0 0.2856728
#> 4  4    0 0.2041096
#> 5  5    0 0.5254897
#> 6  6    0 0.3978214
```
