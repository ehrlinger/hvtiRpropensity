# Fit an ordinary least-squares model bundle

Fits [`stats::lm()`](https://rdrr.io/r/stats/lm.html) to one dataset, or
once per imputation of stacked long-form data, and combines imputations
with Rubin's rules. A single fit is inferred on the t distribution with
its residual degrees of freedom, as `PROC REG` reports it. Pooled
coefficients use Barnard and Rubin's degrees of freedom, described
below.

## Usage

``` r
fit_linear(
  formula,
  data,
  outcome_col = NULL,
  id_col = "id",
  imputation_col = NULL,
  prediction_col = "fitted"
)
```

## Arguments

- formula:

  Model formula. The response may be transformed, e.g. `log(los) ~ age`.
  A `.` on the right-hand side expands to the data's other columns,
  excluding the response, `id_col` and `imputation_col`.

- data:

  A data frame, or stacked imputations when `imputation_col` is set.

- outcome_col:

  Outcome column. By default the one variable in the formula's response;
  if supplied, it must name that variable. The response must read
  exactly one column.

- id_col:

  Patient identifier column. Required only when `imputation_col` is set.

- imputation_col:

  Imputation-index column, or `NULL` for one dataset.

- prediction_col:

  Name of the fitted-value column added to `$data`, averaged across
  imputations.

## Value

An object of class `c("rm_fit", "ps_data")`. `$tables$fit_quality` gives
each fit's residual degrees of freedom, residual standard error and
R-squared; R-squared is `NA` for imputed data, where Rubin's rules
define none.

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

## See also

[`fit_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/fit_logistic.md),
[`bs_continuous()`](https://ehrlinger.github.io/hvtiRpropensity/reference/bs_continuous.md)

## Examples

``` r
dta <- sample_ps_data(n = 150, seed = 42)
fit <- fit_linear(ef ~ age + female, data = dta)
fit$tables$estimates[, c("term", "estimate", "std.error", "p.value")]
#>          term   estimate  std.error      p.value
#> 1 (Intercept) 68.7280984 5.12425021 2.156296e-32
#> 2         age -0.1867446 0.06779147 6.237331e-03
#> 3      female  0.1259720 1.20495525 9.168075e-01
```
