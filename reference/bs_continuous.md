# Estimate a balancing score via linear regression

Fits [`stats::lm()`](https://rdrr.io/r/stats/lm.html) on a (possibly
multiply-imputed) dataset, averages the per-patient **fitted values**
across imputations, and assigns rank-based quintile, decile, and cluster
strata. Mirrors the workflow in `tp.rm.continuous.balncing_score.sas`
(PROC REG + PROC SUMMARY + the `%grp` macro for binary stratum
indicators).

## Usage

``` r
bs_continuous(
  formula,
  data,
  outcome_col = NULL,
  id_col = "id",
  imputation_col = NULL,
  score_col = "bs",
  n_strata = 10L,
  strata_col = "cluster",
  create_binary_strata = TRUE,
  covariates = NULL
)
```

## Arguments

- formula:

  A formula with the continuous outcome / exposure on the LHS, e.g.
  `hct_nadr ~ female + agee + in_bmi + hct_pr + iv_cpb`. A `.` on the
  right-hand side expands to the data's other columns, excluding the
  response, `id_col` and `imputation_col`.

- data:

  A data frame, or a stacked MI data frame when `imputation_col` is set.

- outcome_col:

  Name of the continuous outcome column. `NULL` (default) extracts from
  `formula`, whose response must read exactly one column; a supplied
  value must name that column.

- id_col:

  Patient identifier column. Required, and checked, only when
  `imputation_col` is set. Default `"id"`.

- imputation_col:

  Imputation-index column for stacked MI data. `NULL` (default) assumes
  a single complete dataset.

- score_col:

  Name of the output balancing score column (averaged fitted value).
  Default `"bs"`.

- n_strata:

  Number of strata (clusters). Also controls how many binary indicator
  columns are created when `create_binary_strata = TRUE`. Default `10L`.

- strata_col:

  Name of the cluster/stratum column appended to `$data`. Default
  `"cluster"`.

- create_binary_strata:

  Logical. If `TRUE` (default), binary columns `stra_1` through
  `stra_<n_strata>` are appended, mirroring the `%grp` macro in the SAS
  template.

- covariates:

  Covariate columns for diagnostics (not used in computation; reserved
  for future balance-checking helpers).

## Value

An object of class `c("bs_continuous", "ps_data")` with:

- `$data`:

  Base data frame with `score_col`, `"quintile"`, `"decile"`,
  `strata_col`, and (if `create_binary_strata`) `stra_1` through
  `stra_<n_strata>` appended.

- `$meta`:

  Named list: `formula`, `outcome_col`, `id_col`, `imputation_col`,
  `score_col`, `n_strata`, `strata_col`, `bundle_version`,
  `model_family`, `package_versions`, `method`, `n_imputations`,
  `n_total`.

- `$tables`:

  Named list: `strata_counts` (`stratum`, `n`), `estimates`
  (Rubin-pooled coefficients when imputed, with Barnard and Rubin's
  degrees of freedom), `covariance`, `by_imputation` and `fit_status`,
  laid out as in
  [`fit_linear()`](https://ehrlinger.github.io/hvtiRpropensity/reference/fit_linear.md).

- `$models`:

  The fitted [`stats::lm()`](https://rdrr.io/r/stats/lm.html) models,
  one per imputation, named by imputation index.

## Details

The balancing score is the expected value of the continuous exposure
(e.g. nadir haematocrit) given the full covariate set. Patients in the
same stratum have similar expected exposure values, making
within-stratum comparisons less confounded.

A formula with an aliased term, such as a covariate that duplicates
another, is an error naming the imputation, rather than a score from a
rank-deficient fit.

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

[`bs_count()`](https://ehrlinger.github.io/hvtiRpropensity/reference/bs_count.md),
[`fit_linear()`](https://ehrlinger.github.io/hvtiRpropensity/reference/fit_linear.md),
[`ps_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_logistic.md),
[`sample_ps_data()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sample_ps_data.md)

## Examples

``` r
# Use ef (ejection fraction) as a stand-in for a continuous exposure
dta <- sample_ps_data(n = 300, seed = 42)
obj <- bs_continuous(
  ef ~ age + female + diabetes + hypertension,
  data    = dta,
  id_col  = "id"
)
print(obj)
#> <bs_continuous>
#>   N total     : 600
#>   Outcome     : ef
#>   Score col   : bs
#>   Strata      : 10 clusters (cluster)
#>   Method      : balancing-linear
#>   Tables      : strata_counts, estimates, covariance, by_imputation, fit_status 
table(obj$data$cluster)
#> 
#>  1  2  3  4  5  6  7  8  9 10 
#> 60 60 60 60 60 60 60 60 60 60 
```
