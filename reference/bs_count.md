# Estimate a balancing score via negative-binomial or Poisson regression

Fits a GLM with a count distribution
([`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html) for
negative binomial, or [`stats::glm()`](https://rdrr.io/r/stats/glm.html)
with Poisson family) on a (possibly multiply-imputed) dataset, averages
the per-patient **linear predictors** (log-scale, i.e. xbeta) across
imputations, and assigns rank-based strata. Mirrors
`tp.pm.count.balncing_score.sas` (PROC GENMOD `dist=nb link=log`

- PROC SUMMARY + the `%grp` macro).

## Usage

``` r
bs_count(
  formula,
  data,
  outcome_col = NULL,
  id_col = "id",
  imputation_col = NULL,
  dist = c("negbin", "poisson"),
  score_col = "bs",
  n_strata = 10L,
  strata_col = "cluster",
  create_binary_strata = TRUE,
  covariates = NULL
)
```

## Arguments

- formula:

  A formula with the count outcome, untransformed, on the LHS, e.g.
  `rbc_tot ~ female + agee + in_hct + hx_csurg + iv_cpb`. A `.` on the
  right-hand side expands to the data's other columns, excluding the
  outcome, `id_col` and `imputation_col`.

- data:

  A data frame, or a stacked MI data frame when `imputation_col` is set.

- outcome_col:

  Name of the count outcome column. `NULL` (default) extracts from
  `formula`, whose response must read exactly one column; a supplied
  value must name that column.

- id_col:

  Patient identifier column. Required, and checked, only when
  `imputation_col` is set. Default `"id"`.

- imputation_col:

  Imputation-index column. `NULL` (default) means a single complete
  dataset.

- dist:

  Distribution for the GLM. `"negbin"` (default) calls
  [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html) (negative
  binomial, log link); `"poisson"` calls
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html) with
  `family = poisson(link = "log")`.

- score_col:

  Name of the output balancing score column (averaged linear predictor
  on the log scale). Default `"bs"`.

- n_strata:

  Number of strata. Default `10L`.

- strata_col:

  Name of the cluster/stratum column. Default `"cluster"`.

- create_binary_strata:

  Logical. If `TRUE` (default), binary indicator columns `stra_1`
  through `stra_<n_strata>` are appended.

- covariates:

  Covariate columns for diagnostics.

## Value

An object of class `c("bs_count", "ps_data")` with:

- `$data`:

  `data` (the lowest-numbered imputation when `imputation_col` is set,
  without that column) with `score_col`, `"quintile"`, `"decile"`,
  `strata_col`, and (if `create_binary_strata`) `stra_1` through
  `stra_<n_strata>` appended.

- `$meta`:

  Named list: `formula`, `outcome_col`, `id_col`, `imputation_col`,
  `score_col`, `dist`, `theta` (the negative-binomial dispersion
  estimate, one per imputation, named by imputation index; `NULL` for
  `dist = "poisson"`), `n_strata`, `strata_col`, `bundle_version`,
  `model_family` (`"count"`), `package_versions`, `method`,
  `n_imputations`, `n_total`.

- `$tables`:

  Named list: `strata_counts` (`stratum`, `n`), `estimates`,
  `covariance`, `by_imputation` and `fit_status`, laid out as in
  [`fit_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/fit_logistic.md)
  except that `estimates` has `rate_ratio`, the exponentiated
  coefficient, in place of `odds_ratio`. With `imputation_col`,
  `estimates` and `covariance` are pooled with Rubin's rules, with
  Barnard and Rubin's degrees of freedom.

- `$models`:

  The fitted models, one per imputation, named by imputation index:
  [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html) or
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html) fits.

## Details

The linear predictor (xbeta, log-scale) rather than the fitted count is
used as the balancing score because it is on an unbounded continuous
scale more suitable for rank-based stratification.

## Package dependency

`dist = "negbin"` (the default) requires the MASS package (a recommended
R package, usually pre-installed). Install with
`install.packages("MASS")` if missing, or use `dist = "poisson"`.

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

[`bs_continuous()`](https://ehrlinger.github.io/hvtiRpropensity/reference/bs_continuous.md),
[`ps_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_logistic.md)

## Examples

``` r
# Mirrors the SAS template workflow:
#   PROC GENMOD dist=nb link=log (negative binomial for overdispersed counts)
#   xbeta (log-scale linear predictor) used as the balancing score
#   (not the fitted count -- log scale is more suitable for stratification).
# The %grp macro binary indicators stra_1 ... stra_k are created when
# create_binary_strata = TRUE (the default).
# \donttest{
dta <- sample_ps_data_count(n = 300, seed = 42)

# Poisson fallback (no MASS dependency)
obj <- bs_count(
  rbc_tot ~ female + age + diabetes + hypertension,
  data = dta,
  dist = "poisson"
)
print(obj)
#> <bs_count>
#>   N total     : 300
#>   Outcome     : rbc_tot
#>   Score col   : bs
#>   Distribution: poisson
#>   Strata      : 10 clusters (cluster)
#>   Method      : balancing-poisson
#>   Tables      : strata_counts, estimates, covariance, by_imputation, fit_status 
table(obj$data$cluster)    # 10 equal-ish strata by default
#> 
#>  1  2  3  4  5  6  7  8  9 10 
#> 30 30 30 30 30 30 30 30 30 30 

# Binary stratum indicators mirror the SAS %grp macro output
head(obj$data[, c("id", "bs", "cluster", "stra_1", "stra_2", "stra_10")])
#>   id        bs cluster stra_1 stra_2 stra_10
#> 1  1 0.8957304       9      0      0       0
#> 2  2 0.3222764       3      0      0       0
#> 3  3 0.6961051       8      0      0       0
#> 4  4 0.6770635       7      0      0       0
#> 5  5 0.6753845       7      0      0       0
#> 6  6 0.5995166       6      0      0       0

# Negative binomial (accounts for overdispersion; requires MASS)
if (requireNamespace("MASS", quietly = TRUE)) {
  obj_nb <- bs_count(
    rbc_tot ~ female + age + diabetes + hypertension,
    data = dta,
    dist = "negbin"
  )
  print(obj_nb)
}
#> <bs_count>
#>   N total     : 300
#>   Outcome     : rbc_tot
#>   Score col   : bs
#>   Distribution: negbin
#>   Strata      : 10 clusters (cluster)
#>   Method      : balancing-negbin
#>   Tables      : strata_counts, estimates, covariance, by_imputation, fit_status 
# }
```
