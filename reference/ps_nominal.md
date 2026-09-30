# Estimate propensity scores via nominal (generalised logit) logistic regression

Fits a multinomial logistic regression model via
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) for a
nominal (unordered) treatment variable, equivalent to SAS
`PROC LOGISTIC ... / link=glogit`. One probability column per treatment
level is appended to the dataset.

## Usage

``` r
ps_nominal(
  formula,
  data,
  treatment_col = NULL,
  id_col = "id",
  imputation_col = NULL,
  ref_level = NULL,
  score_col_prefix = "prob",
  trace = FALSE,
  covariates = NULL,
  treatment_levels = NULL,
  smd_pairs = "reference"
)
```

## Arguments

- formula:

  A formula with the nominal treatment on the LHS, e.g.
  `rtyp ~ age + female + ef + diabetes`.

- data:

  A data frame (single or stacked MI).

- treatment_col:

  Name of the nominal treatment column. `NULL` (default) extracts from
  `formula`.

- id_col:

  Patient identifier column. Default `"id"`.

- imputation_col:

  Imputation-index column. `NULL` (default) means a single complete
  dataset.

- ref_level:

  Reference level for the multinomial model. `NULL` (default) uses the
  first factor level, matching `REF=first` in SAS.

- score_col_prefix:

  Prefix for per-level output columns (`prob_<level>` by default).
  Default `"prob"`.

- trace:

  Logical. If `FALSE` (default), suppresses
  [`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html)
  iteration messages.

- covariates:

  Covariate columns for the balance table. `NULL` (default) uses every
  numeric column other than the treatment, score, identifier and strata
  columns.

- treatment_levels:

  Complete nominal treatment levels. `NULL` preserves the levels
  inferred by the historical interface.

- smd_pairs:

  Which pairs of treatment levels the balance table compares.
  `"reference"` (default) compares each level with `ref_level`; `"all"`
  compares every pair. The levels have no order, so there is no
  `"adjacent"` choice.

## Value

An object of class `c("ps_nominal", "ps_data")` with:

- `$data`:

  Base data frame with one `prob_<level>` column per treatment level
  appended.

- `$meta`:

  Named list: `formula`, `treatment_col`, `id_col`, `imputation_col`,
  `score_cols`, `levels`, `ref_level`, `method`, `n_imputations`,
  `n_total`.

- `$tables`:

  Named list: `smd`, `group_counts`. `smd` has one row per covariate and
  pair of levels: `variable`, `level`, `versus`, `smd`.

## Details

When `imputation_col` is supplied, models are fit on each imputed
dataset and per-level probabilities are averaged per patient across
imputations.

## Package dependency

Requires nnet (a recommended R package, usually pre-installed). Install
with `install.packages("nnet")` if missing.

## See also

[`ps_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_logistic.md),
[`ps_ordinal()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_ordinal.md),
[`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md)

## Examples

``` r
# Mirrors the SAS template workflow:
#   PROC LOGISTIC ... / link=glogit  (generalised logit)
# The first factor level is used as the reference category, matching
# SAS REF=first.  Change ref_level to match a different reference.
# \donttest{
if (requireNamespace("nnet", quietly = TRUE)) {
  dta <- sample_ps_data_nominal(n = 300, seed = 42)
  obj <- ps_nominal(
    rtyp ~ age + female + ef + diabetes,
    data = dta
  )
  print(obj)

  # One probability column per treatment level (analogous to p_cos, p_per,
  # p_dev, p_ce from the PROC TRANSPOSE step in the SAS template).
  head(obj$data[, c("id", "rtyp", "prob_COS", "prob_PER", "prob_DEV", "prob_CE")])

  # Explicitly set a different reference level
  obj_ce <- ps_nominal(
    rtyp ~ age + female + ef + diabetes,
    data      = dta,
    ref_level = "CE"    # matches REF=last in SAS
  )
  obj_ce$meta$ref_level
}
#> <ps_nominal>
#>   N total     : 1200
#>   Treatment   : rtyp (4 levels)
#>   Reference   : COS
#>   Score cols  : prob_COS, prob_PER, prob_DEV, prob_CE
#>   Method      : nominal-logistic
#>   Tables      : smd, group_counts, estimates, covariance, by_imputation, fit_status 
#> [1] "CE"
# }
```
