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
  treatment_levels = NULL
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

  Covariate columns for diagnostics.

- treatment_levels:

  Complete nominal treatment levels. `NULL` preserves the levels
  inferred by the historical interface.

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

  Named list: `group_counts`.

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
dta <- sample_ps_data_nominal(n = 300, seed = 42)
obj <- ps_nominal(
  rtyp ~ age + female + ef + diabetes,
  data = dta
)
print(obj)
#> <ps_nominal>
#>   N total     : 1200
#>   Treatment   : rtyp (4 levels)
#>   Reference   : COS
#>   Score cols  : prob_COS, prob_PER, prob_DEV, prob_CE
#>   Method      : nominal-logistic
#>   Tables      : group_counts, estimates, covariance, by_imputation, fit_status 

# One probability column per treatment level (analogous to p_cos, p_per,
# p_dev, p_ce from the PROC TRANSPOSE step in the SAS template).
head(obj$data[, c("id", "rtyp", "prob_COS", "prob_PER", "prob_DEV", "prob_CE")])
#>   id rtyp  prob_COS  prob_PER  prob_DEV    prob_CE
#> 1  1  COS 0.1055278 0.1780586 0.3405289 0.37588468
#> 2  2  COS 0.4016693 0.3427668 0.1668590 0.08870491
#> 3  3  COS 0.3225837 0.3404853 0.2131848 0.12374626
#> 4  4  COS 0.1762023 0.2953493 0.2866436 0.24180480
#> 5  5  COS 0.3849803 0.2885476 0.2067771 0.11969503
#> 6  6  COS 0.4507948 0.2804114 0.1759815 0.09281226

# Explicitly set a different reference level
obj_ce <- ps_nominal(
  rtyp ~ age + female + ef + diabetes,
  data      = dta,
  ref_level = "CE"    # matches REF=last in SAS
)
# }
```
