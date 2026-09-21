# Estimate propensity scores via binary logistic regression

Fits `stats::glm(..., family = binomial())` and appends the estimated
propensity score, its logit, and a matching weight to the dataset.
Supports multiply-imputed data in a stacked "long" format.

## Usage

``` r
ps_logistic(
  formula,
  data,
  treatment_col = NULL,
  id_col = "id",
  imputation_col = NULL,
  score_col = "prob_t",
  logit_col = "logit_t",
  weight_col = "mt_wt",
  covariates = NULL,
  treatment_levels = NULL,
  treated_level = NULL
)
```

## Arguments

- formula:

  A formula with the binary treatment on the left-hand side, e.g.
  `tavr ~ age + female + ef + diabetes`.

- data:

  A data frame. If `imputation_col` is set, `data` must be a stacked
  "long" MI data frame as produced by
  `mice::complete(mids, action = "long")`.

- treatment_col:

  Name of the binary treatment column (0/1 or logical). If `NULL`
  (default), extracted from the LHS of `formula`.

- id_col:

  Patient identifier column. Required when `imputation_col` is set.
  Default `"id"`.

- imputation_col:

  Name of the imputation-index column in a stacked MI data frame. `NULL`
  (default) means a single complete dataset is supplied.

- score_col:

  Output column name for the propensity score (probability of
  treatment). Default `"prob_t"`.

- logit_col:

  Output column name for the logit of the propensity score,
  `log(p / (1 - p))`. Default `"logit_t"`.

- weight_col:

  Output column name for the matching weight. Default `"mt_wt"`.

- covariates:

  Character vector of covariate column names for SMD balance
  diagnostics. If `NULL` (default), all numeric columns other than
  `treatment_col`, `score_col`, `logit_col`, `weight_col`, `id_col`,
  `"quintile"`, and `"decile"` are used.

- treatment_levels:

  Complete binary treatment levels. `NULL` preserves the historical 0/1
  or logical interface.

- treated_level:

  Level whose probability is the propensity score. When
  `treatment_levels` is supplied and this is `NULL`, its last value is
  used.

## Value

An object of class `c("ps_logistic", "ps_data")` with:

- `$data`:

  The base data frame with `score_col`, `logit_col`, `weight_col`,
  `quintile`, and `decile` columns appended.

- `$meta`:

  Named list: `formula`, `treatment_col`, `id_col`, `imputation_col`,
  `score_col`, `logit_col`, `weight_col`, `method`, `n_imputations`,
  `n_total`.

- `$tables`:

  Named list: `smd`, `group_counts`.

## Details

When `imputation_col` is supplied, the model is fit separately on each
imputed dataset; predicted probabilities are averaged per patient across
imputations (Rubin's combination rule for predictions), and the first
imputed dataset's covariate values form the base output data frame. This
mirrors the SAS template workflow of `PROC LOGISTIC ... BY _IMPUTATION_`
followed by `PROC SUMMARY mean=`.

The **matching weight** appended as `weight_col` is
`min(p, 1-p) / (p * trt + (1-p) * (1-trt))`, which equals the `mt_wt`
column produced in the SAS templates (Li & Greene, 2013). It is the ATM
(average treatment effect among the matched) estimand weight and is used
as input to
[`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md)
or
[`ps_weight()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_weight.md).

## See also

[`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md),
[`ps_weight()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_weight.md),
[`ps_ordinal()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_ordinal.md),
[`ps_nominal()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_nominal.md),
[`sample_ps_data()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sample_ps_data.md)

## Examples

``` r
dta <- sample_ps_data(n = 200, seed = 42)
dta$prob_t <- NULL

# --- Single complete dataset (mirrors tp.lm.logistic_propensity_score.nomi.sas)
# Equivalent to: PROC LOGISTIC data=built descending; model tavr = ...;
obj <- ps_logistic(
  tavr ~ age + female + ef + diabetes + hypertension,
  data = dta
)
print(obj)
#> <ps_logistic>
#>   N total     : 400
#>   Treatment   : tavr
#>   PS column   : prob_t
#>   Weight col  : mt_wt
#>   Method      : logistic
#>   Tables      : smd, group_counts, estimates, covariance, by_imputation, fit_status 
summary(obj)
#> Summary of <ps_logistic>
#> 
#> Smd:
#>                  variable     smd
#> age                   age -0.7994
#> female             female  0.3989
#> ef                     ef  0.6139
#> diabetes         diabetes -0.0554
#> hypertension hypertension -0.0310
#> match               match      NA
#> 
#> Group counts:
#>     group   n
#> 1 control 200
#> 2 treated 200
#> 
#> Estimates:
#>           term    estimate  std.error  statistic  df      p.value    conf.low
#> 1  (Intercept)  4.24504327 1.30086827  3.2632384 Inf 1.101468e-03  1.69538831
#> 2          age -0.10638803 0.01575522 -6.7525593 Inf 1.452597e-11 -0.13726769
#> 3       female  0.75688483 0.23165843  3.2672449 Inf 1.085997e-03  0.30284264
#> 4           ef  0.06194245 0.01223150  5.0641747 Inf 4.101735e-07  0.03796915
#> 5     diabetes -0.06529637 0.25459818 -0.2564683 Inf 7.975892e-01 -0.56429962
#> 6 hypertension -0.01084919 0.23719636 -0.0457393 Inf 9.635180e-01 -0.47574552
#>     conf.high odds_ratio pooled
#> 1  6.79469822 69.7587783  FALSE
#> 2 -0.07550838  0.8990757  FALSE
#> 3  1.21092701  2.1316255  FALSE
#> 4  0.08591575  1.0639011  FALSE
#> 5  0.43370689  0.9367898  FALSE
#> 6  0.45404713  0.9892094  FALSE
#> 
#> Covariance:
#>               (Intercept)           age        female            ef
#> (Intercept)   1.692258249 -1.724556e-02  7.125339e-03 -6.768683e-03
#> age          -0.017245556  2.482268e-04 -3.290512e-04 -1.950182e-05
#> female        0.007125339 -3.290512e-04  5.366563e-02 -8.526194e-05
#> ef           -0.006768683 -1.950182e-05 -8.526194e-05  1.496096e-04
#> diabetes     -0.012736931 -1.199707e-04 -1.767623e-03  1.059612e-04
#> hypertension -0.039021927 -6.723910e-06 -1.699813e-03  1.024544e-04
#>                   diabetes  hypertension
#> (Intercept)  -0.0127369307 -3.902193e-02
#> age          -0.0001199707 -6.723910e-06
#> female       -0.0017676234 -1.699813e-03
#> ef            0.0001059612  1.024544e-04
#> diabetes      0.0648202310 -2.589114e-03
#> hypertension -0.0025891136  5.626211e-02
#> 
#> By imputation:
#>              imputation         term    estimate  std.error
#> (Intercept)           1  (Intercept)  4.24504327 1.30086827
#> age                   1          age -0.10638803 0.01575522
#> female                1       female  0.75688483 0.23165843
#> ef                    1           ef  0.06194245 0.01223150
#> diabetes              1     diabetes -0.06529637 0.25459818
#> hypertension          1 hypertension -0.01084919 0.23719636
#> 
#> Fit status:
#>   imputation converged n_input n_analyzed n_excluded
#> 1          1      TRUE     400        400          0
#> 

# The function appends:
#   prob_t   -- propensity score (p-hat; SAS _p_ / _propen_)
#   logit_t  -- log(p/(1-p))    (SAS _logit_)
#   mt_wt    -- matching weight  (SAS mt_wt = min(p,1-p)/(p*trt+(1-p)*(1-trt)))
#   quintile -- rank-based quintile (SAS int(_n_/(nobs/5))+1)
#   decile   -- rank-based decile
head(obj$data[, c("id", "tavr", "prob_t", "logit_t", "mt_wt",
                  "quintile", "decile")])
#>   id tavr    prob_t    logit_t     mt_wt quintile decile
#> 1  1    0 0.1112484 -2.0780524 0.1251738        1      1
#> 2  2    0 0.5270896  0.1084644 1.0000000        3      6
#> 3  3    0 0.3951929 -0.4255358 0.6534196        2      4
#> 4  4    0 0.2836606 -0.9263756 0.3959863        2      3
#> 5  5    0 0.3751620 -0.5101344 0.6004149        2      4
#> 6  6    0 0.3148051 -0.7777498 0.4594387        2      3

# Pass the scored data to ps_match() for downstream matching
matched <- ps_match(obj$data, score_col = obj$meta$score_col, seed = 42)
nrow(matched$data[matched$data$match == 1L, ])
#> [1] 400

# --- Multiply-imputed data (mirrors tp.lm.logistic_propensity_score.sas)
# Equivalent to:
#   PROC LOGISTIC data=built descending; BY _IMPUTATION_; model tavr = ...;
#   PROC SUMMARY data=decile; class ccfid; var _p_; output out=... mean=_propen;
# \donttest{
# Simulate a stacked MI dataset (2 imputations, column "_Imputation_")
dta_mi <- rbind(
  cbind(dta, `_Imputation_` = 1L),
  cbind(dta, `_Imputation_` = 2L)
)
names(dta_mi)[names(dta_mi) == "_Imputation_"] <- "imp"

obj_mi <- ps_logistic(
  tavr ~ age + female + ef + diabetes + hypertension,
  data           = dta_mi,
  imputation_col = "imp",
  id_col         = "id"
)
print(obj_mi)
#> <ps_logistic>
#>   N total     : 400
#>   Treatment   : tavr
#>   PS column   : prob_t
#>   Weight col  : mt_wt
#>   Method      : logistic-MI (2 imputations)
#>   Tables      : smd, group_counts, estimates, covariance, by_imputation, fit_status 
# Per-patient PS is the average across the two imputed-dataset predictions,
# matching the PROC SUMMARY mean= step in the SAS template.
head(obj_mi$data[, c("id", "tavr", "prob_t")])
#>   id tavr    prob_t
#> 1  1    0 0.1112484
#> 2  2    0 0.5270896
#> 3  3    0 0.3951929
#> 4  4    0 0.2836606
#> 5  5    0 0.3751620
#> 6  6    0 0.3148051
# }
```
