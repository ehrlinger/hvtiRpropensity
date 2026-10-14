* `fit_logistic()` accepts a dot formula such as `y ~ .`. It stopped with
  "Required column(s) not found in `data`: ." before. The dot expands to every
  column except the outcome, `id_col` and `imputation_col`, as it already did in
  `bs_continuous()` and `fit_linear()`.
* `bs_count()` checks `outcome_col` the way `bs_continuous()` does: a value that
  is not the formula's response, or a response naming more than one column, is
  refused. A dot formula no longer fits the `id_col` and `imputation_col`
  columns as covariates.
* `fit_logistic()`, `fit_linear()`, `bs_continuous()` and `bs_count()` check
  only the predictors the model keeps. A variable the formula subtracts, as in
  `y ~ . - dropme`, is no longer checked for matching levels across
  imputations, and no longer listed in `$meta$predictors`.
* `fit_linear()`, `bs_continuous()` and `bs_count()` require the response to
  read exactly one column even when `outcome_col` is supplied. `I(count + other)`
  with `outcome_col = "count"` fitted the sum while recording `count` as the
  outcome; `cbind(count, other)` stopped with an unhelpful subscript error.
