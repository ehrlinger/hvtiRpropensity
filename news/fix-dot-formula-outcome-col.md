* `fit_logistic()` accepts a dot formula such as `y ~ .`. It stopped with
  "Required column(s) not found in `data`: ." before. The dot expands to every
  column except the outcome, `id_col` and `imputation_col`, as it already did in
  `bs_continuous()` and `fit_linear()`.
* `bs_count()` checks `outcome_col` the way `bs_continuous()` does: a value that
  is not the formula's response, or a response naming more than one column, is
  refused. A dot formula no longer fits the `id_col` and `imputation_col`
  columns as covariates.
