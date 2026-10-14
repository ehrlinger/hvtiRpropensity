* Pooled estimates from `fit_logistic()`, `fit_linear()`, `bs_continuous()`
  and `bs_count()` use Barnard and Rubin's (1999) small-sample degrees of
  freedom, with complete-data degrees of freedom n - k, as `mice::pool()` does.
  They used Rubin's (1987) large-sample degrees of freedom, which can run into
  the millions when little is missing. The `df` column, and the p-values and
  intervals built on it, change for every pooled model; estimates and standard
  errors do not. A coefficient with no between-imputation variance, which had
  infinite `df` and a normal test, now gets a finite `df`. To reproduce the
  values in SAS, run `PROC MIANALYZE` with `EDF=` set to n - k. The help pages for all four
  functions explain the calculation under "Pooled degrees of freedom".
* Pooling over imputations stops with an error when the first imputation's
  model has no residual degrees of freedom (n - k of zero or less), as a
  saturated model does. A saturated pooled logistic model returned meaningless
  `df` and p-values before.
