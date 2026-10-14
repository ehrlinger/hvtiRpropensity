* Pooled estimates from `fit_logistic()`, `fit_linear()`, `bs_continuous()`
  and `bs_count()` use Barnard and Rubin's (1999) small-sample degrees of
  freedom, with complete-data degrees of freedom n - k, as `mice::pool()` does.
  They used Rubin's (1987) large-sample degrees of freedom, which can run into
  the millions when little is missing. The `df` column, and the p-values and
  intervals built on it, change for every pooled model; estimates and standard
  errors do not. A pooled logistic, ordinal, nominal or count model now reports
  a finite `df` and a t test. To reproduce the values in SAS, run
  `PROC MIANALYZE` with `EDF=` set to n - k. The help pages for all four
  functions explain the calculation under "Pooled degrees of freedom".
* Pooling over imputations stops with an error when a model has no residual
  degrees of freedom, as a saturated model does, instead of returning `NaN`
  p-values and intervals.
