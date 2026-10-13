* `fit_linear()` fits an ordinary least-squares model to one dataset or to
  stacked imputations, with Rubin-pooled inference verified against
  `PROC MIANALYZE`. It returns the same data, models, meta and tables shape
  as `fit_logistic()`, and its `$tables$fit_quality` gives each fit's
  residual standard error and R-squared. A single fit, here and in
  `bs_continuous()`'s new estimates, is inferred on the t distribution with
  residual degrees of freedom, as `PROC REG` reports it. A dot formula such
  as `ef ~ .` expands over the data's columns excluding the identifier and
  imputation columns, here and in `bs_continuous()`, so the patient ID is
  never fitted as a covariate.
* `bs_continuous()` now carries fit status, coefficient estimates, covariance
  and per-imputation estimates, as `bs_count()` does. Its balancing score is
  unchanged. A transformed exposure, such as `log(hct_nadr) ~ ...`, now works.
* `bs_continuous()` now stops, naming the imputation, when a term is aliased
  (for example a covariate that duplicates another), where it used to return
  a score from a rank-deficient fit.
