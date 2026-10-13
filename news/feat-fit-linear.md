* `fit_linear()` fits an ordinary least-squares model to one dataset or to
  stacked imputations, with Rubin-pooled inference verified against
  `PROC MIANALYZE`. It returns the same bundle as `fit_logistic()`.
* `bs_continuous()` now carries fit status, coefficient estimates, covariance
  and per-imputation estimates, as `bs_count()` does. Its balancing score is
  unchanged. A transformed exposure, such as `log(hct_nadr) ~ ...`, now works.
* `bs_continuous()` now stops, naming the imputation, when a term is aliased
  (for example a covariate that duplicates another), where it used to return
  a score from a rank-deficient fit.
