* The pooled-covariance check against `PROC MIANALYZE` in
  `test_lm_sas_oracle.R` now compares like with like. SAS's `TCov` is the
  multivariate total covariance, `(1 + rbar) * Ubar`, not the elementwise Rubin
  total that `fit_logistic()` reports, so the test checks the diagonal against
  the SAS standard errors and rebuilds `TCov` from the per-imputation fits, at
  `1e-5` rather than `1e-4`. No package code changed.
