* Ordinal models in `fit_logistic()` and `ps_ordinal()` report standard
  errors from the expected (Fisher) information, as SAS `PROC LOGISTIC` does,
  instead of `MASS::polr()`'s observed information. The two differ for the
  cumulative logit; on the package's SAS oracle data the standard errors moved
  by up to 0.4%, and the p-values, intervals and pooled results with them.
  `vcov()` on the stored `polr()` model still returns the observed version.
* `polr()` runs at a tighter tolerance (`reltol = 1e-14`), so ordinal
  estimates move slightly (about 2e-5 on the oracle data) and now agree with
  SAS to about 1e-7.
