* The SAS oracle fixtures for `fit_logistic()` and `fit_linear()` now include
  a `PROC MIANALYZE EDF=` run, and the tests check the pooled degrees of
  freedom, p-values and intervals against it directly. No change in behaviour.
