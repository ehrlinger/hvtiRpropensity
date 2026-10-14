* The SAS oracle fixtures for `fit_logistic()` and `bs_count()` are exported at
  full precision, with the `PROC LOGISTIC` fits run at `GCONV=1e-12`. The
  binary, nominal and pooled comparisons, and the `PROC GENMOD` Poisson
  comparison, are tested at tolerances near their measured agreement. No change
  in behaviour.
