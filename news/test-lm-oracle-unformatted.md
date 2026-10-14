* The SAS oracle fixtures for `fit_logistic()` and `bs_count()` carry full
  precision from `PROC LOGISTIC` fits run at `GCONV=1e-12`, and the binary,
  nominal, Poisson and pooled comparisons are tested at tolerances near their
  measured agreement. No change in behaviour.
