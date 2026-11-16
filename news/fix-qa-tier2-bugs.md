* `ps_weight(trim =)` trims only the arm it reweights. For the ATT the control
  weights are winsorised at quantiles of the control weights, and every treated
  weight stays exactly 1, as documented; the ATC is the mirror image. Before,
  the quantiles were pooled across both arms, so with few treated patients the
  treated weights of 1 were capped. The ATE still trims at pooled quantiles,
  since both arms are reweighted. `sa_trim_sweep()` follows the same rule, and
  its warning about infinite weights no longer says they cannot be trimmed.
  This changes trimmed ATT and ATC weights; re-run a study that used them.
* The balance tables of `ps_match()`, `ps_weight()`, `ps_logistic()`,
  `ps_forest()`, `ps_nominal()` and `ps_ordinal()` route each covariate by type.
  A 0/1 or logical column takes the binary standardized difference, a factor or
  character column the categorical one of Yang and Dalton (one row per column),
  and other numeric columns the Gaussian one. Factor, character and logical
  columns were left out of the default tables before, and naming one stopped.
  Binary columns move by a small amount, because they were treated as Gaussian.
  A categorical difference is unsigned whatever number of levels the rows hold,
  so a factor that loses a level after matching is not switched to the signed
  binary formula. The tables keep their `variable` and `smd` columns.
* `fit_logistic()` and `ps_logistic()` stop, with a condition of class
  `hvtiRpropensity_separation`, when a binary fit is separated: its
  coefficients keep growing under further Newton steps, so no maximum
  likelihood estimate exists. A converged fit in which one extreme patient's
  probability rounds to 0 or 1 is not separation; `ps_logistic()` fits it and
  warns, with class `hvtiRpropensity_saturated_score`, naming the patients. A
  term the model cannot estimate is named as aliased instead of reported as a
  non-finite coefficient.
* `ps_rmst()` warns, with class `hvtiRpropensity_bootstrap_failures`, when more
  than 5% of an estimator's bootstrap replicates fail. A replicate whose
  resampled follow-up in an arm ends before `tau` now fails instead of being
  carried flat to `tau`.
* A treatment column holding one arm is refused by every function that needs
  both, instead of failing later or returning empty tables.
* `ps_match()` breaks tied distances by the smallest `id_col`, among the
  controls within the caliper, so the match no longer depends on row order.
* `ps_weight()` names the patients who get a weight of 0 (a score of exactly 0
  for an ATT control, or 1 for an ATC treated patient), instead of the
  `ps_stddiff()` error about non-positive weights.
* The help, README and vignette no longer say `hvtiPlotR::hv_balance()` needs
  reshaped tables: hvtiPlotR 2.8.22 reads a `ps_match()` or `ps_weight()`
  object directly. A `ps_logistic()` or `ps_forest()` object has a single `smd`
  table, and `hv_balance()` stops on it.
