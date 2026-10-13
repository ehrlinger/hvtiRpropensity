* **Behaviour change:** `ps_ordinal()` now stops unless the treatment is an
  ordered factor or `treatment_levels` gives the level order. It used to
  order a character or unordered-factor treatment alphabetically, so
  `"Low"`/`"Medium"`/`"High"` was fitted as High < Low < Medium, with no
  warning and coefficients of the wrong sign. The error says how to declare
  the order (#87).
* `ps_logistic()` with `treated_level` alone now accepts any treatment with
  two values, such as a `"TAVR"`/`"SAVR"` factor, as `ps_forest()` already
  did, and gives the same scores as the 0/1 coding. A treatment with more
  than two values still stops (#87).
* `ps_logistic()`, `ps_ordinal()`, `ps_nominal()`, `fit_logistic()` and
  `fit_linear()` check `id_col` only when `imputation_col` is set, as their
  documentation says. Data without an `id` column no longer stops a
  single-dataset fit (#87).
