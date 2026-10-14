* `fit_linear()`, `bs_continuous()`, `bs_count()` and `ps_forest()` check every
  column the formula names, subtracted ones included, as `fit_logistic()` does.
  A subtracted column missing from `data` reached `lm()` or `glm()` as "object
  not found", or picked up a same-named object from the calling environment;
  in `ps_forest()` a misspelled one was ignored, so the column meant to be
  dropped was fitted.
* `bs_count()` refuses a transformed response such as `log(count + 1)`, which
  it fitted as counts with a warning. `fit_logistic()` refuses one too,
  matching its error message.
* A `.` that leaves no columns once the response, `id_col` and
  `imputation_col` are set aside stops with a clear error.
