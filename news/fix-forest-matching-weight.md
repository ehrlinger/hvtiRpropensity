* **Behaviour change:** `ps_forest()` now writes Li-Greene matching weights
  (ATM; Li & Greene, 2013) to `weight_col` (default `mt_wt`), `min(p, 1 - p) / p`
  for treated patients and `min(p, 1 - p) / (1 - p)` for controls, exactly as
  `ps_logistic()` does. It previously wrote overlap (ATO) weights, `1 - p` and
  `p`, under the same name, so swapping one scoring function for the other
  silently changed the estimand. Anyone who read `mt_wt` from a `ps_forest()`
  result, for weighted balance, `ps_mw_var()` or
  `hvtiPlotR::hv_mirror_hist(weight_col = "mt_wt")`, gets different weights and
  should rerun. Both help pages now state the estimand (#83).
