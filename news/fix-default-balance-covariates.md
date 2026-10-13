* The default covariates of the balance tables in `ps_weight()`, `ps_match()`, `ps_logistic()`, `ps_nominal()` and
  `ps_ordinal()` now leave out the id column, the treatment and score columns, the columns this package appends
  (`prob_t`, `logit_t`, `mt_wt`, `iptw`, `match`, `pair_id`, `quintile`, `decile`, `bs`, `cluster`, and the
  `prob_<level>` scores `ps_nominal()` and `ps_ordinal()` generate) and constant columns. A predictor named in the
  model formula is always kept. Before, `id` and `match` appeared in the `ps_weight()` table of `sample_ps_data()`,
  and `prob_t` in a `ps_match()` table matched on `logit_t`. A `covariates` vector you supply is used unchanged (#85).
* `ps_weight()` gains `id_col = "id"`, as `ps_match()` has, so the identifier stays out of its balance tables.
