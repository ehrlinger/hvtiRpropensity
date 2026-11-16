* Errors follow the hvtiR family rubric. Every error is raised with base
  `stop(call. = FALSE)` and starts with the name of the function you called,
  as `ps_weight(): ...`, also when an internal check raises it. The `rlang`
  conditions are gone. Errors a caller may want to catch carry a class:
  `hvtiRpropensity_missing_column`, `hvtiRpropensity_one_arm`,
  `hvtiRpropensity_meta_conflict`, `hvtiRpropensity_tau`,
  `hvtiRpropensity_weight`, `hvtiRpropensity_aliased`,
  `hvtiRpropensity_convergence` and `hvtiRpropensity_separation`. Warnings
  follow the same pattern.
* `treatment_col` has no default in `ps_weight()`, `ps_match()`,
  `ps_mw_var()`, `ps_stddiff()` and `ps_stddiff_perm()`, nor in
  `sa_overlap()` and `sa_trim_sweep()` for a data frame. `"tavr"` was one
  study's column name. Omitting it still means `"tavr"` for this release, and
  warns once per session.
* `ps_rmst(weights =)` takes the family's upper-case estimand names, `"ATO"`,
  `"ATT"` and `"ATE"`. The `weighting` column, the `estimator` labels (such as
  `"all/ATO"`) and `$meta$stored_weighting` use them too. Lower case still
  works for this release, and warns. The weight columns in `$data` keep their
  lower-case names, such as `w_ato`.
* Renamed arguments, the old names working for this release with a warning:
  `ps_mw_var(n_rep =)` is `n_boot`, and `$meta$n_rep` is `$meta$n_boot` (the
  table column `n_rep`, the replicates that remained, is unchanged).
  `sample_ps_data_count(n_imputations =)` is `m`. The first argument of
  `ps_rmst()`, `ps_support()`, `sa_rosenbaum()`, `sa_overlap()` and
  `sa_trim_sweep()` is `data`, not `x`.
* Arguments that take one or more column names end in `_cols`, the old names
  working for this release with a warning: `covariates` is `covariate_cols` in
  `ps_match()`, `ps_weight()`, `ps_logistic()`, `ps_nominal()`, `ps_ordinal()`,
  `bs_continuous()` and `bs_count()`; `ps_mw_var(outcomes =)` is
  `outcome_cols`; and `ps_stddiff()` and `ps_stddiff_perm()` take
  `gaussian_cols`, `nong_ord_cols`, `binary_cols` and `categorical_cols`.
  The `$meta` field names are unchanged.
* An abbreviated argument name that used to reach a renamed argument by
  partial matching can now match both the old and the new name, and R stops
  with "argument matches multiple formal arguments". For example
  `ps_mw_var(outcome = "age")` matched `outcomes` before and now matches both
  `outcomes` and `outcome_cols`; `covariate =` and `cov =` do the same. Write
  the new name in full.
* `sa_overlap()` and `sa_trim_sweep()` given a `ps_data` object that records no
  `treatment_col` still use `"tavr"` for this release, with the same
  deprecation warning as a data frame.
* Every function that fits a model or draws random numbers records the family
  provenance record in `$meta$provenance`: `package`, `version`, `sha`,
  `method`, `seed`, `rng_kind` and `time`.
