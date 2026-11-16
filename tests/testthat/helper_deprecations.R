# The suite calls many functions without `treatment_col`, which now warns once
# per session (it used to default to "tavr"). Mark those warnings as already
# given, so they do not land inside unrelated expect_warning() and
# expect_no_warning() calls. test_conditions.R resets the record and tests the
# deprecations themselves.
for (fn in c("ps_weight", "ps_match", "ps_mw_var", "ps_stddiff", "ps_stddiff_perm", "sa_overlap", "sa_trim_sweep")) {
  assign(paste0(fn, ":treatment_col"), TRUE, envir = hvtiRpropensity:::.deprecation_seen)
}
rm(fn)
