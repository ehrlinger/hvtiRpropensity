* The vignette's `hvtiPlotR` examples now run. `hv_mirror_hist()` takes column
  names as strings (`score_col`, `group_col`, `match_col`, `weight_col`) and
  returns a data object for `plot()`. `hv_balance()` does not read `$tables`
  directly: it needs one long data frame with `variable`, `group` and
  `std_diff` in percent, and the vignette now shows the reshape. The package
  help and the README no longer say the tables feed it as they are.
  'hvtiPlotR' is now in `Suggests` (#91).
* The vignette now covers `ps_forest()`, `ps_stddiff()`, `ps_stddiff_perm()`,
  `ps_mw_var()`, `ps_rmst()` on matched pairs, `ps_support()`, matching on the
  logit with `caliper_sd = 0.2`, and the `rare` argument `sa_evalue()` needs
  for odds and hazard ratios. Every `ps_match()` call in it gives a seed. It
  no longer says all four `sa_*()` functions measure unmeasured confounding:
  `sa_overlap()` and `sa_trim_sweep()` check positivity and overlap.
* The README quick start runs: it weighted a data frame whose score column it
  had just removed. Its status block now says which parts SAS oracles check
  (since 0.1.10), and its function tables list all 25 exports.
* The package help page lists all 25 exports, says which return `ps_data`
  objects (all but the four data generators and `is_ps_data()`), and points
  to the SAS programs shipped in
  `system.file("sas-templates", package = "hvtiRpropensity")`.
