* `ps_stddiff()` now normalises the weights to mean 1 within each group before
  computing a weighted Gaussian or rank difference, so multiplying a group's
  weights by a constant no longer changes it (#82). As a result
  `ps_weight()$tables$smd_weighted` is the same with `stabilise = TRUE` and
  `stabilise = FALSE`, and `ps_stddiff_perm()` with `weight_col` inherits the
  fix.
* **Weighted SMDs change** for Gaussian and `nong_ord` variables whenever the
  weights do not already average 1 in each group, as matching weights and
  unstabilised IPTW weights do not. They then no longer equal the `%stddiff`
  macro's output, which depends on the scale of the weights; they still equal
  it when the weights average 1 in each group. Binary and categorical
  differences, and every unweighted difference, are unchanged.
