* Tests now assert what their names claim: the caliper, `trim = NULL`,
  estimand and stabilised-weight tests are pinned to hand-computed values, the
  imputation test checks `prob_t` against the mean of the per-imputation fits,
  the AUC test uses `wilcox.test()` as an independent oracle, and the
  `sa_trim_sweep()` tests reach the estimand and machine-epsilon paths they
  name. No package code changed (#93).
