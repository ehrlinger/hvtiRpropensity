* `ps_weight()` now stabilises weights by the marginal treatment probability
  for `estimand = "ATE"` only. For `"ATT"` and `"ATC"`, `stabilise` has no
  effect, so treated weights under the ATT and control weights under the ATC
  are exactly 1, as documented; before, the default `stabilise = TRUE` scaled
  them by the arm's marginal probability. `$meta$stabilised` records whether
  stabilisation was applied, and `sa_trim_sweep()` follows the same rule.
  ATE weights are unchanged (#88).
* `ps_match()` accepts only the number 1 for `ratio`. `ratio = 1.5` and
  `ratio = "1"` were silently accepted and now stop with an error naming the
  value (#88).
* `ps_support()` stops when a `flags` element is named `"trim"` or `"common"`,
  which the built-in rule silently overwrote, and when `trim` is not two
  increasing probabilities in \[0, 1\] (#88).
* An `NA` in `caliper` or `caliper_sd` (`ps_match()`) or `trim`
  (`ps_weight()`, `ps_support()`) now stops with an error naming the argument
  instead of "missing value where TRUE/FALSE needed" (#88).
