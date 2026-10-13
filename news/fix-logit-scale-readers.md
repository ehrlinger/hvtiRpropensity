* `ps_rmst()`, `ps_support()`, `sa_overlap()` and `sa_trim_sweep()` now
  back-transform the score of a `ps_match()` object built with
  `scale = "logit"` with `plogis()`, so its weights, trim rule and overlap
  diagnostics are on the probability scale. Before, they read the logit as a
  probability: `ps_rmst()` clipped it to \[0.001, 0.999\], `ps_support()`
  trimmed logits at 0.1 and 0.9, and the two `sa_*()` functions stopped with a
  range error (#80). `ps_support()` and `ps_rmst()` record `$meta$scale`,
  carried from their input, so a chain such as
  `ps_rmst(ps_support(matched), ...)` still reads the logit correctly.
* `ps_match()` warns when `scale = "logit"` is given a score column whose
  every value lies in \[0, 1\], which is most likely a probability such as the
  default `prob_t`. Its \[0, 1\] range error now suggests `scale = "logit"`
  for a logit score.
