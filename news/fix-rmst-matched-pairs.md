* `ps_rmst()` on a `ps_match()` object now analyses the matched pairs only,
  so `n` is twice the number of matched pairs, and its bootstrap resamples
  whole pairs by `pair_id` instead of rows. A matched input runs
  `"unweighted"` by default and refuses score weighting, which would adjust
  twice. Earlier releases used every row of the cohort and ignored the
  pairing. `$meta$resampling` records `"pair"` or `"row"`.
* `ps_rmst()` on a `ps_weight()` object now uses its stored weight column for
  the scheme that matches its estimand (`"ate"` or `"att"`), so trimming is
  kept, where it used to recompute the weights from the score.
  `$meta$stored_weighting` names that scheme.
* `ps_rmst()` now stops when `tau` is past the last observed time in either
  arm of an analysed subset, counting only rows that enter the weighted curve
  (time and event present, weight above zero), where it used to carry the last survival value
  forward. The help page now also says that failed bootstrap draws are
  counted in `n_failed` and that `refit = FALSE` intervals ignore the
  uncertainty in estimating the score.
