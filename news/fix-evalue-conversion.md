* `sa_evalue()` now applies the odds-ratio and hazard-ratio conversions the
  right way round (VanderWeele & Ding 2017). For a rare outcome the OR or HR
  is used as the risk ratio directly; for a common outcome the OR becomes
  `sqrt(OR)` and the HR becomes `(1 - 0.5^sqrt(HR)) / (1 - 0.5^sqrt(1/HR))`.
  Earlier versions always took the square root, which understated the E-value
  for a rare outcome (an OR of 2.5 gave 2.54 instead of 4.44).
* Behaviour change: `sa_evalue()` gains a `rare` argument, and a call with
  `type = "OR"` or `"HR"` now aborts unless it sets `rare` or supplies `p0`
  (from which `rare` is taken as `p0 <= 0.15`). The returned list records the
  choice as `$rare`. The function also aborts unless
  `ci_lo <= estimate <= ci_hi`.
