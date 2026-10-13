* Breaking change: `sa_rosenbaum()`, `sa_evalue()`, `sa_overlap()` and
  `sa_trim_sweep()` now return `ps_data` objects of classes `sa_rosenbaum`,
  `sa_evalue`, `sa_overlap` and `sa_trim_sweep`, with the four slots `$data`,
  `$meta`, `$tables` and `$models`, a `print()` method each, and the shared
  `summary()`. Code that read the old lists reaches the same values here:
  * `sa_rosenbaum()`: `res$bounds` is `res$tables$bounds`; `res$sensitivity_value`
    and `res$n_pairs` are `res$meta$sensitivity_value` and `res$meta$n_pairs`.
    The pairs used are in the new `res$tables$pairs`.
  * `sa_evalue()`: `evalue_estimate`, `evalue_ci`, `rr_estimate`, `type` and
    `rare` move under `res$meta`, and `res$tables$evalue` holds both E-values
    as a table. `res$data` is an empty data frame.
  * `sa_overlap()`: `summary`, `outside_overlap` and `positivity_flags` move
    under `res$tables`, and `overlap_region` under `res$meta`.
  * `sa_trim_sweep()`: the data frame it returned is `res$tables$sweep`.
* Breaking change: the four `sa_*()` functions follow one rule for arguments a
  `ps_data` input already records. A supplied argument that disagrees with the
  object's `$meta` now aborts, naming both values, where it used to be silently
  overridden: `sa_rosenbaum(m, "ef", treated_level = 0)`,
  `sa_trim_sweep(w_ate, estimand = "ATT")` and
  `sa_overlap(m, score_col = "other")` all abort. An argument not supplied
  takes the metadata, then its default. `sa_rosenbaum()` now also reads
  `pair_id_col` from the metadata.
* `sa_rosenbaum()` corrects its variance for ties, using the sum of squared
  average ranks in place of `n(n+1)(2n+1)/6`. Without tied absolute
  differences the bounds are unchanged; with ties they change slightly. At
  Gamma = 1 the bound now equals the one-sided `wilcox.test()` p-value with
  `exact = FALSE`.
* `sa_rosenbaum()` gains `alternative = c("auto", "greater", "less")`. The
  default `"auto"` keeps the earlier results, picking the direction from the
  sign of the median difference; the documentation now says so, and that
  Rosenbaum's method assumes the direction is fixed in advance.
* `sa_overlap()` excludes patients with a missing treatment as well as those
  with a missing score, and `sa_trim_sweep()` excludes patients with a missing
  score, with a warning giving the count, so the effective sample sizes and
  `pct_trimmed` are computed on the analysed rows (part of #84).
