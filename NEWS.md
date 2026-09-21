# hvtiRpropensity (unreleased)

* `ps_ordinal()` now creates its documented rank-based quintile and decile
  columns from the averaged probability of the highest ordered treatment
  level, including for stacked imputations. Existing columns with either name
  are rejected instead of overwritten.
* `bs_count()` bundles now record their bundle version, count-model family and
  R/package versions, matching the saved-model metadata contract introduced in
  0.1.6.

# hvtiRpropensity 0.1.6

* Adds `fit_logistic()` for binary, proportional-odds ordinal, and nominal
  generalized-logit model bundles, including explicit outcome-level contracts,
  retained per-imputation fits, averaged patient predictions, and Rubin-pooled
  coefficients and covariance.
* Adds `validate_logistic()` for calibration, observed-versus-expected events,
  AUC, and Brier validation of a saved version-1 binary bundle without
  refitting or modifying it.
* Existing binary, ordinal, and nominal propensity functions and count
  balancing scores retain their scored columns and diagnostics while gaining
  fitted models, inference, covariance, and fit-status tables. Negative-binomial
  bundles also retain theta for every imputation.
* Stacked-imputation workflows now require the same patient keys in every
  imputation. They error on missing patients instead of silently averaging a
  patient's score over fewer imputations.

# hvtiRpropensity 0.1.5

* **`ps_mw_var()` estimates a matching-weight treatment effect with a
  bootstrap variance**, ported from the CCF `%mw_var` SAS macro (Rajeswaran
  2014, after Li and Greene 2013). For each outcome it reports the
  `PROC MEANS` weighted mean, SD and sum of weights by group, the difference,
  and from `n_rep` replicates resampled within each group the bootstrap SD,
  `PCTLDEF=1` percentiles (`quantile(type = 4)`), z and p. As in the macro the
  weights are held fixed across replicates, so the SD ignores the
  uncertainty in estimating them; the documentation says so. The macro's
  output labelled the sum of weights `N`, reported here as `sumwgt_0` and
  `sumwgt_1`.

* **`withr` moves from Suggests to Imports.** `ps_stddiff_perm()` and
  `ps_mw_var()` take a `seed` and restore the caller's random number stream
  afterwards, through `withr::local_seed()`. Restoring the stream by hand
  means writing `.Random.seed` in the global environment, which `R CMD check`
  reports as a NOTE.

* **`ps_stddiff_perm()` gives each standardized difference a permutation
  reference**, ported from the CCF `%stddiffci` SAS macro (Artis 2020). It
  reports the observed `ps_stddiff()` value beside the 2.5, 16, 50, 84 and
  97.5 percentiles of the same statistic over `n_perm` shuffles of the group
  labels, using SAS `PCTLDEF=5` (`quantile(type = 2)`), the macro's
  `PROC UNIVARIATE` default. The percentiles describe what label shuffling
  alone produces, not a confidence interval. Weights that depend on the group
  are recomputed for every permutation by a `reweight` function, which is
  required with `weight_col`; the macro instead expected those permuted
  weights to exist as columns already. A `seed` makes the result
  reproducible and leaves the caller's random number stream as it was.

* **`ps_stddiff()` computes standardized differences for every variable
  type**, ported from the CCF `%stddiff` SAS macro (Artis 2019): Gaussian,
  non-Gaussian or ordinal by pooled ranks, binary, and categorical by Yang and
  Dalton's Mahalanobis form, each optionally weighted. It returns a
  `ps_stddiff` object whose `$tables$stddiff` holds one row per variable. The
  denominator averages the two group variances as the macro does. A
  categorical variable whose groups share no levels gets `NA` with a warning.
  Design: hvtiRtemplates
  `dev/specs/2026-09-16-standardized-difference-design.md`; tracked in #34.

* **The SMD tables from `ps_match()`, `ps_logistic()` and `ps_weight()` now
  come from `ps_stddiff()`, and their numbers change.** Every covariate is
  treated as Gaussian, as before, but the denominator is now the macro's
  `sqrt((var1 + var0) / 2)`. The unweighted tables (`smd_before`, `smd_after`,
  `smd`, `smd_unweighted`) used to pool the SD by sample size, which differs
  whenever the groups are unequal, so **`smd_before` is the table most likely
  to move**; with equal groups it does not. The weighted table
  (`smd_weighted`) used to divide its variances by `sum(w)` and now divides by
  `n - 1`, as `PROC MEANS` does. The tables keep their `variable` and `smd`
  columns and 4-place rounding. A covariate named in `covariates` must now be
  numeric; a character column used to yield `NA` with a warning and now stops.

* **`ps_weight()` stops with a clear error when a propensity score of exactly
  0 or 1 gives an infinite weight**, naming the count and pointing to `trim`.
  It used to fail later with "missing value where TRUE/FALSE needed" from
  inside the SMD calculation. With `trim` set, the infinite weight is
  winsorised as before and the call succeeds.

# hvtiRpropensity 0.1.4

* Fixes an empty changelog on the pkgdown site. `NEWS.md` opened with a bare
  `# hvtiRpropensity` title line, so level one was the title level and every
  version heading sat at level two. pkgdown reads the top heading level in the
  file as the version level, found no versions there, and warned "no version
  headings found" on every build. That title line is removed and the four
  version headings are promoted to level one, matching the rest of the family.
  A level-one heading that names no version, such as the `(unreleased)` one
  above, is skipped by pkgdown rather than counted, so it does not reintroduce
  the problem. `utils::news()` was unaffected throughout and still reports the
  same four versions.

# hvtiRpropensity 0.1.3

* Removed the explicit `Maintainer:` field from `DESCRIPTION` and moved the
  maintainer address to `john.ehrlinger@gmail.com` in `Authors@R`, matching the
  rest of the `hvtiR*` family. The two fields had disagreed — `Authors@R` named
  the `ehrlinj@ccf.org` address while `Maintainer:` named the gmail one — which
  `R CMD check --as-cran` reported as a note. `Maintainer:` is now derived from
  the `cre` role at build time, so the two cannot drift apart again.

# hvtiRpropensity 0.1.2

* Qualified the `rnorm()` calls in the `sample_ps_data*()` generators as
  `stats::rnorm()`, matching every other statistics call in the package, and
  declared `stats` under `Imports`. This clears the `R CMD check` note about
  an undefined global function. Generated data is unchanged — the seeds and
  the RNG draw order are identical.

# hvtiRpropensity 0.1.1

* Renamed from `hvtiPropensityScores` into the `hvtiR*` package family
  (`hvtiRtemplates`, `hvtiRutilities`, `hvtiRlifetables`, `hvtiRtables`).
  No user-facing function names or behaviour changed.

# hvtiRpropensity 0.1.0

* Initial development scaffold.
* Added `sample_ps_data()` — reproducible synthetic cardiac-surgery dataset.
* Added `ps_match()` — greedy nearest-neighbour 1:1 propensity score matching
  without replacement, with optional caliper.
* Added `ps_weight()` — IPTW weighting supporting ATE, ATT, and ATC estimands,
  with optional stabilisation and weight winsorisation.
* Added `is_ps_data()` predicate and `print` / `summary` S3 methods for the
  common `ps_data` base class.
