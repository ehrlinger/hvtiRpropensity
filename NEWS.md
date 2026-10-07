# hvtiRpropensity (unreleased)

* `DESCRIPTION` now declares the Quarto command line tool in
  `SystemRequirements`. The vignettes have always needed it to build; the
  field makes that visible to installers and to `R CMD check`.

# hvtiRpropensity 0.1.11

* A seed now makes every random function reproducible whatever the session
  drew before it. `ps_forest()`, `ps_match()`, `ps_rmst()`,
  `ps_stddiff_perm()` and `ps_mw_var()` set the stream to `abs(seed)`
  immediately before their draws and restore the caller's stream afterwards.
  `ps_forest()` now does this as well as passing `seed = -abs(seed)` to
  `randomForestSRC::rfsrc()`; before, it set no R seed. The default stays
  `seed = NULL`, which draws from the caller's stream as before. A positive
  seed gives the same result as in 0.1.10; a negative one now gives the same
  result as its absolute value.
* Each of those five functions records the seed it used in `$meta$seed`, `NA`
  when none was given. `ps_stddiff_perm()` and `ps_mw_var()` used to store
  `NULL` there. A seed that is not one whole number is an error.
* `ps_forest()` documents the thread trade-off: if two runs with the same seed
  disagree, set `options(rf.cores = 1L)`.
* A new test scans every function in the package and fails on a random draw
  that no seeding call precedes.
* `sample_ps_data()`, `sample_ps_data_ordinal()`, `sample_ps_data_nominal()` and
  `sample_ps_data_count()` no longer warn when called with `seed = NULL`. `NULL`
  now draws from the caller's random-number stream, as `ps_match()` and
  `ps_rmst()` do, so `set.seed()` before the call reproduces the data. Before,
  the generators reseeded from the clock and warned that `.Random.seed` was
  `NULL`.
* The generators now check `seed` and use its absolute value, as the five
  functions above do: a negative seed gives the same data as its absolute
  value, and a seed that is not one whole number is an error.
* Now requires R 4.4.0 or newer, up from 4.1.0, to match the rest of the
  HVTI family. `hvtiR::install()` installs the members together, and several
  already required 4.4.0, so on an older R the install failed whatever this
  package declared.

# hvtiRpropensity 0.1.10

* New tests compare `ps_stddiff()`, `ps_stddiff_perm()` and `ps_mw_var()` with
  the CCF macros `%stddiff`, `%stddiffci` and `%mw_var`, run in SAS on a
  synthetic dataset (`tests/testthat/fixtures/stddiff-sas/`). All three agree
  with the macros: standardized differences to within 4e-7, the `%mw_var`
  estimate to 3e-10 (#34).
* New `treated_level` argument on `ps_match()`, `ps_weight()`, `ps_stddiff()`,
  `ps_stddiff_perm()`, `ps_mw_var()`, `sa_overlap()`, `sa_trim_sweep()` and
  `sa_rosenbaum()`. It names the treated value, so the treatment column may
  hold any two values (for example `"surgical"` and `"transcatheter"`), as
  `ps_logistic()` already allowed. Results equal those from the same data
  coded 0/1, the group tables show the study's own level names, and `$data`
  keeps the column as it was given. `ps_forest()` accepts such a column through
  its existing `treated_level`. The default, `NULL`, keeps the 0/1 or logical
  requirement, so existing calls are unchanged (#60).
* `sa_overlap()`, `sa_trim_sweep()` and `sa_rosenbaum()` read the treated level
  from a scored or matched object, so an object scored with
  `treated_level = 0` is no longer analysed with the arms swapped.
* `sa_trim_sweep()` and `sa_rosenbaum()` now check the treatment column. A
  column that is not 0/1 or logical used to be coerced with warnings and
  return a result; it is now an error unless `treated_level` is given.
* `sa_trim_sweep()` excludes patients whose treatment is missing, with a
  warning. It used to keep them, which made every effective sample size `NA`.
* The group tables of `ps_match()`, `ps_weight()`, `ps_forest()` and
  `sa_overlap()` no longer label the groups `control` and `treated`. `group`
  now holds the value in `treatment_col` (`0`/`1`, or `FALSE`/`TRUE`), and a
  new logical `treated` column says which group is treated, as in
  `ps_logistic()`. **Code that matched `group == "treated"` should use the
  `treated` column** (#51).
* `ps_ordinal()` and `ps_nominal()` now return a balance table, `$tables$smd`,
  as `ps_logistic()` does. It has one row per covariate and pair of treatment
  levels (`variable`, `level`, `versus`, `smd`), each computed on the patients
  in those two levels. A new `smd_pairs` argument chooses the pairs:
  `"reference"` (default) compares each level with the reference, `"adjacent"`
  (ordinal only) compares each level with the one below it, and `"all"`
  compares every pair. The choices can be combined (#48).
* `ps_rmst()` now defaults to `seed = NULL`, as `ps_match()`,
  `ps_stddiff_perm()` and `ps_mw_var()` already did. The default was
  `seed = 1024L`, which gave the same bootstrap draws on every call. **A call
  that relied on the default now gives an interval that varies from run to
  run; pass `seed` for a reproducible one.** Point estimates are unchanged.
* `sample_ps_data_count(seed = NULL)` no longer errors and returns a fresh
  dataset, as the other three generators did. The generators keep their fixed
  default, `seed = 42L`, so that examples share one demo dataset (#53).
* `fit_logistic()` and `validate_logistic()` have examples.
* The examples for `ps_ordinal()`, `ps_nominal()` and `ps_rmst()` run only when
  the suggested package they need ('MASS', 'nnet', 'survival') is installed.
* The package help no longer says `sa_rosenbaum()` requires 'rbounds'. It never
  did: the bounds are computed in base R.
* `DESCRIPTION` spells out CORR and SMD, quotes 'hvtiPlotR' and cites
  Rosenbaum and Rubin (1983). Help pages use ASCII punctuation in place of
  em dashes and arrows (#53).

# hvtiRpropensity 0.1.9

* `ps_match()`, `ps_rmst()` and the `sample_ps_data*()` generators restore the
  caller's random number stream when they return. Each called `set.seed()`
  and left the stream reseeded, so everything the caller drew afterwards
  depended on the seed passed here. They now use `withr::local_seed()`, as
  `ps_stddiff_perm()` and `ps_mw_var()` already did. Results for a given seed
  are unchanged.
* `bs_count()` and other count-family models now name the exponentiated
  coefficient column `rate_ratio`. It was labelled `odds_ratio`, but the
  exponent of a log-link Poisson or negative-binomial coefficient is a rate
  ratio. Binary, ordinal and nominal models keep `odds_ratio` (#46).
* `ps_logistic()`'s `$tables$group_counts` names each row by the declared
  treatment level instead of `control`/`treated`, and adds a logical
  `treated` column, so it reads the same as `ps_ordinal()` and `ps_nominal()`.
  Without `treatment_levels` the groups are `0`/`1` (or `FALSE`/`TRUE`) (#47).

# hvtiRpropensity 0.1.8

* New `ps_forest()`: out-of-bag random-forest propensity score with the same
  `ps_data` contract as `ps_logistic()` (needs the suggested
  `randomForestSRC`).
* New `ps_support()`: flags weak-support patients under common-support,
  score-trimming and user-supplied rules (e.g. isolation-forest tails) and
  tabulates pairwise agreement (kappa, Jaccard).
* New `ps_rmst()`: observed-outcome restricted mean survival contrast from
  weighted Kaplan-Meier curves with a stratified bootstrap interval; its
  `$tables$estimates` / `$tables$curves` feed `hvtiPlotR` RMST plots.

# hvtiRpropensity 0.1.7

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
