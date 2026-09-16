# hvtiRpropensity (unreleased)

* **`ps_stddiff()` computes standardized differences for every variable
  type**, ported from the CCF `%stddiff` SAS macro (Artis 2019): Gaussian,
  non-Gaussian or ordinal by pooled ranks, binary, and categorical by Yang and
  Dalton's Mahalanobis form, each optionally weighted. It returns a
  `ps_stddiff` object whose `$tables$stddiff` holds one row per variable. The
  denominator averages the two group variances as the macro does, so it can
  differ from the SMD tables `ps_match()`, `ps_logistic()` and `ps_weight()`
  report today when group sizes are unequal; moving those onto
  `ps_stddiff()` is a separate change. A categorical variable whose groups
  share no levels gets `NA` with a warning. Design:
  hvtiRtemplates `dev/specs/2026-09-16-standardized-difference-design.md`;
  tracked in #34.

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
