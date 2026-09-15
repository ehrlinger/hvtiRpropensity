# Changelog

## hvtiRpropensity 0.1.4

- Fixes an empty changelog on the pkgdown site. `NEWS.md` opened with a
  bare `# hvtiRpropensity` title line, so level one was the title level
  and every version heading sat at level two. pkgdown reads the top
  heading level in the file as the version level, found no versions
  there, and warned “no version headings found” on every build. That
  title line is removed and the four version headings are promoted to
  level one, matching the rest of the family. A level-one heading that
  names no version, such as the `(unreleased)` one above, is skipped by
  pkgdown rather than counted, so it does not reintroduce the problem.
  [`utils::news()`](https://rdrr.io/r/utils/news.html) was unaffected
  throughout and still reports the same four versions.

## hvtiRpropensity 0.1.3

- Removed the explicit `Maintainer:` field from `DESCRIPTION` and moved
  the maintainer address to `john.ehrlinger@gmail.com` in `Authors@R`,
  matching the rest of the `hvtiR*` family. The two fields had disagreed
  — `Authors@R` named the `ehrlinj@ccf.org` address while `Maintainer:`
  named the gmail one — which `R CMD check --as-cran` reported as a
  note. `Maintainer:` is now derived from the `cre` role at build time,
  so the two cannot drift apart again.

## hvtiRpropensity 0.1.2

- Qualified the [`rnorm()`](https://rdrr.io/r/stats/Normal.html) calls
  in the `sample_ps_data*()` generators as
  [`stats::rnorm()`](https://rdrr.io/r/stats/Normal.html), matching
  every other statistics call in the package, and declared `stats` under
  `Imports`. This clears the `R CMD check` note about an undefined
  global function. Generated data is unchanged — the seeds and the RNG
  draw order are identical.

## hvtiRpropensity 0.1.1

- Renamed from `hvtiPropensityScores` into the `hvtiR*` package family
  (`hvtiRtemplates`, `hvtiRutilities`, `hvtiRlifetables`,
  `hvtiRtables`). No user-facing function names or behaviour changed.

## hvtiRpropensity 0.1.0

- Initial development scaffold.
- Added
  [`sample_ps_data()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sample_ps_data.md)
  — reproducible synthetic cardiac-surgery dataset.
- Added
  [`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md)
  — greedy nearest-neighbour 1:1 propensity score matching without
  replacement, with optional caliper.
- Added
  [`ps_weight()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_weight.md)
  — IPTW weighting supporting ATE, ATT, and ATC estimands, with optional
  stabilisation and weight winsorisation.
- Added
  [`is_ps_data()`](https://ehrlinger.github.io/hvtiRpropensity/reference/is_ps_data.md)
  predicate and `print` / `summary` S3 methods for the common `ps_data`
  base class.
