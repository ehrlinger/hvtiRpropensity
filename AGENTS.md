# hvtiRpropensity

Propensity-score analysis for the HVTI CORR group: score estimation
([`ps_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_logistic.md),
[`ps_nominal()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_nominal.md),
[`ps_ordinal()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_ordinal.md),
[`ps_forest()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_forest.md)),
the general model bundle behind them
([`fit_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/fit_logistic.md),
[`validate_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/validate_logistic.md)),
application
([`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md),
[`ps_weight()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_weight.md),
[`ps_support()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_support.md),
[`ps_rmst()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_rmst.md)),
balance statistics
([`bs_continuous()`](https://ehrlinger.github.io/hvtiRpropensity/reference/bs_continuous.md),
[`bs_count()`](https://ehrlinger.github.io/hvtiRpropensity/reference/bs_count.md),
[`ps_stddiff()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_stddiff.md),
[`ps_stddiff_perm()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_stddiff_perm.md),
[`ps_mw_var()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_mw_var.md))
and four sensitivity methods
([`sa_rosenbaum()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_rosenbaum.md),
[`sa_evalue()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_evalue.md),
[`sa_overlap()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_overlap.md),
[`sa_trim_sweep()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_trim_sweep.md)),
plus
[`is_ps_data()`](https://ehrlinger.github.io/hvtiRpropensity/reference/is_ps_data.md)
and four synthetic data generators. Twenty-four exports; `NAMESPACE` is
the authority.

This file is the operational contract and applies in full. It is tool
neutral, so Codex and any other agent read the same rules. Claude Code
affordances live in `CLAUDE.md`, which imports this file.

## Definition of done

- [`devtools::test()`](https://devtools.r-lib.org/reference/test.html)
  passes. The runner is `tests/test-all.R`.
- [`devtools::check()`](https://devtools.r-lib.org/reference/check.html)
  is **0 errors, 0 warnings, 0 notes.** Verified 2026-09-30 at 0.1.10
  under `R CMD check --as-cran` with the manual built, from a clean
  `git archive` export.
- [`devtools::document()`](https://devtools.r-lib.org/reference/document.html)
  has been run and `man/` and `NAMESPACE` are committed with the source
  change.

## The standard of care

**This package is not headed to CRAN, and it is held to the CRAN
standard anyway.** The hvtiR packages are distributed internally to get
feedback faster than a public release allows, not to lower the bar.

- Before a release, audit against every chapter of the [CRAN
  Cookbook](https://contributor.r-project.org/cran-cookbook/), run
  `R CMD check --as-cran` with the manual and vignettes built, check
  reverse dependencies (`hvtiRtemplates` depends on this package) and
  check URLs.
- Every Cookbook or `--as-cran` finding is real work to do. Order
  findings by effort and risk, never by “would CRAN care”, and never
  offer “skip it, it’s internal” as an option.

## The automated gates

| workflow | fails on |
|----|----|
| `R-CMD-check.yaml` | `R CMD check` across platforms |
| `check-manual.yaml` | the PDF manual build |
| `lint.yaml` | [`lintr::lint_package()`](https://lintr.r-lib.org/reference/lint.html) |
| `pkgdown.yaml` | the site build |
| `house-style.yaml` | the composed house style in `.claude/house-style.md` |
| `test-coverage.yaml` | coverage upload |

## Rules for this repo

- **Every public function returns a `ps_data` subclass, built by
  [`new_ps_data()`](https://ehrlinger.github.io/hvtiRpropensity/reference/new_ps_data.md).**
  The object structure is guaranteed across every subclass and code
  depends on it:
  - `$data` — the original data frame with scores or weights appended
  - `$meta` — column names used, method parameters, formula, computed
    statistics
  - `$tables` — diagnostic tables (SMD before/after, group counts,
    effective N); may be an empty list, but the element exists

  Adding a `ps_*()` function means returning that shape through the
  constructor, not inventing a parallel one.
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) methods are defined
  once on the base class and inherited.
- **Lines are 120 characters here**, because the `ps_*()` and `bs_*()`
  constructors take many named arguments and read better whole than
  wrapped. ⚠️ The family runs 80, 100, 120 and 135. Read `.lintr`.
- **`indentation_linter` and `commented_code_linter` are OFF**,
  deliberately — both fire heavily on the aligned-argument style and
  catch no defect here. Everything else is lintr’s default and **is**
  enforced.
- **Test files are `test_*.R` with an underscore, and the runner is
  `tests/test-all.R`.** ⚠️ This matches `hvtiPlotR` and differs from
  `hvtiRutilities`, `hvtiRdatabuild`, `hvtiRtables` and
  `hvtiRbootstrap`, which use `test-*.R` and `tests/testthat.R`.
- **None of the four sensitivity methods needs a package outside
  `Imports`.**
  [`sa_rosenbaum()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_rosenbaum.md)
  computes its bounds in base R and does not use `rbounds`. Keep it that
  way: adding a hard dependency to the cheap methods removes the reason
  they exist.
- **Roxygen markdown is ENABLED** (`Roxygen: list(markdown = TRUE)`). ⚠️
  `hvtiRutilities` and `hvtiRtemplates` have no such field and need Rd
  markup instead.
- **`testthat` edition 3.** `VignetteBuilder` is **quarto**.

## Gotchas

- **`DESCRIPTION`’s `Date:` moves with `Version:`.** Update both in the
  same version-bump commit.
- The package is **0.x**: the API is not frozen, but the `ps_data`
  structure above is what every subclass and both base methods rely on,
  so changing it is a breaking change in practice.
- **Every function that draws random numbers restores the caller’s
  stream**, through
  [`withr::local_seed()`](https://withr.r-lib.org/reference/with_seed.html):
  the sample-data generators,
  [`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md),
  [`ps_rmst()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_rmst.md),
  [`ps_stddiff_perm()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_stddiff_perm.md)
  and
  [`ps_mw_var()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_mw_var.md).
  Qualify RNG calls
  ([`stats::rnorm()`](https://rdrr.io/r/stats/Normal.html)) rather than
  adding a
  [`utils::globalVariables()`](https://rdrr.io/r/utils/globalVariables.html)
  suppression.

## Git and versioning

- **Never push to `main`.** Branch, then open a PR and let the
  maintainer merge.
- **`main` is protected by a GitHub ruleset, and nothing in this repo
  records that.** A clone shows no trace of it, so it is stated here.
  The ruleset is named `protect main`, is identical across the hvtiR
  family repositories, and enforces four rules on the default branch: no
  deletion, no force-push, pull-request-only, and an **automatic Copilot
  code review** on every PR. A rejected push comes from the server, not
  a local hook. ⚠️ It currently requires **zero approvals**.
  `require_code_owner_review` is set but inert because no repository in
  the family has a `CODEOWNERS` file, so a PR can merge unreviewed.
- Versions are **straight three digits** (`0.1.1`). Never a `.9000`
  suffix or a fourth digit.
- **Patch-digit bumps only**, as fixes land. Minor and major are the
  maintainer’s decision.
- **Bump when you name a version, not when you merge.** A pull request
  lands without touching `Version:`. Its entry goes under a
  `# hvtiRpropensity (unreleased)` heading in `NEWS.md`, which you add
  when it is not already there. A separate commit then renames that
  heading to the new version and updates `DESCRIPTION` and its `Date`,
  at most once a day. The heading is gone again after a bump, so the
  next change re-adds it. `.claude/house-style.md` carries the rule and
  the reasoning.
- **A change that ships nothing gets no `NEWS.md` entry and no bump.**
  That is a pull request whose every changed file is left out of the
  tarball `R CMD build` produces, meaning the base branch’s
  `.Rbuildignore` excludes it: here `.github/`, `AGENTS.md` and
  `CLAUDE.md` among others. One shipped file means the change ships, and
  the usual rules apply. No user can observe a change that ships
  nothing, so the pull request and its commit message are the record.
  Read `.Rbuildignore` rather than judging by feel.

## Looking up a dependency’s API

Read the installed documentation before calling a function from a
dependency you have not used in this session. Arguments get renamed,
defaults change and functions are deprecated between releases, so recall
is not a reliable source. The version that matters is the one installed
in the library the gates run against; the minimum in `DESCRIPTION` is
the oldest version the code must still work with.

- Version: `Rscript -e 'packageVersion("pkg")'`.
- Help page as text:
  `Rscript -e 'tools::Rd2txt(utils:::.getHelpFile(help("fn", package = "pkg")))'`.
- Signature and exports: `args(pkg::fn)`, `getNamespaceExports("pkg")`.
- What changed between versions: `news(package = "pkg")` or the
  package’s `NEWS.md`.
- The source, when the help page is thin: `pkg:::fn` prints it.

pkgdown sites and CRAN pages describe the latest release, which may not
be the installed one; confirm the version before relying on them. If the
installed docs contradict what you expected, the docs win. If an
argument you need does not exist, say so instead of substituting a
guess.

## Change discipline

1.  **Think before coding.** Do not assume, ask. If the request is
    ambiguous or a name, path or signature is uncertain, surface the
    confusion rather than running with a guess.
2.  **Simplicity first.** Write the minimum that solves the stated
    problem.
3.  **Surgical changes.** Touch only what the task requires. Raise
    nearby problems separately.
4.  **Goal-driven execution.** State what done looks like before
    starting, and use tests as the criterion.

## Prose

Documentation prose follows the house voice, composed into
`.claude/house-style.md` and checked by `house-style.yaml`. A
propensity-score reader needs to know which estimand a function targets
and what balance it achieved — say both.
