# hvtiRpropensity

Propensity-score analysis for the HVTI CORR group: score estimation
([`ps_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_logistic.md),
[`ps_nominal()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_nominal.md),
[`ps_ordinal()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_ordinal.md)),
application
([`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md),
[`ps_weight()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_weight.md)),
balance statistics
([`bs_continuous()`](https://ehrlinger.github.io/hvtiRpropensity/reference/bs_continuous.md),
[`bs_count()`](https://ehrlinger.github.io/hvtiRpropensity/reference/bs_count.md))
and four sensitivity methods
([`sa_rosenbaum()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_rosenbaum.md),
[`sa_evalue()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_evalue.md),
[`sa_overlap()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_overlap.md),
[`sa_trim_sweep()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_trim_sweep.md)),
plus synthetic data generators. Sixteen exports.

This file is the operational contract and applies in full. It is tool
neutral, so Codex and any other agent read the same rules. Claude Code
affordances live in `CLAUDE.md`, which imports this file.

## Definition of done

- [`devtools::test()`](https://devtools.r-lib.org/reference/test.html)
  passes. The runner is `tests/test-all.R`.
- [`devtools::check()`](https://devtools.r-lib.org/reference/check.html)
  is **0 errors, 0 warnings, and no NEW notes.** ⚠️ **This package is
  not at 0/0/0.** As of 2026-08-20 one NOTE stands: `rnorm` is used in
  the `sample_ps_data*()` generators without
  `importFrom("stats", "rnorm")`. It is a real defect, not an artifact,
  and it is not yours unless you touch those functions — but do not let
  a second note hide behind it.
- [`devtools::document()`](https://devtools.r-lib.org/reference/document.html)
  has been run and `man/` and `NAMESPACE` are committed with the source
  change.

## The automated gates

| workflow             | fails on                                             |
|----------------------|------------------------------------------------------|
| `R-CMD-check.yaml`   | `R CMD check` across platforms                       |
| `check-manual.yaml`  | the PDF manual build                                 |
| `lint.yaml`          | `lintr::lint_package()`                              |
| `pkgdown.yaml`       | the site build                                       |
| `house-style.yaml`   | the composed house style in `.claude/house-style.md` |
| `test-coverage.yaml` | coverage upload                                      |

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
- **[`sa_rosenbaum()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_rosenbaum.md)
  requires `rbounds`; the other three sensitivity methods do not.** Keep
  it that way —
  [`sa_evalue()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_evalue.md)
  is deliberately dependency-free, and adding a hard dependency to the
  cheap methods removes the reason they exist.
- **Roxygen markdown is ENABLED** (`Roxygen: list(markdown = TRUE)`). ⚠️
  `hvtiRutilities` and `hvtiRtemplates` have no such field and need Rd
  markup instead.
- **`testthat` edition 3.** `VignetteBuilder` is **quarto**.

## Gotchas

- **`DESCRIPTION`’s `Date:` is stale** — 2026-04-02 against version
  0.1.1. Refresh it on the next version bump rather than carrying it
  further.
- The package is **0.x**: the API is not frozen, but the `ps_data`
  structure above is what every subclass and both base methods rely on,
  so changing it is a breaking change in practice.
- The sample-data generators are the only functions using the RNG, and
  they are the source of the outstanding NOTE. If you touch them, fix
  the import rather than adding a
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
