# hvtiRpropensity 0.2.2

* The README's Installation section now leads with `hvtiR::install()`, which
  installs the whole HVTI R family in one call, and then gives
  `pak::pak("ehrlinger/hvtiRpropensity")` for this package alone.

# hvtiRpropensity 0.2.1

0.2.1 fixes formula handling across the model functions and the findings of a
rollup review of the 0.2.0 series. Two changes move published numbers; check
them before upgrading a study:

* Pooled degrees of freedom are Barnard and Rubin's, as `mice::pool()` reports
  them, so pooled p-values and intervals change.
* Ordinal standard errors come from the expected information, as SAS reports
  them, so ordinal p-values and intervals change, by a fraction of a percent
  on the package's test data.

* `fit_logistic()` accepts a dot formula such as `y ~ .`. It stopped with
  "Required column(s) not found in `data`: ." before. The dot expands to every
  column except the outcome, `id_col` and `imputation_col`, as it already did in
  `bs_continuous()` and `fit_linear()`.
* `bs_count()` checks `outcome_col` the way `bs_continuous()` does: a value that
  is not the formula's response, or a response naming more than one column, is
  refused. A dot formula no longer fits the `id_col` and `imputation_col`
  columns as covariates.
* `fit_logistic()`, `fit_linear()`, `bs_continuous()` and `bs_count()` check
  only the predictors the model keeps. A variable the formula subtracts, as in
  `y ~ . - dropme`, is no longer checked for matching levels across
  imputations, and `fit_logistic()` no longer lists it in `$meta$predictors`.
* `fit_linear()`, `bs_continuous()` and `bs_count()` require the response to
  read exactly one column even when `outcome_col` is supplied.
  `I(count + other)` with `outcome_col = "count"` fitted the sum while
  recording `count` as the outcome; `cbind(count, other)` stopped with an
  unhelpful subscript error.

* Pooled estimates from `fit_logistic()`, `fit_linear()`, `bs_continuous()`
  and `bs_count()` use Barnard and Rubin's (1999) small-sample degrees of
  freedom, with complete-data degrees of freedom n - k, as `mice::pool()` does.
  They used Rubin's (1987) large-sample degrees of freedom, which can run into
  the millions when little is missing. The `df` column, and the p-values and
  intervals built on it, change for every pooled model; estimates and standard
  errors do not. A coefficient with no between-imputation variance, which had
  infinite `df` and a normal test, now gets a finite `df`. To reproduce the
  values in SAS, run `PROC MIANALYZE` with `EDF=` set to n - k. The help pages
  for all four functions explain the calculation under "Pooled degrees of
  freedom".
* Pooling over imputations stops with an error when the first imputation's
  model has no residual degrees of freedom (n - k of zero or less), as a
  saturated model does. A saturated pooled logistic model returned meaningless
  `df` and p-values before.

* The SAS oracle fixtures for `fit_logistic()` and `fit_linear()` now include
  a `PROC MIANALYZE EDF=` run, and the tests check the pooled degrees of
  freedom, p-values and intervals against it directly. No change in behaviour.

* `ps_forest()` accepts a dot formula such as `tavr ~ .`. It stopped with
  "Required column(s) not found in `data`: ." before. The dot expands to every
  column except the treatment and `id_col`, as it does in `ps_logistic()`.
  A column the formula subtracts, as in `. - x3`, is no longer fitted, checked
  for missing values or reported in the SMD table.
* `ps_logistic()`, `ps_ordinal()` and `ps_nominal()` leave a constant column
  the formula subtracts out of the SMD table.

* The SAS oracle fixtures for `fit_logistic()` and `bs_count()` are exported at
  full precision, with the `PROC LOGISTIC` fits run at `GCONV=1e-12`. The
  binary, nominal and pooled comparisons, and the `PROC GENMOD` Poisson
  comparison, are tested at tolerances near their measured agreement. No change
  in behaviour.

* `fit_linear()`, `bs_continuous()`, `bs_count()` and `ps_forest()` check every
  column the formula names, subtracted ones included, as `fit_logistic()` does.
  A subtracted column missing from `data` reached `lm()` or `glm()` as "object
  not found", or picked up a same-named object from the calling environment;
  in `ps_forest()` a misspelled one was ignored, so the column meant to be
  dropped was fitted.
* `bs_count()` refuses a transformed response such as `log(count + 1)`, which
  it fitted as counts with a warning. `fit_logistic()` refuses one too,
  matching its error message.
* A `.` that leaves no columns once the response, `id_col` and
  `imputation_col` are set aside stops with a clear error.

* Ordinal models in `fit_logistic()` and `ps_ordinal()` report standard
  errors from the expected (Fisher) information, as SAS `PROC LOGISTIC` does,
  instead of `MASS::polr()`'s observed information. The two differ for the
  cumulative logit; on the package's SAS oracle data the standard errors moved
  by up to 0.4%, and the p-values, intervals and pooled results with them.
  `vcov()` on the stored `polr()` model still returns the observed version.
* `polr()` runs at a tighter tolerance (`reltol = 1e-14`), so ordinal
  estimates move slightly (about 2e-5 on the oracle data) and now agree with
  SAS to about 1e-7.

# hvtiRpropensity 0.2.0

0.2.0 consolidates the work since 0.1.0 and the fixes from a rollup review of
that series (#79 to #93). Before upgrading a study, check these changes, each
described in full below:

* The four `sa_*()` functions return `ps_data` objects; results move into
  `$meta` and `$tables`.
* `sa_evalue()` needs `rare = TRUE` or `FALSE` for an odds ratio or a hazard
  ratio, and its earlier conversions were the wrong way round.
* `ps_rmst()` on a `ps_match()` object analyses the matched pairs and
  resamples them by pair; score weighting of a matched object is refused.
* Weighted SMDs normalise the weights within each arm, so they no longer
  depend on `stabilise` and can differ from SAS `%stddiff` when the weights do
  not average 1.
* `ps_forest()` writes matching weights, not overlap weights, to `mt_wt`.
* `ps_ordinal()` needs an ordered factor or `treatment_levels`.
* `ps_weight()` stabilises ATE weights only.
* Rows with a missing treatment or score are excluded and counted, not
  crashed on or silently kept.

* `DESCRIPTION` now declares the Quarto command line tool in
  `SystemRequirements`. The vignettes have always needed it to build; the
  field makes that visible to installers and to `R CMD check`.

* Articles on the pkgdown site put the table of contents on the left and use
  the full width of the window, through `pkgdown/extra.css`. The installed
  vignettes are unchanged: the Quarto vignette engine renders them in its own
  minimal format, which has no sidebar layout.

* `ps_match()` gains `scale = c("probability", "logit")` and `caliper_sd`, so
  pairs can be matched on the logit of the propensity score with a caliper in
  pooled standard-deviation units, `sqrt((var_t + var_c) / 2)`. `scale = "logit", caliper_sd = 0.2` is the
  Austin (2011) recommendation. The default stays `"probability"`, which keeps
  the \[0, 1\] check on `score_col`. `$meta` now records `scale` and
  `caliper_sd`, and `$meta$caliper` is the absolute caliper applied.

* The SAS templates in `inst/sas-templates/` no longer carry the study they
  were taken from: each study path is now `/studies/<study>`, and the study
  title and cohort line read `<study description>` and
  `CCF, 19xx to 20xx, n=xxxx`, as the `%let STUDYDESC` and `%let STUDYPOP`
  lines already did. The analysis code is unchanged.

* `fit_linear()` fits an ordinary least-squares model to one dataset or to
  stacked imputations, with Rubin-pooled inference verified against
  `PROC MIANALYZE`. It returns the same data, models, meta and tables shape
  as `fit_logistic()`, and its `$tables$fit_quality` gives each fit's
  residual standard error and R-squared. A single fit, here and in
  `bs_continuous()`'s new estimates, is inferred on the t distribution with
  residual degrees of freedom, as `PROC REG` reports it. A dot formula such
  as `ef ~ .` expands over the data's columns excluding the identifier and
  imputation columns, here and in `bs_continuous()`, so the patient ID is
  never fitted as a covariate.
* `bs_continuous()` now carries fit status, coefficient estimates, covariance
  and per-imputation estimates, as `bs_count()` does. Its balancing score is
  unchanged. A transformed exposure, such as `log(hct_nadr) ~ ...`, now works.
* `bs_continuous()` now stops, naming the imputation, when a term is aliased
  (for example a covariate that duplicates another), where it used to return
  a score from a rank-deficient fit.

* The pooled-covariance check against `PROC MIANALYZE` in
  `test_lm_sas_oracle.R` now compares like with like. SAS's `TCov` is the
  multivariate total covariance, `(1 + rbar) * Ubar`, not the elementwise Rubin
  total that `fit_logistic()` reports, so the test checks the diagonal against
  the SAS standard errors and rebuilds `TCov` from the per-imputation fits, at
  `1e-5` rather than `1e-4`. No package code changed.

* The help pages for `ps_match()`, `ps_weight()` and the package now name
  `hvtiPlotR::hv_mirror_hist()` and `hvtiPlotR::hv_balance()` as code rather
  than linking to them. 'hvtiPlotR' is not a declared dependency, so
  `R CMD check` could not resolve the links on a machine without it and
  reported "Unknown package 'hvtiPlotR' in Rd xrefs".

* **Behaviour change:** `ps_forest()` now writes Li-Greene matching weights
  (ATM; Li & Greene, 2013) to `weight_col` (default `mt_wt`), `min(p, 1 - p) / p`
  for treated patients and `min(p, 1 - p) / (1 - p)` for controls, exactly as
  `ps_logistic()` does. It previously wrote overlap (ATO) weights, `1 - p` and
  `p`, under the same name, so swapping one scoring function for the other
  silently changed the estimand. Anyone who read `mt_wt` from a `ps_forest()`
  result, for weighted balance, `ps_mw_var()` or
  `hvtiPlotR::hv_mirror_hist(weight_col = "mt_wt")`, gets different weights and
  should rerun. Both help pages now state the estimand (#83).

* **Behaviour change:** `ps_ordinal()` now stops unless the treatment is an
  ordered factor or `treatment_levels` gives the level order. It used to
  order a character or unordered-factor treatment alphabetically, so
  `"Low"`/`"Medium"`/`"High"` was fitted as High < Low < Medium, with no
  warning and coefficients of the wrong sign. The error says how to declare
  the order (#87).
* `ps_logistic()` with `treated_level` alone now accepts any treatment with
  two values, such as a `"TAVR"`/`"SAVR"` factor, as `ps_forest()` already
  did, and gives the same scores as the 0/1 coding. A treatment with more
  than two values still stops (#87).
* `ps_logistic()`, `ps_ordinal()`, `ps_nominal()`, `fit_logistic()` and
  `fit_linear()` check `id_col` only when `imputation_col` is set, as their
  documentation says. Data without an `id` column no longer stops a
  single-dataset fit (#87).

* `ps_stddiff()` now normalises the weights to mean 1 within each group before
  computing a weighted Gaussian or rank difference, so multiplying a group's
  weights by a constant no longer changes it (#82). As a result
  `ps_weight()$tables$smd_weighted` is the same with `stabilise = TRUE` and
  `stabilise = FALSE`, and `ps_stddiff_perm()` with `weight_col` inherits the
  fix.
* **Weighted SMDs change** for Gaussian and `nong_ord` variables whenever the
  weights do not already average 1 in each group, as matching weights and
  unstabilised IPTW weights do not. They then no longer equal the `%stddiff`
  macro's output, which depends on the scale of the weights; they still equal
  it when the weights average 1 in each group. Binary and categorical
  differences, and every unweighted difference, are unchanged.

* `sa_evalue()` now applies the odds-ratio and hazard-ratio conversions the
  right way round (VanderWeele & Ding 2017). For a rare outcome the OR or HR
  is used as the risk ratio directly; for a common outcome the OR becomes
  `sqrt(OR)` and the HR becomes `(1 - 0.5^sqrt(HR)) / (1 - 0.5^sqrt(1/HR))`.
  Earlier versions always took the square root, which understated the E-value
  for a rare outcome (an OR of 2.5 gave 2.54 instead of 4.44).
* Behaviour change: `sa_evalue()` gains a `rare` argument, and a call with
  `type = "OR"` or `"HR"` now aborts unless it sets `rare = TRUE` or
  `rare = FALSE`. Rarity is the outcome's prevalence in the study population
  (under about 15% at the end of follow-up) and is not inferred from `p0`,
  which is used only for `type = "RD"`. The returned list records the choice
  as `$rare`. The function also aborts unless `ci_lo <= estimate <= ci_hi`.

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

* `NEWS.md` corrections: the 0.1.0 entry now lists all sixteen initial exports, and the 0.1.8 `ps_forest()`, 0.1.9 `group_counts` and 0.1.11 seed entries no longer claim more than the code did. Each correction is dated in place (#92).

* The default covariates of the balance tables in `ps_weight()`, `ps_match()`, `ps_logistic()`, `ps_nominal()` and
  `ps_ordinal()` now leave out the id column, the treatment and score columns, the columns this package appends
  (`prob_t`, `logit_t`, `mt_wt`, `iptw`, `match`, `pair_id`, `quintile`, `decile`, `bs`, `cluster`, and the
  `prob_<level>` scores `ps_nominal()` and `ps_ordinal()` generate) and constant columns. A predictor named in the
  model formula is always kept. Before, `id` and `match` appeared in the `ps_weight()` table of `sample_ps_data()`,
  and `prob_t` in a `ps_match()` table matched on `logit_t`. A `covariates` vector you supply is used unchanged (#85).
* `ps_weight()` gains `id_col = "id"`, as `ps_match()` has, so the identifier stays out of its balance tables.

* `fit_logistic()` checks that 'MASS' (ordinal) or 'nnet' (nominal) is
  installed before it fits, and the tests that use either package skip when
  it is absent (#89).
* The test suite and the `ps_forest()` example cap 'randomForestSRC' at two
  cores (#89).

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

* `ps_weight()` now stabilises weights by the marginal treatment probability
  for `estimand = "ATE"` only. For `"ATT"` and `"ATC"`, `stabilise` has no
  effect, so treated weights under the ATT and control weights under the ATC
  are exactly 1, as documented; before, the default `stabilise = TRUE` scaled
  them by the arm's marginal probability. `$meta$stabilised` records whether
  stabilisation was applied, and `sa_trim_sweep()` follows the same rule.
  ATE weights are unchanged (#88).
* `ps_match()` accepts only the number 1 for `ratio`. `ratio = 1.5` and
  `ratio = "1"` were silently accepted and now stop with an error naming the
  value (#88).
* `ps_support()` stops when a `flags` element is named `"trim"` or `"common"`,
  which the built-in rule silently overwrote, and when `trim` is not two
  increasing probabilities in \[0, 1\] (#88).
* An `NA` in `caliper` or `caliper_sd` (`ps_match()`) or `trim`
  (`ps_weight()`, `ps_support()`) now stops with an error naming the argument
  instead of "missing value where TRUE/FALSE needed" (#88).

- Rows with a missing treatment or a missing score are now excluded and counted rather than crashing or misleading (#84). `ps_match()` leaves them out of the treated and control pools, `ps_weight()` gives them a missing weight and leaves them out of the stabilising proportion and the balance tables, and `ps_support()` computes the common-support interval on the remaining rows and gives them missing flags. Each states the count in one warning and records it in `$meta$n_excluded`. `ps_support()` also works with a single rule (`trim = NULL` or `common = FALSE`).
- `ps_logistic()` counts a patient with a missing treatment out of `group_counts`, the strata and the balance tables, keeps their predicted score, says so in one accurate warning, and records the count in `$meta$n_missing_treatment`.
- `fit_status` now reports the rows a model dropped for every family. `ps_nominal()` and `fit_logistic(family = "nominal")` reported none, because `nnet` has no `nobs()` method.
- `bs_continuous()` and `bs_count()` give a missing balancing score missing strata instead of ranking it into the top stratum, and need `id_col` only when `imputation_col` is set.

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

* Tests now assert what their names claim: the caliper, `trim = NULL`,
  estimand and stabilised-weight tests are pinned to hand-computed values, the
  imputation test checks `prob_t` against the mean of the per-imputation fits,
  the AUC test uses `wilcox.test()` as an independent oracle, and the
  `sa_trim_sweep()` tests reach the estimand and machine-epsilon paths they
  name. No package code changed (#93).

- The help pages now list everything each function returns (#90). `fit_logistic()`, `validate_logistic()`, `ps_logistic()`, `ps_ordinal()`, `ps_nominal()` and `bs_count()` name every `$data` column, `$meta` field, `$tables` entry and their `$models`; `ps_match()` and `ps_weight()` add `treated_level`, and `ps_match()` says that `n_unmatched` counts unmatched treated patients only. `fit_logistic()` documents the choices for `family` and the role of `id_col`.
- `ps_logistic()`, `ps_ordinal()`, `ps_nominal()` and `fit_logistic()` describe their handling of multiple imputation correctly: per-patient probabilities are averaged across imputations on the probability scale, the "across" approach of Leyrat et al. (2019), and not Rubin's rules, which pool only the coefficients. `$data` and the balance table use the lowest-numbered imputation's covariates. `ps_logistic()` gives the "within" recipe that Leyrat et al. recommend.
- The help pages name the estimand: `ps_match()` the ATT, narrowed by a caliper to the treated patients who could be matched; `ps_weight()` the ATE, ATT or ATC; `ps_mw_var()` with the matching weight the ATM; and `ps_support()` says that trimming changes it.
- Corrections: `ps_logistic()`'s matching weight `mt_wt` feeds neither `ps_match()` nor `ps_weight()`; `ps_stddiff()` no longer calls `ps_weight()`'s weights matching weights; `ps_forest()` lists its own `$meta`; `ps_mw_var()` cites Li and Greene (2013) with its DOI; `sample_ps_data_count(n_imputations > 1)` says that it draws independent patients per imputation; and a stray backtick no longer breaks the `treatment_col` entry of `sa_overlap()` and `sa_trim_sweep()`.

* The vignette's `hvtiPlotR` examples now run. `hv_mirror_hist()` takes column
  names as strings (`score_col`, `group_col`, `match_col`, `weight_col`) and
  returns a data object for `plot()`. `hv_balance()` does not read `$tables`
  directly: it needs one long data frame with `variable`, `group` and
  `std_diff` in percent, and the vignette now shows the reshape. The package
  help and the README no longer say the tables feed it as they are.
  'hvtiPlotR' is now in `Suggests` (#91).
* The vignette now covers `ps_forest()`, `ps_stddiff()`, `ps_stddiff_perm()`,
  `ps_mw_var()`, `ps_rmst()` on matched pairs, `ps_support()`, matching on the
  logit with `caliper_sd = 0.2`, and the `rare` argument `sa_evalue()` needs
  for odds and hazard ratios. Every `ps_match()` call in it gives a seed. It
  no longer says all four `sa_*()` functions measure unmeasured confounding:
  `sa_overlap()` and `sa_trim_sweep()` check positivity and overlap.
* The README quick start runs: it weighted a data frame whose score column it
  had just removed. Its status block now says which parts SAS oracles check
  (since 0.1.10), and its function tables list all 25 exports.
* The package help page lists all 25 exports, says which return `ps_data`
  objects (all but the four data generators and `is_ps_data()`), and points
  to the SAS programs shipped in
  `system.file("sas-templates", package = "hvtiRpropensity")`.

# hvtiRpropensity 0.1.11

* A seed now makes every random function reproducible whatever the session
  drew before it. `ps_forest()`, `ps_match()`, `ps_rmst()`,
  `ps_stddiff_perm()` and `ps_mw_var()` set the stream to `abs(seed)`
  immediately before their draws and restore the caller's stream afterwards.
  `ps_forest()` now does this as well as passing `seed = -abs(seed)` to
  `randomForestSRC::rfsrc()`; before, it set no R seed. The default stays
  `seed = NULL`, which draws from the caller's stream as before. A positive
  seed gives the same result as in 0.1.10, except in `ps_forest()`, whose
  0.1.10 result depended on the session's earlier draws; a negative one now
  gives the same result as its absolute value. (Corrected 2026-10-13: the
  original entry did not exempt `ps_forest()`.)
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
  `treated` column. (Corrected 2026-10-13: the original entry said this made
  it read the same as `ps_ordinal()` and `ps_nominal()`, whose
  `group_counts` have `group` and `n` only.)
  Without `treatment_levels` the groups are `0`/`1` (or `FALSE`/`TRUE`) (#47).

# hvtiRpropensity 0.1.8

* New `ps_forest()`: out-of-bag random-forest propensity score in a
  `ps_data` object (needs the suggested `randomForestSRC`). (Corrected
  2026-10-13: the original entry said it had the same contract as
  `ps_logistic()`, but its `mt_wt` held overlap weights rather than matching
  weights, and its `$meta` and `$tables` lack the model-bundle fields.)
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
  rest of the `hvtiR*` family. The two fields had disagreed: `Authors@R` named
  the `ehrlinj@ccf.org` address while `Maintainer:` named the gmail one, which
  `R CMD check --as-cran` reported as a note. `Maintainer:` is now derived from
  the `cre` role at build time, so the two cannot drift apart again.

# hvtiRpropensity 0.1.2

* Qualified the `rnorm()` calls in the `sample_ps_data*()` generators as
  `stats::rnorm()`, matching every other statistics call in the package, and
  declared `stats` under `Imports`. This clears the `R CMD check` note about
  an undefined global function. Generated data is unchanged: the seeds and
  the RNG draw order are identical.

# hvtiRpropensity 0.1.1

* Renamed from `hvtiPropensityScores` into the `hvtiR*` package family
  (`hvtiRtemplates`, `hvtiRutilities`, `hvtiRlifetables`, `hvtiRtables`).
  No user-facing function names or behaviour changed.

# hvtiRpropensity 0.1.0

* Initial development scaffold.
* (Added 2026-10-13: the original 0.1.0 entry listed four exports. These were
  exported too: `ps_logistic()`, `ps_ordinal()`, `ps_nominal()`,
  `bs_continuous()`, `bs_count()`, `sa_rosenbaum()`, `sa_evalue()`,
  `sa_overlap()`, `sa_trim_sweep()`, `sample_ps_data_ordinal()`,
  `sample_ps_data_nominal()` and `sample_ps_data_count()`.)
* Added `sample_ps_data()`: reproducible synthetic cardiac-surgery dataset.
* Added `ps_match()`: greedy nearest-neighbour 1:1 propensity score matching
  without replacement, with optional caliper.
* Added `ps_weight()`: IPTW weighting supporting ATE, ATT, and ATC estimands,
  with optional stabilisation and weight winsorisation.
* Added `is_ps_data()` predicate and `print` / `summary` S3 methods for the
  common `ps_data` base class.
