###############################################################################
## test_conditions.R -- the hvtiR family error, deprecation and provenance
## conventions (rubric adopted 2026-11-16).
###############################################################################

dta <- sample_ps_data(n = 120, seed = 4)

reset_deprecations <- function(env = parent.frame()) {
  seen <- hvtiRpropensity:::.deprecation_seen
  saved <- mget(ls(seen), envir = seen)
  rm(list = ls(seen), envir = seen)
  withr::defer({
    rm(list = ls(seen), envir = seen)
    list2env(saved, envir = seen)
  }, envir = env)
}

# ---- Errors -----------------------------------------------------------------

test_that("no rlang::abort() or rlang::warn() call is left in R/", {
  src <- list.files(test_path("..", "..", "R"), pattern = "[.]R$", full.names = TRUE)
  skip_if(length(src) == 0L, "package sources not available")
  hits <- unlist(lapply(src, function(f) grep("rlang::(abort|warn)\\(", readLines(f), value = TRUE)))
  expect_length(hits, 0L)
})

test_that("errors name the exported function and carry no call", {
  e <- expect_error(ps_weight(dta, treatment_col = "tavr", trim = 2))
  expect_match(conditionMessage(e), "^ps_weight\\(\\): ")
  expect_null(conditionCall(e))
})

test_that("an internal helper's error names the exported function the user called", {
  e <- expect_error(ps_match(dta, treatment_col = "nope"), class = "hvtiRpropensity_missing_column")
  expect_match(conditionMessage(e), "^ps_match\\(\\): Required column")
  e <- expect_error(hvtiRpropensity::ps_weight(dta, treatment_col = "nope"), class = "hvtiRpropensity_missing_column")
  expect_match(conditionMessage(e), "^ps_weight\\(\\): ")
})

test_that("conditions a caller may catch are classed", {
  one_arm <- transform(dta, tavr = 1L)
  expect_error(ps_weight(one_arm, treatment_col = "tavr"), class = "hvtiRpropensity_one_arm")
  w <- ps_weight(dta, treatment_col = "tavr", estimand = "ATT")
  expect_error(sa_trim_sweep(w, estimand = "ATE"), class = "hvtiRpropensity_meta_conflict")
})

# ---- Deprecations ------------------------------------------------------------

test_that("an omitted treatment_col warns once per session and means \"tavr\"", {
  reset_deprecations()
  expect_warning(a <- ps_weight(dta), class = "hvtiRpropensity_deprecated")
  expect_no_warning(b <- ps_weight(dta))
  expect_identical(a$meta$treatment_col, "tavr")
  expect_identical(a$data, b$data)
  expect_warning(sa_overlap(dta), "treatment_col")
  # A ps_data input records its treatment column, so nothing is deprecated.
  expect_no_warning(sa_overlap(ps_weight(dta, treatment_col = "tavr")))
  expect_warning(m <- ps_match(dta, seed = 1), "treatment_col")
  expect_identical(m$meta$treatment_col, "tavr")
  expect_warning(ps_stddiff(dta, gaussian = "age"), "treatment_col")
  expect_warning(ps_stddiff_perm(dta, gaussian = "age", n_perm = 3, seed = 1), "treatment_col")
  expect_warning(ps_mw_var(dta, outcomes = "age", weight_col = "prob_t", n_boot = 3, seed = 1), "treatment_col")
  expect_warning(sa_trim_sweep(dta), "treatment_col")
})

test_that("an error raised through lapply() still names the function", {
  e <- expect_error(lapply(list(dta), ps_match, treatment_col = "nope"), class = "hvtiRpropensity_missing_column")
  expect_match(conditionMessage(e), "^ps_match\\(\\): ")
})

test_that("a ps_data object without a recorded treatment column is refused by name", {
  obj <- ps_weight(dta, treatment_col = "tavr")
  obj$meta$treatment_col <- NULL
  expect_error(sa_overlap(obj), "records no `treatment_col`")
})

test_that("renamed arguments keep working for one release", {
  reset_deprecations()
  expect_warning(old <- ps_mw_var(dta, treatment_col = "tavr", outcomes = "age", weight_col = "prob_t",
                                  n_rep = 5, seed = 1), "`n_rep` is renamed `n_boot`")
  new <- ps_mw_var(dta, treatment_col = "tavr", outcomes = "age", weight_col = "prob_t", n_boot = 5, seed = 1)
  expect_identical(old$tables, new$tables)
  expect_error(ps_mw_var(dta, treatment_col = "tavr", outcomes = "age", weight_col = "prob_t", n_rep = 5,
                         n_boot = 5), "give `n_boot` only")

  expect_warning(old <- sample_ps_data_count(n = 30, seed = 1, n_imputations = 2), "renamed `m`")
  expect_identical(old, sample_ps_data_count(n = 30, seed = 1, m = 2))

  w <- ps_weight(dta, treatment_col = "tavr")
  expect_warning(old <- sa_overlap(x = w), "`x` is renamed `data`")
  expect_identical(old$tables, sa_overlap(w)$tables)
})

test_that("ps_rmst() takes upper-case estimands and deprecates lower case", {
  skip_if_not_installed("survival")
  reset_deprecations()
  d <- dta
  d$t <- rep(c(1, 2, 3, 4, 5), length.out = nrow(d))
  d$e <- rep(c(1L, 0L), length.out = nrow(d))
  fit <- ps_logistic(tavr ~ age + ef, d[setdiff(names(d), "prob_t")])
  up <- ps_rmst(fit, "t", "e", tau = 4, weights = "ATO", n_boot = 2, seed = 1)
  expect_identical(up$tables$estimates$weighting, "ATO")
  expect_identical(up$tables$estimates$estimator, "all/ATO")
  expect_true("w_ato" %in% names(up$data))
  expect_warning(low <- ps_rmst(fit, "t", "e", tau = 4, weights = "ato", n_boot = 2, seed = 1),
                 class = "hvtiRpropensity_deprecated")
  expect_identical(low$tables, up$tables)
  both <- ps_rmst(fit, "t", "e", tau = 4, weights = c("ATO", "ato"), n_boot = 2, seed = 1)
  expect_identical(nrow(both$tables$estimates), 1L)
})

# ---- Provenance -----------------------------------------------------------

conforms <- function(p) {
  fields <- c("package", "version", "sha", "method", "seed", "rng_kind", "time")
  testthat::expect_true(all(fields %in% names(p)))
  testthat::expect_true(is.character(p$package) && length(p$package) == 1L)
  testthat::expect_true(is.character(p$version) && length(p$version) == 1L)
  testthat::expect_true(is.character(p$sha) && length(p$sha) == 1L)
  testthat::expect_true(is.character(p$method) && length(p$method) >= 1L)
  testthat::expect_true(is.integer(p$seed) && length(p$seed) == 1L)
  testthat::expect_true(is.character(p$rng_kind) && length(p$rng_kind) == 3L)
  testthat::expect_true(inherits(p$time, "POSIXct") && length(p$time) == 1L)
}

test_that("the provenance record has the family's seven fields", {
  conforms(.provenance("x"))
  conforms(.provenance("x", 5L))
  expect_identical(.provenance("x", 5L)$seed, 5L)
  expect_identical(.provenance("x")$seed, NA_integer_)
})

test_that("every fit and random function records provenance", {
  fit <- ps_logistic(tavr ~ age + ef, dta[setdiff(names(dta), "prob_t")])
  m <- ps_match(dta, treatment_col = "tavr", seed = 3)
  objs <- list(
    fit,
    fit_logistic(tavr ~ age + ef, dta, outcome_levels = c(0, 1), event_level = 1),
    fit_linear(ef ~ age + female, dta),
    m,
    ps_stddiff_perm(dta, treatment_col = "tavr", gaussian = "age", n_perm = 5, seed = 1),
    ps_mw_var(dta, treatment_col = "tavr", outcomes = "age", weight_col = "prob_t", n_boot = 5, seed = 1)
  )
  for (o in objs) conforms(o$meta$provenance)
  expect_identical(m$meta$provenance$seed, 3L)
  expect_identical(m$meta$provenance$package, "hvtiRpropensity")
})
