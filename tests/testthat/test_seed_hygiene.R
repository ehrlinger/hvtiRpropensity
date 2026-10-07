###############################################################################
## test_seed_hygiene.R
##
## With a seed, a function that draws random numbers must give the same answer
## whatever the session drew before it. Two parts:
##   1. Each random function, run with one seed after two different set.seed()
##      states, returns identical results and records the seed in $meta.
##   2. An enforcement scan of every function in the namespace fails on a
##      random draw that no seeding call precedes, with a self-test showing the
##      scan catches one.
###############################################################################

# Same seed, different session histories: the result must not move.
expect_seed_stable <- function(run) {
  set.seed(1)
  first <- run()
  set.seed(2026)
  stats::runif(17)
  second <- run()
  testthat::expect_identical(first, second)
}

surv_data <- function() {
  d <- sample_ps_data(n = 120, seed = 9)[, c("id", "tavr", "age", "ef")]
  d$t <- rep_len(c(1, 2.5, 4, 5), nrow(d))
  d$e <- rep_len(c(1L, 0L, 1L), nrow(d))
  d
}

mw_data <- function() {
  d <- sample_ps_data(n = 120, seed = 42)
  d$mt_wt <- pmin(d$prob_t, 1 - d$prob_t) / ifelse(d$tavr == 1, d$prob_t, 1 - d$prob_t)
  d
}

test_that("ps_forest() with a seed ignores the session's history and records the seed", {
  skip_if_not_installed("randomForestSRC")
  d <- sample_ps_data(n = 120, seed = 42)[, c("id", "tavr", "age", "ef")]
  fit <- function(seed) ps_forest(tavr ~ age + ef, d, ntree = 40, seed = seed)
  expect_seed_stable(function() fit(5)$data)
  expect_identical(fit(-5)$data, fit(5)$data)
  expect_false(identical(fit(5)$data$prob_t, fit(6)$data$prob_t))
  expect_identical(fit(-5)$meta$seed, 5L)
  expect_identical(ps_forest(tavr ~ age + ef, d, ntree = 10)$meta$seed, NA_integer_)
  set.seed(2026)
  before <- .Random.seed
  fit(5)
  expect_identical(.Random.seed, before)
})

test_that("ps_match() with a seed ignores the session's history and records the seed", {
  d <- sample_ps_data(n = 150, seed = 42)
  expect_seed_stable(function() ps_match(d, seed = 7)$data)
  expect_identical(ps_match(d, seed = -7)$data, ps_match(d, seed = 7)$data)
  expect_identical(ps_match(d, seed = 7)$meta$seed, 7L)
  expect_identical(ps_match(d)$meta$seed, NA_integer_)
})

test_that("ps_rmst() with a seed ignores the session's history and records the seed", {
  skip_if_not_installed("survival")
  fit <- ps_logistic(tavr ~ age + ef, surv_data())
  run <- function(seed) ps_rmst(fit, "t", "e", tau = 4, weights = "unweighted", n_boot = 10, seed = seed)
  expect_seed_stable(function() run(3)$tables)
  expect_identical(run(-3)$tables, run(3)$tables)
  expect_identical(run(3)$meta$seed, 3L)
})

test_that("ps_stddiff_perm() with a seed ignores the session's history and records the seed", {
  d <- sample_ps_data(n = 120, seed = 42)
  run <- function(seed) ps_stddiff_perm(d, gaussian = "age", n_perm = 20, seed = seed)
  expect_seed_stable(function() run(4)$tables)
  expect_identical(run(-4)$tables, run(4)$tables)
  expect_identical(run(4)$meta$seed, 4L)
  expect_identical(ps_stddiff_perm(d, gaussian = "age", n_perm = 5)$meta$seed, NA_integer_)
})

test_that("ps_mw_var() with a seed ignores the session's history and records the seed", {
  d <- mw_data()
  run <- function(seed) ps_mw_var(d, outcomes = "age", weight_col = "mt_wt", n_rep = 20, seed = seed)
  expect_seed_stable(function() run(8)$tables)
  expect_identical(run(-8)$tables, run(8)$tables)
  expect_identical(run(8)$meta$seed, 8L)
})

test_that("a seed that is not one whole number is rejected", {
  d <- sample_ps_data(n = 60, seed = 42)
  expect_error(ps_match(d, seed = 1.5), "whole number")
  expect_error(ps_match(d, seed = c(1, 2)), "whole number")
  expect_error(ps_match(d, seed = "1"), "whole number")
  expect_error(ps_stddiff_perm(d, gaussian = "age", n_perm = 5, seed = NA), "whole number")
})


# ---------------------------------------------------------------------------
# Enforcement: no random draw without a seeding call ahead of it
# ---------------------------------------------------------------------------

# Calls that draw from R's stream, or (randomForestSRC, varPro) seed
# themselves from it when given no seed.
random_calls <- c("sample", "sample.int", "runif", "rnorm", "rbinom", "rpois", "rexp", "rgamma", "rbeta",
                  "rlnorm", "rweibull", "rlogis", "rcauchy", "rchisq", "rgeom", "rhyper", "rnbinom",
                  "rmultinom", "rsignrank", "rwilcox", "rfsrc", "rfsrc.fast", "varpro")
seeding_calls <- c("set.seed", "local_seed", "with_seed")

# Every symbol in a function body, depth first and left to right, which is
# source order. `pkg::fn` contributes `fn`, and a function passed by name
# (do.call(randomForestSRC::rfsrc, args)) counts as a call to it.
body_symbols <- function(e) {
  if (is.name(e)) return(as.character(e))
  if (is.call(e) || is.pairlist(e) || is.expression(e)) {
    return(unlist(lapply(as.list(e), function(x) if (missing(x)) character() else body_symbols(x))))
  }
  character()
}

# Names of functions, among `fns`, with a random draw that no seeding call
# precedes, as "name: first unseeded draw".
seed_violations <- function(fns) {
  bad <- vapply(names(fns), function(nm) {
    syms <- body_symbols(body(fns[[nm]]))
    draws <- which(syms %in% random_calls)
    if (!length(draws)) return(NA_character_)
    seeds <- which(syms %in% seeding_calls)
    if (length(seeds) && min(seeds) < min(draws)) return(NA_character_)
    paste0(nm, ": ", syms[min(draws)])
  }, character(1L), USE.NAMES = FALSE)
  bad[!is.na(bad)]
}

namespace_functions <- function() {
  ns <- asNamespace("hvtiRpropensity")
  nms <- ls(ns, all.names = TRUE)
  fns <- mget(nms, envir = ns)
  fns[vapply(fns, is.function, logical(1L))]
}

test_that("the scan catches an unseeded draw and passes a seeded one", {
  fns <- list(
    unseeded = function(x) sample(x),
    late     = function(x, seed = NULL) {
      y <- sample(x)
      if (!is.null(seed)) withr::local_seed(seed)
      y
    },
    forest   = function(d) do.call(randomForestSRC::rfsrc, list(data = d)),
    seeded   = function(x, seed = NULL) {
      if (!is.null(seed)) withr::local_seed(seed)
      stats::rnorm(x)
    },
    no_draw  = function(x) mean(x)
  )
  expect_identical(seed_violations(fns), c("unseeded: sample", "late: sample", "forest: rfsrc"))
})

test_that("every random draw in the package follows a seeding call", {
  fns <- namespace_functions()
  # Guard against a vacuous pass: the scan must see the known random functions.
  drawing <- names(fns)[vapply(fns, function(f) any(body_symbols(body(f)) %in% random_calls), logical(1L))]
  expect_true(all(c("ps_forest", "ps_match", "ps_rmst", "ps_stddiff_perm", "ps_mw_var",
                    "sample_ps_data") %in% drawing))
  expect_identical(seed_violations(fns), character())
})
