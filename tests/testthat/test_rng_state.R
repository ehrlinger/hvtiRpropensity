###############################################################################
## test_rng_state.R
##
## A function that takes a seed must leave the caller's random number stream
## as it found it: a set.seed() inside a function silently reseeds everything
## the caller draws afterwards.
###############################################################################

rng_state <- function() get(".Random.seed", envir = globalenv())

expect_rng_untouched <- function(call) {
  set.seed(2026)
  before <- rng_state()
  force(call)
  testthat::expect_identical(rng_state(), before)
}

test_that("the sample-data generators restore the caller's stream", {
  expect_rng_untouched(sample_ps_data(n = 50, seed = 1))
  expect_rng_untouched(sample_ps_data_ordinal(n = 50, seed = 1))
  expect_rng_untouched(sample_ps_data_nominal(n = 50, seed = 1))
  expect_rng_untouched(sample_ps_data_count(n = 50, seed = 1))
  expect_rng_untouched(sample_ps_data_count(n = 50, seed = 1, n_imputations = 2))
})

test_that("a seed still makes the sample data reproducible", {
  expect_identical(sample_ps_data(n = 50, seed = 3), sample_ps_data(n = 50, seed = 3))
  expect_identical(sample_ps_data_count(n = 50, seed = 3, n_imputations = 2),
                   sample_ps_data_count(n = 50, seed = 3, n_imputations = 2))
})

test_that("ps_match() with a seed restores the caller's stream", {
  dta <- sample_ps_data(n = 100, seed = 42)
  expect_rng_untouched(ps_match(dta, seed = 7))
  expect_identical(ps_match(dta, seed = 7)$data, ps_match(dta, seed = 7)$data)
})

test_that("ps_rmst() with a seed restores the caller's stream", {
  d <- sample_ps_data(n = 120, seed = 9)[, c("id", "tavr", "age", "ef")]
  d$t <- rep_len(c(1, 2.5, 4, 5), nrow(d))
  d$e <- rep_len(c(1L, 0L, 1L), nrow(d))
  fit <- ps_logistic(tavr ~ age + ef, d)
  expect_rng_untouched(ps_rmst(fit, "t", "e", tau = 4, weights = "unweighted", n_boot = 5, seed = 7))
})

test_that("ps_rmst() without a seed draws from the caller's stream", {
  d <- sample_ps_data(n = 120, seed = 9)[, c("id", "tavr", "age", "ef")]
  d$t <- rep_len(c(1, 2.5, 4, 5), nrow(d))
  d$e <- rep_len(c(1L, 0L, 1L), nrow(d))
  fit <- ps_logistic(tavr ~ age + ef, d)
  run <- function() ps_rmst(fit, "t", "e", tau = 4, weights = "unweighted", n_boot = 20)$tables$estimates
  set.seed(11)
  first <- run()
  second <- run()
  expect_false(identical(first$lo_days, second$lo_days))
  set.seed(11)
  expect_identical(run(), first)
})

test_that("the sample-data generators accept seed = NULL and draw fresh data", {
  set.seed(12)
  expect_false(identical(sample_ps_data(n = 50, seed = NULL), sample_ps_data(n = 50, seed = NULL)))
  count <- sample_ps_data_count(n = 50, seed = NULL, n_imputations = 2)
  expect_equal(nrow(count), 100L)
  expect_false(identical(count$age[1:50], count$age[51:100]))
})
