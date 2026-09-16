###############################################################################
## test_ps_stddiff_perm.R
##
## Tests for ps_stddiff_perm(), a port of stddiffci.sas (Artis 2020): the
## observed standardized difference beside percentiles of its permutation
## distribution.
###############################################################################

dta <- sample_ps_data(n = 80, seed = 7)

# ---- Percentile definition --------------------------------------------------

test_that("percentiles use PROC UNIVARIATE's default, PCTLDEF=5 (R type 2)", {
  # 1..10: np = .25 -> x[1]; 1.6 -> x[2]; 5 exactly -> mean(x[5], x[6]); 8.4 -> x[9]; 9.75 -> x[10].
  # type 7 would give 2.44 at .16, and type 4 (PCTLDEF=1) would give 1.6.
  expect_equal(hvtiRpropensity:::.perm_percentiles(c(4, 1, 10, 7, 2, 9, 3, 8, 5, 6)), c(1, 2, 5.5, 9, 10))
})

test_that("percentiles skip missing permutations, as PROC UNIVARIATE does", {
  expect_equal(hvtiRpropensity:::.perm_percentiles(c(1:10, NA, NA)), c(1, 2, 5.5, 9, 10))
})

# ---- Class and structure ----------------------------------------------------

test_that("ps_stddiff_perm() returns a ps_stddiff_perm / ps_data object", {
  obj <- ps_stddiff_perm(dta, gaussian = "age", binary = "female", n_perm = 20, seed = 1)
  expect_s3_class(obj, "ps_stddiff_perm")
  expect_s3_class(obj, "ps_data")
  expect_identical(obj$data, dta)
  expect_named(obj$tables$stddiff_perm,
               c("variable", "label", "type", "observed", "p2_5", "p16", "p50", "p84", "p97_5"))
  expect_identical(obj$meta$n_perm, 20L)
})

test_that("the observed column is ps_stddiff() on the unpermuted data", {
  obj <- ps_stddiff_perm(dta, gaussian = c("age", "ef"), binary = "female", n_perm = 20, seed = 1)
  ref <- ps_stddiff(dta, gaussian = c("age", "ef"), binary = "female")$tables$stddiff
  expect_identical(obj$tables$stddiff_perm$observed, ref$stddiff)
  expect_identical(obj$tables$stddiff_perm$variable, ref$variable)
})

# ---- Permutation ------------------------------------------------------------

test_that("the percentiles come from ps_stddiff() on permuted groups", {
  # Rebuild the permutations by hand with the same seed: the function draws
  # one sample() of the group column per permutation, in order.
  obj <- ps_stddiff_perm(dta, gaussian = "age", n_perm = 30, seed = 11)
  set.seed(11)
  sims <- vapply(seq_len(30), function(i) {
    d <- dta
    d$tavr <- sample(d$tavr)
    ps_stddiff(d, gaussian = "age")$tables$stddiff$stddiff
  }, numeric(1))
  expect_equal(unlist(obj$tables$stddiff_perm[1, c("p2_5", "p16", "p50", "p84", "p97_5")], use.names = FALSE),
               quantile(sims, c(.025, .16, .5, .84, .975), type = 2, names = FALSE))
})

test_that("rows with a missing group are not permuted into the groups", {
  d <- dta
  d$tavr[1:5] <- NA
  seen <- integer(0)
  obj <- ps_stddiff_perm(d, gaussian = "age", weight_col = "w", n_perm = 5, seed = 2,
                         reweight = function(p) {
                           seen <<- c(seen, sum(is.na(p$tavr[1:5])))
                           rep(1, nrow(p))
                         })
  expect_identical(seen, rep(5L, 6))  # the observed call, then 5 permutations
})

test_that("a seed makes the result reproducible, and different seeds differ", {
  a <- ps_stddiff_perm(dta, gaussian = "age", n_perm = 25, seed = 3)$tables$stddiff_perm
  b <- ps_stddiff_perm(dta, gaussian = "age", n_perm = 25, seed = 3)$tables$stddiff_perm
  c <- ps_stddiff_perm(dta, gaussian = "age", n_perm = 25, seed = 4)$tables$stddiff_perm
  expect_identical(a, b)
  expect_false(identical(a$p16, c$p16))
})

test_that("the caller's random number stream is left as it was", {
  set.seed(99)
  before <- runif(3)
  set.seed(99)
  ps_stddiff_perm(dta, gaussian = "age", n_perm = 10, seed = 5)
  expect_identical(runif(3), before)
})

# ---- Weights ----------------------------------------------------------------

test_that("with weights, reweight is called on each permuted data set and its weights are used", {
  calls <- 0L
  obj <- ps_stddiff_perm(dta, gaussian = "age", weight_col = "w", n_perm = 12, seed = 6,
                         reweight = function(p) {
                           calls <<- calls + 1L
                           ifelse(p$tavr == 1, 2, 1)
                         })
  # Once for the unpermuted data, so the observed row is weighted by the same
  # rule as every permutation, then once per permutation.
  expect_identical(calls, 13L)
  ref <- dta
  ref$w <- ifelse(ref$tavr == 1, 2, 1)
  expect_identical(obj$tables$stddiff_perm$observed,
                   ps_stddiff(ref, gaussian = "age", weight_col = "w")$tables$stddiff$stddiff)
})

test_that("a weight column without reweight stops, and so does reweight without a weight column", {
  d <- dta
  d$w <- 1
  expect_error(ps_stddiff_perm(d, gaussian = "age", weight_col = "w", n_perm = 5), "reweight")
  expect_error(ps_stddiff_perm(d, gaussian = "age", reweight = function(p) rep(1, nrow(p)), n_perm = 5),
               "weight_col")
})

test_that("reweight must return one finite weight per row", {
  expect_error(ps_stddiff_perm(dta, gaussian = "age", weight_col = "w", n_perm = 5,
                               reweight = function(p) 1), "one weight per row")
})

test_that("n_perm must be a positive whole number", {
  expect_error(ps_stddiff_perm(dta, gaussian = "age", n_perm = 0), "n_perm")
  expect_error(ps_stddiff_perm(dta, gaussian = "age", n_perm = 2.5), "n_perm")
})

test_that("n_perm rejects Inf rather than failing later", {
  expect_error(ps_stddiff_perm(dta, gaussian = "age", n_perm = Inf), "n_perm")
})

test_that("a seed also covers randomness inside reweight, including for the observed row", {
  rw <- function(p) runif(nrow(p), 1, 2)
  a <- ps_stddiff_perm(dta, gaussian = "age", weight_col = "w", reweight = rw, n_perm = 5, seed = 8)
  b <- ps_stddiff_perm(dta, gaussian = "age", weight_col = "w", reweight = rw, n_perm = 5, seed = 8)
  expect_identical(a$tables$stddiff_perm, b$tables$stddiff_perm)
})

test_that("warnings from reweight on a permuted data set are not hidden", {
  # Warn only after the first (observed) call, so the warning must come from a permutation.
  calls <- 0L
  expect_warning(
    ps_stddiff_perm(dta, gaussian = "age", weight_col = "w", n_perm = 3, seed = 1,
                    reweight = function(p) {
                      calls <<- calls + 1L
                      if (calls > 1L) warning("model did not converge")
                      rep(1, nrow(p))
                    }),
    "did not converge"
  )
})
