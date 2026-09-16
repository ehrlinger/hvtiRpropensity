###############################################################################
## test_ps_stddiff.R
##
## Tests for ps_stddiff(), a port of stddiff.sas (Artis 2019).
##
## Expected values are worked by hand from the SAS formulas on small vectors,
## not produced by the function under test, so a formula error cannot agree
## with itself. The arithmetic is shown beside each one.
###############################################################################

sd_of <- function(obj, v) obj$tables$stddiff$stddiff[obj$tables$stddiff$variable == v]

# ---- Class and structure ----------------------------------------------------

test_that("ps_stddiff() returns a ps_stddiff / ps_data object", {
  d <- data.frame(g = c(1, 1, 1, 0, 0, 0), x = c(1, 2, 3, 2, 4, 6))
  obj <- ps_stddiff(d, treatment_col = "g", gaussian = "x")
  expect_s3_class(obj, "ps_stddiff")
  expect_s3_class(obj, "ps_data")
  expect_identical(obj$data, d)
  expect_named(obj$tables$stddiff, c("variable", "label", "type", "stddiff"))
})

# ---- Gaussian ---------------------------------------------------------------

test_that("gaussian uses sqrt((var1 + var0) / 2), not a size-weighted pool", {
  # g1: 1 2 3, mean 2, var 1.   g0: 2 4 6 8, mean 5, var 20/3.
  # (2 - 5) / sqrt((1 + 20/3) / 2) = -3 / sqrt(23/6) = -1.532262
  # A sample-size-weighted pool gives -3 / sqrt((2*1 + 3*20/3) / 5) = -1.391365.
  d <- data.frame(g = c(1, 1, 1, 0, 0, 0, 0), x = c(1, 2, 3, 2, 4, 6, 8))
  obj <- ps_stddiff(d, treatment_col = "g", gaussian = "x")
  expect_equal(sd_of(obj, "x"), -3 / sqrt(23 / 6))
})

test_that("weighted gaussian uses PROC MEANS weighted variance, n - 1", {
  # g1: x 1 2 3, w 1 1 2. mean 9/4. sum w(x-m)^2 = 1.5625 + .0625 + 1.125 = 2.75, /2 = 1.375
  # g0: x 2 4 6, w 2 1 1. mean 14/4. sum w(x-m)^2 = 4.5 + .25 + 6.25 = 11, /2 = 5.5
  # (2.25 - 3.5) / sqrt((1.375 + 5.5) / 2)
  d <- data.frame(g = c(1, 1, 1, 0, 0, 0), x = c(1, 2, 3, 2, 4, 6), w = c(1, 1, 2, 2, 1, 1))
  obj <- ps_stddiff(d, treatment_col = "g", gaussian = "x", weight_col = "w")
  expect_equal(sd_of(obj, "x"), -1.25 / sqrt(3.4375))
})

test_that("missing values are dropped per variable, missing group per row", {
  d <- data.frame(g = c(1, 1, 1, 0, 0, 0, NA), x = c(1, 2, 3, 2, 4, 6, 100), y = c(1, NA, 3, 2, 4, 6, 0))
  obj <- ps_stddiff(d, treatment_col = "g", gaussian = c("x", "y"))
  expect_equal(sd_of(obj, "x"), -2 / sqrt(2.5))
  # y, g1: 1 3 -> mean 2, var 2.   g0: 2 4 6 -> mean 4, var 4.
  expect_equal(sd_of(obj, "y"), -2 / sqrt(3))
})

# ---- Non-Gaussian / ordinal -------------------------------------------------

test_that("nong_ord uses pooled ranks with ties averaged", {
  # pooled 3 4 5 | 1 2 3 -> ranks 3.5 5 6 | 1 2 3.5
  # g1 mean 29/6, var 19/12.  g0 mean 13/6, var 19/12.
  # (29/6 - 13/6) / sqrt(19/12)
  d <- data.frame(g = c(1, 1, 1, 0, 0, 0), x = c(3, 4, 5, 1, 2, 3))
  obj <- ps_stddiff(d, treatment_col = "g", nong_ord = "x")
  expect_equal(sd_of(obj, "x"), (16 / 6) / sqrt(19 / 12))
})

test_that("weighted nong_ord weights the rank moments, not the ranks", {
  # ranks as above: g1 3.5 5 6, w 1 1 2 -> mean 20.5/4 = 5.125,
  #   sum w(r-m)^2 = 2.640625 + .015625 + 1.53125 = 4.1875, /2 = 2.09375
  # g0 1 2 3.5, w 2 1 1 -> mean 7.5/4 = 1.875, the same sum, /2 = 2.09375
  d <- data.frame(g = c(1, 1, 1, 0, 0, 0), x = c(3, 4, 5, 1, 2, 3), w = c(1, 1, 2, 2, 1, 1))
  obj <- ps_stddiff(d, treatment_col = "g", nong_ord = "x", weight_col = "w")
  expect_equal(sd_of(obj, "x"), 3.25 / sqrt(2.09375))
})

# ---- Binary and categorical -------------------------------------------------

test_that("binary uses the proportion difference over the average Bernoulli variance", {
  # p1 = .5, p0 = .25: .25 / sqrt((.25 + .1875) / 2)
  d <- data.frame(g = c(1, 1, 1, 1, 0, 0, 0, 0), x = c(1, 1, 0, 0, 1, 0, 0, 0))
  obj <- ps_stddiff(d, treatment_col = "g", binary = "x")
  expect_equal(sd_of(obj, "x"), 0.25 / sqrt(0.21875))
  expect_identical(obj$tables$stddiff$type, "binary")
})

test_that("weighted binary proportions come from the weights", {
  # g1: x 1 0, w 3 1 -> p1 = .75.  g0: x 1 0, w 1 1 -> p0 = .5.
  # .25 / sqrt((.1875 + .25) / 2)
  d <- data.frame(g = c(1, 1, 0, 0), x = c(1, 0, 1, 0), w = c(3, 1, 1, 1))
  obj <- ps_stddiff(d, treatment_col = "g", binary = "x", weight_col = "w")
  expect_equal(sd_of(obj, "x"), 0.25 / sqrt(0.21875))
})

test_that("categorical uses Yang and Dalton's Mahalanobis form, last level dropped", {
  # g1: a a b c -> .5 .25 .25.   g0: a b b c -> .25 .5 .25.   drop c.
  # t - c = (.25, -.25).  S = [.21875 -.125; -.125 .21875], det = .0322265625
  # d^2 = 2 * .25 * (.21875 - .125) * .25 / det
  d <- data.frame(g = c(1, 1, 1, 1, 0, 0, 0, 0), x = c("a", "a", "b", "c", "a", "b", "b", "c"))
  obj <- ps_stddiff(d, treatment_col = "g", categorical = "x")
  expect_equal(sd_of(obj, "x"), sqrt(2 * 0.25 * 0.09375 * 0.25 / 0.0322265625))
  expect_identical(obj$tables$stddiff$type, "categorical")
})

test_that("categorical is routed by observed level count, as the macro does", {
  d <- data.frame(g = c(1, 1, 0, 0), one = c(1, 1, 1, 1), two = c("m", "f", "f", "f"),
                  three = c("a", "b", "c", "a"))
  obj <- ps_stddiff(d, treatment_col = "g", binary = "one", categorical = c("two", "three"))
  tbl <- obj$tables$stddiff
  expect_identical(tbl$type[tbl$variable == "one"], "onelevel")
  expect_identical(sd_of(obj, "one"), 0)
  expect_identical(tbl$type[tbl$variable == "three"], "categorical")
  expect_identical(tbl$type[tbl$variable == "two"], "binary")
})

test_that("a binary variable must be coded 0/1, rather than silently becoming categorical", {
  d <- data.frame(g = c(1, 1, 0, 0), x = c(1, 2, 3, 1))
  expect_error(ps_stddiff(d, treatment_col = "g", binary = "x"), "0/1")
  d <- data.frame(g = c(1, 1, 0, 0), x = c(TRUE, FALSE, FALSE, FALSE))
  expect_identical(ps_stddiff(d, treatment_col = "g", binary = "x")$tables$stddiff$type, "binary")
})

test_that("an unordered factor is accepted as categorical", {
  d <- data.frame(g = c(1, 1, 1, 1, 0, 0, 0, 0),
                  x = factor(c("a", "a", "b", "c", "a", "b", "b", "c")))
  obj <- ps_stddiff(d, treatment_col = "g", categorical = "x")
  expect_equal(sd_of(obj, "x"), sqrt(2 * 0.25 * 0.09375 * 0.25 / 0.0322265625))
})

test_that("weighted categorical proportions come from the weights", {
  # Unweighted every level is 1/3 in both groups, so the difference would be 0.
  # g1: a b c, w 2 1 1 -> .5 .25 .25.   g0: a b c, w 1 2 1 -> .25 .5 .25.
  # Same proportions as the unweighted categorical case above.
  d <- data.frame(g = c(1, 1, 1, 0, 0, 0), x = c("a", "b", "c", "a", "b", "c"), w = c(2, 1, 1, 1, 2, 1))
  obj <- ps_stddiff(d, treatment_col = "g", categorical = "x", weight_col = "w")
  expect_equal(sd_of(obj, "x"), sqrt(2 * 0.25 * 0.09375 * 0.25 / 0.0322265625))
})

test_that("a variable with no non-missing values is NA and does not stop the table", {
  d <- data.frame(g = c(1, 1, 0, 0), x = c(NA, NA, NA, NA), y = c(1, 2, 3, 5))
  obj <- ps_stddiff(d, treatment_col = "g", categorical = "x", gaussian = "y")
  expect_identical(sd_of(obj, "x"), NA_real_)
  expect_true(is.finite(sd_of(obj, "y")))
})

# ---- Output -----------------------------------------------------------------

test_that("rows keep the order variables were given, and labels are carried", {
  d <- data.frame(g = c(1, 1, 1, 0, 0, 0), x = c(1, 2, 3, 2, 4, 6), b = c(1, 0, 1, 0, 0, 1))
  attr(d$x, "label") <- "Age at surgery"
  obj <- ps_stddiff(d, treatment_col = "g", binary = "b", gaussian = "x")
  expect_identical(obj$tables$stddiff$variable, c("x", "b"))
  expect_identical(obj$tables$stddiff$label, c("Age at surgery", NA_character_))
})

test_that("a categorical variable whose groups share no levels is NA with a warning, and does not stop the table", {
  # S is singular here: group 1 holds only a and b, group 0 only c.
  d <- data.frame(g = c(1, 1, 0, 0), x = c("a", "b", "c", "c"), y = c(1, 2, 3, 4))
  expect_warning(obj <- ps_stddiff(d, treatment_col = "g", categorical = "x", gaussian = "y"), "`x`")
  expect_true(is.finite(sd_of(obj, "y")))
  expect_identical(sd_of(obj, "x"), NA_real_)
})

test_that("a variable with zero pooled variance gets NA, not Inf", {
  d <- data.frame(g = c(1, 1, 0, 0), x = c(5, 5, 5, 5))
  obj <- ps_stddiff(d, treatment_col = "g", gaussian = "x")
  expect_identical(sd_of(obj, "x"), NA_real_)
})

# ---- Validation -------------------------------------------------------------

test_that("the group must be 0/1 with both groups present", {
  d <- data.frame(g = c(1, 2, 1, 2), x = 1:4)
  expect_error(ps_stddiff(d, treatment_col = "g", gaussian = "x"), "binary")
  d <- data.frame(g = c(1, 1, 1, 1), x = 1:4)
  expect_error(ps_stddiff(d, treatment_col = "g", gaussian = "x"), "both")
})

test_that("at least one variable is required, and each only once", {
  d <- data.frame(g = c(1, 1, 0, 0), x = 1:4)
  expect_error(ps_stddiff(d, treatment_col = "g"), "at least one")
  expect_error(ps_stddiff(d, treatment_col = "g", gaussian = "x", nong_ord = "x"), "more than one")
})

test_that("gaussian and nong_ord variables must be numeric", {
  d <- data.frame(g = c(1, 1, 0, 0), x = c("a", "b", "a", "b"))
  expect_error(ps_stddiff(d, treatment_col = "g", gaussian = "x"), "numeric")
})

test_that("weights must be positive and not missing", {
  d <- data.frame(g = c(1, 1, 0, 0), x = 1:4, w = c(1, 0, 1, 1))
  expect_error(ps_stddiff(d, treatment_col = "g", gaussian = "x", weight_col = "w"), "positive")
  d$w <- c(1, Inf, 1, 1)
  expect_error(ps_stddiff(d, treatment_col = "g", gaussian = "x", weight_col = "w"), "finite")
})
