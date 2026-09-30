###############################################################################
## test_ps_weight.R
##
## Tests for ps_weight().
###############################################################################

dta <- sample_ps_data(n = 100, seed = 42)

# ---- Class and structure ----------------------------------------------------

test_that("ps_weight() returns a ps_weight / ps_data object", {
  obj <- ps_weight(dta)
  expect_s3_class(obj, "ps_weight")
  expect_s3_class(obj, "ps_data")
})

test_that("ps_weight() object has $data, $meta, $tables slots", {
  obj <- ps_weight(dta)
  expect_true(!is.null(obj$data))
  expect_true(!is.null(obj$meta))
  expect_true(!is.null(obj$tables))
})

# ---- $data ------------------------------------------------------------------

test_that("$data has the same number of rows as the input", {
  obj <- ps_weight(dta)
  expect_equal(nrow(obj$data), nrow(dta))
})

test_that("$data gains an 'iptw' column of positive numerics", {
  obj <- ps_weight(dta)
  expect_true("iptw" %in% names(obj$data))
  expect_true(is.numeric(obj$data$iptw))
  expect_true(all(obj$data$iptw > 0))
})

# ---- Estimands --------------------------------------------------------------

test_that("ATE, ATT, and ATC estimands all produce valid weights", {
  for (est in c("ATE", "ATT", "ATC")) {
    obj <- ps_weight(dta, estimand = est)
    expect_true(all(obj$data$iptw > 0),
                info = paste("estimand =", est))
    expect_equal(obj$meta$estimand, est)
  }
})

# ---- Stabilisation ----------------------------------------------------------

test_that("stabilised weights have smaller variance than unstabilised", {
  obj_stab   <- ps_weight(dta, stabilise = TRUE)
  obj_unstab <- ps_weight(dta, stabilise = FALSE)
  expect_lte(stats::var(obj_stab$data$iptw),
             stats::var(obj_unstab$data$iptw))
})

# ---- Trimming ---------------------------------------------------------------

test_that("trim winsorises weights to the expected quantile range", {
  obj <- ps_weight(dta, trim = 0.05)
  w   <- obj$data$iptw
  q_lo <- stats::quantile(ps_weight(dta, trim = NULL)$data$iptw, 0.05)
  q_hi <- stats::quantile(ps_weight(dta, trim = NULL)$data$iptw, 0.95)
  expect_true(all(w >= min(q_lo) - .Machine$double.eps * 10))
  expect_true(all(w <= max(q_hi) + .Machine$double.eps * 10))
})

test_that("trim = NULL does not trim anything", {
  obj_notrim <- ps_weight(dta, trim = NULL, stabilise = FALSE)
  obj_trim   <- ps_weight(dta, trim = 0.01, stabilise = FALSE)
  # range with trimming should be <= range without
  expect_lte(diff(range(obj_trim$data$iptw)),
             diff(range(obj_notrim$data$iptw)) + .Machine$double.eps)
})

test_that("invalid trim values raise an error", {
  expect_error(ps_weight(dta, trim = 0))
  expect_error(ps_weight(dta, trim = 0.5))
  expect_error(ps_weight(dta, trim = -0.1))
  expect_error(ps_weight(dta, trim = "a"))
})

# ---- $meta ------------------------------------------------------------------

test_that("$meta contains required keys", {
  obj <- ps_weight(dta)
  expect_true(all(c("treatment_col", "score_col", "weight_col", "estimand",
                    "stabilised", "trim", "method", "n_total") %in%
                    names(obj$meta)))
})

test_that("n_total equals nrow(data)", {
  obj <- ps_weight(dta)
  expect_equal(obj$meta$n_total, nrow(dta))
})

# ---- $tables ----------------------------------------------------------------

test_that("$tables contains SMD and effective N tables", {
  obj <- ps_weight(dta)
  expect_true(all(c("smd_unweighted", "smd_weighted",
                    "group_counts", "effective_n") %in%
                    names(obj$tables)))
})

test_that("effective_n has control and treated rows", {
  obj <- ps_weight(dta)
  expect_equal(nrow(obj$tables$effective_n), 2L)
  expect_true(all(c("group", "n_effective") %in%
                    names(obj$tables$effective_n)))
  expect_true(all(obj$tables$effective_n$n_effective > 0))
})

# ---- Custom weight_col ------------------------------------------------------

test_that("custom weight_col name is respected", {
  obj <- ps_weight(dta, weight_col = "mt_wt")
  expect_true("mt_wt" %in% names(obj$data))
})

# ---- Input validation -------------------------------------------------------

test_that("ps_weight() errors on non-data-frame input", {
  expect_error(ps_weight(list(tavr = 1, prob_t = 0.5)))
})

test_that("ps_weight() errors on missing columns", {
  expect_error(ps_weight(dta, treatment_col = "not_here"))
  expect_error(ps_weight(dta, score_col = "not_here"))
})

test_that("ps_weight() errors on non-binary treatment column", {
  bad <- dta
  bad$tavr <- 2L
  expect_error(ps_weight(bad))
})

test_that("ps_weight() errors on out-of-range propensity score", {
  bad <- dta
  bad$prob_t[1] <- -0.1
  expect_error(ps_weight(bad))
})

test_that("group tables name groups by their treatment value", {
  obj <- ps_weight(dta)
  expect_named(obj$tables$group_counts, c("group", "treated", "n"))
  expect_named(obj$tables$effective_n, c("group", "treated", "n_effective"))
  expect_identical(obj$tables$group_counts$group, c("0", "1"))
  expect_identical(obj$tables$effective_n$treated, c(FALSE, TRUE))

  logical_dta <- dta
  logical_dta$tavr <- logical_dta$tavr == 1
  expect_identical(ps_weight(logical_dta)$tables$group_counts$group, c("FALSE", "TRUE"))
})
