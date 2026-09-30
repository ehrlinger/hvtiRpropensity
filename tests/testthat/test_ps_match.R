###############################################################################
## test_ps_match.R
##
## Tests for ps_match().
###############################################################################

dta <- sample_ps_data(n = 100, seed = 42)

# ---- Class and structure ----------------------------------------------------

test_that("ps_match() returns a ps_match / ps_data object", {
  obj <- ps_match(dta)
  expect_s3_class(obj, "ps_match")
  expect_s3_class(obj, "ps_data")
})

test_that("ps_match() object has $data, $meta, $tables slots", {
  obj <- ps_match(dta)
  expect_true(!is.null(obj$data))
  expect_true(!is.null(obj$meta))
  expect_true(!is.null(obj$tables))
})

# ---- $data ------------------------------------------------------------------

test_that("$data has the same number of rows as the input", {
  obj <- ps_match(dta)
  expect_equal(nrow(obj$data), nrow(dta))
})

test_that("$data gains a 'match' column with values 0 or 1", {
  obj <- ps_match(dta)
  expect_true("match" %in% names(obj$data))
  expect_true(all(obj$data$match %in% c(0L, 1L)))
})

test_that("matched subset has equal group sizes", {
  obj <- ps_match(dta)
  matched <- obj$data[obj$data$match == 1L, ]
  expect_equal(sum(matched$tavr == 1L), sum(matched$tavr == 0L))
})

# ---- $meta ------------------------------------------------------------------

test_that("$meta contains required keys", {
  obj <- ps_match(dta)
  expect_true(all(c("treatment_col", "score_col", "id_col", "match_col",
                    "method", "n_total", "n_matched", "n_unmatched") %in%
                    names(obj$meta)))
})

test_that("n_total equals nrow(data)", {
  obj <- ps_match(dta)
  expect_equal(obj$meta$n_total, nrow(dta))
})

test_that("n_matched + n_unmatched equals number of treated patients", {
  obj <- ps_match(dta)
  n_trt <- sum(dta$tavr == 1L)
  expect_equal(obj$meta$n_matched + obj$meta$n_unmatched, n_trt)
})

# ---- $tables ----------------------------------------------------------------

test_that("$tables contains smd_before and smd_after", {
  obj <- ps_match(dta)
  expect_true(all(c("smd_before", "smd_after",
                    "group_counts_before", "group_counts_after") %in%
                    names(obj$tables)))
})

test_that("SMD tables have 'variable' and 'smd' columns", {
  obj <- ps_match(dta)
  expect_true(all(c("variable", "smd") %in% names(obj$tables$smd_before)))
  expect_true(all(c("variable", "smd") %in% names(obj$tables$smd_after)))
})

# ---- Caliper ----------------------------------------------------------------

test_that("caliper reduces or maintains number of matched pairs", {
  obj_no_cal <- ps_match(dta)
  obj_cal    <- ps_match(dta, caliper = 0.01)
  expect_lte(obj_cal$meta$n_matched, obj_no_cal$meta$n_matched)
})

test_that("very tight caliper can match zero pairs without error", {
  expect_no_error(obj <- ps_match(dta, caliper = 1e-10))
  expect_gte(obj$meta$n_matched, 0L)
})

# ---- Custom column names ----------------------------------------------------

test_that("custom match_col name is respected", {
  obj <- ps_match(dta, match_col = "is_matched")
  expect_true("is_matched" %in% names(obj$data))
})

# ---- Input validation -------------------------------------------------------

test_that("ps_match() errors on non-data-frame input", {
  expect_error(ps_match(list(tavr = 1, prob_t = 0.5, id = 1)))
})

test_that("ps_match() errors on missing treatment column", {
  expect_error(ps_match(dta, treatment_col = "missing_col"))
})

test_that("ps_match() errors on non-binary treatment column", {
  bad <- dta
  bad$tavr <- 2L
  expect_error(ps_match(bad))
})

test_that("ps_match() errors on out-of-range propensity score", {
  bad <- dta
  bad$prob_t[1] <- 1.5
  expect_error(ps_match(bad))
})

test_that("ps_match() errors on invalid caliper", {
  expect_error(ps_match(dta, caliper = -0.1))
  expect_error(ps_match(dta, caliper = "a"))
})

test_that("ps_match() errors when ratio != 1", {
  expect_error(ps_match(dta, ratio = 2))
})

test_that("group count tables name groups by their treatment value", {
  obj <- ps_match(dta)
  for (tbl in obj$tables[c("group_counts_before", "group_counts_after")]) {
    expect_named(tbl, c("group", "treated", "n"))
    expect_identical(tbl$group, c("0", "1"))
    expect_identical(tbl$treated, c(FALSE, TRUE))
  }
  expect_equal(obj$tables$group_counts_before$n, c(sum(dta$tavr == 0), sum(dta$tavr == 1)))
})
