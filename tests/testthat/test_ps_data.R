###############################################################################
## test_ps_data.R
##
## Tests for the ps_data base class (new_ps_data, is_ps_data,
## print.ps_data, summary.ps_data).
###############################################################################

test_that("new_ps_data() returns correct class vector", {
  obj <- hvtiRpropensity:::new_ps_data(
    data     = data.frame(x = 1),
    meta     = list(method = "test"),
    tables   = list(),
    subclass = "ps_test"
  )
  expect_s3_class(obj, "ps_test")
  expect_s3_class(obj, "ps_data")
  expect_identical(class(obj), c("ps_test", "ps_data"))
})

test_that("new_ps_data() stores the four stable slots", {
  df  <- data.frame(a = 1:3)
  meta <- list(n_total = 3L, method = "unit_test")
  tbl  <- list(counts = data.frame(n = 3L))
  models <- list(one = stats::lm(a ~ 1, data = df))

  obj <- hvtiRpropensity:::new_ps_data(
    data     = df,
    meta     = meta,
    tables   = tbl,
    models   = models,
    subclass = "ps_test"
  )
  expect_named(obj, c("data", "meta", "tables", "models"))
  expect_identical(obj$data,   df)
  expect_identical(obj$meta,   meta)
  expect_identical(obj$tables, tbl)
  expect_identical(obj$models, models)
})

test_that("new_ps_data() defaults to an empty model list", {
  obj <- hvtiRpropensity:::new_ps_data(
    data = data.frame(x = 1), meta = list(), subclass = "ps_test"
  )
  expect_identical(obj$models, list())
})

test_that("new_ps_data() errors on invalid inputs", {
  expect_error(
    hvtiRpropensity:::new_ps_data(
      data = list(x = 1),   # not a data.frame
      meta = list(), tables = list(), subclass = "ps_test"
    )
  )
  expect_error(
    hvtiRpropensity:::new_ps_data(
      data = data.frame(x = 1), meta = list(), tables = list(),
      subclass = character(0)  # zero-length
    )
  )
  expect_error(
    hvtiRpropensity:::new_ps_data(
      data = data.frame(x = 1), meta = list(), models = 1,
      subclass = "ps_test"
    )
  )
})

test_that("is_ps_data() returns TRUE for ps_data objects and FALSE otherwise", {
  dta <- sample_ps_data(n = 50, seed = 1)
  obj <- ps_match(dta)
  expect_true(is_ps_data(obj))
  expect_false(is_ps_data(list()))
  expect_false(is_ps_data(NULL))
  expect_false(is_ps_data(42))
})

test_that("print.ps_data() runs without error and returns invisibly", {
  dta <- sample_ps_data(n = 50, seed = 1)
  obj <- ps_match(dta)
  expect_output(out <- print(obj))
  expect_identical(out, obj)
})

test_that("summary.ps_data() runs without error and returns tables invisibly", {
  dta <- sample_ps_data(n = 50, seed = 1)
  obj <- ps_match(dta)
  expect_no_error(tbls <- summary(obj))
  expect_identical(tbls, obj$tables)
})
