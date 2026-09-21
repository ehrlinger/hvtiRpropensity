###############################################################################
## test_lm_validate.R
##
## Validation of saved binary model bundles without refitting.
###############################################################################

validation_fixture <- function() {
  d <- utils::read.csv(testthat::test_path("fixtures", "lm-sas", "input.csv"))
  d[d$imputation == 1L, ]
}

validation_model <- function() {
  fit_logistic(
    binary ~ x1_imp + x2,
    validation_fixture(),
    outcome_levels = c(0, 1),
    event_level = 1
  )
}

test_that("validate_logistic scores a binary bundle without refitting", {
  model <- validation_model()
  original <- model
  d <- validation_fixture()
  testthat::local_mocked_bindings(
    .fit_imputations = function(...) stop("refit attempted"),
    .model_inference = function(...) stop("inference attempted")
  )

  result <- validate_logistic(model, d, groups = 5L)

  expect_s3_class(result, "lm_validation")
  expect_named(result, c("data", "meta", "tables", "models"))
  expect_named(result$tables, c("calibration", "performance"))
  expect_named(result$tables$calibration,
               c("group", "n", "observed", "expected"))
  expect_named(result$tables$performance,
               c("n", "observed", "expected", "oe_ratio", "auc", "brier"))
  expect_equal(result$data$predicted,
               unname(stats::predict(model$models[[1L]], newdata = d,
                                     type = "response")))
  expect_equal(sum(result$tables$calibration$n), nrow(d))
  expect_identical(model, original)
})

test_that("validate_logistic averages every stored imputation model", {
  d <- utils::read.csv(testthat::test_path("fixtures", "lm-sas", "input.csv"))
  model <- fit_logistic(
    binary ~ x1_imp + x2,
    d,
    id_col = "id",
    imputation_col = "imputation",
    outcome_levels = c(0, 1),
    event_level = 1
  )
  validation <- validation_fixture()
  result <- validate_logistic(model, validation)
  expected <- rowMeans(vapply(
    model$models,
    stats::predict,
    numeric(nrow(validation)),
    newdata = validation,
    type = "response"
  ))

  expect_equal(result$data$predicted, unname(expected))
  expect_length(result$models, 3L)
})

test_that("validate_logistic checks bundle and cohort compatibility", {
  model <- validation_model()
  d <- validation_fixture()

  expect_error(validate_logistic(model, d[, setdiff(names(d), "x2")]),
               "missing.*x2")

  ordinal <- fit_logistic(
    ordinal ~ x1_imp + x2,
    d,
    family = "ordinal",
    outcome_levels = c("low", "mid", "high")
  )
  expect_error(validate_logistic(ordinal, d), "binary")

  changed <- d
  changed$binary[[1L]] <- 2L
  expect_error(validate_logistic(model, changed), "undeclared.*2")

  wrong_version <- model
  wrong_version$meta$bundle_version <- 2L
  expect_error(validate_logistic(wrong_version, d), "bundle version.*1")

  expect_error(validate_logistic(list(), d), "lm_fit")
})

test_that("validate_logistic checks grouping and output names", {
  model <- validation_model()
  d <- validation_fixture()

  expect_error(validate_logistic(model, d, groups = 1L), "groups.*integer.*2")
  expect_error(validate_logistic(model, d, groups = 2.5), "groups.*integer.*2")
  expect_error(validate_logistic(model, d, groups = nrow(d) + 1L),
               "groups.*observations")
  d$predicted <- 0
  expect_error(validate_logistic(model, d), "predicted.*already exists")
})

test_that("validate_logistic computes standard binary performance measures", {
  model <- validation_model()
  d <- validation_fixture()
  result <- validate_logistic(model, d, groups = 6L)
  y <- as.integer(d$binary == 1L)
  p <- result$data$predicted
  ranks <- rank(p, ties.method = "average")
  auc <- (sum(ranks[y == 1L]) - sum(y) * (sum(y) + 1) / 2) /
    (sum(y) * sum(y == 0L))

  expect_equal(result$tables$performance$n, nrow(d))
  expect_equal(result$tables$performance$observed, sum(y))
  expect_equal(result$tables$performance$expected, sum(p))
  expect_equal(result$tables$performance$oe_ratio, sum(y) / sum(p))
  expect_equal(result$tables$performance$auc, auc)
  expect_equal(result$tables$performance$brier, mean((y - p)^2))
})
