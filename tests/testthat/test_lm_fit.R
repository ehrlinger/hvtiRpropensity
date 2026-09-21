###############################################################################
## test_lm_fit.R
##
## General logistic model bundles.
###############################################################################

lm_fixture <- function() {
  utils::read.csv(testthat::test_path("fixtures", "lm-sas", "input.csv"))
}

test_that("fit_logistic models the declared binary event", {
  d <- lm_fixture()
  one <- d[d$imputation == 1L, ]
  obj <- fit_logistic(
    binary ~ x1_imp + x2,
    one,
    family = "binary",
    outcome_levels = c(0, 1),
    event_level = 1
  )

  expect_s3_class(obj, "lm_fit")
  expect_named(obj, c("data", "meta", "tables", "models"))
  expect_identical(obj$meta$bundle_version, 1L)
  expect_identical(obj$meta$model_family, "binary")
  expect_identical(obj$meta$outcome_levels, c("0", "1"))
  expect_identical(obj$meta$event_level, "1")
  expect_identical(obj$meta$predictors, c("x1_imp", "x2"))
  expect_named(obj$meta$package_versions,
               c("R", "hvtiRpropensity", "stats"))
  expect_named(
    obj$tables,
    c("estimates", "covariance", "by_imputation", "fit_status")
  )
  expect_equal(
    obj$data$prob,
    unname(stats::predict(obj$models[[1L]], newdata = one, type = "response")),
    tolerance = 1e-10
  )
})

test_that("fit_logistic ignores the input factor order", {
  d <- lm_fixture()
  d <- d[d$imputation == 1L, ]
  d$binary <- factor(d$binary, levels = c(1, 0))

  obj <- fit_logistic(
    binary ~ x1_imp + x2,
    d,
    family = "binary",
    outcome_levels = c(0, 1),
    event_level = 1
  )
  direct <- stats::glm(
    factor(as.character(d$binary), levels = c("0", "1")) ~ x1_imp + x2,
    data = d,
    family = stats::binomial(),
    na.action = stats::na.exclude
  )

  expect_equal(stats::coef(obj$models[[1L]]), stats::coef(direct))
  expect_equal(obj$data$prob,
               unname(stats::predict(direct, newdata = d, type = "response")))
})

test_that("fit_logistic enforces the declared binary levels", {
  d <- lm_fixture()
  d <- d[d$imputation == 1L, ]
  extra <- d
  extra$binary[[1L]] <- 2L

  expect_error(
    fit_logistic(binary ~ x1_imp, extra, outcome_levels = c(0, 1),
                 event_level = 1),
    "undeclared.*2"
  )
  expect_error(
    fit_logistic(binary ~ x1_imp, d[d$binary == 0L, ],
                 outcome_levels = c(0, 1), event_level = 1),
    "declared.*1.*absent"
  )
  expect_error(
    fit_logistic(binary ~ x1_imp, d, outcome_levels = c(0, 0),
                 event_level = 0),
    "exactly two unique"
  )
  expect_error(
    fit_logistic(binary ~ x1_imp, d, outcome_levels = c(0, 1),
                 event_level = 2),
    "event_level.*declared"
  )
})

test_that("fit_logistic refuses a non-converged constituent fit", {
  d <- lm_fixture()
  d <- d[d$imputation == 1L, ]
  testthat::local_mocked_bindings(
    .model_converged = function(fit) FALSE
  )

  expect_error(
    fit_logistic(binary ~ x1_imp, d, outcome_levels = c(0, 1),
                 event_level = 1),
    "failed to converge.*1"
  )
})

test_that("fit_logistic pools stacked binary imputations", {
  d <- lm_fixture()
  obj <- fit_logistic(
    binary ~ x1_imp + x2,
    d,
    family = "binary",
    id_col = "id",
    imputation_col = "imputation",
    outcome_levels = c(0, 1),
    event_level = 1
  )

  expect_length(obj$models, 3L)
  expect_identical(names(obj$models), c("1", "2", "3"))
  expect_true(all(obj$tables$estimates$pooled))
  expect_identical(unique(obj$tables$by_imputation$imputation), c("1", "2", "3"))
  expect_identical(obj$meta$n_imputations, 3L)
  expect_equal(nrow(obj$data), 90L)
  expect_equal(obj$data$prob, unname(rowMeans(vapply(
    seq_along(obj$models),
    function(i) {
      partition <- d[d$imputation == i, ]
      stats::predict(obj$models[[i]], newdata = partition, type = "response")
    },
    numeric(90L)
  ))))
})

test_that("fit_logistic fits the declared ordinal outcome", {
  d <- lm_fixture()
  d <- d[d$imputation == 1L, ]
  obj <- fit_logistic(
    ordinal ~ x1_imp + x2,
    d,
    family = "ordinal",
    outcome_levels = c("low", "mid", "high")
  )

  expect_true(is.ordered(obj$data$ordinal))
  expect_identical(levels(obj$data$ordinal), c("low", "mid", "high"))
  expect_identical(
    obj$meta$cumulative_direction,
    "P(Y <= level) = logistic(threshold - linear predictor)"
  )
  expect_identical(
    obj$tables$estimates$term,
    c("x1_imp", "x2", "threshold:low|mid", "threshold:mid|high")
  )
  expect_named(obj$data, c(names(d), "prob_low", "prob_mid", "prob_high"))
  expect_equal(unname(rowSums(obj$data[c("prob_low", "prob_mid", "prob_high")])),
               rep(1, nrow(d)), tolerance = 1e-10)

  direct_data <- d
  direct_data$ordinal <- ordered(
    as.character(direct_data$ordinal),
    levels = c("low", "mid", "high")
  )
  direct <- MASS::polr(ordinal ~ x1_imp + x2, direct_data, Hess = TRUE,
                       na.action = stats::na.exclude)
  expect_equal(obj$tables$estimates$estimate,
               unname(c(stats::coef(direct), direct$zeta)), tolerance = 1e-10)
  expect_equal(unname(obj$tables$covariance), unname(stats::vcov(direct)),
               tolerance = 1e-10)
})

test_that("ordinal fitting ignores the input factor order", {
  d <- lm_fixture()
  d <- d[d$imputation == 1L, ]
  d$ordinal <- ordered(d$ordinal, levels = c("high", "mid", "low"))

  obj <- fit_logistic(
    ordinal ~ x1_imp + x2,
    d,
    family = "ordinal",
    outcome_levels = c("low", "mid", "high")
  )
  expected_data <- d
  expected_data$ordinal <- ordered(
    as.character(expected_data$ordinal),
    levels = c("low", "mid", "high")
  )
  expected <- MASS::polr(ordinal ~ x1_imp + x2, expected_data, Hess = TRUE,
                         na.action = stats::na.exclude)

  expect_equal(stats::coef(obj$models[[1L]]), stats::coef(expected))
  expect_equal(obj$models[[1L]]$zeta, expected$zeta)
})

test_that("multicategory fitting refuses to overwrite a probability column", {
  d <- lm_fixture()
  d <- d[d$imputation == 1L, ]
  d$prob_low <- -1

  expect_error(
    fit_logistic(
      ordinal ~ x1_imp + x2,
      d,
      family = "ordinal",
      outcome_levels = c("low", "mid", "high")
    ),
    "Column `prob_low` already exists"
  )
})

test_that("fit_logistic fits the declared nominal outcome and reference", {
  d <- lm_fixture()
  d <- d[d$imputation == 1L, ]
  obj <- fit_logistic(
    nominal ~ x1_imp + x2,
    d,
    family = "nominal",
    outcome_levels = c("a", "b", "c"),
    reference_level = "a"
  )

  expect_identical(obj$meta$reference_level, "a")
  expect_identical(
    obj$tables$estimates$term,
    c("b:(Intercept)", "b:x1_imp", "b:x2",
      "c:(Intercept)", "c:x1_imp", "c:x2")
  )
  expect_named(obj$data, c(names(d), "prob_a", "prob_b", "prob_c"))
  expect_equal(unname(rowSums(obj$data[c("prob_a", "prob_b", "prob_c")])),
               rep(1, nrow(d)), tolerance = 1e-10)
  expect_identical(
    dimnames(obj$tables$covariance),
    list(obj$tables$estimates$term, obj$tables$estimates$term)
  )

  direct_data <- d
  direct_data$nominal <- factor(
    as.character(direct_data$nominal),
    levels = c("a", "b", "c")
  )
  direct <- nnet::multinom(nominal ~ x1_imp + x2, direct_data,
                           Hess = TRUE, trace = FALSE,
                           na.action = stats::na.exclude)
  expect_equal(obj$tables$estimates$estimate,
               unname(as.vector(t(stats::coef(direct)))), tolerance = 1e-10)
  expect_equal(unname(obj$tables$covariance), unname(stats::vcov(direct)),
               tolerance = 1e-10)
})

test_that("nominal fitting reports a level missing from one imputation", {
  d <- lm_fixture()
  d$nominal[d$imputation == 2L & d$nominal == "c"] <- "b"

  expect_error(
    fit_logistic(
      nominal ~ x1_imp + x2,
      d,
      family = "nominal",
      id_col = "id",
      imputation_col = "imputation",
      outcome_levels = c("a", "b", "c"),
      reference_level = "a"
    ),
    "declared.*c.*absent.*imputation 2"
  )
})

test_that("multicategory fits pool stacked imputations", {
  d <- lm_fixture()
  ordinal <- fit_logistic(
    ordinal ~ x1_imp + x2,
    d,
    family = "ordinal",
    id_col = "id",
    imputation_col = "imputation",
    outcome_levels = c("low", "mid", "high")
  )
  nominal <- fit_logistic(
    nominal ~ x1_imp + x2,
    d,
    family = "nominal",
    id_col = "id",
    imputation_col = "imputation",
    outcome_levels = c("a", "b", "c"),
    reference_level = "a"
  )

  expect_length(ordinal$models, 3L)
  expect_true(all(ordinal$tables$estimates$pooled))
  expect_length(nominal$models, 3L)
  expect_true(all(nominal$tables$estimates$pooled))
})

test_that("stacked fitting names an imputation with a missing predictor level", {
  d <- lm_fixture()
  d$category <- c("first", "second", "third")[(d$id %% 3L) + 1L]
  d$category[d$imputation == 2L & d$category == "third"] <- "second"

  expect_error(
    fit_logistic(
      binary ~ category + x1_imp,
      d,
      id_col = "id",
      imputation_col = "imputation",
      outcome_levels = c(0, 1),
      event_level = 1
    ),
    "Predictor `category`.*imputation 2"
  )
})

test_that("stacked fitting names an imputation with an aliased term", {
  d <- lm_fixture()
  d$x_alias <- d$x1_imp + ifelse(d$imputation == 2L, 0, d$id / 1000)

  expect_error(
    fit_logistic(
      binary ~ x1_imp + x_alias + x2,
      d,
      id_col = "id",
      imputation_col = "imputation",
      outcome_levels = c(0, 1),
      event_level = 1
    ),
    "non-finite coefficient.*imputation 2"
  )
})
