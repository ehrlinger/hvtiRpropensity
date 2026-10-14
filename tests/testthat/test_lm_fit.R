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
  skip_if_not_installed("MASS")
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
  skip_if_not_installed("MASS")
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
  skip_if_not_installed("nnet")
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
  skip_if_not_installed("MASS")
  skip_if_not_installed("nnet")
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

test_that("fit_linear matches lm() on one dataset and reports R-squared", {
  d <- sample_ps_data(n = 150, seed = 42)
  fit <- fit_linear(ef ~ age + female, d)
  ref <- stats::lm(ef ~ age + female, d)
  expect_s3_class(fit, c("rm_fit", "ps_data"))
  expect_equal(fit$tables$estimates$estimate, unname(stats::coef(ref)))
  expect_equal(fit$tables$estimates$std.error, unname(sqrt(diag(stats::vcov(ref)))))
  expect_equal(fit$data$fitted, unname(stats::fitted(ref)))
  expect_identical(fit$meta$model_family, "linear")
  expect_identical(fit$meta$method, "ols")
  expect_equal(fit$tables$fit_quality$r_squared, summary(ref)$r.squared)
  expect_equal(fit$tables$fit_quality$residual_se, stats::sigma(ref))
  expect_true(all(c("R", "hvtiRpropensity", "stats") %in% names(fit$meta$package_versions)))
})

test_that("fit_linear pools imputations and withholds R-squared", {
  d <- sample_ps_data(n = 150, seed = 42)
  second <- d
  second$age <- second$age + stats::rnorm(nrow(d), 0, 0.5)
  stacked <- rbind(transform(d, imp = 1L), transform(second, imp = 2L))
  fit <- fit_linear(ef ~ age + female, stacked, imputation_col = "imp")
  expect_identical(fit$meta$method, "ols-MI")
  expect_identical(fit$meta$n_imputations, 2L)
  expect_true(all(fit$tables$estimates$pooled))
  expect_true(all(is.na(fit$tables$fit_quality$r_squared)))
  expect_identical(nrow(fit$tables$fit_quality), 2L)
  expect_identical(nrow(fit$data), nrow(d))
})

test_that("fit_linear reads the outcome inside a transformed response", {
  d <- sample_ps_data(n = 150, seed = 42)
  fit <- fit_linear(log(ef) ~ age, d)
  expect_identical(fit$meta$outcome_col, "ef")
})

test_that("fit_linear stops on a non-numeric outcome and on an aliased term", {
  d <- sample_ps_data(n = 150, seed = 42)
  d$label <- as.character(d$female)
  expect_error(fit_linear(label ~ age, d), "must be numeric")
  d$age_copy <- d$age
  expect_error(fit_linear(ef ~ age + age_copy, d), "aliased in imputation 1: age_copy")
})

test_that("fit_linear refuses to overwrite a data column with its predictions", {
  d <- sample_ps_data(n = 150, seed = 42)
  d$fitted <- 1
  expect_error(fit_linear(ef ~ age, d), "already exist")
  d$fitted <- NULL
  expect_error(fit_linear(ef ~ age, d, prediction_col = "ef"), "already exist")
  expect_error(fit_linear(ef ~ age, d, prediction_col = "id"), "already exist")
})

test_that("fit_linear requires prediction_col to be one non-empty string", {
  d <- sample_ps_data(n = 150, seed = 42)
  before <- d
  expect_error(fit_linear(ef ~ age, d, prediction_col = 2), "`prediction_col` must be")
  expect_error(fit_linear(ef ~ age, d, prediction_col = NA_character_), "`prediction_col` must be")
  expect_error(fit_linear(ef ~ age, d, prediction_col = ""), "`prediction_col` must be")
  expect_error(fit_linear(ef ~ age, d, prediction_col = c("a", "b")), "`prediction_col` must be")
  expect_identical(d, before)
})

test_that("fit_linear requires outcome_col to name the formula's response", {
  d <- sample_ps_data(n = 150, seed = 42)
  second <- d
  second$age <- second$age + 1
  stacked <- rbind(transform(d, imp = 1L), transform(second, imp = 2L))
  expect_error(
    fit_linear(ef ~ female, stacked, outcome_col = "age", imputation_col = "imp"),
    "`outcome_col` (\"age\") is not the formula's response (ef)",
    fixed = TRUE
  )
  expect_identical(fit_linear(ef ~ age, d, outcome_col = "ef")$meta$outcome_col, "ef")
  expect_identical(fit_linear(log(ef) ~ age, d, outcome_col = "ef")$meta$outcome_col, "ef")
})

test_that("fit_linear expands a dot formula without the id or imputation column", {
  d <- sample_ps_data(n = 150, seed = 42)
  small <- d[c("ef", "age", "female", "id")]
  fit <- fit_linear(ef ~ ., small)
  expect_identical(names(stats::coef(fit$models[[1L]])), c("(Intercept)", "age", "female"))
  expect_identical(deparse(fit$meta$formula), "ef ~ age + female")
  expect_equal(unname(stats::coef(fit$models[[1L]])),
               unname(stats::coef(stats::lm(ef ~ age + female, small))))

  second <- small
  second$age <- second$age + 1
  stacked <- rbind(transform(small, imp = 1L), transform(second, imp = 2L))
  pooled <- fit_linear(ef ~ ., stacked, imputation_col = "imp")
  expect_setequal(pooled$tables$estimates$term, c("(Intercept)", "age", "female"))
  for (m in pooled$models) {
    expect_identical(names(stats::coef(m)), c("(Intercept)", "age", "female"))
  }
  expect_identical(deparse(pooled$meta$formula), "ef ~ age + female")
})

test_that("fit_logistic expands a dot formula without the id or imputation column", {
  d <- lm_fixture()
  small <- d[c("binary", "x1_imp", "x2", "id", "imputation")]
  obj <- fit_logistic(binary ~ ., small, family = "binary",
                      imputation_col = "imputation",
                      outcome_levels = c(0, 1), event_level = 1)
  expect_identical(obj$meta$predictors, c("x1_imp", "x2"))
  expect_identical(names(stats::coef(obj$models[[1L]])), c("(Intercept)", "x1_imp", "x2"))
})

# Stacked imputations whose `dropme` column has levels a and b in imputation 1
# but a, b and c in imputation 2. Only a model that keeps `dropme` should care.
mi_with_dropme <- function() {
  withr::local_seed(11)
  n <- 60
  base <- data.frame(id = seq_len(n), binary = rep(0:1, n / 2), count = stats::rpois(n, 3),
                     other = stats::rpois(n, 2))
  one <- function(i, levels) {
    cbind(base, imp = i, x = stats::rnorm(n), dropme = factor(rep_len(levels, n)))
  }
  rbind(one(1L, c("a", "b")), one(2L, c("a", "b", "c")))
}

test_that(".retained_predictors drops subtracted variables and keeps offsets", {
  expect_identical(.retained_predictors(y ~ x + d + log(z) + offset(log(t)) - d + `odd name`),
                   c("x", "z", "odd name", "t"))
  expect_identical(.retained_predictors(y ~ 1), character())
})

test_that("predictors a formula subtracts are not checked across imputations", {
  d <- mi_with_dropme()
  dot <- fit_logistic(binary ~ . - dropme - count - other, d, imputation_col = "imp",
                      outcome_levels = c(0, 1), event_level = 1)
  explicit <- fit_logistic(binary ~ x, d, imputation_col = "imp",
                           outcome_levels = c(0, 1), event_level = 1)
  expect_identical(dot$meta$predictors, "x")
  expect_equal(dot$tables$estimates, explicit$tables$estimates)
  expect_error(fit_logistic(binary ~ x + dropme, d, imputation_col = "imp",
                            outcome_levels = c(0, 1), event_level = 1),
               "`dropme` has different observed levels")

  lin <- fit_linear(count ~ . - dropme - binary - other, d, imputation_col = "imp")
  expect_equal(lin$tables$estimates, fit_linear(count ~ x, d, imputation_col = "imp")$tables$estimates)
})

test_that("fit_linear requires a one-column response even with outcome_col", {
  d <- mi_with_dropme()[1:60, ]
  expect_error(fit_linear(I(count + other) ~ x, d, outcome_col = "count"), "exactly one column")
  expect_error(fit_linear(cbind(count, other) ~ x, d, outcome_col = "count"), "exactly one column")
  expect_error(fit_linear(I(count + other) ~ x, d), "exactly one column")
  expect_identical(fit_linear(log1p(count) ~ x, d, outcome_col = "count")$meta$outcome_col, "count")
})

test_that("a subtracted column missing from data is caught by the package, not lm()", {
  d <- mi_with_dropme()[1:60, ]
  zz <- 1 # a same-named object in the calling environment must not be used
  expect_error(fit_linear(count ~ x - zz, d), "Required column(s) not found in `data`: zz", fixed = TRUE)
  expect_error(fit_linear(count ~ . - zz, d), "Required column(s) not found in `data`: zz", fixed = TRUE)
  expect_error(fit_logistic(binary ~ . - zz, d, outcome_levels = c(0, 1)),
               "Required column(s) not found in `data`: zz", fixed = TRUE)
})

test_that("fit_logistic refuses a transformed outcome", {
  d <- mi_with_dropme()[1:60, ]
  expect_error(fit_logistic(log(binary + 1) ~ x, d, outcome_levels = c(0, 1)), "untransformed outcome")
})

test_that("a dot that leaves no columns gets a clear error", {
  d <- mi_with_dropme()[1:60, c("id", "count")]
  expect_error(fit_linear(count ~ ., d), "leaves no columns")
})
