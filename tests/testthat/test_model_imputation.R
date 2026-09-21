###############################################################################
## test_model_imputation.R
##
## Strict single-dataset and stacked-imputation fitting.
###############################################################################

stacked_probe <- function() {
  base <- data.frame(
    id = 1:4,
    response = c(0, 1, 0, 1),
    x = c(-1, 0, 1, 2)
  )
  do.call(rbind, lapply(1:3, function(imp) {
    partition <- base
    partition$imp <- imp
    partition$x <- partition$x + imp / 10
    partition
  }))
}

probe_fit <- function(d) stats::lm(response ~ x, data = d, na.action = stats::na.exclude)
probe_predict <- function(fit, d) as.numeric(stats::predict(fit, newdata = d))

test_that("single-dataset fitting retains the fit and row accounting", {
  d <- stacked_probe()[1:4, c("id", "response", "x")]
  result <- .fit_imputations(d, "response", "id", NULL, probe_fit, probe_predict)

  expect_named(result, c("data", "predictions", "models", "status", "imputations"))
  expect_named(result$models, "1")
  expect_equal(result$data, d)
  expect_equal(result$predictions, probe_predict(result$models[[1L]], d))
  expect_equal(result$status$n_input, 4L)
  expect_equal(result$status$n_analyzed, 4L)
  expect_equal(result$status$n_excluded, 0L)
  expect_true(result$status$converged)
})

test_that("stacked fitting aligns patients and averages complete predictions", {
  d <- stacked_probe()
  d <- d[order(d$imp, -d$id), ]
  result <- .fit_imputations(d, "response", "id", "imp", probe_fit, probe_predict)

  expect_identical(result$imputations, 1:3)
  expect_identical(result$data$id, 4:1)
  expect_false("imp" %in% names(result$data))
  expect_length(result$models, 3L)
  want <- rowMeans(vapply(1:3, function(imp) {
    partition <- d[d$imp == imp, , drop = FALSE]
    partition <- partition[match(result$data$id, partition$id), , drop = FALSE]
    probe_predict(probe_fit(partition), partition)
  }, numeric(4L)))
  expect_equal(result$predictions, want)
})

test_that("stacked fitting refuses an incomplete patient panel", {
  d <- stacked_probe()
  d <- d[!(d$id == 4 & d$imp == 2), ]
  expect_error(
    .fit_imputations(d, "response", "id", "imp", probe_fit, probe_predict),
    "Imputation 2.*missing 1 patient"
  )
})

test_that("stacked fitting refuses duplicate patient keys", {
  d <- rbind(stacked_probe(), stacked_probe()[1L, ])
  expect_error(
    .fit_imputations(d, "response", "id", "imp", probe_fit, probe_predict),
    "Imputation 1.*duplicate patient"
  )
})

test_that("stacked fitting refuses missing patient keys", {
  d <- stacked_probe()
  d$id[d$id == 3L] <- NA_integer_
  expect_error(
    .fit_imputations(d, "response", "id", "imp", probe_fit, probe_predict),
    "Imputation 1.*missing patient key"
  )
})

test_that("stacked fitting requires at least two imputations", {
  d <- stacked_probe()
  d <- d[d$imp == 1L, ]
  expect_error(
    .fit_imputations(d, "response", "id", "imp", probe_fit, probe_predict),
    "at least 2 distinct"
  )
})

test_that("stacked fitting refuses a response that changes by patient", {
  d <- stacked_probe()
  d$response[d$id == 2L & d$imp == 3L] <- 0
  expect_error(
    .fit_imputations(d, "response", "id", "imp", probe_fit, probe_predict),
    "Imputation 3.*response differs for 1 patient"
  )
})

test_that("stacked fitting refuses a prediction with the wrong row count", {
  d <- stacked_probe()
  short_predict <- function(fit, newdata) probe_predict(fit, newdata)[-1L]
  expect_error(
    .fit_imputations(d, "response", "id", "imp", probe_fit, short_predict),
    "Imputation 1.*3 predictions.*4 rows"
  )
})

test_that("stacked fitting refuses changed prediction columns", {
  d <- stacked_probe()
  matrix_predict <- function(fit, newdata) {
    out <- cbind(first = probe_predict(fit, newdata), second = 0)
    if (newdata$x[[1L]] > -0.8) colnames(out) <- c("first", "changed")
    out
  }
  expect_error(
    .fit_imputations(d, "response", "id", "imp", probe_fit, matrix_predict),
    "Imputation 3.*prediction columns"
  )
})
