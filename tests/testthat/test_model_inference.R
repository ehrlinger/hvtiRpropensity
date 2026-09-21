###############################################################################
## test_model_inference.R
##
## Normalized single-fit inference and Rubin pooling.
###############################################################################

rubin_inputs <- function() {
  q <- list(c(a = 1, b = 2), c(a = 3, b = 4))
  u <- list(diag(c(a = 1, b = 4)), diag(c(a = 3, b = 2)))
  u <- lapply(u, function(x) {
    dimnames(x) <- list(c("a", "b"), c("a", "b"))
    x
  })
  list(q = q, u = u)
}

test_that("Rubin pooling uses within and between-imputation covariance", {
  input <- rubin_inputs()
  result <- .rubin_pool(input$q, input$u)
  expected_covariance <- diag(c(2, 3)) + 1.5 * matrix(2, 2, 2)
  dimnames(expected_covariance) <- list(c("a", "b"), c("a", "b"))
  expected_df <- c(
    a = (1 + 2 / 3)^2,
    b = (1 + 3 / 3)^2
  )

  expect_equal(result$estimates$estimate, c(2, 3))
  expect_equal(result$covariance, expected_covariance)
  expect_equal(result$estimates$std.error, unname(sqrt(diag(expected_covariance))))
  expect_equal(result$estimates$df, unname(expected_df))
  expect_true(all(result$estimates$pooled))
  expect_equal(result$estimates$odds_ratio, exp(c(2, 3)))
  expect_equal(nrow(result$by_imputation), 4L)
})

test_that("Rubin pooling uses infinite df when between variance is zero", {
  estimates <- list(c(a = 1), c(a = 1), c(a = 1))
  covariance <- matrix(2, 1, 1, dimnames = list("a", "a"))
  result <- .rubin_pool(estimates, rep(list(covariance), 3L))
  expect_identical(result$estimates$df, Inf)
  expect_equal(result$covariance, covariance)
})

test_that("model inference preserves imputation labels", {
  d <- data.frame(y = rep(c(0, 1), 5), x = seq(-2, 2, length.out = 10))
  fit <- stats::glm(y ~ x, data = d, family = stats::binomial())
  result <- .model_inference(list(original = fit, perturbed = fit), "binary")

  expect_identical(unique(result$by_imputation$imputation),
                   c("original", "perturbed"))
})

test_that("Rubin pooling refuses changed parameter names", {
  input <- rubin_inputs()
  changed <- input$q
  names(changed[[2L]]) <- c("a", "changed")
  expect_error(.rubin_pool(changed, input$u), "coefficient names.*imputation 2")

  reordered <- input$q
  reordered[[2L]] <- reordered[[2L]][2:1]
  expect_error(.rubin_pool(reordered, input$u), "coefficient names.*imputation 2")
})

test_that("Rubin pooling refuses invalid coefficient values", {
  input <- rubin_inputs()
  input$q[[2L]][[1L]] <- Inf
  expect_error(.rubin_pool(input$q, input$u), "non-finite coefficient.*imputation 2")
})

test_that("Rubin pooling errors use the declared imputation labels", {
  input <- rubin_inputs()
  names(input$q) <- c("10", "20")
  names(input$u) <- c("10", "20")
  input$q[[2L]][[1L]] <- Inf

  expect_error(.rubin_pool(input$q, input$u),
               "non-finite coefficient.*imputation 20")
})

test_that("Rubin pooling refuses invalid covariance shape and names", {
  input <- rubin_inputs()
  nonsquare <- input$u
  nonsquare[[2L]] <- matrix(1, 2, 3)
  expect_error(.rubin_pool(input$q, nonsquare), "covariance.*square.*imputation 2")

  renamed <- input$u
  dimnames(renamed[[2L]]) <- list(c("a", "changed"), c("a", "changed"))
  expect_error(.rubin_pool(input$q, renamed), "covariance names.*imputation 2")

  nonfinite <- input$u
  nonfinite[[2L]][1, 1] <- NA_real_
  expect_error(.rubin_pool(input$q, nonfinite), "non-finite covariance.*imputation 2")
})

test_that("binary model inference exposes coefficients and covariance", {
  d <- data.frame(y = rep(c(0, 1), 5), x = seq(-2, 2, length.out = 10))
  fit <- stats::glm(y ~ x, data = d, family = stats::binomial())
  result <- .model_inference(list("1" = fit), "binary")

  expect_identical(result$estimates$term, names(stats::coef(fit)))
  expect_equal(result$estimates$estimate, unname(stats::coef(fit)))
  expect_equal(result$covariance, stats::vcov(fit))
  expect_false(any(result$estimates$pooled))
  expect_true(all(is.infinite(result$estimates$df)))
})

test_that("ordinal model inference distinguishes thresholds", {
  skip_if_not_installed("MASS")
  d <- data.frame(
    y = ordered(rep(c("low", "mid", "high"), each = 8),
                levels = c("low", "mid", "high")),
    x = rep(seq(-1, 1, length.out = 8), 3) + rep(c(-0.2, 0, 0.2), each = 8)
  )
  fit <- MASS::polr(y ~ x, data = d, Hess = TRUE)
  result <- .model_inference(list("1" = fit), "ordinal")

  expect_identical(result$estimates$term,
                   c("x", "threshold:low|mid", "threshold:mid|high"))
  expect_true(all(is.na(result$estimates$odds_ratio[2:3])))
  expect_equal(dimnames(result$covariance),
               list(result$estimates$term, result$estimates$term))
})

test_that("nominal model inference uses outcome-qualified parameter names", {
  skip_if_not_installed("nnet")
  d <- data.frame(
    y = factor(rep(c("a", "b", "c"), each = 12), levels = c("a", "b", "c")),
    x = rep(seq(-1, 1, length.out = 12), 3) + rep(c(-0.3, 0, 0.3), each = 12)
  )
  fit <- nnet::multinom(y ~ x, data = d, Hess = TRUE, trace = FALSE)
  result <- .model_inference(list("1" = fit), "nominal")

  expect_identical(result$estimates$term,
                   c("b:(Intercept)", "b:x", "c:(Intercept)", "c:x"))
  expect_equal(dimnames(result$covariance),
               list(result$estimates$term, result$estimates$term))
})

test_that("model inference names a non-converged imputation", {
  d <- data.frame(y = rep(c(0, 1), 5), x = seq(-2, 2, length.out = 10))
  good <- stats::glm(y ~ x, data = d, family = stats::binomial())
  failed <- good
  failed$converged <- FALSE

  expect_error(
    .model_inference(list("1" = good, "2" = failed), "binary"),
    "failed to converge.*2"
  )
})
