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

test_that(".model_inference() labels exponentiated count coefficients rate_ratio", {
  dta <- data.frame(y = c(0, 1, 2, 4, 3, 1, 0, 5), x = c(1, 2, 3, 4, 5, 6, 7, 8))
  fit <- stats::glm(y ~ x, family = stats::poisson(link = "log"), data = dta)
  result <- .model_inference(list("1" = fit), "count")
  expect_false("odds_ratio" %in% names(result$estimates))
  expect_equal(result$estimates$rate_ratio, unname(exp(stats::coef(fit))))
  binary <- .model_inference(list("1" = stats::glm(I(y > 1) ~ x, family = stats::binomial(), data = dta)),
                             "binary")
  expect_true("odds_ratio" %in% names(binary$estimates))
  expect_false("rate_ratio" %in% names(binary$estimates))
})

test_that("a single linear fit is inferred on the t with residual degrees of freedom", {
  d <- sample_ps_data(n = 120, seed = 7)
  fit <- stats::lm(ef ~ age + female, data = d)
  out <- .model_inference(list("1" = fit), "linear")
  expected <- summary(fit)$coefficients
  expect_equal(out$estimates$estimate, unname(expected[, "Estimate"]))
  expect_equal(out$estimates$std.error, unname(expected[, "Std. Error"]))
  expect_equal(out$estimates$p.value, unname(expected[, "Pr(>|t|)"]))
  expect_equal(out$estimates$df, rep(fit$df.residual, 3L))
  expect_equal(out$estimates$conf.low, unname(stats::confint(fit)[, 1L]))
  expect_equal(out$estimates$conf.high, unname(stats::confint(fit)[, 2L]))
  expect_false("odds_ratio" %in% names(out$estimates))
  expect_false(any(out$estimates$pooled))
})

test_that("linear fits over imputations are Rubin pooled without an odds ratio", {
  d <- sample_ps_data(n = 120, seed = 7)
  second <- d
  second$age <- second$age + 0.5
  fits <- list("1" = stats::lm(ef ~ age, d), "2" = stats::lm(ef ~ age, second))
  out <- .model_inference(fits, "linear")
  pooled <- .rubin_pool(lapply(fits, stats::coef), lapply(fits, stats::vcov),
                        df_com = stats::df.residual(fits[[1L]]))
  expect_equal(out$estimates$estimate, pooled$estimates$estimate)
  expect_equal(out$estimates$df, pooled$estimates$df)
  expect_true(all(out$estimates$pooled))
  expect_false("odds_ratio" %in% names(out$estimates))
})

test_that("Barnard-Rubin df matches mice's formula and caps at the complete-data df", {
  input <- rubin_inputs()
  df_com <- 20
  result <- .rubin_pool(input$q, input$u, df_com = df_com)
  lambda <- c(a = 1.5 * 2 / 5, b = 1.5 * 2 / 6)
  df_old <- (2 - 1) / lambda^2
  df_obs <- (df_com + 1) / (df_com + 3) * df_com * (1 - lambda)
  expect_equal(result$estimates$df, unname(df_old * df_obs / (df_old + df_obs)))
  expect_true(all(result$estimates$df < df_com))

  # No between-imputation variance leaves the observed-data df alone.
  estimates <- list(c(a = 1), c(a = 1), c(a = 1))
  covariance <- matrix(2, 1, 1, dimnames = list("a", "a"))
  flat <- .rubin_pool(estimates, rep(list(covariance), 3L), df_com = df_com)
  expect_equal(flat$estimates$df, df_com * (df_com + 1) / (df_com + 3))
})

test_that("complete-data df is n minus the estimated parameters for every engine", {
  d <- sample_ps_data(n = 120, seed = 7)
  expect_identical(.complete_df(stats::lm(ef ~ age, d)), as.numeric(nrow(d) - 2))
  expect_identical(.complete_df(stats::glm(female ~ age, binomial, d)), as.numeric(nrow(d) - 2))
  expect_identical(.complete_df(list()), Inf)
  skip_if_not_installed("nnet")
  d$grp <- factor(rep_len(c("a", "b", "c"), nrow(d)))
  fit <- nnet::multinom(grp ~ age, d, trace = FALSE)
  expect_identical(.complete_df(fit), as.numeric(nrow(d) - 4))
})

test_that("pooling refuses zero or negative complete-data df", {
  input <- rubin_inputs()
  for (bad in list(0, -2, -0.5, NA_real_, c(10, 20), "20")) {
    expect_error(.rubin_pool(input$q, input$u, df_com = bad),
                 "positive complete-data degrees of freedom")
  }
  expect_no_error(.rubin_pool(input$q, input$u, df_com = Inf))

  # A saturated Poisson fit converges with finite covariance but no residual df.
  d <- data.frame(id = rep(1:3, 2), imp = rep(1:2, each = 3), x = rep(1:3, 2), y = rep(1:3, 2))
  expect_error(
    suppressWarnings(bs_count(y ~ x + I(x^2), d, imputation_col = "imp", dist = "poisson")),
    "the model has 0. A saturated model leaves none"
  )
})
