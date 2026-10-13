###############################################################################
## test_rm_sas_oracle.R
##
## Independent SAS 9.4M8 PROC REG and PROC MIANALYZE oracles for fit_linear()
## and bs_continuous().
###############################################################################

rm_oracle <- function(file, ...) {
  utils::read.csv(testthat::test_path("fixtures", "rm-sas", file),
                  check.names = FALSE, ...)
}

rm_input <- function() {
  d <- rm_oracle("input.csv")
  names(d)[names(d) == "_Imputation_"] <- "imp"
  d
}

rm_close <- function(actual, expected, tolerance = 1e-8, relative = FALSE) {
  testthat::expect_identical(length(actual), length(expected))
  err <- abs(as.numeric(actual) - as.numeric(expected))
  if (relative) err <- err / abs(as.numeric(expected))
  testthat::expect_lt(max(err), tolerance)
}

rm_terms <- c("(Intercept)", "age", "female")

rm_matrix <- function(file) {
  raw <- rm_oracle(file)
  keep <- c("Intercept", "age", "female")
  m <- as.matrix(raw[keep])
  storage.mode(m) <- "double"
  unname(m)
}

test_that("the rm SAS oracle run is complete and clean", {
  fixture_dir <- testthat::test_path("fixtures", "rm-sas")
  expected <- c(
    "single-estimates.csv", "single-covariance.csv", "single-anova.csv",
    "single-scores.csv", "pooled-estimates.csv", "pooled-covariance.csv",
    "mi-scores.csv", "sas-environment.csv", "sas-run-summary.txt"
  )
  expect_true(all(file.exists(file.path(fixture_dir, expected))))
  summary <- readLines(file.path(fixture_dir, "sas-run-summary.txt"), warn = FALSE)
  expect_true("errors=0" %in% summary)
  expect_true("warnings=0" %in% summary)
  expect_true("completion_marker=RM_ORACLE_COMPLETE expected_outputs=8" %in% summary)
})

test_that("fit_linear reproduces PROC REG on one dataset", {
  d <- rm_input()
  fit <- fit_linear(y ~ age + female, d[d$imp == 1L, ])
  est <- fit$tables$estimates
  sas <- rm_oracle("single-estimates.csv")
  anova <- rm_oracle("single-anova.csv", na.strings = "_")
  residual_df <- anova$DF[anova$Source == "Error"]
  expect_identical(length(residual_df), 1L)
  expect_identical(est$term, rm_terms)

  rm_close(est$estimate, sas$Estimate)
  rm_close(est$std.error, sas$StdErr)
  rm_close(est$statistic, sas$tValue)
  rm_close(est$p.value, sas$Probt, relative = TRUE)
  expect_true(all(est$df == residual_df))

  scores <- rm_oracle("single-scores.csv")
  rm_close(fit$data$fitted, scores$fitted[match(fit$data$id, scores$id)])
  rm_close(unname(fit$tables$covariance), rm_matrix("single-covariance.csv"))
})

test_that("fit_linear's Rubin pool reproduces PROC MIANALYZE", {
  fit <- fit_linear(y ~ age + female, rm_input(), imputation_col = "imp")
  est <- fit$tables$estimates
  sas <- rm_oracle("pooled-estimates.csv")

  expect_identical(est$term, rm_terms)
  rm_close(est$statistic, sas$tValue, relative = TRUE)
  rm_close(est$estimate, sas$Estimate)
  rm_close(est$std.error, sas$StdErr)
  rm_close(est$conf.low, sas$LCLMean)
  rm_close(est$conf.high, sas$UCLMean)
  rm_close(est$df, sas$DF, relative = TRUE)
  rm_close(est$p.value, sas$Probt, relative = TRUE)

  # SAS's TCov is the multivariate total covariance under MIANALYZE's
  # proportionality assumption, (1 + rbar) * Ubar, not the elementwise Rubin
  # total T = Ubar + (1 + 1/m) * B that R reports. So compare R's diagonal to
  # the SAS standard errors, and rebuild TCov from the per-imputation fits.
  rm_close(diag(fit$tables$covariance), sas$StdErr^2)

  vcovs <- lapply(fit$models, stats::vcov)
  coefs <- do.call(rbind, lapply(fit$models, stats::coef))
  m <- length(vcovs)
  k <- ncol(coefs)
  expect_identical(colnames(coefs), rm_terms)
  ubar <- Reduce(`+`, vcovs) / m
  between <- stats::cov(coefs)
  rbar <- (1 + 1 / m) * sum(diag(between %*% solve(ubar))) / k
  rm_close(unname((1 + rbar) * ubar), rm_matrix("pooled-covariance.csv"))
})

test_that("bs_continuous's score reproduces the SAS averaged fitted value", {
  obj <- bs_continuous(y ~ age + female, rm_input(), imputation_col = "imp")
  sas <- rm_oracle("mi-scores.csv")
  rm_close(obj$data$bs, sas$fitted[match(obj$data$id, sas$id)])
})
