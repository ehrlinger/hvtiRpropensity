###############################################################################
## test_lm_sas_oracle.R
##
## Independent SAS 9.4M8 model and Rubin-pooling oracles.
###############################################################################

sas_lm_oracle <- function(file, check_names = TRUE) {
  utils::read.csv(
    testthat::test_path("fixtures", "lm-sas", file),
    check.names = check_names
  )
}

sas_lm_covariance <- function(file, row_name) {
  raw <- sas_lm_oracle(file, check_names = FALSE)
  values <- as.matrix(raw[setdiff(names(raw), row_name)])
  storage.mode(values) <- "double"
  rownames(values) <- raw[[row_name]]
  values
}

lm_sas_fixture <- function() {
  sas_lm_oracle("input.csv")
}

expect_sas_close <- function(actual, expected, tolerance) {
  testthat::expect_identical(length(actual), length(expected))
  testthat::expect_lt(max(abs(as.numeric(actual) - as.numeric(expected))), tolerance)
}

test_that("the controlled SAS oracle run is complete and clean", {
  expected <- c(
    paste0(rep(c("binary", "ordinal", "nominal", "count"), each = 3L),
           rep(c("-estimates.csv", "-covariance.csv", "-scores.csv"), 4L)),
    "mi-estimates.csv", "mi-covariance.csv", "mi-outest.csv", "mi-scores.csv",
    "pooled-estimates.csv", "pooled-covariance.csv", "sas-environment.csv"
  )
  fixture_dir <- testthat::test_path("fixtures", "lm-sas")
  expect_true(all(file.exists(file.path(fixture_dir, expected))))

  summary <- readLines(file.path(fixture_dir, "sas-run-summary.txt"), warn = FALSE)
  expect_true("errors=0" %in% summary)
  expect_true("warnings=0" %in% summary)
  expect_true("completion_marker=LM_ORACLE_COMPLETE expected_outputs=19" %in% summary)

  environment <- sas_lm_oracle("sas-environment.csv")
  expect_identical(environment$product_version, "9.04.01M8P02222023")
  expect_identical(environment$platform, "Linux")
})

test_that("binary logistic output agrees with the SAS oracle", {
  d <- lm_sas_fixture()
  one <- d[d$imputation == 1L, ]
  fit <- fit_logistic(
    binary ~ x1_imp + x2,
    one,
    family = "binary",
    outcome_levels = c(0, 1),
    event_level = 1
  )
  estimates <- sas_lm_oracle("binary-estimates.csv")
  terms <- sub("^Intercept$", "(Intercept)", estimates$Variable)
  covariance <- sas_lm_covariance("binary-covariance.csv", "Parameter")
  dimnames(covariance) <- lapply(dimnames(covariance), function(x) {
    sub("^Intercept$", "(Intercept)", x)
  })
  scores <- sas_lm_oracle("binary-scores.csv")
  scores <- scores[match(fit$data$id, scores$id), ]

  expect_sas_close(unname(stats::coef(fit$models[[1L]])[terms]),
                   estimates$Estimate, tolerance = 1e-4)
  expect_sas_close(fit$tables$covariance, covariance, tolerance = 1e-4)
  expect_sas_close(fit$data$prob, scores$probability, tolerance = 5e-5)
})

test_that("ordinal logistic output agrees with the SAS cumulative-logit convention", {
  d <- lm_sas_fixture()
  one <- d[d$imputation == 1L, ]
  fit <- fit_logistic(
    ordinal ~ x1_imp + x2,
    one,
    family = "ordinal",
    outcome_levels = c("low", "mid", "high")
  )
  estimates <- sas_lm_oracle("ordinal-estimates.csv")
  sas_estimate <- c(
    x1_imp = -estimates$Estimate[estimates$Variable == "x1_imp"],
    x2 = -estimates$Estimate[estimates$Variable == "x2"],
    `threshold:low|mid` = estimates$Estimate[
      estimates$Variable == "Intercept" & estimates$ClassVal0 == "low"
    ],
    `threshold:mid|high` = estimates$Estimate[
      estimates$Variable == "Intercept" & estimates$ClassVal0 == "mid"
    ]
  )
  covariance <- sas_lm_covariance("ordinal-covariance.csv", "Parameter")
  order <- c("x1_imp", "x2", "Intercept_low", "Intercept_mid")
  signs <- c(-1, -1, 1, 1)
  covariance <- covariance[order, order] * outer(signs, signs)
  scores <- sas_lm_oracle("ordinal-scores.csv")
  scores <- scores[match(fit$data$id, scores$id), ]

  expect_sas_close(fit$tables$estimates$estimate, sas_estimate, tolerance = 1e-4)
  expect_sas_close(fit$tables$covariance, covariance, tolerance = 4e-3)
  expect_sas_close(
    unname(as.matrix(fit$data[c("prob_low", "prob_mid", "prob_high")])),
    unname(as.matrix(scores[c("IP_low", "IP_mid", "IP_high")])),
    tolerance = 5e-5
  )
})

test_that("nominal generalized-logit output agrees with the SAS oracle", {
  d <- lm_sas_fixture()
  one <- d[d$imputation == 1L, ]
  fit <- fit_logistic(
    nominal ~ x1_imp + x2,
    one,
    family = "nominal",
    outcome_levels = c("a", "b", "c"),
    reference_level = "a"
  )
  estimates <- sas_lm_oracle("nominal-estimates.csv")
  terms <- paste0(
    estimates$Response,
    ":",
    ifelse(estimates$Variable == "Intercept", "(Intercept)", estimates$Variable)
  )
  sas_estimate <- stats::setNames(estimates$Estimate, terms)
  covariance <- sas_lm_covariance("nominal-covariance.csv", "Parameter")
  sas_covariance_terms <- c(
    "b:(Intercept)", "c:(Intercept)", "b:x1_imp",
    "c:x1_imp", "b:x2", "c:x2"
  )
  dimnames(covariance) <- list(sas_covariance_terms, sas_covariance_terms)
  r_terms <- fit$tables$estimates$term
  covariance <- covariance[r_terms, r_terms]
  scores <- sas_lm_oracle("nominal-scores.csv")
  scores <- scores[match(fit$data$id, scores$id), ]

  expect_sas_close(fit$tables$estimates$estimate,
                   sas_estimate[r_terms], tolerance = 1e-4)
  expect_sas_close(fit$tables$covariance, covariance, tolerance = 1e-4)
  expect_sas_close(
    unname(as.matrix(fit$data[c("prob_a", "prob_b", "prob_c")])),
    unname(as.matrix(scores[c("IP_a", "IP_b", "IP_c")])),
    tolerance = 5e-5
  )
})

test_that("Poisson balancing output agrees with the SAS oracle", {
  d <- lm_sas_fixture()
  one <- d[d$imputation == 1L, ]
  fit <- bs_count(count ~ x1_imp + x2, one, dist = "poisson")
  estimates <- sas_lm_oracle("count-estimates.csv")
  estimates <- estimates[estimates$Parameter != "Scale", ]
  covariance <- sas_lm_covariance("count-covariance.csv", "RowName")
  scores <- sas_lm_oracle("count-scores.csv")
  scores <- scores[match(fit$data$id, scores$id), ]

  expect_sas_close(fit$tables$estimates$estimate, estimates$Estimate,
                   tolerance = 1e-4)
  expect_sas_close(fit$tables$covariance, covariance, tolerance = 1e-4)
  expect_sas_close(fit$data$bs, scores$linear_predictor, tolerance = 5e-5)
})

test_that("Rubin pooled inference agrees with PROC MIANALYZE", {
  d <- lm_sas_fixture()
  fit <- fit_logistic(
    binary ~ x1_imp + x2,
    d,
    family = "binary",
    id_col = "id",
    imputation_col = "imputation",
    outcome_levels = c(0, 1),
    event_level = 1
  )
  estimates <- sas_lm_oracle("pooled-estimates.csv")
  covariance <- sas_lm_covariance("pooled-covariance.csv", "Variable")
  scores <- sas_lm_oracle("mi-scores.csv")
  averaged <- stats::aggregate(probability ~ id, scores, mean)
  averaged <- averaged[match(fit$data$id, averaged$id), ]

  expect_sas_close(fit$tables$estimates$estimate, estimates$Estimate,
                   tolerance = 1e-4)
  expect_sas_close(fit$tables$estimates$std.error, estimates$StdErr,
                   tolerance = 1e-4)
  expect_lt(max(abs(fit$tables$estimates$df - estimates$DF) / estimates$DF), 0.01)
  expect_sas_close(fit$tables$estimates$p.value, estimates$Probt,
                   tolerance = 1e-4)
  expect_sas_close(fit$tables$estimates$conf.low, estimates$LCLMean,
                   tolerance = 1e-4)
  expect_sas_close(fit$tables$estimates$conf.high, estimates$UCLMean,
                   tolerance = 1e-4)
  expect_sas_close(fit$tables$covariance, covariance, tolerance = 1e-4)
  expect_sas_close(fit$data$prob, averaged$probability, tolerance = 5e-5)
})
