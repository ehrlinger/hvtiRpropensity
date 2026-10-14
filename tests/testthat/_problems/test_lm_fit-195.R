# Extracted from test_lm_fit.R:195

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "hvtiRpropensity", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
lm_fixture <- function() {
  utils::read.csv(testthat::test_path("fixtures", "lm-sas", "input.csv"))
}

# test -------------------------------------------------------------------------
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
