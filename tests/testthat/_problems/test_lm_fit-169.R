# Extracted from test_lm_fit.R:169

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
