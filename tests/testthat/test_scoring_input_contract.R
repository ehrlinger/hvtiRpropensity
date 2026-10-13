###############################################################################
## test_scoring_input_contract.R
##
## What the scoring functions accept as a treatment and an identifier (#87):
## ps_ordinal() needs a declared level order, ps_logistic() takes any two
## values when `treated_level` names one, and `id_col` is needed only with
## `imputation_col`.
###############################################################################

binary_unscored <- function(...) {
  d <- sample_ps_data(...)
  d$prob_t <- NULL
  d$match <- NULL
  d
}

# ---------------------------------------------------------------------------
# ps_ordinal(): the level order must be declared
# ---------------------------------------------------------------------------

test_that("ps_ordinal() refuses a character treatment without treatment_levels", {
  skip_if_not_installed("MASS")
  dta <- sample_ps_data_ordinal(n = 60, seed = 31)
  dta$nyha_grp <- as.character(dta$nyha_grp)
  expect_error(ps_ordinal(nyha_grp ~ age + ef, dta), "must be an ordered factor")
  expect_error(ps_ordinal(nyha_grp ~ age + ef, dta), "treatment_levels")
})

test_that("ps_ordinal() refuses an unordered factor without treatment_levels", {
  skip_if_not_installed("MASS")
  dta <- sample_ps_data_ordinal(n = 60, seed = 32)
  dta$nyha_grp <- factor(as.character(dta$nyha_grp))
  expect_error(ps_ordinal(nyha_grp ~ age + ef, dta), "must be an ordered factor")
})

test_that("ps_ordinal() takes the order of an ordered factor", {
  skip_if_not_installed("MASS")
  dta <- sample_ps_data_ordinal(n = 60, seed = 33)
  dta$nyha_grp <- ordered(as.character(dta$nyha_grp), levels = c("III", "II", "I"))
  obj <- ps_ordinal(nyha_grp ~ age + ef, dta)
  expect_identical(obj$meta$levels, c("III", "II", "I"))
})

test_that("ps_ordinal() fits a character treatment in the treatment_levels order", {
  skip_if_not_installed("MASS")
  dta <- sample_ps_data_ordinal(n = 60, seed = 34)
  ref <- ps_ordinal(nyha_grp ~ age + ef, dta)
  dta$nyha_grp <- as.character(dta$nyha_grp)
  obj <- ps_ordinal(nyha_grp ~ age + ef, dta, treatment_levels = c("I", "II", "III"))
  expect_identical(obj$meta$levels, c("I", "II", "III"))
  expect_equal(obj$data[c("prob_I", "prob_II", "prob_III")], ref$data[c("prob_I", "prob_II", "prob_III")])
})

# ---------------------------------------------------------------------------
# ps_logistic(): treated_level alone names the treated value
# ---------------------------------------------------------------------------

test_that("ps_logistic() infers two factor levels from treated_level", {
  binary <- binary_unscored(n = 150, seed = 35)
  labelled <- binary
  labelled$tavr <- factor(ifelse(binary$tavr == 1L, "TAVR", "SAVR"))
  ref <- ps_logistic(tavr ~ age + ef + diabetes, binary)
  obj <- ps_logistic(tavr ~ age + ef + diabetes, labelled, treated_level = "TAVR")
  expect_equal(obj$data$prob_t, ref$data$prob_t)
  expect_equal(obj$data$mt_wt, ref$data$mt_wt)
  expect_identical(obj$meta$treated_level, "TAVR")
  expect_setequal(obj$meta$treatment_levels, c("SAVR", "TAVR"))
})

test_that("ps_logistic() infers two character values from treated_level", {
  binary <- binary_unscored(n = 150, seed = 36)
  labelled <- binary
  labelled$tavr <- ifelse(binary$tavr == 1L, "TAVR", "SAVR")
  ref <- ps_logistic(tavr ~ age + ef, binary)
  obj <- ps_logistic(tavr ~ age + ef, labelled, treated_level = "TAVR")
  expect_equal(obj$data$prob_t, ref$data$prob_t)
})

test_that("ps_logistic() with treated_level still refuses more than two values", {
  dta <- binary_unscored(n = 150, seed = 37)
  dta$tavr <- rep(c("TAVR", "SAVR", "medical"), length.out = nrow(dta))
  expect_error(ps_logistic(tavr ~ age + ef, dta, treated_level = "TAVR"), "two treatment values")
})

test_that("ps_logistic() without treated_level still requires 0/1 or logical", {
  dta <- binary_unscored(n = 150, seed = 38)
  dta$tavr <- ifelse(dta$tavr == 1L, "TAVR", "SAVR")
  expect_error(ps_logistic(tavr ~ age + ef, dta), "must be binary")
})

# ---------------------------------------------------------------------------
# id_col: required only with imputation_col
# ---------------------------------------------------------------------------

test_that("the scoring functions run on data without an id column", {
  dta <- binary_unscored(n = 120, seed = 39)
  dta$id <- NULL
  expect_s3_class(ps_logistic(tavr ~ age + ef, dta), "ps_logistic")
  expect_s3_class(
    fit_logistic(tavr ~ age + ef, dta, family = "binary", outcome_levels = c(0, 1), event_level = 1),
    "ps_data"
  )

  nominal <- sample_ps_data_nominal(n = 60, seed = 40)
  nominal$id <- NULL
  expect_s3_class(ps_nominal(rtyp ~ age + ef, nominal), "ps_nominal")

  skip_if_not_installed("MASS")
  ordinal <- sample_ps_data_ordinal(n = 60, seed = 41)
  ordinal$id <- NULL
  expect_s3_class(ps_ordinal(nyha_grp ~ age + ef, ordinal), "ps_ordinal")
})

test_that("an id column is still required with imputation_col", {
  base <- binary_unscored(n = 60, seed = 42)
  stacked <- do.call(rbind, lapply(1:2, function(imputation) {
    partition <- base
    partition$imputation <- imputation
    partition
  }))
  stacked$id <- NULL
  expect_error(ps_logistic(tavr ~ age + ef, stacked, imputation_col = "imputation"), "id")
  expect_error(
    fit_logistic(tavr ~ age + ef, stacked, family = "binary", outcome_levels = c(0, 1), event_level = 1,
                 imputation_col = "imputation"),
    "id"
  )
  expect_error(
    ps_logistic(tavr ~ age + ef, stacked, id_col = NULL, imputation_col = "imputation"),
    "`id_col` is required"
  )
})
