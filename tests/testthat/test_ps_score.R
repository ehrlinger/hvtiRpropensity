###############################################################################
## test_ps_score.R
##
## Tests for ps_logistic(), ps_ordinal(), and ps_nominal().
###############################################################################

unscored_ps_data <- function(...) {
  d <- sample_ps_data(...)
  d$prob_t <- NULL
  d
}

# ---------------------------------------------------------------------------
# ps_logistic — binary treatment, single dataset
# ---------------------------------------------------------------------------

test_that("ps_logistic() returns correct class and slots", {
  dta <- unscored_ps_data(n = 100, seed = 1)
  obj <- ps_logistic(
    tavr ~ age + female + ef + diabetes + hypertension,
    data = dta
  )
  expect_s3_class(obj, "ps_logistic")
  expect_s3_class(obj, "ps_data")
  expect_true(is_ps_data(obj))
  expect_named(obj, c("data", "meta", "tables", "models"))
  expect_length(obj$models, 1L)
  expect_true(all(c("estimates", "covariance", "by_imputation", "fit_status") %in%
                    names(obj$tables)))
})

test_that("ps_logistic() appends score, logit, weight, quintile, decile", {
  dta <- unscored_ps_data(n = 100, seed = 2)
  obj <- ps_logistic(
    tavr ~ age + female + ef,
    data       = dta,
    score_col  = "ps",
    logit_col  = "logit_ps",
    weight_col = "wt"
  )
  expect_true(all(c("ps", "logit_ps", "wt", "quintile", "decile") %in%
                    names(obj$data)))
  expect_true(all(obj$data$ps >= 0 & obj$data$ps <= 1))
  expect_equal(range(obj$data$quintile), c(1L, 5L))
  expect_equal(range(obj$data$decile),   c(1L, 10L))
})

test_that("ps_logistic() matching weights are non-negative", {
  dta <- unscored_ps_data(n = 150, seed = 3)
  obj <- ps_logistic(tavr ~ age + ef, data = dta)
  expect_true(all(obj$data$mt_wt >= 0, na.rm = TRUE))
})

test_that("ps_logistic() meta fields are correctly populated", {
  dta <- unscored_ps_data(n = 80, seed = 4)
  obj <- ps_logistic(tavr ~ age + ef, data = dta)
  expect_equal(obj$meta$treatment_col,  "tavr")
  expect_equal(obj$meta$score_col,      "prob_t")
  expect_equal(obj$meta$logit_col,      "logit_t")
  expect_equal(obj$meta$weight_col,     "mt_wt")
  expect_equal(obj$meta$method,         "logistic")
  expect_equal(obj$meta$n_imputations,  1L)
  expect_equal(obj$meta$n_total,        nrow(dta))
})

test_that("ps_logistic() tables contain smd and group_counts", {
  dta <- unscored_ps_data(n = 100, seed = 5)
  obj <- ps_logistic(tavr ~ age + ef, data = dta)
  expect_named(obj$tables$smd, c("variable", "smd"))
  expect_named(obj$tables$group_counts, c("group", "n"))
  expect_equal(sum(obj$tables$group_counts$n), nrow(dta))
})

test_that("ps_logistic() treatment_col can be inferred from formula", {
  dta <- unscored_ps_data(n = 80, seed = 6)
  obj <- ps_logistic(tavr ~ age, data = dta)
  expect_equal(obj$meta$treatment_col, "tavr")
})

test_that("ps_logistic() errors on non-formula", {
  dta <- sample_ps_data(n = 50, seed = 7)
  expect_error(ps_logistic("tavr ~ age", data = dta))
})

test_that("ps_logistic() errors when treatment column is missing", {
  dta <- sample_ps_data(n = 50, seed = 8)
  expect_error(ps_logistic(no_such_col ~ age, data = dta))
})

# ---------------------------------------------------------------------------
# ps_logistic — multiply-imputed (stacked) data
# ---------------------------------------------------------------------------

test_that("ps_logistic() works with stacked MI data", {
  dta_mi <- sample_ps_data_count(n = 80, seed = 10, n_imputations = 3)
  # Build a binary treatment from diabetes and stack MI-style
  set.seed(99)
  base <- unscored_ps_data(n = 80, seed = 10)
  stacked <- do.call(rbind, lapply(1:3, function(m) {
    d                    <- base
    d$age                <- d$age + rnorm(nrow(d), 0, 0.5)  # jitter
    d[["_IMPUTATION_"]] <- m
    d
  }))

  obj <- ps_logistic(
    tavr ~ age + female + ef,
    data           = stacked,
    imputation_col = "_IMPUTATION_"
  )
  expect_s3_class(obj, "ps_logistic")
  expect_equal(obj$meta$n_imputations, 3L)
  expect_equal(obj$meta$method,        "logistic-MI")
  expect_equal(nrow(obj$data),         nrow(base))
  expect_length(obj$models, 3L)
  expect_true(all(obj$tables$estimates$pooled))
  # Imputation column should be gone from output
  expect_false("_IMPUTATION_" %in% names(obj$data))

  general <- fit_logistic(
    tavr ~ age + female + ef,
    stacked,
    id_col = "id",
    imputation_col = "_IMPUTATION_",
    outcome_levels = c(0, 1),
    event_level = 1,
    prediction_prefix = ".comparison"
  )
  expect_equal(obj$tables$estimates, general$tables$estimates)
  expect_equal(obj$tables$covariance, general$tables$covariance)
})

test_that("ps_logistic() honors explicit labelled treatment levels", {
  dta <- unscored_ps_data(n = 80, seed = 12)
  dta$treatment <- factor(
    ifelse(dta$tavr == 1L, "treated", "control"),
    levels = c("treated", "control")
  )
  obj <- ps_logistic(
    treatment ~ age + ef,
    dta,
    treatment_levels = c("control", "treated"),
    treated_level = "treated"
  )

  expect_identical(obj$meta$treatment_levels, c("control", "treated"))
  expect_identical(obj$meta$treated_level, "treated")
  expect_equal(obj$data$prob_t,
               unname(stats::predict(obj$models[[1L]], newdata = dta,
                                     type = "response")))
  expect_equal(obj$tables$group_counts$n, c(80L, 80L))
})

# ---------------------------------------------------------------------------
# ps_logistic — print
# ---------------------------------------------------------------------------

test_that("print.ps_logistic() produces output and returns invisibly", {
  dta <- unscored_ps_data(n = 60, seed = 11)
  obj <- ps_logistic(tavr ~ age + ef, data = dta)
  expect_output(out <- print(obj), "ps_logistic")
  expect_identical(out, obj)
})

# ---------------------------------------------------------------------------
# ps_ordinal — single dataset
# ---------------------------------------------------------------------------

test_that("ps_ordinal() returns correct class and slots", {
  skip_if_not_installed("MASS")
  dta <- sample_ps_data_ordinal(n = 80, seed = 20)
  obj <- ps_ordinal(
    nyha_grp ~ age + female + ef + diabetes,
    data = dta
  )
  expect_s3_class(obj, "ps_ordinal")
  expect_s3_class(obj, "ps_data")
  expect_named(obj, c("data", "meta", "tables", "models"))
  expect_length(obj$models, 1L)
  expect_true(all(c("estimates", "covariance", "by_imputation", "fit_status") %in%
                    names(obj$tables)))
})

test_that("ps_ordinal() appends one prob column per level", {
  skip_if_not_installed("MASS")
  dta <- sample_ps_data_ordinal(n = 80, seed = 21)
  obj <- ps_ordinal(nyha_grp ~ age + female + ef, data = dta)
  score_cols <- obj$meta$score_cols
  expect_length(score_cols, 3L)   # I, II, III
  expect_true(all(score_cols %in% names(obj$data)))
  # Probabilities sum to 1 (within floating-point tolerance)
  row_sums <- rowSums(obj$data[, score_cols])
  expect_true(all(abs(row_sums - 1) < 1e-6))
  # Each probability is in [0, 1]
  for (sc in score_cols)
    expect_true(all(obj$data[[sc]] >= 0 & obj$data[[sc]] <= 1))
})

test_that("ps_ordinal() meta levels match factor levels", {
  skip_if_not_installed("MASS")
  dta <- sample_ps_data_ordinal(n = 80, seed = 22)
  obj <- ps_ordinal(nyha_grp ~ age + ef, data = dta)
  expect_equal(obj$meta$levels, c("I", "II", "III"))
  expect_equal(obj$meta$treatment_col, "nyha_grp")
})

test_that("ps_ordinal() honors explicit levels over raw factor order", {
  skip_if_not_installed("MASS")
  dta <- sample_ps_data_ordinal(n = 80, seed = 24)
  dta$nyha_grp <- ordered(dta$nyha_grp, levels = c("III", "II", "I"))
  obj <- ps_ordinal(
    nyha_grp ~ age + ef,
    dta,
    treatment_levels = c("I", "II", "III")
  )

  expect_identical(obj$meta$levels, c("I", "II", "III"))
  expect_identical(levels(obj$data$nyha_grp), c("I", "II", "III"))
})

test_that("ps_ordinal() preserves its two-level interface", {
  skip_if_not_installed("MASS")
  dta <- unscored_ps_data(n = 120, seed = 25)
  dta$group <- ifelse(dta$tavr == 1L, "high", "low")
  obj <- ps_ordinal(
    group ~ age + ef,
    dta,
    treatment_levels = c("low", "high")
  )

  expect_identical(obj$meta$levels, c("low", "high"))
  expect_equal(rowSums(obj$data[c("prob_low", "prob_high")]), rep(1, nrow(dta)))
  expect_true(all(c("threshold:low|high", "age", "ef") %in% obj$tables$estimates$term))
})

test_that("print.ps_ordinal() runs without error", {
  skip_if_not_installed("MASS")
  dta <- sample_ps_data_ordinal(n = 80, seed = 23)
  obj <- ps_ordinal(nyha_grp ~ age + ef, data = dta)
  expect_output(out <- print(obj), "ps_ordinal")
  expect_identical(out, obj)
})

# ---------------------------------------------------------------------------
# ps_nominal — single dataset
# ---------------------------------------------------------------------------

test_that("ps_nominal() returns correct class and slots", {
  skip_if_not_installed("nnet")
  dta <- sample_ps_data_nominal(n = 60, seed = 30)
  obj <- ps_nominal(
    rtyp ~ age + female + ef + diabetes,
    data  = dta,
    trace = FALSE
  )
  expect_s3_class(obj, "ps_nominal")
  expect_s3_class(obj, "ps_data")
  expect_named(obj, c("data", "meta", "tables", "models"))
  expect_length(obj$models, 1L)
  expect_true(all(c("estimates", "covariance", "by_imputation", "fit_status") %in%
                    names(obj$tables)))
})

test_that("ps_nominal() appends one prob column per level", {
  skip_if_not_installed("nnet")
  dta <- sample_ps_data_nominal(n = 60, seed = 31)
  obj <- ps_nominal(rtyp ~ age + female + ef, data = dta, trace = FALSE)
  score_cols <- obj$meta$score_cols
  expect_length(score_cols, 4L)   # COS, PER, DEV, CE
  expect_true(all(score_cols %in% names(obj$data)))
  row_sums <- rowSums(obj$data[, score_cols])
  expect_true(all(abs(row_sums - 1) < 1e-5))
})

test_that("ps_nominal() ref_level is the first level by default", {
  skip_if_not_installed("nnet")
  dta <- sample_ps_data_nominal(n = 60, seed = 32)
  obj <- ps_nominal(rtyp ~ age + ef, data = dta, trace = FALSE)
  expect_equal(obj$meta$ref_level, "COS")
})

test_that("ps_nominal() honors explicit levels and reference", {
  skip_if_not_installed("nnet")
  dta <- sample_ps_data_nominal(n = 60, seed = 34)
  dta$rtyp <- factor(dta$rtyp, levels = rev(levels(dta$rtyp)))
  obj <- ps_nominal(
    rtyp ~ age + ef,
    dta,
    treatment_levels = c("COS", "PER", "DEV", "CE"),
    ref_level = "COS",
    trace = FALSE
  )

  expect_identical(obj$meta$levels, c("COS", "PER", "DEV", "CE"))
  expect_identical(obj$meta$ref_level, "COS")
  expect_identical(levels(obj$data$rtyp), c("COS", "PER", "DEV", "CE"))
})

test_that("ps_nominal() preserves its two-level interface", {
  skip_if_not_installed("nnet")
  dta <- unscored_ps_data(n = 120, seed = 35)
  dta$group <- ifelse(dta$tavr == 1L, "treated", "control")
  obj <- ps_nominal(
    group ~ age + ef,
    dta,
    treatment_levels = c("control", "treated"),
    ref_level = "control",
    trace = FALSE
  )

  expect_identical(obj$meta$levels, c("control", "treated"))
  expect_equal(rowSums(obj$data[c("prob_control", "prob_treated")]), rep(1, nrow(dta)))
})

test_that("propensity functions refuse to overwrite output columns", {
  binary <- sample_ps_data(n = 80, seed = 36)
  binary$prob_t <- 999
  expect_error(ps_logistic(tavr ~ age + ef, binary), "already exist.*prob_t")

  ordinal <- sample_ps_data_ordinal(n = 80, seed = 37)
  ordinal$prob_I <- 999
  expect_error(ps_ordinal(nyha_grp ~ age + ef, ordinal), "already exist.*prob_I")

  nominal <- sample_ps_data_nominal(n = 80, seed = 38)
  nominal$prob_COS <- 999
  expect_error(ps_nominal(rtyp ~ age + ef, nominal), "already exist.*prob_COS")
})

test_that("print.ps_nominal() runs without error", {
  skip_if_not_installed("nnet")
  dta <- sample_ps_data_nominal(n = 60, seed = 33)
  obj <- ps_nominal(rtyp ~ age + ef, data = dta, trace = FALSE)
  expect_output(out <- print(obj), "ps_nominal")
  expect_identical(out, obj)
})

# ---------------------------------------------------------------------------
# sample data generators
# ---------------------------------------------------------------------------

test_that("sample_ps_data_ordinal() returns expected structure", {
  dta <- sample_ps_data_ordinal(n = 50, seed = 1)
  expect_equal(nrow(dta), 150L)
  expect_true(is.ordered(dta$nyha_grp))
  expect_equal(levels(dta$nyha_grp), c("I", "II", "III"))
})

test_that("sample_ps_data_nominal() returns expected structure", {
  dta <- sample_ps_data_nominal(n = 50, seed = 1)
  expect_equal(nrow(dta), 200L)
  expect_true(is.factor(dta$rtyp) && !is.ordered(dta$rtyp))
  expect_equal(levels(dta$rtyp), c("COS", "PER", "DEV", "CE"))
})
