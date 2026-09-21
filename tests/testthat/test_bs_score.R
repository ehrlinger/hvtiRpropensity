###############################################################################
## test_bs_score.R
##
## Tests for bs_continuous() and bs_count().
###############################################################################

# ---------------------------------------------------------------------------
# bs_continuous — single dataset
# ---------------------------------------------------------------------------

test_that("bs_continuous() returns correct class and slots", {
  dta <- sample_ps_data(n = 150, seed = 1)
  obj <- bs_continuous(
    ef ~ age + female + diabetes + hypertension,
    data = dta
  )
  expect_s3_class(obj, "bs_continuous")
  expect_s3_class(obj, "ps_data")
  expect_true(is_ps_data(obj))
  expect_named(obj, c("data", "meta", "tables", "models"))
  expect_length(obj$models, 1L)
})

test_that("bs_continuous() appends score, quintile, decile, cluster", {
  dta <- sample_ps_data(n = 150, seed = 2)
  obj <- bs_continuous(ef ~ age + female + diabetes, data = dta)
  expect_true(all(c("bs", "quintile", "decile", "cluster") %in%
                    names(obj$data)))
  expect_equal(range(obj$data$quintile), c(1L, 5L))
  expect_equal(range(obj$data$decile),   c(1L, 10L))
  expect_equal(range(obj$data$cluster),  c(1L, 10L))
})

test_that("bs_continuous() binary strata columns are created by default", {
  dta <- sample_ps_data(n = 120, seed = 3)
  obj <- bs_continuous(ef ~ age + female, data = dta, n_strata = 5L)
  stra_cols <- paste0("stra_", 1:5)
  expect_true(all(stra_cols %in% names(obj$data)))
  # Each row's stra columns sum to 1
  row_sums <- rowSums(obj$data[, stra_cols])
  expect_true(all(row_sums == 1L))
})

test_that("bs_continuous() binary strata suppressed with create_binary_strata = FALSE", {
  dta <- sample_ps_data(n = 100, seed = 4)
  obj <- bs_continuous(ef ~ age + female,
                       data = dta, create_binary_strata = FALSE)
  expect_false(any(grepl("^stra_", names(obj$data))))
})

test_that("bs_continuous() meta fields are correct", {
  dta <- sample_ps_data(n = 100, seed = 5)
  obj <- bs_continuous(ef ~ age + female, data = dta)
  expect_equal(obj$meta$outcome_col,   "ef")
  expect_equal(obj$meta$score_col,     "bs")
  expect_equal(obj$meta$n_strata,      10L)
  expect_equal(obj$meta$strata_col,    "cluster")
  expect_equal(obj$meta$method,        "balancing-linear")
  expect_equal(obj$meta$n_imputations, 1L)
  expect_equal(obj$meta$n_total,       nrow(dta))
})

test_that("bs_continuous() strata_counts table is present", {
  dta <- sample_ps_data(n = 100, seed = 6)
  obj <- bs_continuous(ef ~ age + female, data = dta, n_strata = 5L)
  sc  <- obj$tables$strata_counts
  expect_named(sc, c("stratum", "n"))
  expect_equal(nrow(sc), 5L)
  expect_equal(sum(sc$n), nrow(dta))
})

test_that("bs_continuous() errors on bad n_strata", {
  dta <- sample_ps_data(n = 100, seed = 7)
  expect_error(bs_continuous(ef ~ age, data = dta, n_strata = 1L))
  expect_error(bs_continuous(ef ~ age, data = dta, n_strata = "ten"))
})

# ---------------------------------------------------------------------------
# bs_continuous — stacked MI data
# ---------------------------------------------------------------------------

test_that("bs_continuous() works with stacked MI data", {
  base <- sample_ps_data(n = 80, seed = 10)
  stacked <- do.call(rbind, lapply(1:3, function(m) {
    d                    <- base
    d$ef                 <- d$ef + rnorm(nrow(d), 0, 0.3)
    d[["_IMPUTATION_"]] <- m
    d
  }))
  obj <- bs_continuous(
    ef ~ age + female + diabetes,
    data           = stacked,
    imputation_col = "_IMPUTATION_"
  )
  expect_s3_class(obj, "bs_continuous")
  expect_equal(obj$meta$n_imputations, 3L)
  expect_equal(obj$meta$method,        "balancing-linear-MI")
  expect_equal(nrow(obj$data),         nrow(base))
  expect_false("_IMPUTATION_" %in% names(obj$data))
  expect_length(obj$models, 3L)
})

# ---------------------------------------------------------------------------
# bs_continuous — print
# ---------------------------------------------------------------------------

test_that("print.bs_continuous() produces output and returns invisibly", {
  dta <- sample_ps_data(n = 80, seed = 11)
  obj <- bs_continuous(ef ~ age, data = dta)
  expect_output(out <- print(obj), "bs_continuous")
  expect_identical(out, obj)
})

# ---------------------------------------------------------------------------
# bs_count — single dataset (Poisson, no MASS required)
# ---------------------------------------------------------------------------

test_that("bs_count() with dist='poisson' returns correct class", {
  dta <- sample_ps_data_count(n = 150, seed = 40)
  obj <- bs_count(
    rbc_tot ~ age + female + diabetes + hct_pr,
    data = dta,
    dist = "poisson"
  )
  expect_s3_class(obj, "bs_count")
  expect_s3_class(obj, "ps_data")
  expect_true(is_ps_data(obj))
  expect_length(obj$models, 1L)
  expect_named(
    obj$tables$estimates,
    c("term", "estimate", "std.error", "statistic", "df", "p.value",
      "conf.low", "conf.high", "odds_ratio", "pooled")
  )
  expect_named(obj$tables,
               c("strata_counts", "estimates", "covariance",
                 "by_imputation", "fit_status"))
})

test_that("bs_count() appends score, quintile, decile, cluster", {
  dta <- sample_ps_data_count(n = 150, seed = 41)
  obj <- bs_count(rbc_tot ~ age + female, data = dta, dist = "poisson")
  expect_true(all(c("bs", "quintile", "decile", "cluster") %in%
                    names(obj$data)))
  expect_equal(range(obj$data$quintile), c(1L, 5L))
  expect_equal(range(obj$data$cluster),  c(1L, 10L))
})

test_that("bs_count() binary strata sum to 1 per row", {
  dta <- sample_ps_data_count(n = 100, seed = 42)
  obj <- bs_count(rbc_tot ~ age + female, data = dta,
                  dist = "poisson", n_strata = 5L)
  stra_cols <- paste0("stra_", 1:5)
  row_sums <- rowSums(obj$data[, stra_cols])
  expect_true(all(row_sums == 1L))
})

test_that("bs_count() meta fields are correct for Poisson", {
  dta <- sample_ps_data_count(n = 100, seed = 43)
  obj <- bs_count(rbc_tot ~ age + female, data = dta, dist = "poisson")
  expect_equal(obj$meta$outcome_col,   "rbc_tot")
  expect_equal(obj$meta$score_col,     "bs")
  expect_equal(obj$meta$dist,          "poisson")
  expect_equal(obj$meta$n_strata,      10L)
  expect_equal(obj$meta$method,        "balancing-poisson")
  expect_equal(obj$meta$n_total,       nrow(dta))
})

test_that("bs_count() strata_counts sums to n", {
  dta <- sample_ps_data_count(n = 100, seed = 44)
  obj <- bs_count(rbc_tot ~ age + female, data = dta, dist = "poisson")
  expect_equal(sum(obj$tables$strata_counts$n), nrow(dta))
})

# ---------------------------------------------------------------------------
# bs_count — negbin (requires MASS)
# ---------------------------------------------------------------------------

test_that("bs_count() with dist='negbin' returns correct method string", {
  skip_if_not_installed("MASS")
  dta <- sample_ps_data_count(n = 150, seed = 50)
  obj <- bs_count(
    rbc_tot ~ age + female + diabetes + hct_pr,
    data = dta,
    dist = "negbin"
  )
  expect_s3_class(obj, "bs_count")
  expect_equal(obj$meta$dist,   "negbin")
  expect_equal(obj$meta$method, "balancing-negbin")
  expect_length(obj$meta$theta, 1L)
  expect_false(any(grepl("theta", obj$tables$estimates$term,
                         ignore.case = TRUE)))
  expect_false(any(grepl("theta", rownames(obj$tables$covariance),
                         ignore.case = TRUE)))
})

# ---------------------------------------------------------------------------
# bs_count — stacked MI data
# ---------------------------------------------------------------------------

test_that("bs_count() works with stacked MI data", {
  dta_mi <- sample_ps_data_count(n = 80, seed = 60, n_imputations = 3)
  obj <- bs_count(
    rbc_tot ~ age + female + diabetes,
    data           = dta_mi,
    imputation_col = "_IMPUTATION_",
    dist           = "poisson"
  )
  expect_equal(obj$meta$n_imputations, 3L)
  expect_equal(obj$meta$method,        "balancing-poisson-MI")
  expect_equal(nrow(obj$data),         80L)
  expect_false("_IMPUTATION_" %in% names(obj$data))
  expect_length(obj$models, 3L)
  expect_true(all(obj$tables$estimates$pooled))
  expect_identical(unique(obj$tables$by_imputation$imputation),
                   c("1", "2", "3"))
})

test_that("bs_count() matches direct Poisson inference and link scores", {
  dta <- sample_ps_data_count(n = 100, seed = 62)
  obj <- bs_count(rbc_tot ~ age + female, dta, dist = "poisson")
  direct <- stats::glm(
    rbc_tot ~ age + female,
    dta,
    family = stats::poisson(link = "log"),
    na.action = stats::na.exclude
  )

  expect_equal(obj$tables$estimates$estimate,
               unname(stats::coef(direct)), tolerance = 1e-10)
  expect_equal(obj$tables$covariance, stats::vcov(direct), tolerance = 1e-10)
  expect_equal(obj$data$bs,
               unname(stats::predict(direct, newdata = dta, type = "link")),
               tolerance = 1e-10)
})

test_that("bs_count() rejects an unknown distribution", {
  dta <- sample_ps_data_count(n = 50, seed = 63)
  expect_error(bs_count(rbc_tot ~ age, dta, dist = "auto"),
               "arg.*negbin.*poisson")
})

# ---------------------------------------------------------------------------
# bs_count — print
# ---------------------------------------------------------------------------

test_that("print.bs_count() produces output and returns invisibly", {
  dta <- sample_ps_data_count(n = 80, seed = 61)
  obj <- bs_count(rbc_tot ~ age + female, data = dta, dist = "poisson")
  expect_output(out <- print(obj), "bs_count")
  expect_identical(out, obj)
})

# ---------------------------------------------------------------------------
# sample_ps_data_count
# ---------------------------------------------------------------------------

test_that("sample_ps_data_count() returns expected structure (single)", {
  dta <- sample_ps_data_count(n = 100, seed = 1)
  expect_equal(nrow(dta), 100L)
  expect_true(all(c("id", "rbc_tot", "age", "female", "diabetes",
                    "hypertension", "hct_pr", "cpb_time") %in% names(dta)))
  expect_true(all(dta$rbc_tot >= 0L))
  expect_false("_IMPUTATION_" %in% names(dta))
})

test_that("sample_ps_data_count() returns stacked MI format", {
  dta <- sample_ps_data_count(n = 100, seed = 1, n_imputations = 4)
  expect_equal(nrow(dta), 400L)
  expect_true("_IMPUTATION_" %in% names(dta))
  expect_equal(sort(unique(dta[["_IMPUTATION_"]])), 1:4)
})
