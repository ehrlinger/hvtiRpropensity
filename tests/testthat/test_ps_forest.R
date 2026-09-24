###############################################################################
## test_ps_forest.R -- ps_forest(), ps_support(), ps_rmst()
###############################################################################

skip_if_not_installed("randomForestSRC")
skip_if_not_installed("survival")

raw_data <- function(n, seed) sample_ps_data(n = n, seed = seed)[, c("id", "tavr", "age", "female", "ef")]
dta <- raw_data(150, 42)
fo  <- ps_forest(tavr ~ age + ef, dta, ntree = 60, seed = 7)

test_that("ps_forest() follows the ps_data contract", {
  expect_s3_class(fo, "ps_forest")
  expect_s3_class(fo, "ps_data")
  expect_true(all(c("prob_t", "logit_t", "mt_wt", "quintile", "decile") %in% names(fo$data)))
  expect_true(all(fo$data$prob_t > 0 & fo$data$prob_t < 1))
  expect_equal(fo$meta$method, "forest-oob")
  expect_equal(fo$meta$score_col, "prob_t")
  expect_equal(nrow(fo$data), nrow(dta))
})

test_that("ps_forest() feeds ps_weight() and sa_overlap()", {
  w <- ps_weight(fo$data, treatment_col = "tavr", score_col = "prob_t")
  expect_s3_class(w, "ps_weight")
  expect_no_error(sa_overlap(fo))
})

test_that("ps_forest() is reproducible with a seed and rejects missing values", {
  again <- ps_forest(tavr ~ age + ef, dta, ntree = 60, seed = 7)
  expect_equal(again$data$prob_t, fo$data$prob_t)
  bad <- dta
  bad$age[1] <- NA
  expect_error(ps_forest(tavr ~ age + ef, bad, ntree = 20), "missing")
})

test_that("ps_support() flags rules and measures agreement", {
  sup <- ps_support(fo, flags = list(iso = seq_len(nrow(dta)) <= 15))
  expect_s3_class(sup, "ps_support")
  expect_true(all(c("weak_common", "weak_trim", "weak_iso") %in% names(sup$data)))
  expect_equal(nrow(sup$tables$agreement), 3L)
  ag <- sup$tables$agreement
  expect_true(all(ag$n_both + ag$n_neither <= ag$n))
  expect_error(ps_support(fo, flags = list(bad = c(TRUE, FALSE))), "one value per row")
  expect_error(ps_support(fo, flags = list(TRUE)), "named")
})

test_that(".support_agreement() gives known kappa and jaccard", {
  a <- c(rep(TRUE, 10), rep(TRUE, 5), rep(FALSE, 5), rep(FALSE, 80))
  b <- c(rep(TRUE, 10), rep(FALSE, 5), rep(TRUE, 5), rep(FALSE, 80))
  ag <- .support_agreement(a, b)
  pe <- 0.15 * 0.15 + 0.85 * 0.85
  expect_equal(ag$kappa, (0.9 - pe) / (1 - pe))
  expect_equal(ag$jaccard, 0.5)
  expect_true(is.na(.support_agreement(!a & FALSE, !a & FALSE)$kappa))
})

test_that(".rmst_step() matches a hand-computed curve", {
  # events at 1 (x2) and 2 (x1) of 4 people, tau = 3: 1 + 0.5 + 0.25
  km <- .weighted_km(c(1, 1, 2, 5), c(1, 1, 1, 0), rep(1, 4), 3)
  expect_equal(km$rmst, 1.75)
  expect_equal(.weighted_km(c(1, 1, 2, 5), c(1, 1, 1, 0), rep(2, 4), 3)$rmst, 1.75)
  expect_lt(.weighted_km(c(1, 1, 2, 5), c(1, 1, 1, 0), c(5, 1, 1, 1), 3)$rmst, 1.75)
})

surv_dta <- local({
  set.seed(9)
  d <- raw_data(200, 9)
  d$t <- stats::rexp(nrow(d), 0.4 * exp(-1.5 * d$tavr))
  d$e <- as.integer(d$t < 5)
  d$t <- pmin(d$t, 5)
  d
})

test_that("ps_rmst() returns estimates and curves with the documented columns", {
  res <- ps_rmst(ps_logistic(tavr ~ age + ef, surv_dta), "t", "e", tau = 4,
                 weights = c("unweighted", "ato"), subsets = list(half = surv_dta$age > 60), n_boot = 25)
  expect_s3_class(res, "ps_rmst")
  est <- res$tables$estimates
  expect_setequal(est$estimator, c("all/unweighted", "all/ato", "half/unweighted", "half/ato"))
  expect_true(all(c("diff_days", "lo_days", "hi_days", "ess_treated", "ess_control", "n_failed") %in% names(est)))
  expect_equal(est$diff_days, est$diff * 365.2425)
  expect_true(all(est$lo_days <= est$hi_days))
  expect_true(all(c("estimator", "arm", "time", "surv") %in% names(res$tables$curves)))
  # treatment lowers hazard in the simulation, so unweighted contrast favours treated
  expect_gt(est$diff_days[est$estimator == "all/unweighted"], 0)
  # deterministic given seed
  again <- ps_rmst(ps_logistic(tavr ~ age + ef, surv_dta), "t", "e", tau = 4,
                   weights = c("unweighted", "ato"), subsets = list(half = surv_dta$age > 60), n_boot = 25)
  expect_equal(again$tables$estimates, est)
})

test_that("ps_rmst() supports refit and validates input", {
  x <- ps_logistic(tavr ~ age + ef, surv_dta)
  r <- ps_rmst(x, "t", "e", tau = 4, weights = "ato", n_boot = 8, refit = TRUE)
  expect_equal(r$tables$estimates$n_failed, 0L)
  rf <- ps_rmst(fo, "age", "tavr", tau = 50, weights = "unweighted", n_boot = 3)
  expect_s3_class(rf, "ps_rmst")
  expect_error(ps_rmst(x, "t", "e", tau = -1), "positive")
  expect_error(ps_rmst(x, "nope", "e", tau = 4), "nope|not found|missing", ignore.case = TRUE)
  expect_error(ps_rmst(x, "t", "e", tau = 4, subsets = list(a = TRUE)), "one value per row")
  expect_error(ps_rmst(dta, "t", "e", tau = 4), "ps_data")
})
