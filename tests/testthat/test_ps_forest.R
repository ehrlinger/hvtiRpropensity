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

test_that("ps_forest() writes Li-Greene matching weights to mt_wt", {
  withr::local_options(rf.cores = 2L)
  fit <- ps_forest(tavr ~ age + ef, dta, ntree = 60, seed = 7)
  p   <- fit$data$prob_t
  trt <- fit$data$tavr == 1
  expect_true(any(trt) && any(!trt))
  expect_equal(fit$data$mt_wt[trt],  pmin(p[trt], 1 - p[trt]) / p[trt])
  expect_equal(fit$data$mt_wt[!trt], pmin(p[!trt], 1 - p[!trt]) / (1 - p[!trt]))
  expect_true(all(fit$data$mt_wt > 0 & fit$data$mt_wt <= 1))
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
                 weights = c("unweighted", "ato"), subsets = list(half = surv_dta$age > 60), n_boot = 25,
                 seed = 1024)
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
                   weights = c("unweighted", "ato"), subsets = list(half = surv_dta$age > 60), n_boot = 25,
                   seed = 1024)
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

test_that("ps_forest() requires treatment_col to match the formula response", {
  expect_no_error(ps_forest(tavr ~ age + ef, dta, treatment_col = "tavr", ntree = 20, seed = 1))
  expect_error(ps_forest(tavr ~ age + ef, dta, treatment_col = "female", ntree = 20), "must match")
})

test_that("degenerate inputs give NA, not errors or NaN", {
  expect_true(is.na(.weighted_km(numeric(0), numeric(0), numeric(0), 3)$rmst))
  ag <- .support_agreement(c(NA, NA), c(TRUE, NA))
  expect_equal(ag$n, 0L)
  expect_true(is.na(ag$pct_agree) && is.na(ag$kappa) && is.na(ag$jaccard))
  agr <- ps_support(fo)$tables$agreement
  expect_false("stringsAsFactors" %in% names(agr))
  expect_true(all(c("rule_a", "rule_b") %in% names(agr)))
})

test_that("ps_rmst() reports NA intervals when a subset has only one arm", {
  x <- ps_logistic(tavr ~ age + ef, surv_dta)
  only_ctl <- surv_dta$tavr == 0
  r <- expect_no_warning(ps_rmst(x, "t", "e", tau = 4, weights = "unweighted",
                                 subsets = list(ctl = only_ctl), n_boot = 5))
  row <- r$tables$estimates[r$tables$estimates$estimator == "ctl/unweighted", ]
  expect_true(is.na(row$diff_days) && is.na(row$lo_days) && is.na(row$hi_days))
  expect_equal(row$n_failed, 5L)
})

test_that("ps_forest() group_counts names groups by their treatment value", {
  counts <- fo$tables$group_counts
  expect_named(counts, c("group", "treated", "n"))
  expect_identical(counts$group, c("0", "1"))
  expect_identical(counts$treated, c(FALSE, TRUE))
  expect_equal(counts$n, c(sum(dta$tavr == 0), sum(dta$tavr == 1)))
})

test_that("ps_forest expands a dot formula without id and drops subtracted terms", {
  d <- sample_ps_data(n = 100, seed = 1)[c("id", "tavr", "age", "female", "ef")]
  explicit <- ps_forest(tavr ~ age + female + ef, d, ntree = 50, seed = 1)
  dot <- ps_forest(tavr ~ ., d, ntree = 50, seed = 1)
  expect_identical(dot$meta$formula, tavr ~ age + female + ef, ignore_formula_env = TRUE)
  expect_equal(dot$data$prob_t, explicit$data$prob_t)

  # A subtracted column is neither fitted nor checked for missing values.
  d$dropme <- c(NA, rep(1, nrow(d) - 1L))
  minus <- ps_forest(tavr ~ . - dropme, d, ntree = 50, seed = 1)
  expect_equal(minus$data$prob_t, explicit$data$prob_t)
  expect_identical(minus$tables$smd$variable, c("age", "female", "ef"))
})

test_that("ps_forest keeps offset() terms, which rfsrc fits as predictors", {
  d <- sample_ps_data(n = 100, seed = 1)[c("id", "tavr", "age", "female", "ef")]
  mixed <- ps_forest(tavr ~ age + offset(ef), d, ntree = 50, seed = 1)
  expect_identical(mixed$meta$formula, tavr ~ age + offset(ef), ignore_formula_env = TRUE)
  expect_identical(mixed$tables$smd$variable, c("age", "ef"))

  only <- ps_forest(tavr ~ offset(ef), d, ntree = 50, seed = 1)
  expect_identical(only$meta$formula, tavr ~ offset(ef), ignore_formula_env = TRUE)
  expect_identical(only$tables$smd$variable, "ef")

  dot <- ps_forest(tavr ~ . - female + offset(log(ef)), d, ntree = 50, seed = 1)
  expect_identical(dot$meta$formula, tavr ~ age + ef + offset(log(ef)), ignore_formula_env = TRUE)
})

test_that("ps_forest refuses a subtracted column that is not in data", {
  d <- sample_ps_data(n = 100, seed = 1)[c("id", "tavr", "age", "female", "ef")]
  expect_error(ps_forest(tavr ~ . - eff, d, ntree = 50, seed = 1),
               "Required column(s) not found in `data`: eff", fixed = TRUE)
  expect_error(ps_forest(tavr ~ age - eff, d, ntree = 50, seed = 1),
               "Required column(s) not found in `data`: eff", fixed = TRUE)
})
