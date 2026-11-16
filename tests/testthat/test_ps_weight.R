###############################################################################
## test_ps_weight.R
##
## Tests for ps_weight().
###############################################################################

dta <- sample_ps_data(n = 100, seed = 42)

# ---- Class and structure ----------------------------------------------------

test_that("ps_weight() returns a ps_weight / ps_data object", {
  obj <- ps_weight(dta)
  expect_s3_class(obj, "ps_weight")
  expect_s3_class(obj, "ps_data")
})

test_that("ps_weight() object has $data, $meta, $tables slots", {
  obj <- ps_weight(dta)
  expect_true(!is.null(obj$data))
  expect_true(!is.null(obj$meta))
  expect_true(!is.null(obj$tables))
})

# ---- $data ------------------------------------------------------------------

test_that("$data has the same number of rows as the input", {
  obj <- ps_weight(dta)
  expect_equal(nrow(obj$data), nrow(dta))
})

test_that("$data gains an 'iptw' column of positive numerics", {
  obj <- ps_weight(dta)
  expect_true("iptw" %in% names(obj$data))
  expect_true(is.numeric(obj$data$iptw))
  expect_true(all(obj$data$iptw > 0))
})

# ---- Estimands --------------------------------------------------------------

test_that("ATE, ATT, and ATC weights match the inverse-probability formulas", {
  ps  <- dta$prob_t
  trt <- dta$tavr == 1L
  expected <- list(
    ATE = ifelse(trt, 1 / ps, 1 / (1 - ps)),
    ATT = ifelse(trt, 1, ps / (1 - ps)),
    ATC = ifelse(trt, (1 - ps) / ps, 1)
  )
  for (est in names(expected)) {
    obj <- ps_weight(dta, estimand = est, stabilise = FALSE)
    expect_equal(obj$data$iptw, expected[[est]], info = paste("estimand =", est))
    expect_equal(obj$meta$estimand, est)
  }
})

# ---- Stabilisation ----------------------------------------------------------

test_that("stabilising the weights does not change the weighted SMDs", {
  # Stabilising multiplies each group's weights by a constant, which cannot
  # change balance (issue #82).
  expect_equal(ps_weight(dta, stabilise = TRUE)$tables$smd_weighted,
               ps_weight(dta, stabilise = FALSE)$tables$smd_weighted)
})

test_that("stabilised ATE weights are the raw weights times each arm's share", {
  p1   <- mean(dta$tavr == 1L)
  raw  <- ps_weight(dta, estimand = "ATE", stabilise = FALSE)$data$iptw
  stab <- ps_weight(dta, estimand = "ATE", stabilise = TRUE)$data$iptw
  expect_equal(stab, raw * ifelse(dta$tavr == 1L, p1, 1 - p1))
})

# ---- Trimming ---------------------------------------------------------------

test_that("trim winsorises weights to the expected quantile range", {
  obj <- ps_weight(dta, trim = 0.05)
  w   <- obj$data$iptw
  q_lo <- stats::quantile(ps_weight(dta, trim = NULL)$data$iptw, 0.05)
  q_hi <- stats::quantile(ps_weight(dta, trim = NULL)$data$iptw, 0.95)
  expect_true(all(w >= min(q_lo) - .Machine$double.eps * 10))
  expect_true(all(w <= max(q_hi) + .Machine$double.eps * 10))
})

test_that("trim = NULL does not trim anything", {
  p1 <- mean(dta$tavr == 1L)
  untrimmed <- ifelse(dta$tavr == 1L, 1 / dta$prob_t, 1 / (1 - dta$prob_t))
  expect_identical(ps_weight(dta, trim = NULL, stabilise = FALSE)$data$iptw, untrimmed)
  expect_identical(ps_weight(dta, trim = NULL, stabilise = TRUE)$data$iptw,
                   untrimmed * ifelse(dta$tavr == 1L, p1, 1 - p1))
  expect_null(ps_weight(dta, trim = NULL)$meta$trim)
})

test_that("invalid trim values raise an error", {
  expect_error(ps_weight(dta, trim = 0))
  expect_error(ps_weight(dta, trim = 0.5))
  expect_error(ps_weight(dta, trim = -0.1))
  expect_error(ps_weight(dta, trim = "a"))
})

# ---- $meta ------------------------------------------------------------------

test_that("$meta contains required keys", {
  obj <- ps_weight(dta)
  expect_true(all(c("treatment_col", "score_col", "weight_col", "estimand",
                    "stabilised", "trim", "method", "n_total") %in%
                    names(obj$meta)))
})

test_that("n_total equals nrow(data)", {
  obj <- ps_weight(dta)
  expect_equal(obj$meta$n_total, nrow(dta))
})

# ---- $tables ----------------------------------------------------------------

test_that("$tables contains SMD and effective N tables", {
  obj <- ps_weight(dta)
  expect_true(all(c("smd_unweighted", "smd_weighted",
                    "group_counts", "effective_n") %in%
                    names(obj$tables)))
})

test_that("effective_n has control and treated rows", {
  obj <- ps_weight(dta)
  expect_equal(nrow(obj$tables$effective_n), 2L)
  expect_true(all(c("group", "n_effective") %in%
                    names(obj$tables$effective_n)))
  expect_true(all(obj$tables$effective_n$n_effective > 0))
})

# ---- Custom weight_col ------------------------------------------------------

test_that("custom weight_col name is respected", {
  obj <- ps_weight(dta, weight_col = "mt_wt")
  expect_true("mt_wt" %in% names(obj$data))
})

# ---- Input validation -------------------------------------------------------

test_that("ps_weight() errors on non-data-frame input", {
  expect_error(ps_weight(list(tavr = 1, prob_t = 0.5)))
})

test_that("ps_weight() errors on missing columns", {
  expect_error(ps_weight(dta, treatment_col = "not_here"))
  expect_error(ps_weight(dta, score_col = "not_here"))
})

test_that("ps_weight() errors on non-binary treatment column", {
  bad <- dta
  bad$tavr <- 2L
  expect_error(ps_weight(bad))
})

test_that("ps_weight() errors on out-of-range propensity score", {
  bad <- dta
  bad$prob_t[1] <- -0.1
  expect_error(ps_weight(bad))
})

test_that("group tables name groups by their treatment value", {
  obj <- ps_weight(dta)
  expect_named(obj$tables$group_counts, c("group", "treated", "n"))
  expect_named(obj$tables$effective_n, c("group", "treated", "n_effective"))
  expect_identical(obj$tables$group_counts$group, c("0", "1"))
  expect_identical(obj$tables$effective_n$treated, c(FALSE, TRUE))

  logical_dta <- dta
  logical_dta$tavr <- logical_dta$tavr == 1
  expect_identical(ps_weight(logical_dta)$tables$group_counts$group, c("FALSE", "TRUE"))
})

test_that("default balance covariates leave out identifiers and package outputs", {
  d <- sample_ps_data(300, seed = 1)
  obj <- ps_weight(d)
  for (tbl in obj$tables[c("smd_unweighted", "smd_weighted")]) {
    expect_false(any(c("id", "match", "tavr", "prob_t", "iptw") %in% tbl$variable))
    expect_true(all(c("age", "female", "ef", "diabetes", "hypertension") %in% tbl$variable))
  }

  fit <- ps_logistic(tavr ~ age + female + ef + diabetes + hypertension,
                     d[setdiff(names(d), c("prob_t", "match"))])
  vars <- ps_weight(fit$data)$tables$smd_weighted$variable
  expect_false(any(c("id", "logit_t", "mt_wt", "quintile", "decile", "prob_t") %in% vars))

  d$site <- 7L
  expect_false("site" %in% ps_weight(d)$tables$smd_weighted$variable)

  d$patient <- d$id
  d$id <- NULL
  expect_false("patient" %in% ps_weight(d, id_col = "patient")$tables$smd_weighted$variable)
})

test_that("explicit covariates are used as given", {
  d <- sample_ps_data(300, seed = 1)
  obj <- ps_weight(d, covariates = c("id", "age"))
  expect_identical(obj$tables$smd_weighted$variable, c("id", "age"))
})

test_that("a covariate that only looks like a score column stays in the default set", {
  d <- sample_ps_data(300, seed = 1)
  d <- d[setdiff(names(d), c("prob_t", "match"))]
  d$prob_baseline_risk <- d$age / 100
  fit <- ps_logistic(tavr ~ prob_baseline_risk, d)
  expect_identical(fit$tables$smd$variable[fit$tables$smd$variable == "prob_baseline_risk"], "prob_baseline_risk")
  expect_true("prob_baseline_risk" %in% ps_weight(fit$data)$tables$smd_weighted$variable)
})

# ---- Trimming the reweighted arm only (QA tier 2, 2026-11-16) --------------

test_that("ATT trimming leaves every treated weight at exactly 1", {
  # About 10% treated: pooled quantiles put the treated weights of 1 in the
  # trimmed tail, which the arm-specific rule must not touch.
  d <- sample_ps_data(n = 400, seed = 11)
  set.seed(2)
  d$tavr <- as.integer(seq_len(nrow(d)) %in% sample(nrow(d), 40))
  d$prob_t <- stats::plogis(-2.2 + 0.03 * (d$age - 70) + stats::rnorm(nrow(d), 0, 0.5))
  obj <- ps_weight(d, estimand = "ATT", trim = 0.1)
  expect_identical(obj$data$iptw[d$tavr == 1L], rep(1, 40))
  ctl <- d$prob_t[d$tavr == 0L] / (1 - d$prob_t[d$tavr == 0L])
  q <- stats::quantile(ctl, c(0.1, 0.9), names = FALSE)
  expect_equal(obj$data$iptw[d$tavr == 0L], pmin(pmax(ctl, q[1]), q[2]))
})

test_that("ATC trimming leaves every control weight at exactly 1", {
  obj <- ps_weight(dta, estimand = "ATC", trim = 0.1)
  expect_identical(obj$data$iptw[dta$tavr == 0L], rep(1, sum(dta$tavr == 0L)))
})

test_that("a zero ATT weight is reported by patient", {
  d <- sample_ps_data(n = 60, seed = 2)
  ctl <- which(d$tavr == 0L)[1:2]
  d$prob_t[ctl] <- 0
  expect_error(ps_weight(d, estimand = "ATT"), paste0("`id` ", d$id[ctl[1]], ", ", d$id[ctl[2]]))
})

# ---- Non-numeric covariates in the balance tables ---------------------------

test_that("factor, character and logical covariates enter the default balance tables", {
  d <- dta
  d$site <- factor(rep(c("a", "b", "c"), length.out = nrow(d)))
  d$nyha <- rep(c("I", "II"), length.out = nrow(d))
  d$smoker <- rep(c(TRUE, FALSE, FALSE), length.out = nrow(d))
  d$name <- paste0("pt", seq_len(nrow(d)))  # an identifier, never a covariate
  obj <- ps_weight(d)
  expect_true(all(c("site", "nyha", "smoker") %in% obj$tables$smd_unweighted$variable))
  expect_false("name" %in% obj$tables$smd_unweighted$variable)
  expect_identical(obj$tables$smd_weighted$variable, obj$tables$smd_unweighted$variable)
  ref <- ps_stddiff(d, treatment_col = "tavr", categorical = "site", binary = "female")$tables$stddiff
  smd <- obj$tables$smd_unweighted
  expect_equal(smd$smd[smd$variable == "site"], round(ref$stddiff[ref$variable == "site"], 4))
  expect_equal(smd$smd[smd$variable == "female"], round(ref$stddiff[ref$variable == "female"], 4))
  expect_true(all(is.finite(smd$smd)))
})

test_that("a one-arm treatment column is refused", {
  d <- dta
  d$tavr <- 1L
  expect_error(ps_weight(d), "both treatment values")
  expect_error(ps_match(d), "both treatment values")
  expect_error(ps_logistic(tavr ~ age + ef, data = d), "both treatment values")
})
