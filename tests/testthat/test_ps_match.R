###############################################################################
## test_ps_match.R
##
## Tests for ps_match().
###############################################################################

dta <- sample_ps_data(n = 100, seed = 42)

# ---- Class and structure ----------------------------------------------------

test_that("ps_match() returns a ps_match / ps_data object", {
  obj <- ps_match(dta)
  expect_s3_class(obj, "ps_match")
  expect_s3_class(obj, "ps_data")
})

test_that("ps_match() object has $data, $meta, $tables slots", {
  obj <- ps_match(dta)
  expect_true(!is.null(obj$data))
  expect_true(!is.null(obj$meta))
  expect_true(!is.null(obj$tables))
})

# ---- $data ------------------------------------------------------------------

test_that("$data has the same number of rows as the input", {
  obj <- ps_match(dta)
  expect_equal(nrow(obj$data), nrow(dta))
})

test_that("$data gains a 'match' column with values 0 or 1", {
  obj <- ps_match(dta)
  expect_true("match" %in% names(obj$data))
  expect_true(all(obj$data$match %in% c(0L, 1L)))
})

test_that("matched subset has equal group sizes", {
  obj <- ps_match(dta)
  matched <- obj$data[obj$data$match == 1L, ]
  expect_equal(sum(matched$tavr == 1L), sum(matched$tavr == 0L))
})

# ---- $meta ------------------------------------------------------------------

test_that("$meta contains required keys", {
  obj <- ps_match(dta)
  expect_true(all(c("treatment_col", "score_col", "id_col", "match_col",
                    "method", "n_total", "n_matched", "n_unmatched") %in%
                    names(obj$meta)))
})

test_that("n_total equals nrow(data)", {
  obj <- ps_match(dta)
  expect_equal(obj$meta$n_total, nrow(dta))
})

test_that("n_matched + n_unmatched equals number of treated patients", {
  obj <- ps_match(dta)
  n_trt <- sum(dta$tavr == 1L)
  expect_equal(obj$meta$n_matched + obj$meta$n_unmatched, n_trt)
})

# ---- $tables ----------------------------------------------------------------

test_that("$tables contains smd_before and smd_after", {
  obj <- ps_match(dta)
  expect_true(all(c("smd_before", "smd_after",
                    "group_counts_before", "group_counts_after") %in%
                    names(obj$tables)))
})

test_that("SMD tables have 'variable' and 'smd' columns", {
  obj <- ps_match(dta)
  expect_true(all(c("variable", "smd") %in% names(obj$tables$smd_before)))
  expect_true(all(c("variable", "smd") %in% names(obj$tables$smd_after)))
})

# ---- Caliper ----------------------------------------------------------------

test_that("caliper reduces or maintains number of matched pairs", {
  obj_no_cal <- ps_match(dta)
  obj_cal    <- ps_match(dta, caliper = 0.01)
  expect_lte(obj_cal$meta$n_matched, obj_no_cal$meta$n_matched)
})

test_that("very tight caliper can match zero pairs without error", {
  expect_no_error(obj <- ps_match(dta, caliper = 1e-10))
  expect_gte(obj$meta$n_matched, 0L)
})

# ---- Custom column names ----------------------------------------------------

test_that("custom match_col name is respected", {
  obj <- ps_match(dta, match_col = "is_matched")
  expect_true("is_matched" %in% names(obj$data))
})

# ---- Input validation -------------------------------------------------------

test_that("ps_match() errors on non-data-frame input", {
  expect_error(ps_match(list(tavr = 1, prob_t = 0.5, id = 1)))
})

test_that("ps_match() errors on missing treatment column", {
  expect_error(ps_match(dta, treatment_col = "missing_col"))
})

test_that("ps_match() errors on non-binary treatment column", {
  bad <- dta
  bad$tavr <- 2L
  expect_error(ps_match(bad))
})

test_that("ps_match() errors on out-of-range propensity score", {
  bad <- dta
  bad$prob_t[1] <- 1.5
  expect_error(ps_match(bad))
})

test_that("ps_match() errors on invalid caliper", {
  expect_error(ps_match(dta, caliper = -0.1))
  expect_error(ps_match(dta, caliper = "a"))
})

test_that("ps_match() errors when ratio != 1", {
  expect_error(ps_match(dta, ratio = 2))
})

test_that("group count tables name groups by their treatment value", {
  obj <- ps_match(dta)
  for (tbl in obj$tables[c("group_counts_before", "group_counts_after")]) {
    expect_named(tbl, c("group", "treated", "n"))
    expect_identical(tbl$group, c("0", "1"))
    expect_identical(tbl$treated, c(FALSE, TRUE))
  }
  expect_equal(obj$tables$group_counts_before$n, c(sum(dta$tavr == 0), sum(dta$tavr == 1)))
})

# ---- Logit scale ------------------------------------------------------------

test_that("scale = 'logit' gives the same pairs as the linear-rescale workaround", {
  d <- sample_ps_data(n = 300, seed = 42)
  d$logit_t <- stats::qlogis(d$prob_t)
  rng <- diff(range(d$logit_t))
  d$logit_01 <- (d$logit_t - min(d$logit_t)) / rng
  cal <- 0.2 * stats::sd(d$logit_t)

  obj_logit  <- ps_match(d, score_col = "logit_t", scale = "logit", caliper = cal, seed = 101)
  obj_manual <- ps_match(d, score_col = "logit_01", caliper = cal / rng, seed = 101)

  expect_identical(obj_logit$data$match, obj_manual$data$match)
  expect_identical(obj_logit$data$pair_id, obj_manual$data$pair_id)
  expect_lt(sum(obj_logit$data$match), nrow(d))  # the caliper actually binds
})

test_that("caliper_sd is pooled SD units of score_col on the chosen scale", {
  d <- sample_ps_data(n = 300, seed = 42)
  d$logit_t <- stats::qlogis(d$prob_t)
  pooled_sd <- sqrt((stats::var(d$logit_t[d$tavr == 1L]) + stats::var(d$logit_t[d$tavr == 0L])) / 2)
  obj_sd  <- ps_match(d, score_col = "logit_t", scale = "logit", caliper_sd = 0.2, seed = 101)
  obj_abs <- ps_match(d, score_col = "logit_t", scale = "logit", caliper = 0.2 * pooled_sd, seed = 101)
  expect_identical(obj_sd$data$pair_id, obj_abs$data$pair_id)
  expect_equal(obj_sd$meta$caliper, 0.2 * pooled_sd)
  expect_identical(obj_sd$meta$caliper_sd, 0.2)
  expect_identical(obj_sd$meta$scale, "logit")
})

test_that("matching on logit_t leaves prob_t and identifiers out of the default balance table", {
  d <- sample_ps_data(n = 300, seed = 42)
  d$logit_t <- stats::qlogis(d$prob_t)
  obj <- ps_match(d, score_col = "logit_t", scale = "logit", caliper_sd = 0.2, seed = 101)
  for (tbl in obj$tables[c("smd_before", "smd_after")]) {
    expect_false(any(c("prob_t", "logit_t", "id", "match", "pair_id") %in% tbl$variable))
  }
  explicit <- ps_match(d, score_col = "logit_t", scale = "logit", caliper_sd = 0.2, seed = 101,
                       covariates = c("prob_t", "age"))
  expect_identical(explicit$tables$smd_before$variable, c("prob_t", "age"))
})

test_that("scale defaults to probability and still rejects values outside [0, 1]", {
  d <- dta
  d$logit_t <- stats::qlogis(d$prob_t)
  expect_identical(ps_match(d, seed = 1)$meta$scale, "probability")
  expect_error(ps_match(d, score_col = "logit_t"), "\\[0, 1\\]")
  expect_no_error(ps_match(d, score_col = "logit_t", scale = "logit"))
})

test_that("logit scale and caliper_sd reject bad input", {
  d <- dta
  d$logit_t <- stats::qlogis(d$prob_t)
  d$logit_t[1] <- Inf
  expect_error(ps_match(d, score_col = "logit_t", scale = "logit"), "finite")
  expect_error(ps_match(dta, scale = "odds"))
  expect_error(ps_match(dta, caliper_sd = -1), "caliper_sd")
  expect_error(ps_match(dta, caliper_sd = 0.2, caliper = 0.01), "at most one")
})

test_that("caliper_sd ignores patients with missing treatment", {
  d <- sample_ps_data(n = 300, seed = 42)
  d$logit_t <- stats::qlogis(d$prob_t)
  extra <- d[1L, ]
  extra$id <- max(d$id) + 1L
  extra$tavr <- NA
  extra$logit_t <- 50
  base <- ps_match(d, score_col = "logit_t", scale = "logit", caliper_sd = 0.2, seed = 101)
  with_na <- suppressWarnings(
    ps_match(rbind(d, extra), score_col = "logit_t", scale = "logit", caliper_sd = 0.2, seed = 101)
  )
  expect_equal(with_na$meta$caliper, base$meta$caliper)
})
