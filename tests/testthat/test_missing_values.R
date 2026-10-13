###############################################################################
## test_missing_values.R
##
## Rows with a missing treatment or a missing score are excluded, counted and
## stated in one warning (issue #84).
###############################################################################

missing_rows <- function(seed = 84) {
  dta <- sample_ps_data(n = 120, seed = seed)
  dta$tavr[c(3, 40)] <- NA
  dta$prob_t[c(7, 41, 90)] <- NA
  dta
}

# ---------------------------------------------------------------------------
# ps_match
# ---------------------------------------------------------------------------

test_that("ps_match() excludes rows with a missing score or treatment, with one warning", {
  dta <- missing_rows()
  expect_warning(m <- ps_match(dta, seed = 1), "5 row\\(s\\)")
  bad <- is.na(dta$tavr) | is.na(dta$prob_t)
  expect_true(all(m$data$match[bad] == 0L))
  expect_true(all(is.na(m$data$pair_id[bad])))
  expect_identical(m$meta$n_excluded, 5L)
  expect_equal(sum(m$tables$group_counts_before$n), nrow(dta) - 5)
  expect_false(anyNA(m$tables$smd_before$smd))
})

test_that("ps_match() with a caliper does not crash on a missing score", {
  dta <- missing_rows()
  expect_warning(m <- ps_match(dta, caliper = 0.05, seed = 1), "5 row\\(s\\)")
  expect_true(all(m$data$match[is.na(dta$prob_t)] == 0L))
})

test_that("ps_match() warns exactly once when both treatment and score are missing", {
  dta <- missing_rows()
  warnings <- testthat::capture_warnings(ps_match(dta, seed = 1))
  expect_length(warnings, 1L)
})

# ---------------------------------------------------------------------------
# ps_weight
# ---------------------------------------------------------------------------

test_that("ps_weight() gives excluded rows a missing weight and keeps the rest finite", {
  dta <- missing_rows()
  warnings <- testthat::capture_warnings(w <- ps_weight(dta))
  expect_length(warnings, 1L)
  expect_match(warnings, "5 row\\(s\\)")
  bad <- is.na(dta$tavr) | is.na(dta$prob_t)
  expect_true(all(is.na(w$data$iptw[bad])))
  expect_true(all(is.finite(w$data$iptw[!bad])))
  expect_identical(w$meta$n_excluded, 5L)
  expect_equal(sum(w$tables$group_counts$n), nrow(dta) - 5)
  expect_false(anyNA(w$tables$effective_n$n_effective))
  expect_false(anyNA(w$tables$smd_weighted$smd))
})

test_that("ps_weight() stabilises with the treated share of the analysable rows", {
  dta <- missing_rows()
  keep <- !(is.na(dta$tavr) | is.na(dta$prob_t))
  full <- ps_weight(dta[keep, ])
  part <- suppressWarnings(ps_weight(dta))
  expect_equal(part$data$iptw[keep], full$data$iptw)
})

# ---------------------------------------------------------------------------
# ps_support
# ---------------------------------------------------------------------------

test_that("ps_support() computes bounds on analysable rows and flags excluded rows NA", {
  dta <- missing_rows()
  m <- suppressWarnings(ps_match(dta, seed = 1))
  expect_warning(s <- ps_support(m), "5 row\\(s\\)")
  bad <- is.na(dta$tavr) | is.na(dta$prob_t)
  expect_false(anyNA(s$meta$common_support))
  expect_true(all(is.na(s$data$weak_common[bad])))
  expect_true(all(is.na(s$data$weak_trim[bad])))
  expect_false(anyNA(s$data$weak_common[!bad]))
  expect_equal(s$tables$counts$n, rep(nrow(dta) - 5L, 2L))
  expect_equal(s$tables$counts$n_excluded, c(5L, 5L))
  expect_gt(s$tables$counts$n_weak[s$tables$counts$rule == "common"], 0)
})

test_that("ps_support() works with a single rule", {
  dta <- sample_ps_data(n = 120, seed = 5)
  m <- ps_match(dta, seed = 1)
  s1 <- ps_support(m, trim = NULL)
  expect_identical(s1$meta$rules, "common")
  expect_equal(nrow(s1$tables$agreement), 0L)
  s2 <- ps_support(m, common = FALSE)
  expect_identical(s2$meta$rules, "trim")
})

# ---------------------------------------------------------------------------
# ps_logistic
# ---------------------------------------------------------------------------

test_that("ps_logistic() counts groups and strata on patients with a treatment", {
  dta <- sample_ps_data(n = 120, seed = 6)
  dta$prob_t <- NULL
  dta$tavr[c(2, 30)] <- NA
  warnings <- testthat::capture_warnings(
    obj <- ps_logistic(tavr ~ age + female + ef, data = dta)
  )
  expect_length(warnings, 1L)
  expect_match(warnings, "missing for 2 patient")
  expect_false(anyNA(obj$tables$group_counts$n))
  expect_equal(sum(obj$tables$group_counts$n), nrow(dta) - 2)
  expect_identical(obj$meta$n_missing_treatment, 2L)
  expect_false(anyNA(obj$data$prob_t[c(2, 30)]))
  expect_true(all(is.na(obj$data$quintile[c(2, 30)])))
  expect_true(all(is.na(obj$data$decile[c(2, 30)])))
  expect_true(all(is.na(obj$data$mt_wt[c(2, 30)])))
  expect_equal(as.vector(table(obj$data$quintile)), rep((nrow(dta) - 2) / 5, 5), tolerance = 0.05)
  expect_identical(obj$tables$fit_status$n_excluded, 2L)
})

# ---------------------------------------------------------------------------
# fit_status for every family
# ---------------------------------------------------------------------------

test_that("fit_status counts rows dropped by a nominal model", {
  dta <- sample_ps_data_nominal(n = 100, seed = 7)
  dta$age[c(1, 5, 9)] <- NA
  obj <- ps_nominal(rtyp ~ age + female + ef, data = dta)
  expect_identical(obj$tables$fit_status$n_analyzed, nrow(dta) - 3L)
  expect_identical(obj$tables$fit_status$n_excluded, 3L)
  fit <- fit_logistic(rtyp ~ age + female + ef, data = dta, family = "nominal",
                      outcome_levels = levels(dta$rtyp), reference_level = levels(dta$rtyp)[1L])
  expect_identical(fit$tables$fit_status$n_excluded, 3L)
})

test_that("fit_status counts rows dropped by binary, ordinal and linear models", {
  dta <- sample_ps_data(n = 100, seed = 8)
  dta$age[c(2, 4)] <- NA
  bin <- fit_logistic(tavr ~ age + ef, data = dta, family = "binary",
                      outcome_levels = c(0, 1), event_level = 1)
  expect_identical(bin$tables$fit_status$n_excluded, 2L)
  lin <- fit_linear(ef ~ age + female, data = dta)
  expect_identical(lin$tables$fit_status$n_excluded, 2L)
  ord <- sample_ps_data_ordinal(n = 100, seed = 8)
  ord$age[c(2, 4)] <- NA
  ofit <- ps_ordinal(stats::reformulate(c("age", "female"), names(ord)[2L]), data = ord)
  expect_identical(ofit$tables$fit_status$n_excluded, 2L)
})

# ---------------------------------------------------------------------------
# bs_continuous / bs_count
# ---------------------------------------------------------------------------

test_that("bs_continuous() gives a missing balancing score missing strata", {
  dta <- sample_ps_data(n = 100, seed = 9)
  dta$age[c(1, 2)] <- NA
  obj <- bs_continuous(ef ~ age + female + diabetes, data = dta)
  expect_true(all(is.na(obj$data$bs[1:2])))
  expect_true(all(is.na(obj$data$quintile[1:2])))
  expect_true(all(is.na(obj$data$cluster[1:2])))
  expect_equal(sum(obj$tables$strata_counts$n), nrow(dta) - 2)
  expect_equal(range(obj$data$decile, na.rm = TRUE), c(1L, 10L))
})

test_that("bs_count() gives a missing balancing score missing strata", {
  dta <- sample_ps_data_count(n = 100, seed = 9)
  dta$age[c(1, 2)] <- NA
  obj <- bs_count(rbc_tot ~ age + female + hct_pr, data = dta, dist = "poisson")
  expect_true(all(is.na(obj$data$cluster[1:2])))
  expect_equal(sum(obj$tables$strata_counts$n), nrow(dta) - 2)
})

test_that("bs_*() need id_col only with imputations", {
  dta <- sample_ps_data(n = 100, seed = 10)
  dta$id <- NULL
  expect_s3_class(bs_continuous(ef ~ age + female, data = dta), "bs_continuous")
  cnt <- sample_ps_data_count(n = 100, seed = 10)
  cnt$id <- NULL
  expect_s3_class(bs_count(rbc_tot ~ age + female, data = cnt, dist = "poisson"), "bs_count")
})
