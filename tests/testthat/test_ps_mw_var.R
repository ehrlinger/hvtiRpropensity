###############################################################################
## test_ps_mw_var.R
##
## Tests for ps_mw_var(), a port of %mw_var (wt_mtch_boot_sd.sas, Rajeswaran
## 2014): a matching-weighted difference in group means with a bootstrap SD,
## percentiles, z and p.
###############################################################################

dta <- sample_ps_data(n = 80, seed = 12)
dta$mt_wt <- ifelse(dta$tavr == 1, 1 - dta$prob_t, dta$prob_t) / pmax(dta$prob_t, 1 - dta$prob_t)

# ---- Point estimate ---------------------------------------------------------

test_that("group summaries are PROC MEANS weighted mean, SD and sum of weights", {
  # g0: y 1 2 3, w 1 1 2 -> mean 9/4, var (1.5625 + .0625 + 1.125) / 2 = 1.375, sumwgt 4
  # g1: y 2 4 6, w 2 1 1 -> mean 14/4, var (4.5 + .25 + 6.25) / 2 = 5.5, sumwgt 4
  d <- data.frame(g = c(0, 0, 0, 1, 1, 1), y = c(1, 2, 3, 2, 4, 6), w = c(1, 1, 2, 2, 1, 1))
  tbl <- ps_mw_var(d, treatment_col = "g", outcomes = "y", weight_col = "w", n_rep = 20, seed = 1)$tables$mw_var
  expect_equal(tbl$mean_0, 2.25)
  expect_equal(tbl$sd_0, sqrt(1.375))
  expect_equal(tbl$sumwgt_0, 4)
  expect_equal(tbl$mean_1, 3.5)
  expect_equal(tbl$sd_1, sqrt(5.5))
  expect_equal(tbl$sumwgt_1, 4)
  expect_equal(tbl$estimate, 1.25)
})

# ---- Bootstrap --------------------------------------------------------------

test_that("the bootstrap resamples within each group, with the weights held fixed", {
  # Rebuild by hand: per replicate, draw group 0's rows then group 1's, each
  # with replacement to its own size, and keep each row's original weight.
  obj <- ps_mw_var(dta, outcomes = c("age", "ef"), weight_col = "mt_wt", n_rep = 40, seed = 21)
  i0 <- which(dta$tavr == 0)
  i1 <- which(dta$tavr == 1)
  wmean <- function(r, v) sum(dta$mt_wt[r] * dta[[v]][r]) / sum(dta$mt_wt[r])
  set.seed(21)
  reps <- t(vapply(seq_len(40), function(i) {
    r0 <- i0[sample.int(length(i0), replace = TRUE)]
    r1 <- i1[sample.int(length(i1), replace = TRUE)]
    c(wmean(r1, "age") - wmean(r0, "age"), wmean(r1, "ef") - wmean(r0, "ef"))
  }, numeric(2)))
  tbl <- obj$tables$mw_var
  expect_equal(tbl$sd, c(sd(reps[, 1]), sd(reps[, 2])))
  expect_equal(unlist(tbl[1, c("p2_5", "p16", "p50", "p84", "p97_5")], use.names = FALSE),
               quantile(reps[, 1], c(.025, .16, .5, .84, .975), type = 4, names = FALSE))
  expect_identical(tbl$n_rep, c(40L, 40L))
  expect_equal(tbl$z, tbl$estimate / tbl$sd)
  expect_equal(tbl$p, 2 * (1 - pnorm(abs(tbl$z))))
})

test_that("bootstrap percentiles use PROC STDIZE PCTLDEF=1 (R type 4)", {
  # 1..10 at type 4: .025 -> 1, .16 -> 1.6, .5 -> 5, .84 -> 8.4, .975 -> 9.75.
  # PCTLDEF=5 (type 2), which stddiffci uses, would give 1, 2, 5.5, 9, 10.
  expect_equal(hvtiRpropensity:::.perm_percentiles(c(4, 1, 10, 7, 2, 9, 3, 8, 5, 6), type = 4),
               c(1, 1.6, 5, 8.4, 9.75))
})

test_that("a missing outcome is dropped for that outcome only", {
  d <- dta
  d$age[1:4] <- NA
  obj <- ps_mw_var(d, outcomes = c("age", "ef"), weight_col = "mt_wt", n_rep = 10, seed = 2)
  ok <- !is.na(d$age)
  m1 <- with(d[ok & d$tavr == 1, ], sum(mt_wt * age) / sum(mt_wt))
  m0 <- with(d[ok & d$tavr == 0, ], sum(mt_wt * age) / sum(mt_wt))
  expect_equal(obj$tables$mw_var$estimate[1], m1 - m0)
  expect_false(anyNA(obj$tables$mw_var$sd))
})

test_that("a seed is reproducible and leaves the caller's random number stream as it was", {
  a <- ps_mw_var(dta, outcomes = "age", weight_col = "mt_wt", n_rep = 15, seed = 5)
  b <- ps_mw_var(dta, outcomes = "age", weight_col = "mt_wt", n_rep = 15, seed = 5)
  expect_identical(a$tables, b$tables)
  set.seed(77)
  before <- runif(3)
  set.seed(77)
  ps_mw_var(dta, outcomes = "age", weight_col = "mt_wt", n_rep = 5, seed = 9)
  expect_identical(runif(3), before)
})

# ---- Structure and validation -----------------------------------------------

test_that("ps_mw_var() returns a ps_mw_var / ps_data object", {
  obj <- ps_mw_var(dta, outcomes = "age", weight_col = "mt_wt", n_rep = 5, seed = 1)
  expect_s3_class(obj, "ps_mw_var")
  expect_s3_class(obj, "ps_data")
  expect_identical(obj$data, dta)
  expect_named(obj$tables$mw_var,
               c("outcome", "label", "sumwgt_0", "mean_0", "sd_0", "sumwgt_1", "mean_1", "sd_1",
                 "estimate", "n_rep", "sd", "p2_5", "p16", "p50", "p84", "p97_5", "z", "p"))
})

test_that("outcomes must be numeric and the weight positive and finite", {
  d <- dta
  d$site <- "a"
  expect_error(ps_mw_var(d, outcomes = "site", weight_col = "mt_wt", n_rep = 5), "numeric")
  d$mt_wt[1] <- 0
  expect_error(ps_mw_var(d, outcomes = "age", weight_col = "mt_wt", n_rep = 5), "positive")
})

test_that("n_rep must be a positive whole number", {
  expect_error(ps_mw_var(dta, outcomes = "age", weight_col = "mt_wt", n_rep = 0), "n_rep")
  expect_error(ps_mw_var(dta, outcomes = "age", weight_col = "mt_wt", n_rep = Inf), "n_rep")
})
