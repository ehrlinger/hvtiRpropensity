###############################################################################
## test_smd_tables.R
##
## The SMD tables ps_match(), ps_logistic() and ps_weight() report come from
## ps_stddiff(), so the package has one standardized-difference formula.
## Until 2026-09-16 the unweighted tables pooled the SD by sample size and the
## weighted table divided its variance by sum(w); neither reproduced the SAS
## %stddiff macro.
###############################################################################

test_that("ps_match() smd_before uses the SAS denominator when groups are unequal", {
  # tavr 1: age 1 2 3, mean 2, var 1.   tavr 0: age 2 4 6 8, mean 5, var 20/3.
  # SAS: -3 / sqrt((1 + 20/3) / 2) = -1.5323.
  # The sample-size-weighted pool this replaced gave -3 / sqrt(22/5) = -1.4302.
  d <- data.frame(id = 1:7, tavr = c(1, 1, 1, 0, 0, 0, 0), age = c(1, 2, 3, 2, 4, 6, 8),
                  prob_t = c(.6, .5, .4, .55, .45, .35, .3))
  obj <- ps_match(d, covariate_cols = "age", seed = 1)
  expect_equal(obj$tables$smd_before$smd, round(-3 / sqrt(23 / 6), 4L))
})

test_that("ps_logistic() smd matches ps_stddiff(), binary columns taking the binary formula", {
  dta <- sample_ps_data(n = 150, seed = 3)
  dta$prob_t <- NULL
  # sample_ps_data() splits the arms evenly, where the two denominators agree;
  # dropping treated patients makes them differ.
  dta <- dta[!(dta$tavr == 1 & seq_len(nrow(dta)) %% 3 != 0), ]
  obj <- ps_logistic(tavr ~ age + female + ef, data = dta, covariate_cols = c("age", "female", "ef"))
  ref <- ps_stddiff(dta, treatment_col = "tavr", gaussian_cols = c("age", "ef"), binary_cols = "female")$tables$stddiff
  expect_equal(obj$tables$smd$smd, round(ref$stddiff[match(c("age", "female", "ef"), ref$variable)], 4L))
})

test_that("ps_weight() weighted smd uses the PROC MEANS n - 1 weighted variance", {
  dta <- sample_ps_data(n = 150, seed = 4)
  obj <- ps_weight(dta, covariate_cols = c("age", "ef"))
  ref <- ps_stddiff(obj$data, treatment_col = "tavr", gaussian_cols = c("age", "ef"),
                    weight_col = "iptw")$tables$stddiff
  expect_equal(obj$tables$smd_weighted$smd, round(ref$stddiff, 4L))
})

test_that("a match with no pairs still reports an all-NA smd_after", {
  dta <- sample_ps_data(n = 100, seed = 42)
  obj <- ps_match(dta, caliper = 1e-12)
  expect_identical(obj$meta$n_matched, 0L)
  expect_true(all(is.na(obj$tables$smd_after$smd)))
  expect_identical(obj$tables$smd_after$variable, obj$tables$smd_before$variable)
})

test_that("an unsupported covariate type stops even when a match finds no pairs", {
  dta <- sample_ps_data(n = 100, seed = 42)
  dta$op_date <- as.Date("2020-01-01") + seq_len(nrow(dta))
  expect_error(ps_match(dta, covariate_cols = c("age", "op_date"), caliper = 1e-12), "op_date")
})

test_that("a character covariate gets a balance row even when a match finds no pairs", {
  dta <- sample_ps_data(n = 100, seed = 42)
  dta$site <- rep(c("a", "b"), length.out = nrow(dta))
  obj <- ps_match(dta, covariate_cols = c("age", "site"), caliper = 1e-12)
  expect_identical(obj$tables$smd_after$variable, c("age", "site"))
})

test_that("ps_weight() stops clearly when a score of 0 or 1 gives an infinite weight", {
  dta <- sample_ps_data(n = 60, seed = 2)
  dta$prob_t[which(dta$tavr == 1)[1]] <- 0
  expect_error(ps_weight(dta), "trim")
  # Winsorising clips the infinite weight, so trim is the documented way through.
  expect_s3_class(ps_weight(dta, trim = 0.05), "ps_weight")
})

test_that("a covariate takes the same formula before and after a match", {
  # 0/1/2 in the full data, only 0/1 in the subset: the subset table keeps the
  # Gaussian formula the full data chose.
  d <- data.frame(tavr = rep(0:1, each = 6), x = c(0, 1, 2, 0, 1, 1, 1, 0, 1, 1, 0, 1))
  types <- .covariate_types(d, "x")
  expect_identical(unname(types), "gaussian")
  sub <- d[d$x != 2, ]
  ref <- ps_stddiff(sub, treatment_col = "tavr", gaussian_cols = "x")$tables$stddiff$stddiff
  expect_equal(.smd_table(sub, "tavr", "x", types = types)$smd, round(ref, 4))
  expect_false(isTRUE(all.equal(.smd_table(sub, "tavr", "x")$smd, round(ref, 4))))
})

test_that("a three-level factor keeps the unsigned multilevel distance when a level drops out", {
  # PR #122 review: after matching, only two levels remain.
  d <- data.frame(tavr = rep(0:1, each = 9),
                  site = factor(c("a", "a", "b", "b", "b", "c", "c", "c", "c",
                                  "a", "a", "a", "a", "b", "b", "c", "c", "c")))
  types <- .covariate_types(d, "site")
  sub <- d[d$site != "c", ]
  sub_smd <- .smd_table(sub, "tavr", "site", types = types)$smd
  ref <- ps_stddiff(sub, treatment_col = "tavr", categorical_cols = "site")$tables$stddiff$stddiff
  expect_true(ref < 0)
  expect_equal(sub_smd, round(abs(ref), 4))
  expect_gte(.smd_table(d, "tavr", "site", types = types)$smd, 0)
})
