###############################################################################
## test_treated_level.R
##
## A treatment column holding any two values, with `treated_level` naming the
## treated one, must give the same numbers as the same data coded 0/1.
###############################################################################

binary <- sample_ps_data(n = 150, seed = 42)
labelled <- binary
labelled$tavr <- ifelse(binary$tavr == 1L, "transcatheter", "surgical")
lvl <- "transcatheter"

test_that("ps_match() gives the 0/1 result and keeps the caller's column", {
  ref <- ps_match(binary, seed = 3)
  obj <- ps_match(labelled, seed = 3, treated_level = lvl)
  expect_identical(obj$data$tavr, labelled$tavr)
  expect_identical(obj$data$pair_id, ref$data$pair_id)
  expect_equal(obj$tables$smd_after, ref$tables$smd_after)
  expect_identical(obj$tables$group_counts_before$group, c("surgical", "transcatheter"))
  expect_identical(obj$tables$group_counts_before$treated, c(FALSE, TRUE))
  expect_equal(obj$tables$group_counts_before$n, ref$tables$group_counts_before$n)
  expect_identical(obj$meta$treated_level, lvl)
})

test_that("ps_weight() gives the 0/1 result and keeps the caller's column", {
  ref <- ps_weight(binary)
  obj <- ps_weight(labelled, treated_level = lvl)
  expect_identical(obj$data$tavr, labelled$tavr)
  expect_equal(obj$data$iptw, ref$data$iptw)
  expect_equal(obj$tables$smd_weighted, ref$tables$smd_weighted)
  expect_identical(obj$tables$effective_n$group, c("surgical", "transcatheter"))
  expect_equal(obj$tables$effective_n$n_effective, ref$tables$effective_n$n_effective)
})

test_that("treated_level can name 0 as the treated value", {
  flipped <- binary
  flipped$prob_t <- 1 - binary$prob_t
  obj <- ps_weight(flipped, treated_level = 0)
  ref <- ps_weight(transform(flipped, tavr = 1L - tavr))
  expect_equal(obj$data$iptw, ref$data$iptw)
  expect_identical(obj$tables$group_counts$group, c("1", "0"))
  expect_identical(obj$data$tavr, flipped$tavr)
})

test_that("ps_stddiff(), ps_stddiff_perm() and ps_mw_var() give the 0/1 result", {
  expect_equal(
    ps_stddiff(labelled, gaussian = c("age", "ef"), binary = "female", treated_level = lvl)$tables$stddiff,
    ps_stddiff(binary, gaussian = c("age", "ef"), binary = "female")$tables$stddiff
  )
  expect_equal(
    ps_stddiff_perm(labelled, gaussian = "age", n_perm = 20, seed = 5, treated_level = lvl)$tables,
    ps_stddiff_perm(binary, gaussian = "age", n_perm = 20, seed = 5)$tables
  )
  weighted <- ps_weight(binary)$data
  weighted_labelled <- weighted
  weighted_labelled$tavr <- labelled$tavr
  expect_equal(
    ps_mw_var(weighted_labelled, outcomes = "ef", weight_col = "iptw", n_rep = 20, seed = 5,
              treated_level = lvl)$tables,
    ps_mw_var(weighted, outcomes = "ef", weight_col = "iptw", n_rep = 20, seed = 5)$tables
  )
})

test_that("sa_overlap() and sa_trim_sweep() give the 0/1 result", {
  ref <- sa_overlap(binary)
  obj <- sa_overlap(labelled, treated_level = lvl)
  expect_identical(obj$summary$group, c("surgical", "transcatheter"))
  expect_equal(obj$summary[, -1], ref$summary[, -1])
  expect_equal(obj$positivity_flags[, -1], ref$positivity_flags[, -1])
  expect_equal(sa_trim_sweep(labelled, treated_level = lvl), sa_trim_sweep(binary))
})

test_that("a scored object carries its treated level downstream", {
  unscored <- labelled
  unscored$prob_t <- NULL
  fit <- ps_logistic(tavr ~ age + ef, unscored,
                     treatment_levels = c("surgical", "transcatheter"), treated_level = lvl)
  overlap <- sa_overlap(fit)
  expect_identical(overlap$summary$group, c("surgical", "transcatheter"))
  expect_equal(overlap$summary$n, c(sum(binary$tavr == 0), sum(binary$tavr == 1)))
  expect_no_warning(sa_trim_sweep(fit))

  matched <- ps_match(fit$data, seed = 3, treated_level = lvl)
  matched_ref <- ps_match(transform(fit$data, tavr = as.integer(tavr == lvl)), seed = 3)
  expect_equal(sa_rosenbaum(matched, outcome_col = "ef", gamma_max = 2),
               sa_rosenbaum(matched_ref, outcome_col = "ef", gamma_max = 2))
})

test_that("ps_forest() accepts a labelled treatment", {
  skip_if_not_installed("randomForestSRC")
  unscored <- labelled[, c("id", "tavr", "age", "ef")]
  obj <- ps_forest(tavr ~ age + ef, unscored, ntree = 60, seed = 7, treated_level = lvl)
  ref <- ps_forest(tavr ~ age + ef, transform(unscored, tavr = as.integer(tavr == lvl)), ntree = 60, seed = 7)
  expect_equal(obj$data$prob_t, ref$data$prob_t)
  expect_identical(obj$data$tavr, unscored$tavr)
  expect_identical(obj$tables$group_counts$group, c("surgical", "transcatheter"))
  expect_identical(obj$meta$treatment_levels, c("surgical", "transcatheter"))
})

test_that("a labelled column without treated_level is refused, not coerced", {
  expect_error(ps_match(labelled), "binary")
  expect_error(ps_weight(labelled), "binary")
  expect_error(sa_overlap(labelled), "binary")
  expect_error(sa_trim_sweep(labelled), "binary")
  expect_error(ps_stddiff(labelled, gaussian = "age"), "binary")
})

test_that("treated_level is validated", {
  expect_error(ps_weight(labelled, treated_level = "tavr"), "is not a value of column")
  three <- labelled
  three$tavr[1:5] <- "hybrid"
  expect_error(ps_weight(three, treated_level = lvl), "must hold two treatment values; it holds 3")
  expect_error(ps_weight(labelled, treated_level = c("a", "b")), "single non-missing value")
  expect_error(ps_weight(labelled, treated_level = NA), "single non-missing value")
})

test_that("treated_level needs both treatment values present", {
  one <- labelled[labelled$tavr == lvl, ]
  expect_error(ps_weight(one, treated_level = lvl), "must hold two treatment values; it holds 1")
  none <- labelled
  none$tavr <- NA_character_
  expect_error(ps_weight(none, treated_level = lvl), "must hold two treatment values; it holds 0")
})

test_that("sa_trim_sweep() excludes patients with a missing treatment", {
  missing <- binary
  missing$tavr[c(1, 200)] <- NA
  expect_warning(swept <- sa_trim_sweep(missing), "2 NA value")
  expect_equal(swept, sa_trim_sweep(binary[-c(1, 200), ]))
  expect_false(anyNA(swept))
})
