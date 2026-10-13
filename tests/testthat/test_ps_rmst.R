###############################################################################
## test_ps_rmst.R -- ps_rmst() on ps_match and ps_weight inputs, and tau
###############################################################################

skip_if_not_installed("survival")

rmst_data <- local({
  d <- sample_ps_data(n = 200, seed = 9)
  u <- (seq_len(nrow(d)) * 0.618034) %% 1
  d$t <- pmin(-log(u + 1e-3) / (0.4 * exp(-1.5 * d$tavr)), 5)
  d$e <- as.integer(d$t < 5)
  d
})

# Area under the weighted Kaplan-Meier curve on [0, tau], computed directly.
manual_rmst <- function(time, event, w, tau) {
  km <- survival::survfit(survival::Surv(time, event) ~ 1, weights = w)
  tt <- c(0, km$time)
  ss <- c(1, km$surv)
  keep <- tt <= tau
  sum(diff(c(tt[keep], tau)) * ss[keep])
}

test_that("a ps_match input is analysed on its matched rows and resampled by pair", {
  m <- ps_match(rmst_data, seed = 1, caliper = 0.05)
  r <- ps_rmst(m, "t", "e", tau = 4, n_boot = 20, seed = 3)
  est <- r$tables$estimates
  expect_identical(est$weighting, "unweighted")
  expect_equal(est$n, 2L * m$meta$n_matched)
  expect_lt(est$n, nrow(rmst_data))
  expect_identical(r$meta$resampling, "pair")

  md <- m$data[m$data$match == 1L, ]
  a <- md$tavr == 1
  expect_equal(est$rmst_treated, manual_rmst(md$t[a], md$e[a], rep(1, sum(a)), 4))
  expect_equal(est$rmst_control, manual_rmst(md$t[!a], md$e[!a], rep(1, sum(!a)), 4))

  expect_identical(ps_rmst(m, "t", "e", tau = 4, n_boot = 20, seed = 3)$tables, r$tables)
  set.seed(21)
  before <- .Random.seed
  ps_rmst(m, "t", "e", tau = 4, n_boot = 5, seed = 3)
  expect_identical(.Random.seed, before)
})

test_that("the pair bootstrap keeps both members of every drawn pair", {
  # Each pair's two members share a time and event, so any draw of whole pairs gives identical arms
  # and a contrast of exactly zero. A row bootstrap would not.
  k <- 12
  d <- data.frame(id = seq_len(2 * k), tavr = rep(c(1L, 0L), each = k),
                  s = rep(seq(0.2, 0.75, length.out = k), 2),
                  t = rep(seq(1, 5, length.out = k), 2), e = rep(rep(c(1L, 0L), length.out = k), 2))
  m <- ps_match(d, score_col = "s", seed = 1)
  expect_equal(m$meta$n_matched, k)
  expect_true(all(tapply(m$data$t, m$data$pair_id, function(v) length(unique(v)) == 1L)))
  r <- ps_rmst(m, "t", "e", tau = 4, n_boot = 50, seed = 2)
  est <- r$tables$estimates
  expect_equal(est$n, 2L * k)
  expect_equal(est$diff, 0)
  expect_equal(c(est$lo_days, est$hi_days), c(0, 0))
  expect_equal(est$n_failed, 0L)

  row_boot <- ps_rmst(ps_logistic(tavr ~ s, d), "t", "e", tau = 4, weights = "unweighted", n_boot = 50,
                      seed = 2)
  expect_identical(row_boot$meta$resampling, "row")
  expect_gt(row_boot$tables$estimates$hi_days, row_boot$tables$estimates$lo_days)
})

test_that("a ps_match input refuses score weighting", {
  m <- ps_match(rmst_data, seed = 1, caliper = 0.05)
  expect_error(ps_rmst(m, "t", "e", tau = 4, weights = "ato", n_boot = 2), "unweighted")
  expect_error(ps_rmst(m, "t", "e", tau = 4, weights = c("unweighted", "att"), n_boot = 2), "unweighted")
})

test_that("a ps_weight input uses its stored weight column", {
  w <- ps_weight(rmst_data, estimand = "ATT", trim = 0.05)
  r <- ps_rmst(w, "t", "e", tau = 4, weights = c("unweighted", "att"), n_boot = 5, seed = 1)
  expect_identical(r$meta$stored_weighting, "att")
  expect_equal(r$data$w_att, w$data$iptw)
  est <- r$tables$estimates[r$tables$estimates$weighting == "att", ]
  a <- rmst_data$tavr == 1
  expect_equal(est$rmst_treated, manual_rmst(rmst_data$t[a], rmst_data$e[a], w$data$iptw[a], 4))
  expect_equal(est$rmst_control, manual_rmst(rmst_data$t[!a], rmst_data$e[!a], w$data$iptw[!a], 4))
  wc <- w$data$iptw[!a]
  expect_equal(est$ess_control, sum(wc)^2 / sum(wc^2))

  # A scheme that does not match the estimand is still computed from the score.
  r_ato <- ps_rmst(w, "t", "e", tau = 4, weights = "ato", n_boot = 2, seed = 1)
  expect_true(is.na(r_ato$meta$stored_weighting))
  expect_equal(r_ato$data$w_ato, ifelse(a, 1 - rmst_data$prob_t, rmst_data$prob_t))
})

test_that("tau past the follow-up of either arm stops", {
  fit <- ps_logistic(tavr ~ age + ef, rmst_data[setdiff(names(rmst_data), "prob_t")])
  expect_error(ps_rmst(fit, "t", "e", tau = 6, n_boot = 2), "past the last observed time")
  short <- rmst_data$tavr == 0 | rmst_data$t < 3
  expect_error(ps_rmst(fit, "t", "e", tau = 4, subsets = list(short = short), n_boot = 2),
               "treated arm of estimator \"short/unweighted\"")
  expect_no_error(ps_rmst(fit, "t", "e", tau = 5, n_boot = 2))
})

test_that("stored ps_weight weights survive a weight_col named like an output column", {
  base <- ps_weight(rmst_data, estimand = "ATE")
  clash <- ps_weight(rmst_data, estimand = "ATE", weight_col = "w_unweighted")
  r_base <- ps_rmst(base, "t", "e", tau = 4, weights = c("unweighted", "ate"), n_boot = 2, seed = 1)
  r_clash <- ps_rmst(clash, "t", "e", tau = 4, weights = c("unweighted", "ate"), n_boot = 2, seed = 1)
  ate <- function(r) r$tables$estimates[r$tables$estimates$weighting == "ate", ]
  expect_equal(ate(r_clash)$diff, ate(r_base)$diff)
  expect_false(isTRUE(all.equal(ate(r_clash)$diff, 0)))
  expect_equal(r_clash$data$w_ate, base$data$iptw)
  expect_true(all(r_clash$data$w_unweighted == 1))
})

short_control <- function() {
  # Controls are followed to 2 at most and treated to 7, plus one control row at 8.
  n <- 40
  d <- data.frame(id = seq_len(n + 1), tavr = c(rep(c(1L, 0L), n / 2), 0L),
                  x = c(seq(-1, 1, length.out = n), 0))
  d$t <- ifelse(d$tavr == 1, seq(0.5, 7, length.out = n + 1), seq(0.2, 2, length.out = n + 1))
  d$t[n + 1] <- 8
  d$e <- rep_len(c(1L, 0L, 1L), n + 1)
  d
}

test_that("the horizon guard ignores a zero-weight row past follow-up", {
  d <- short_control()
  d$prob_t <- 0.5
  w <- ps_weight(d, estimand = "ATT", stabilise = FALSE)
  w$data$iptw[nrow(d)] <- 0
  expect_error(ps_rmst(w, "t", "e", tau = 6, weights = "att", n_boot = 2), "control arm")
  expect_no_error(ps_rmst(w, "t", "e", tau = 6, weights = "unweighted", n_boot = 2))
})

test_that("the horizon guard ignores a row with a missing event past follow-up", {
  d <- short_control()
  d$e[nrow(d)] <- NA
  fit <- ps_logistic(tavr ~ x, d)
  expect_error(ps_rmst(fit, "t", "e", tau = 6, weights = "unweighted", n_boot = 2), "control arm")
})
