###############################################################################
## test_input_validation.R
##
## Argument validation in ps_weight(), ps_match() and ps_support() (#88).
###############################################################################

dta <- sample_ps_data(n = 100, seed = 42)

# ---- ps_weight(): stabilisation applies to the ATE only --------------------

test_that("ATE weights are unchanged: stabilised = raw weight times own-arm probability", {
  obj <- ps_weight(dta, estimand = "ATE", stabilise = TRUE)
  ps  <- dta$prob_t
  trt <- dta$tavr
  p1  <- mean(trt == 1L)
  expected <- ifelse(trt == 1L, p1 / ps, (1 - p1) / (1 - ps))
  expect_equal(obj$data$iptw, expected)
  expect_true(obj$meta$stabilised)
  expect_false(ps_weight(dta, estimand = "ATE", stabilise = FALSE)$meta$stabilised)
})

test_that("ATT treated weights are exactly 1 with stabilise = TRUE", {
  obj <- ps_weight(dta, estimand = "ATT", stabilise = TRUE)
  expect_identical(unique(obj$data$iptw[dta$tavr == 1L]), 1)
  expect_equal(obj$data$iptw, ps_weight(dta, estimand = "ATT", stabilise = FALSE)$data$iptw)
  expect_false(obj$meta$stabilised)
})

test_that("ATC control weights are exactly 1 with stabilise = TRUE", {
  obj <- ps_weight(dta, estimand = "ATC", stabilise = TRUE)
  expect_identical(unique(obj$data$iptw[dta$tavr == 0L]), 1)
  expect_equal(obj$data$iptw, ps_weight(dta, estimand = "ATC", stabilise = FALSE)$data$iptw)
  expect_false(obj$meta$stabilised)
})

test_that("sa_trim_sweep() does not stabilise ATT weights either", {
  sw_t <- sa_trim_sweep(dta, trim_seq = 0, estimand = "ATT", stabilise = TRUE)
  sw_f <- sa_trim_sweep(dta, trim_seq = 0, estimand = "ATT", stabilise = FALSE)
  expect_equal(sw_t, sw_f)
})

test_that("ps_weight() names an NA trim", {
  expect_error(ps_weight(dta, trim = NA_real_), "`trim`")
  expect_error(ps_weight(dta, trim = NA), "`trim`")
})

# ---- ps_match(): ratio must be the number 1 ---------------------------------

test_that("ps_match() refuses a non-integer, character, NA or vector ratio", {
  expect_error(ps_match(dta, ratio = 1.5), "ratio = 1.5")
  expect_error(ps_match(dta, ratio = "1"), "ratio = \"1\"", fixed = TRUE)
  expect_error(ps_match(dta, ratio = NA_real_), "ratio = NA")
  expect_error(ps_match(dta, ratio = c(1, 1)), "ratio = c(1, 1)", fixed = TRUE)
  expect_no_error(ps_match(dta, ratio = 1))
  expect_no_error(ps_match(dta, ratio = 1L))
})

test_that("ps_match() names an NA caliper or caliper_sd", {
  expect_error(ps_match(dta, caliper = NA_real_), "`caliper`")
  expect_error(ps_match(dta, caliper_sd = NA_real_), "`caliper_sd`")
})

# ---- ps_support(): flag names and trim --------------------------------------

sup_in <- ps_logistic(tavr ~ age + ef, dta[, c("id", "tavr", "age", "ef")])

test_that("ps_support() refuses a flag named after a built-in rule", {
  flag <- rep(FALSE, nrow(sup_in$data))
  expect_error(ps_support(sup_in, flags = list(trim = flag)), "\"trim\"")
  expect_error(ps_support(sup_in, flags = list(common = flag)), "\"common\"")
  expect_no_error(ps_support(sup_in, flags = list(other = flag)))
})

test_that("ps_support() refuses an out-of-range or NA trim", {
  expect_error(ps_support(sup_in, trim = c(-0.1, 0.9)), "`trim`")
  expect_error(ps_support(sup_in, trim = c(0.1, 1.5)), "`trim`")
  expect_error(ps_support(sup_in, trim = c(NA, 0.9)), "`trim`")
  expect_error(ps_support(sup_in, trim = c(0.9, 0.1)), "`trim`")
  expect_no_error(ps_support(sup_in, trim = c(0.05, 0.95)))
})
