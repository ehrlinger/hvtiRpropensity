###############################################################################
## test_sa_sensitivity.R
##
## Tests for sa_rosenbaum(), sa_evalue(), sa_overlap(), and sa_trim_sweep().
###############################################################################


# ---------------------------------------------------------------------------
# Shared fixture
# ---------------------------------------------------------------------------

make_matched <- function(n = 200, seed = 1) {
  dta <- sample_ps_data(n = n, seed = seed)
  ps_match(dta)
}

make_weighted <- function(n = 200, seed = 1, estimand = "ATE") {
  dta <- sample_ps_data(n = n, seed = seed)
  ps_weight(dta, estimand = estimand)
}


# ---------------------------------------------------------------------------
# ps_match: pair_id column (prerequisite for sa_rosenbaum)
# ---------------------------------------------------------------------------

test_that("ps_match() adds pair_id column to $data", {
  m <- make_matched()
  expect_true("pair_id" %in% names(m$data))
})

test_that("ps_match() pair_id is NA for unmatched, integer for matched", {
  m <- make_matched()
  matched_rows   <- m$data[m$data$match == 1L, ]
  unmatched_rows <- m$data[m$data$match == 0L, ]
  expect_true(all(!is.na(matched_rows$pair_id)))
  expect_true(all(is.na(unmatched_rows$pair_id)))
})

test_that("ps_match() pair_id links exactly one treated to one control per pair", {
  m <- make_matched()
  matched <- m$data[m$data$match == 1L, ]
  pair_counts <- table(matched$pair_id)
  expect_true(all(pair_counts == 2L))
})

test_that("ps_match() meta contains pair_id_col", {
  m <- make_matched()
  expect_equal(m$meta$pair_id_col, "pair_id")
})


# ---------------------------------------------------------------------------
# sa_evalue
# ---------------------------------------------------------------------------

test_that("sa_evalue() returns an sa_evalue ps_data object", {
  res <- sa_evalue(2.5, ci_lo = 1.5, ci_hi = 3.5)
  expect_s3_class(res, c("sa_evalue", "ps_data"))
  expect_named(res, c("data", "meta", "tables", "models"))
  expect_identical(nrow(res$data), 0L)
  expect_true(all(c("evalue_estimate", "evalue_ci", "rr_estimate", "type", "rare") %in% names(res$meta)))
  expect_identical(res$tables$evalue$quantity, c("estimate", "ci"))
  expect_equal(res$tables$evalue$evalue, c(res$meta$evalue_estimate, res$meta$evalue_ci))
  expect_equal(res$tables$evalue$value[2], 1.5)
  expect_output(print(res), "E-value")
  expect_output(summary(res), "Evalue")
})

test_that("sa_evalue() E-value for RR = 1 is 1", {
  res <- sa_evalue(1.0)
  expect_equal(res$meta$evalue_estimate, 1)
})

test_that("sa_evalue() E-value for RR = 2 is 2 + sqrt(2)", {
  res <- sa_evalue(2.0)
  expect_equal(res$meta$evalue_estimate, 2 + sqrt(2 * (2 - 1)), tolerance = 1e-8)
})

test_that("sa_evalue() is symmetric: RR and 1/RR give the same E-value", {
  e1 <- sa_evalue(2.0)$meta$evalue_estimate
  e2 <- sa_evalue(0.5)$meta$evalue_estimate
  expect_equal(e1, e2, tolerance = 1e-8)
})

test_that("sa_evalue() evalue_ci is NA when no CI is supplied", {
  res <- sa_evalue(2.5)
  expect_true(is.na(res$meta$evalue_ci))
})

test_that("sa_evalue() evalue_ci is 1 when CI crosses null", {
  res <- sa_evalue(2.5, ci_lo = 0.9, ci_hi = 3.5)
  expect_equal(res$meta$evalue_ci, 1)
})

test_that("sa_evalue() evalue_ci > 1 when CI does not cross null (RR > 1)", {
  res <- sa_evalue(2.5, ci_lo = 1.5, ci_hi = 3.5)
  expect_gt(res$meta$evalue_ci, 1)
  # CI E-value should be based on ci_lo (closer to null for RR > 1)
  rr <- 1.5
  expected <- rr + sqrt(rr * (rr - 1))
  expect_equal(res$meta$evalue_ci, expected, tolerance = 1e-8)
})

# Hand-computed references: E(RR) = RR + sqrt(RR * (RR - 1)) for RR >= 1.
.ev_ref <- function(rr) {
  if (rr < 1) rr <- 1 / rr
  rr + sqrt(rr * (rr - 1))
}

test_that("sa_evalue() uses the OR directly for a rare outcome", {
  res <- sa_evalue(2.5, type = "OR", rare = TRUE)
  expect_equal(res$meta$rr_estimate, 2.5)
  expect_equal(res$meta$evalue_estimate, 2.5 + sqrt(2.5 * 1.5), tolerance = 1e-8)
  expect_equal(res$meta$evalue_estimate, 4.4365, tolerance = 1e-4)
  expect_true(res$meta$rare)
  expect_equal(res$meta$type, "OR")
})

test_that("sa_evalue() takes sqrt(OR) for a common outcome", {
  res <- sa_evalue(2.5, type = "OR", rare = FALSE)
  expect_equal(res$meta$rr_estimate, sqrt(2.5), tolerance = 1e-8)
  expect_equal(res$meta$evalue_estimate, .ev_ref(sqrt(2.5)), tolerance = 1e-8)
  expect_false(res$meta$rare)
})

test_that("sa_evalue() uses the HR directly for a rare outcome", {
  res <- sa_evalue(2.5, type = "HR", rare = TRUE)
  expect_equal(res$meta$rr_estimate, 2.5)
  expect_equal(res$meta$evalue_estimate, 4.4365, tolerance = 1e-4)
})

test_that("sa_evalue() applies the common-outcome HR conversion", {
  rr <- (1 - 0.5^sqrt(2.5)) / (1 - 0.5^sqrt(1 / 2.5))
  res <- sa_evalue(2.5, type = "HR", rare = FALSE)
  expect_equal(res$meta$rr_estimate, rr, tolerance = 1e-8)
  expect_equal(res$meta$rr_estimate, 1.8759, tolerance = 1e-4)
  expect_equal(res$meta$evalue_estimate, .ev_ref(rr), tolerance = 1e-8)
  expect_equal(res$meta$evalue_estimate, 3.1576, tolerance = 1e-4)
})

test_that("sa_evalue() handles protective OR and HR estimates", {
  expect_equal(sa_evalue(0.4, type = "OR", rare = TRUE)$meta$evalue_estimate,
               .ev_ref(1 / 0.4), tolerance = 1e-8)
  expect_equal(sa_evalue(0.4, type = "OR", rare = FALSE)$meta$evalue_estimate,
               .ev_ref(sqrt(0.4)), tolerance = 1e-8)
  rr <- (1 - 0.5^sqrt(0.4)) / (1 - 0.5^sqrt(1 / 0.4))
  res <- sa_evalue(0.4, type = "HR", rare = FALSE)
  expect_lt(res$meta$rr_estimate, 1)
  expect_equal(res$meta$evalue_estimate, .ev_ref(rr), tolerance = 1e-8)
  # The common-outcome HR conversion is symmetric under HR -> 1/HR
  expect_equal(res$meta$evalue_estimate,
               sa_evalue(2.5, type = "HR", rare = FALSE)$meta$evalue_estimate,
               tolerance = 1e-8)
})

test_that("sa_evalue() converts the CI bound on the same scale", {
  res <- sa_evalue(2.5, ci_lo = 1.6, ci_hi = 3.9, type = "OR", rare = FALSE)
  expect_equal(res$meta$evalue_ci, .ev_ref(sqrt(1.6)), tolerance = 1e-8)
  res <- sa_evalue(0.5, ci_lo = 0.3, ci_hi = 0.8, type = "HR", rare = FALSE)
  rr_hi <- (1 - 0.5^sqrt(0.8)) / (1 - 0.5^sqrt(1 / 0.8))
  expect_equal(res$meta$evalue_ci, .ev_ref(rr_hi), tolerance = 1e-8)
})

test_that("sa_evalue() does not infer rare from p0", {
  # Reviewer case: OR = 4, p0 = 0.10 gives exposed risk ~0.31; an equal-arm
  # cohort has ~20% prevalence, so p0 alone cannot establish rarity.
  expect_error(sa_evalue(4, type = "OR", p0 = 0.10), "rare = TRUE")
  expect_error(sa_evalue(2.5, type = "HR", p0 = 0.05), "cannot be inferred")
  # With rare given, p0 has no effect on OR or HR
  expect_equal(sa_evalue(4, type = "OR", p0 = 0.10, rare = FALSE)$meta$rr_estimate, 2)
  expect_equal(sa_evalue(4, type = "OR", p0 = 0.30, rare = TRUE)$meta$rr_estimate, 4)
  expect_no_warning(sa_evalue(4, type = "OR", p0 = 0.30, rare = TRUE))
})

test_that("sa_evalue() aborts for OR and HR when rarity is unknown", {
  expect_error(sa_evalue(2.5, type = "OR"), "rare")
  expect_error(sa_evalue(2.5, type = "HR"), "rare")
  expect_error(sa_evalue(2.5, type = "OR", rare = NA), "rare")
})

test_that("sa_evalue() ignores rare for RR, IRR and RD", {
  expect_true(is.na(sa_evalue(2.5)$meta$rare))
  expect_equal(sa_evalue(2.5, rare = FALSE)$meta$rr_estimate, 2.5)
})

test_that("sa_evalue() aborts on an inverted or excluding CI", {
  expect_error(sa_evalue(2, ci_lo = 3, ci_hi = 1.5), "ci_lo <= estimate <= ci_hi")
  expect_error(sa_evalue(2, ci_lo = 2.5, ci_hi = 3), "ci_lo <= estimate <= ci_hi")
  expect_error(sa_evalue(0.1, ci_lo = 0.12, ci_hi = 0.2, type = "RD", p0 = 0.2),
               "ci_lo <= estimate <= ci_hi")
})

test_that("sa_evalue() works for RD type with p0", {
  res <- sa_evalue(0.10, type = "RD", p0 = 0.20)
  expect_equal(res$meta$rr_estimate, (0.10 + 0.20) / 0.20, tolerance = 1e-8)
})

test_that("sa_evalue() errors on non-numeric estimate", {
  expect_error(sa_evalue("2.5"))
})

test_that("sa_evalue() errors when p0 missing for RD", {
  expect_error(sa_evalue(0.1, type = "RD"))
})

test_that("sa_evalue() errors when only one CI bound is supplied", {
  expect_error(sa_evalue(2.5, ci_lo = 1.5))
  expect_error(sa_evalue(2.5, ci_hi = 3.5))
})

test_that("sa_evalue() errors on non-positive OR", {
  expect_error(sa_evalue(-0.5, type = "OR"))
})


# ---------------------------------------------------------------------------
# sa_overlap
# ---------------------------------------------------------------------------

test_that("sa_overlap() returns list with expected names (plain data frame)", {
  dta <- sample_ps_data(n = 300, seed = 10)
  res <- sa_overlap(dta, score_col = "prob_t", treatment_col = "tavr")
  expect_s3_class(res, c("sa_overlap", "ps_data"))
  expect_identical(res$data, dta)
  expect_named(res$tables, c("summary", "outside_overlap", "positivity_flags"))
  expect_output(print(res), "Overlap")
  expect_output(summary(res), "Positivity flags")
})

test_that("sa_overlap() summary has two rows and correct columns", {
  dta <- sample_ps_data(n = 300, seed = 11)
  res <- sa_overlap(dta)
  expect_equal(nrow(res$tables$summary), 2L)
  expect_true(all(c("group", "n", "ps_min", "ps_max") %in% names(res$tables$summary)))
})

test_that("sa_overlap() overlap_region is a named numeric vector of length 2", {
  dta <- sample_ps_data(n = 300, seed = 12)
  res <- sa_overlap(dta)
  expect_named(res$meta$overlap_region, c("lower", "upper"))
  expect_true(res$meta$overlap_region["lower"] <= res$meta$overlap_region["upper"])
})

test_that("sa_overlap() pct_outside is in [0, 100]", {
  dta <- sample_ps_data(n = 300, seed = 13)
  res <- sa_overlap(dta)
  expect_true(all(res$tables$outside_overlap$pct_outside >= 0))
  expect_true(all(res$tables$outside_overlap$pct_outside <= 100))
})

test_that("sa_overlap() positivity_flags has correct structure", {
  dta <- sample_ps_data(n = 300, seed = 14)
  res <- sa_overlap(dta)
  pf  <- res$tables$positivity_flags
  expect_named(pf, c("group", "treated", "n_near_zero", "n_near_one",
                     "pct_near_zero", "pct_near_one"))
  expect_equal(nrow(pf), 2L)
})

test_that("sa_overlap() accepts a ps_data object", {
  m   <- make_matched(n = 200, seed = 15)
  res <- sa_overlap(m)
  expect_s3_class(res, "sa_overlap")
  expect_identical(res$data, m$data)
  expect_identical(res$meta$score_col, m$meta$score_col)
})

test_that("sa_overlap() summary n values sum to nrow(data)", {
  dta <- sample_ps_data(n = 300, seed = 16)
  res <- sa_overlap(dta)
  expect_equal(sum(res$tables$summary$n), nrow(dta))
})

test_that("sa_overlap() errors on bad trim_threshold", {
  dta <- sample_ps_data(n = 100, seed = 17)
  expect_error(sa_overlap(dta, trim_threshold = 0.6))
  expect_error(sa_overlap(dta, trim_threshold = -0.1))
})


# ---------------------------------------------------------------------------
# sa_trim_sweep
# ---------------------------------------------------------------------------

test_that("sa_trim_sweep() returns a data frame with correct columns", {
  w   <- make_weighted()
  obj <- sa_trim_sweep(w)
  expect_s3_class(obj, c("sa_trim_sweep", "ps_data"))
  expect_identical(obj$data, w$data)
  expect_identical(obj$meta$estimand, "ATE")
  expect_output(print(obj), "IPTW-ATE")
  expect_output(summary(obj), "Sweep")
  res <- obj$tables$sweep
  expect_s3_class(res, "data.frame")
  expect_true(all(c("trim", "ess_control", "ess_treated",
                    "max_weight", "sd_weight", "pct_trimmed") %in% names(res)))
})

test_that("sa_trim_sweep() has one row per element in trim_seq", {
  w    <- make_weighted()
  tseq <- seq(0, 0.05, by = 0.01)
  res  <- sa_trim_sweep(w, trim_seq = tseq)$tables$sweep
  expect_equal(nrow(res), length(tseq))
  expect_equal(res$trim, sort(tseq))
})

test_that("sa_trim_sweep() pct_trimmed is 0 at trim = 0", {
  w   <- make_weighted()
  res <- sa_trim_sweep(w, trim_seq = c(0, 0.05))$tables$sweep
  expect_equal(res$pct_trimmed[res$trim == 0], 0)
})

test_that("sa_trim_sweep() max_weight decreases (weakly) as trim increases", {
  w   <- make_weighted()
  res <- sa_trim_sweep(w, trim_seq = seq(0, 0.10, by = 0.01))$tables$sweep
  # max_weight should be non-increasing with trim
  expect_true(all(diff(res$max_weight) <= 1e-6))
})

test_that("sa_trim_sweep() accepts a plain data frame", {
  dta <- sample_ps_data(n = 300, seed = 20)
  res <- sa_trim_sweep(dta, trim_seq = c(0, 0.05),
                       score_col = "prob_t", treatment_col = "tavr")
  expect_equal(nrow(res$tables$sweep), 2L)
})

test_that("sa_trim_sweep() inherits estimand from ps_weight object", {
  w_att <- make_weighted(estimand = "ATT")
  res   <- sa_trim_sweep(w_att, trim_seq = c(0, 0.05))
  expect_identical(res$meta$estimand, "ATT")
  # The untrimmed row summarises the ATT weights themselves, not the ATE ones.
  untrimmed <- res$tables$sweep[1L, ]
  expect_equal(untrimmed$max_weight, round(max(w_att$data$iptw), 4))
  expect_equal(untrimmed$ess_treated, sum(w_att$data$tavr == 1L))
  expect_equal(res$tables$sweep,
               sa_trim_sweep(w_att$data, trim_seq = c(0, 0.05), estimand = "ATT")$tables$sweep)
})

test_that("sa_trim_sweep() errors on bad trim_seq values", {
  w <- make_weighted()
  expect_error(sa_trim_sweep(w, trim_seq = c(0, 0.6)))
  expect_error(sa_trim_sweep(w, trim_seq = c(-0.01, 0.05)))
})


# ---------------------------------------------------------------------------
# sa_rosenbaum (direct implementation; no external package required)
# ---------------------------------------------------------------------------

test_that("sa_rosenbaum() returns list with expected names", {
  m   <- make_matched(n = 300, seed = 42)
  res <- sa_rosenbaum(m, outcome_col = "ef", gamma_max = 2, gamma_inc = 0.5)
  expect_s3_class(res, c("sa_rosenbaum", "ps_data"))
  expect_named(res, c("data", "meta", "tables", "models"))
  expect_identical(res$data, m$data)
  expect_named(res$tables, c("bounds", "pairs"))
  expect_equal(nrow(res$tables$pairs), res$meta$n_pairs)
  expect_equal(res$meta$gamma, c(1, 1.5, 2))
  expect_output(print(res), "Sensitivity")
  expect_output(summary(res), "Bounds")
})

test_that("sa_rosenbaum() bounds data frame has correct columns and rows", {
  m   <- make_matched(n = 300, seed = 43)
  res <- sa_rosenbaum(m, outcome_col = "ef", gamma_max = 2, gamma_inc = 0.5)
  expect_named(res$tables$bounds, c("gamma", "p_upper", "p_lower", "reject_upper"))
  # gamma values: 1.0, 1.5, 2.0 -> 3 rows
  expect_equal(nrow(res$tables$bounds), 3L)
  expect_equal(res$tables$bounds$gamma[1L], 1.0)
  expect_equal(res$tables$bounds$gamma[3L], 2.0)
})

test_that("sa_rosenbaum() n_pairs is a positive integer", {
  m   <- make_matched(n = 300, seed = 44)
  res <- sa_rosenbaum(m, outcome_col = "ef", gamma_max = 1.5, gamma_inc = 0.5)
  expect_true(res$meta$n_pairs > 0L)
  expect_true(is.numeric(res$meta$n_pairs))
})

test_that("sa_rosenbaum() p_upper >= p_lower for all Gamma (monotone bounds)", {
  m   <- make_matched(n = 300, seed = 49)
  res <- sa_rosenbaum(m, outcome_col = "ef", gamma_max = 3, gamma_inc = 0.5)
  expect_true(all(res$tables$bounds$p_upper >= res$tables$bounds$p_lower - 1e-10))
})

test_that("sa_rosenbaum() p_upper is non-decreasing in Gamma", {
  m   <- make_matched(n = 300, seed = 50)
  res <- sa_rosenbaum(m, outcome_col = "ef", gamma_max = 3, gamma_inc = 0.5)
  expect_true(all(diff(res$tables$bounds$p_upper) >= -1e-10))
})

test_that("sa_rosenbaum() sensitivity_value is NA or a Gamma in the tested range", {
  m   <- make_matched(n = 300, seed = 48)
  res <- sa_rosenbaum(m, outcome_col = "ef", gamma_max = 2, gamma_inc = 0.5,
                      alpha = 0.05)
  sv  <- res$meta$sensitivity_value
  expect_true(is.na(sv) || (is.numeric(sv) && sv >= 1 && sv <= 2))
})

test_that("sa_rosenbaum() accepts a plain data frame with treatment_col", {
  m       <- make_matched(n = 300, seed = 46)
  matched <- m$data[!is.na(m$data$pair_id), ]
  res     <- sa_rosenbaum(matched, outcome_col = "ef",
                          treatment_col = "tavr", gamma_max = 1.5,
                          gamma_inc = 0.5)
  expect_s3_class(res, "sa_rosenbaum")
  expect_named(res$tables$bounds, c("gamma", "p_upper", "p_lower", "reject_upper"))
})

test_that("sa_rosenbaum() errors when treatment_col missing for plain data frame", {
  m       <- make_matched(n = 200, seed = 47)
  matched <- m$data[!is.na(m$data$pair_id), ]
  expect_error(sa_rosenbaum(matched, outcome_col = "ef"))
})


# ---------------------------------------------------------------------------
# NEW GUARD TESTS: sa_rosenbaum()
# ---------------------------------------------------------------------------

test_that("sa_rosenbaum() errors when gamma_max < 1", {
  m <- make_matched(n = 200, seed = 51)
  expect_error(sa_rosenbaum(m, outcome_col = "ef", gamma_max = 0.5),
               regexp = "gamma_max")
})

test_that("sa_rosenbaum() errors when gamma_inc <= 0", {
  m <- make_matched(n = 200, seed = 52)
  expect_error(sa_rosenbaum(m, outcome_col = "ef", gamma_inc = -0.1),
               regexp = "gamma_inc")
  expect_error(sa_rosenbaum(m, outcome_col = "ef", gamma_inc = 0),
               regexp = "gamma_inc")
})

test_that("sa_rosenbaum() errors when outcome_col is a factor", {
  m    <- make_matched(n = 200, seed = 53)
  dta  <- m$data
  dta$ef_factor <- factor(round(dta$ef))
  m$data <- dta
  expect_error(sa_rosenbaum(m, outcome_col = "ef_factor"),
               regexp = "numeric")
})

test_that("sa_rosenbaum() warns when NA pairs are dropped", {
  m   <- make_matched(n = 300, seed = 54)
  dta <- m$data
  # Introduce NAs into the outcome for a few matched patients
  matched_idx <- which(!is.na(dta$pair_id))[1:4]
  dta$ef[matched_idx] <- NA
  m$data <- dta
  expect_warning(sa_rosenbaum(m, outcome_col = "ef", gamma_max = 1.5),
                 regexp = "dropped")
})

test_that("sa_rosenbaum() seq gamma values are exact (float accumulation check)", {
  m   <- make_matched(n = 300, seed = 55)
  res <- sa_rosenbaum(m, outcome_col = "ef", gamma_max = 3, gamma_inc = 0.25)
  # All gamma values should be representable to at least 8 decimal places
  expect_equal(res$tables$bounds$gamma[nrow(res$tables$bounds)], 3.0, tolerance = 1e-8)
})

test_that("sa_rosenbaum() warns and re-orients when median difference is negative", {
  m    <- make_matched(n = 300, seed = 56)
  dta  <- m$data
  # Flip the outcome so treated < control (negative effect)
  dta$ef_neg <- -dta$ef
  m$data <- dta
  expect_warning(sa_rosenbaum(m, outcome_col = "ef_neg", gamma_max = 2),
                 regexp = "negative")
})

test_that("ps_match() seed makes the pair_id that sa_rosenbaum() reads reproducible", {
  dta <- sample_ps_data(n = 300, seed = 99)
  m1  <- ps_match(dta, seed = 42)
  m2  <- ps_match(dta, seed = 42)
  m3  <- ps_match(dta, seed = 7)
  expect_equal(m1$data$pair_id, m2$data$pair_id)
  # Different seeds should (almost certainly) produce different matches
  # for data with ties in the PS distribution
  expect_false(identical(m1$data$pair_id, m3$data$pair_id))
})

test_that("ps_match() warns when pair_id column already exists in input data", {
  dta          <- sample_ps_data(n = 200, seed = 60)
  dta$pair_id  <- seq_len(nrow(dta))   # pre-existing column
  expect_warning(ps_match(dta, seed = 1), regexp = "pair_id")
})


# ---------------------------------------------------------------------------
# NEW GUARD TESTS: sa_evalue()
# ---------------------------------------------------------------------------

test_that("sa_evalue() errors when RD + p0 <= 0", {
  # RD = -0.25, p0 = 0.20 -> implied p1 = -0.05 (impossible)
  expect_error(sa_evalue(-0.25, type = "RD", p0 = 0.20),
               regexp = "non-positive")
})

test_that("sa_evalue() errors when RD + p0 > 1", {
  # RD = 0.85, p0 = 0.20 -> implied p1 = 1.05 (impossible)
  expect_error(sa_evalue(0.85, type = "RD", p0 = 0.20),
               regexp = "exceeds 1")
})

test_that("sa_evalue() errors on non-positive CI bound for RR", {
  expect_error(sa_evalue(2.0, ci_lo = 0, ci_hi = 3.0, type = "RR"),
               regexp = "positive")
})

test_that("sa_evalue() does NOT warn for OR with rare given", {
  expect_no_warning(sa_evalue(2.0, type = "OR", rare = FALSE))
  expect_no_warning(sa_evalue(2.0, type = "OR", rare = TRUE))
})

test_that("sa_evalue() supports type = 'IRR' with the same formula as RR", {
  rr_res  <- sa_evalue(2.5, type = "RR")
  irr_res <- sa_evalue(2.5, type = "IRR")
  expect_equal(irr_res$meta$evalue_estimate, rr_res$meta$evalue_estimate)
  expect_equal(irr_res$meta$rr_estimate,     rr_res$meta$rr_estimate)
  expect_equal(irr_res$meta$type, "IRR")
})

test_that("sa_evalue() null-crossing check works for near-null estimate (float safe)", {
  # estimate very close to 1 should be treated as null (evalue = 1)
  res <- sa_evalue(1 + 1e-10, ci_lo = 0.8, ci_hi = 1.2)
  expect_equal(res$meta$evalue_ci, 1)
})


# ---------------------------------------------------------------------------
# NEW GUARD TESTS: sa_overlap()
# ---------------------------------------------------------------------------

test_that("sa_overlap() errors on non-binary treatment column", {
  dta          <- sample_ps_data(n = 200, seed = 70)
  dta$bad_trt  <- dta$tavr + 1L   # values 1 and 2, not 0/1
  expect_error(sa_overlap(dta, treatment_col = "bad_trt"),
               regexp = "binary")
})

test_that("sa_overlap() warns and excludes patients with NA propensity score", {
  dta          <- sample_ps_data(n = 300, seed = 71)
  dta$prob_t[1:5] <- NA
  expect_warning(sa_overlap(dta), regexp = "missing")
  # After exclusion the n counts should sum to nrow(dta) - 5
  res <- suppressWarnings(sa_overlap(dta))
  expect_equal(sum(res$tables$summary$n), nrow(dta) - 5L)
})

test_that("sa_overlap() warns when PS distributions do not overlap", {
  dta           <- sample_ps_data(n = 200, seed = 72)
  # Force complete separation: control gets PS near 0, treated gets PS near 1
  dta$prob_t[dta$tavr == 0L] <- runif(sum(dta$tavr == 0L), 0.01, 0.10)
  dta$prob_t[dta$tavr == 1L] <- runif(sum(dta$tavr == 1L), 0.90, 0.99)
  expect_warning(sa_overlap(dta), regexp = "common support")
})

test_that("sa_overlap() warns when a group has fewer than 2 patients", {
  dta        <- sample_ps_data(n = 200, seed = 73)
  # Keep only 1 treated patient + all controls
  one_treated <- which(dta$tavr == 1L)[1L]
  all_control  <- which(dta$tavr == 0L)
  small_dta   <- dta[c(one_treated, all_control), ]
  expect_warning(sa_overlap(small_dta), regexp = "fewer than 2")
})


# ---------------------------------------------------------------------------
# NEW GUARD TESTS: sa_trim_sweep()
# ---------------------------------------------------------------------------

test_that("sa_trim_sweep() warns when ps = 0 or ps = 1 produces Inf weights", {
  dta            <- sample_ps_data(n = 200, seed = 80)
  # Inject a ps = 0 for a treated patient -> Inf ATE weight
  dta$prob_t[dta$tavr == 1L][1L] <- 0
  expect_warning(sa_trim_sweep(dta, score_col = "prob_t",
                               treatment_col = "tavr",
                               trim_seq = c(0, 0.05)),
                 regexp = "infinite")
})

test_that("sa_trim_sweep() treats a trim below machine epsilon as no trim", {
  # 2e-16 is below .Machine$double.eps, so the guard applies. Without it the
  # quantile() bounds at 2e-16 and 1 - 2e-16 would still clip two weights.
  w   <- make_weighted(seed = 81)
  tr  <- 2e-16
  expect_lt(tr, .Machine$double.eps)
  raw <- w$data$iptw
  expect_gt(sum(raw < stats::quantile(raw, tr) | raw > stats::quantile(raw, 1 - tr)), 0L)
  res <- sa_trim_sweep(w, trim_seq = c(0, tr))$tables$sweep
  expect_equal(res$pct_trimmed[2L], 0)
  expect_equal(res[2L, names(res) != "trim"], res[1L, names(res) != "trim"], ignore_attr = TRUE)
})

test_that("sa_overlap() tables name groups by their treatment value", {
  dta <- sample_ps_data(n = 300, seed = 15)
  res <- sa_overlap(dta)
  for (tbl in res$tables[c("summary", "outside_overlap", "positivity_flags")]) {
    expect_identical(tbl$group, c("0", "1"))
    expect_identical(tbl$treated, c(FALSE, TRUE))
  }
  expect_equal(res$tables$summary$n, c(sum(dta$tavr == 0), sum(dta$tavr == 1)))
})

test_that("sa_overlap() follows the treated level a scored object declares", {
  dta <- sample_ps_data(n = 200, seed = 16)
  dta$prob_t <- NULL
  fit <- ps_logistic(tavr ~ age + ef, dta, treatment_levels = c(1, 0), treated_level = 0)
  res <- sa_overlap(fit)
  expect_identical(res$tables$summary$group, c("1", "0"))
  expect_identical(res$tables$summary$treated, c(FALSE, TRUE))
  expect_equal(res$tables$summary$n, c(sum(dta$tavr == 1), sum(dta$tavr == 0)))
  # The treated row summarises the score of the patients whose value is 0.
  expect_equal(res$tables$summary$ps_mean[2], round(mean(fit$data$prob_t[fit$data$tavr == 0]), 4))
})


# ---------------------------------------------------------------------------
# sa_rosenbaum(): reference values, ties and direction
# ---------------------------------------------------------------------------

# One pair per element of d, with control outcome 10 and treated outcome 10 + d.
rosenbaum_pairs <- function(d) {
  data.frame(pair_id = rep(seq_along(d), each = 2L),
             tavr    = rep(c(1L, 0L), length(d)),
             ef      = as.vector(rbind(10 + d, 10)))
}

test_that("sa_rosenbaum() at Gamma = 1 equals the one-sided wilcox.test() p-value, with ties", {
  d   <- c(1, 1, -2, 3, 3, 3, 4, -5)
  dta <- rosenbaum_pairs(d)
  res <- sa_rosenbaum(dta, "ef", treatment_col = "tavr", gamma_max = 1, alternative = "greater")
  ref <- stats::wilcox.test(10 + d, rep(10, length(d)), paired = TRUE, alternative = "greater",
                            exact = FALSE, correct = TRUE)$p.value
  expect_equal(res$tables$bounds$p_upper, ref, tolerance = 1e-12)
})

test_that("sa_rosenbaum() uses the tie-corrected variance sum(r^2)", {
  # |d| ranks 1.5 1.5 3 5 5 5 7 8; T+ = 25, sum(r) = 36, sum(r^2) = 201.5.
  dta <- rosenbaum_pairs(c(1, 1, -2, 3, 3, 3, 4, -5))
  res <- sa_rosenbaum(dta, "ef", treatment_col = "tavr", gamma_max = 2, gamma_inc = 1,
                      alternative = "greater")
  sd2 <- sqrt(2 / 9 * 201.5)
  expect_equal(res$meta$t_plus, 25)
  expect_equal(res$tables$bounds$p_upper[2], 1 - stats::pnorm((24.5 - 24) / sd2), tolerance = 1e-12)
  expect_equal(res$tables$bounds$p_lower[2], 1 - stats::pnorm((24.5 - 12) / sd2), tolerance = 1e-12)
})

test_that("sa_rosenbaum() without ties matches the n(n+1)(2n+1)/6 variance", {
  # T+ = 1 + 2 + 4 + 5 + 6 + 8 = 26; n = 8. At Gamma = 1.5, p = 0.6:
  # mean 0.6 * 36 = 21.6, variance 1.5 / 2.5^2 * 8 * 9 * 17 / 6 = 48.96.
  dta <- rosenbaum_pairs(c(1, 2, -3, 4, 5, 6, -7, 8))
  res <- sa_rosenbaum(dta, "ef", treatment_col = "tavr", gamma_max = 1.5, gamma_inc = 0.5)
  expect_equal(res$tables$bounds$p_upper[2], 1 - stats::pnorm((25.5 - 21.6) / sqrt(48.96)),
               tolerance = 1e-12)
})

test_that("sa_rosenbaum() alternative fixes the direction in advance", {
  dta <- rosenbaum_pairs(-c(1, 2, -3, 4, 5, 6, -7, 8))
  expect_warning(auto <- sa_rosenbaum(dta, "ef", treatment_col = "tavr", gamma_max = 2), "negative")
  expect_no_warning(less <- sa_rosenbaum(dta, "ef", treatment_col = "tavr", gamma_max = 2,
                                         alternative = "less"))
  expect_equal(auto$tables$bounds, less$tables$bounds)
  expect_true(auto$meta$direction_from_data)
  expect_false(less$meta$direction_from_data)
  expect_identical(less$meta$alternative, "less")
  expect_no_warning(greater <- sa_rosenbaum(dta, "ef", treatment_col = "tavr", gamma_max = 2,
                                            alternative = "greater"))
  expect_gt(greater$tables$bounds$p_lower[1], 0.5)
})


# ---------------------------------------------------------------------------
# One argument-precedence rule for every sa_*() function
# ---------------------------------------------------------------------------

test_that("an explicit argument that conflicts with the object's metadata aborts", {
  m <- make_matched(n = 300, seed = 42)
  w <- make_weighted(estimand = "ATE")
  expect_error(sa_rosenbaum(m, "ef", treated_level = 0), "conflicts with the object's metadata")
  expect_error(sa_rosenbaum(m, "ef", treatment_col = "female"), "treatment_col")
  expect_error(sa_rosenbaum(m, "ef", pair_id_col = "id"), "pair_id_col")
  expect_error(sa_trim_sweep(w, estimand = "ATT"), "`estimand = \"ATT\"`.*\\$meta\\$estimand = \"ATE\"")
  expect_error(sa_trim_sweep(w, stabilise = !w$meta$stabilised), "meta\\$stabilised")
  expect_error(sa_overlap(m, score_col = "other"), "score_col")
  expect_error(sa_overlap(m, treated_level = 0), "treated_level")
})

test_that("an explicit argument that agrees with the metadata gives the default result", {
  m <- make_matched(n = 300, seed = 42)
  w <- make_weighted(estimand = "ATE")
  expect_equal(sa_rosenbaum(m, "ef", treated_level = 1, treatment_col = "tavr"),
               sa_rosenbaum(m, "ef"))
  expect_equal(sa_trim_sweep(w, estimand = "ATE", stabilise = w$meta$stabilised), sa_trim_sweep(w))
  expect_equal(sa_overlap(m, score_col = m$meta$score_col), sa_overlap(m))
})

test_that("unsupplied arguments take the object's metadata", {
  w_att <- make_weighted(estimand = "ATT")
  expect_identical(sa_trim_sweep(w_att)$meta$estimand, "ATT")
  expect_identical(sa_trim_sweep(w_att)$meta$stabilised, w_att$meta$stabilised)
})


# ---------------------------------------------------------------------------
# Missing treatment or score (issue #84)
# ---------------------------------------------------------------------------

test_that("sa_overlap() excludes patients with a missing treatment", {
  d <- sample_ps_data(n = 200, seed = 90)
  d$tavr[1:3] <- NA
  res <- suppressWarnings(sa_overlap(d))
  expect_equal(sum(res$tables$summary$n), nrow(d) - 3L)
  expect_equal(res$meta$n_used, nrow(d) - 3L)
  expect_false(anyNA(res$tables$summary))
})

test_that("sa_trim_sweep() excludes patients with a missing score, with a warning", {
  d <- sample_ps_data(n = 200, seed = 91)
  d$prob_t[1:2] <- NA
  expect_warning(res <- sa_trim_sweep(d), "2 patient\\(s\\) with missing 'prob_t'")
  sweep <- res$tables$sweep
  expect_true(all(is.finite(sweep$ess_control)))
  expect_true(all(is.finite(sweep$ess_treated)))
  expect_equal(res$meta$n_used, nrow(d) - 2L)
  expect_equal(sweep, sa_trim_sweep(d[-(1:2), ])$tables$sweep)
})
