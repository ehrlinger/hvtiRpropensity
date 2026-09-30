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
# ps_match — pair_id column (prerequisite for sa_rosenbaum)
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

test_that("sa_evalue() returns list with expected names for RR", {
  res <- sa_evalue(2.5)
  expect_named(res, c("evalue_estimate", "evalue_ci", "rr_estimate", "type"))
})

test_that("sa_evalue() E-value for RR = 1 is 1", {
  res <- sa_evalue(1.0)
  expect_equal(res$evalue_estimate, 1)
})

test_that("sa_evalue() E-value for RR = 2 is 2 + sqrt(2)", {
  res <- sa_evalue(2.0)
  expect_equal(res$evalue_estimate, 2 + sqrt(2 * (2 - 1)), tolerance = 1e-8)
})

test_that("sa_evalue() is symmetric: RR and 1/RR give the same E-value", {
  e1 <- sa_evalue(2.0)$evalue_estimate
  e2 <- sa_evalue(0.5)$evalue_estimate
  expect_equal(e1, e2, tolerance = 1e-8)
})

test_that("sa_evalue() evalue_ci is NA when no CI is supplied", {
  res <- sa_evalue(2.5)
  expect_true(is.na(res$evalue_ci))
})

test_that("sa_evalue() evalue_ci is 1 when CI crosses null", {
  res <- sa_evalue(2.5, ci_lo = 0.9, ci_hi = 3.5)
  expect_equal(res$evalue_ci, 1)
})

test_that("sa_evalue() evalue_ci > 1 when CI does not cross null (RR > 1)", {
  res <- sa_evalue(2.5, ci_lo = 1.5, ci_hi = 3.5)
  expect_gt(res$evalue_ci, 1)
  # CI E-value should be based on ci_lo (closer to null for RR > 1)
  rr <- 1.5
  expected <- rr + sqrt(rr * (rr - 1))
  expect_equal(res$evalue_ci, expected, tolerance = 1e-8)
})

test_that("sa_evalue() works for OR type with rare outcome approximation", {
  res <- sa_evalue(4.0, type = "OR")
  expect_equal(res$rr_estimate, sqrt(4), tolerance = 1e-8)
  expect_equal(res$type, "OR")
})

test_that("sa_evalue() works for HR type", {
  res <- sa_evalue(3.0, type = "HR")
  expect_equal(res$rr_estimate, sqrt(3), tolerance = 1e-8)
})

test_that("sa_evalue() works for RD type with p0", {
  res <- sa_evalue(0.10, type = "RD", p0 = 0.20)
  expect_equal(res$rr_estimate, (0.10 + 0.20) / 0.20, tolerance = 1e-8)
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
  expect_named(res, c("summary", "overlap_region", "outside_overlap",
                      "positivity_flags"))
})

test_that("sa_overlap() summary has two rows and correct columns", {
  dta <- sample_ps_data(n = 300, seed = 11)
  res <- sa_overlap(dta)
  expect_equal(nrow(res$summary), 2L)
  expect_true(all(c("group", "n", "ps_min", "ps_max") %in% names(res$summary)))
})

test_that("sa_overlap() overlap_region is a named numeric vector of length 2", {
  dta <- sample_ps_data(n = 300, seed = 12)
  res <- sa_overlap(dta)
  expect_named(res$overlap_region, c("lower", "upper"))
  expect_true(res$overlap_region["lower"] <= res$overlap_region["upper"])
})

test_that("sa_overlap() pct_outside is in [0, 100]", {
  dta <- sample_ps_data(n = 300, seed = 13)
  res <- sa_overlap(dta)
  expect_true(all(res$outside_overlap$pct_outside >= 0))
  expect_true(all(res$outside_overlap$pct_outside <= 100))
})

test_that("sa_overlap() positivity_flags has correct structure", {
  dta <- sample_ps_data(n = 300, seed = 14)
  res <- sa_overlap(dta)
  pf  <- res$positivity_flags
  expect_named(pf, c("group", "treated", "n_near_zero", "n_near_one",
                     "pct_near_zero", "pct_near_one"))
  expect_equal(nrow(pf), 2L)
})

test_that("sa_overlap() accepts a ps_data object", {
  m   <- make_matched(n = 200, seed = 15)
  res <- sa_overlap(m)
  expect_named(res, c("summary", "overlap_region", "outside_overlap",
                      "positivity_flags"))
})

test_that("sa_overlap() summary n values sum to nrow(data)", {
  dta <- sample_ps_data(n = 300, seed = 16)
  res <- sa_overlap(dta)
  expect_equal(sum(res$summary$n), nrow(dta))
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
  res <- sa_trim_sweep(w)
  expect_s3_class(res, "data.frame")
  expect_true(all(c("trim", "ess_control", "ess_treated",
                    "max_weight", "sd_weight", "pct_trimmed") %in% names(res)))
})

test_that("sa_trim_sweep() has one row per element in trim_seq", {
  w    <- make_weighted()
  tseq <- seq(0, 0.05, by = 0.01)
  res  <- sa_trim_sweep(w, trim_seq = tseq)
  expect_equal(nrow(res), length(tseq))
  expect_equal(res$trim, sort(tseq))
})

test_that("sa_trim_sweep() pct_trimmed is 0 at trim = 0", {
  w   <- make_weighted()
  res <- sa_trim_sweep(w, trim_seq = c(0, 0.05))
  expect_equal(res$pct_trimmed[res$trim == 0], 0)
})

test_that("sa_trim_sweep() max_weight decreases (weakly) as trim increases", {
  w   <- make_weighted()
  res <- sa_trim_sweep(w, trim_seq = seq(0, 0.10, by = 0.01))
  # max_weight should be non-increasing with trim
  expect_true(all(diff(res$max_weight) <= 1e-6))
})

test_that("sa_trim_sweep() accepts a plain data frame", {
  dta <- sample_ps_data(n = 300, seed = 20)
  res <- sa_trim_sweep(dta, trim_seq = c(0, 0.05),
                       score_col = "prob_t", treatment_col = "tavr")
  expect_equal(nrow(res), 2L)
})

test_that("sa_trim_sweep() inherits estimand from ps_weight object", {
  w_att <- make_weighted(estimand = "ATT")
  # Should run without error; the ATT estimand is picked up from meta
  res   <- sa_trim_sweep(w_att, trim_seq = c(0, 0.05))
  expect_equal(nrow(res), 2L)
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
  expect_named(res, c("bounds", "sensitivity_value", "n_pairs"))
})

test_that("sa_rosenbaum() bounds data frame has correct columns and rows", {
  m   <- make_matched(n = 300, seed = 43)
  res <- sa_rosenbaum(m, outcome_col = "ef", gamma_max = 2, gamma_inc = 0.5)
  expect_named(res$bounds, c("gamma", "p_upper", "p_lower", "reject_upper"))
  # gamma values: 1.0, 1.5, 2.0 -> 3 rows
  expect_equal(nrow(res$bounds), 3L)
  expect_equal(res$bounds$gamma[1L], 1.0)
  expect_equal(res$bounds$gamma[3L], 2.0)
})

test_that("sa_rosenbaum() n_pairs is a positive integer", {
  m   <- make_matched(n = 300, seed = 44)
  res <- sa_rosenbaum(m, outcome_col = "ef", gamma_max = 1.5, gamma_inc = 0.5)
  expect_true(res$n_pairs > 0L)
  expect_true(is.numeric(res$n_pairs))
})

test_that("sa_rosenbaum() p_upper and p_lower at Gamma=1 are equal and valid", {
  # At Gamma=1 (no hidden bias), upper and lower bounds coincide with the
  # unadjusted normal approximation to the Wilcoxon signed-rank p-value.
  m   <- make_matched(n = 300, seed = 45)
  res <- sa_rosenbaum(m, outcome_col = "ef", gamma_max = 1, gamma_inc = 1)
  p_up <- res$bounds$p_upper[1L]
  p_lo <- res$bounds$p_lower[1L]
  expect_equal(p_up, p_lo, tolerance = 1e-8)
  expect_true(p_up >= 0 && p_up <= 1)
})

test_that("sa_rosenbaum() p_upper >= p_lower for all Gamma (monotone bounds)", {
  m   <- make_matched(n = 300, seed = 49)
  res <- sa_rosenbaum(m, outcome_col = "ef", gamma_max = 3, gamma_inc = 0.5)
  expect_true(all(res$bounds$p_upper >= res$bounds$p_lower - 1e-10))
})

test_that("sa_rosenbaum() p_upper is non-decreasing in Gamma", {
  m   <- make_matched(n = 300, seed = 50)
  res <- sa_rosenbaum(m, outcome_col = "ef", gamma_max = 3, gamma_inc = 0.5)
  expect_true(all(diff(res$bounds$p_upper) >= -1e-10))
})

test_that("sa_rosenbaum() sensitivity_value is NA or a Gamma in the tested range", {
  m   <- make_matched(n = 300, seed = 48)
  res <- sa_rosenbaum(m, outcome_col = "ef", gamma_max = 2, gamma_inc = 0.5,
                      alpha = 0.05)
  sv  <- res$sensitivity_value
  expect_true(is.na(sv) || (is.numeric(sv) && sv >= 1 && sv <= 2))
})

test_that("sa_rosenbaum() accepts a plain data frame with treatment_col", {
  m       <- make_matched(n = 300, seed = 46)
  matched <- m$data[!is.na(m$data$pair_id), ]
  res     <- sa_rosenbaum(matched, outcome_col = "ef",
                          treatment_col = "tavr", gamma_max = 1.5,
                          gamma_inc = 0.5)
  expect_named(res, c("bounds", "sensitivity_value", "n_pairs"))
})

test_that("sa_rosenbaum() errors when treatment_col missing for plain data frame", {
  m       <- make_matched(n = 200, seed = 47)
  matched <- m$data[!is.na(m$data$pair_id), ]
  expect_error(sa_rosenbaum(matched, outcome_col = "ef"))
})


# ---------------------------------------------------------------------------
# NEW GUARD TESTS — sa_rosenbaum()
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
  expect_equal(res$bounds$gamma[nrow(res$bounds)], 3.0, tolerance = 1e-8)
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

test_that("sa_rosenbaum() seed argument makes ps_match result reproducible", {
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
# NEW GUARD TESTS — sa_evalue()
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

test_that("sa_evalue() warns on OR with high prevalence outcome (p0 > 0.15)", {
  expect_warning(sa_evalue(2.0, type = "OR", p0 = 0.30),
                 regexp = "rare-outcome")
})

test_that("sa_evalue() does NOT warn for OR when p0 <= 0.15", {
  expect_no_warning(sa_evalue(2.0, type = "OR", p0 = 0.10))
})

test_that("sa_evalue() supports type = 'IRR' with the same formula as RR", {
  rr_res  <- sa_evalue(2.5, type = "RR")
  irr_res <- sa_evalue(2.5, type = "IRR")
  expect_equal(irr_res$evalue_estimate, rr_res$evalue_estimate)
  expect_equal(irr_res$rr_estimate,     rr_res$rr_estimate)
  expect_equal(irr_res$type, "IRR")
})

test_that("sa_evalue() null-crossing check works for near-null estimate (float safe)", {
  # estimate very close to 1 should be treated as null (evalue = 1)
  res <- sa_evalue(1 + 1e-10, ci_lo = 0.8, ci_hi = 1.2)
  expect_equal(res$evalue_ci, 1)
})


# ---------------------------------------------------------------------------
# NEW GUARD TESTS — sa_overlap()
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
  expect_equal(sum(res$summary$n), nrow(dta) - 5L)
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
# NEW GUARD TESTS — sa_trim_sweep()
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

test_that("sa_trim_sweep() pct_trimmed is 0 for all near-zero trim values", {
  # seq() may produce very small non-zero values instead of exact 0
  w   <- make_weighted(seed = 81)
  res <- sa_trim_sweep(w, trim_seq = seq(0, 0.10, by = 0.01))
  expect_equal(res$pct_trimmed[1L], 0)
})

test_that("sa_overlap() tables name groups by their treatment value", {
  dta <- sample_ps_data(n = 300, seed = 15)
  res <- sa_overlap(dta)
  for (tbl in res[c("summary", "outside_overlap", "positivity_flags")]) {
    expect_identical(tbl$group, c("0", "1"))
    expect_identical(tbl$treated, c(FALSE, TRUE))
  }
  expect_equal(res$summary$n, c(sum(dta$tavr == 0), sum(dta$tavr == 1)))
})

test_that("sa_overlap() follows the treated level a scored object declares", {
  dta <- sample_ps_data(n = 200, seed = 16)
  dta$prob_t <- NULL
  fit <- ps_logistic(tavr ~ age + ef, dta, treatment_levels = c(1, 0), treated_level = 0)
  res <- sa_overlap(fit)
  expect_identical(res$summary$group, c("1", "0"))
  expect_identical(res$summary$treated, c(FALSE, TRUE))
  expect_equal(res$summary$n, c(sum(dta$tavr == 1), sum(dta$tavr == 0)))
  # The treated row summarises the score of the patients whose value is 0.
  expect_equal(res$summary$ps_mean[2], round(mean(fit$data$prob_t[fit$data$tavr == 0]), 4))
})
