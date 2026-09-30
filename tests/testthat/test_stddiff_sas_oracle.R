###############################################################################
## test_stddiff_sas_oracle.R
##
## ps_stddiff(), ps_stddiff_perm() and ps_mw_var() against the CCF macros
## %stddiff, %stddiffci and %mw_var, run by fixtures/stddiff-sas/
## generate-fixtures.sas on the synthetic fixtures/stddiff-sas/input.csv.
## The SAS output is the oracle; see fixtures/stddiff-sas/README.md. Issue #34.
###############################################################################

oracle_dir <- test_path("fixtures", "stddiff-sas")
oracle_files <- c("sd-unweighted.csv", "sd-weighted.csv", "ci-unweighted.csv", "ci-weighted.csv",
                  "mw-summary.csv", "mw-replicates.csv")

skip_without_oracle <- function(files = oracle_files) {
  missing <- files[!file.exists(file.path(oracle_dir, files))]
  if (length(missing) > 0L) {
    testthat::skip(paste("SAS oracle output not yet returned:", paste(missing, collapse = ", ")))
  }
}

read_oracle <- function(file) {
  d <- utils::read.csv(file.path(oracle_dir, file), check.names = FALSE, na.strings = c("", "."),
                       stringsAsFactors = FALSE)
  names(d) <- tolower(names(d))
  d
}

fixture <- utils::read.csv(file.path(oracle_dir, "input.csv"), na.strings = "", stringsAsFactors = FALSE)
variables <- c(gaussian = "x_gauss", nong_ord = "x_ord", binary = "x_bin", categorical = "x_cat")
n_perm <- 20L

r_stddiff <- function(group = "grp", weight = NULL) {
  tbl <- ps_stddiff(fixture, treatment_col = group, gaussian = variables[["gaussian"]],
                    nong_ord = variables[["nong_ord"]], binary = variables[["binary"]],
                    categorical = variables[["categorical"]], weight_col = weight)$tables$stddiff
  stats::setNames(tbl$stddiff, tbl$variable)
}

sas_by_variable <- function(d, value) {
  stats::setNames(d[[value]], tolower(d$varname))[unname(variables)]
}

# ---- %stddiff ---------------------------------------------------------------

test_that("ps_stddiff() matches %stddiff to 4 decimals, unweighted", {
  skip_without_oracle("sd-unweighted.csv")
  sas <- sas_by_variable(read_oracle("sd-unweighted.csv"), "stddiff")
  expect_equal(r_stddiff()[unname(variables)], sas, tolerance = 5e-5, ignore_attr = TRUE)
})

test_that("ps_stddiff() matches %stddiff to 4 decimals, weighted", {
  skip_without_oracle("sd-weighted.csv")
  sas <- sas_by_variable(read_oracle("sd-weighted.csv"), "stddiff")
  expect_equal(r_stddiff(weight = "w")[unname(variables)], sas, tolerance = 5e-5, ignore_attr = TRUE)
})

# ---- %stddiffci -------------------------------------------------------------
# %stddiffci takes its permutations from the input (fgrp_1..fgrp_20), so both
# sides score the same permuted groups and the percentiles can be compared
# exactly. ps_stddiff_perm() draws its own permutations, so its per-permutation
# step is reproduced here on the shared ones: put the permuted group in the
# treatment column, re-derive the weights with `reweight`, score, then take
# .perm_percentiles(type = 2).

matching_weight <- function(p, g) round(pmin(p, 1 - p) / (p * g + (1 - p) * (1 - g)), 6)
reweight <- function(d) matching_weight(d$p_score, d$grp)

permuted_stddiff <- function(k, weighted) {
  d <- fixture
  d$grp <- d[[paste0("fgrp_", k)]]
  if (weighted) d$w <- reweight(d)
  tbl <- ps_stddiff(d, treatment_col = "grp", gaussian = variables[["gaussian"]], nong_ord = variables[["nong_ord"]],
                    binary = variables[["binary"]], categorical = variables[["categorical"]],
                    weight_col = if (weighted) "w")$tables$stddiff
  stats::setNames(tbl$stddiff, tbl$variable)[unname(variables)]
}

check_permutation_reference <- function(file, weighted) {
  sas <- read_oracle(file)
  observed <- r_stddiff(weight = if (weighted) "w")
  perms <- vapply(seq_len(n_perm), permuted_stddiff, numeric(length(variables)), weighted = weighted)
  pct <- t(apply(perms, 1L, .perm_percentiles, type = 2))
  colnames(pct) <- c("p_2_5", "p_16", "p_50", "p_84", "p_97_5")
  for (v in unname(variables)) {
    row <- sas[tolower(sas$varname) == v, ]
    testthat::expect_equal(unname(observed[[v]]), row$stddiff, tolerance = 5e-5, label = paste(v, "observed"))
    testthat::expect_equal(unname(pct[v, ]), unname(unlist(row[colnames(pct)])), tolerance = 5e-5,
                           label = paste(v, "percentiles"))
  }
}

test_that("the w_k columns SAS reads are the weights reweight() derives", {
  for (k in seq_len(n_perm)) {
    d <- fixture
    d$grp <- d[[paste0("fgrp_", k)]]
    expect_equal(d[[paste0("w_", k)]], reweight(d), label = paste0("w_", k))
  }
})

test_that("the permutation reference matches %stddiffci, unweighted", {
  skip_without_oracle("ci-unweighted.csv")
  check_permutation_reference("ci-unweighted.csv", weighted = FALSE)
})

test_that("the permutation reference matches %stddiffci, weighted", {
  skip_without_oracle("ci-weighted.csv")
  check_permutation_reference("ci-weighted.csv", weighted = TRUE)
})

# ---- %mw_var ----------------------------------------------------------------
# SURVEYSELECT and R draw different bootstrap samples, so the replicates cannot
# be compared. The observed difference and group summaries can, and so can
# ps_mw_var()'s SD and percentile definitions, applied to SAS's own replicates.

test_that("ps_mw_var() matches %mw_var on the observed difference and group summaries", {
  skip_without_oracle("mw-summary.csv")
  sas <- read_oracle("mw-summary.csv")
  r <- ps_mw_var(fixture, treatment_col = "grp", outcomes = c("y_cont", "y_bin"), weight_col = "w",
                 n_rep = 10, seed = 1)$tables$mw_var
  for (v in c("y_cont", "y_bin")) {
    s <- sas[tolower(sas$`_outcome`) == v, ]
    rr <- r[r$outcome == v, ]
    expect_equal(rr$estimate, s$`_mw_trt`, tolerance = 1e-6, label = paste(v, "estimate"))
    expect_equal(c(rr$mean_0, rr$mean_1), c(s$`_mean_0`, s$`_mean_1`), tolerance = 1e-6,
                 label = paste(v, "group means"))
    expect_equal(c(rr$sd_0, rr$sd_1), c(s$`_std_0`, s$`_std_1`), tolerance = 1e-6, label = paste(v, "group SDs"))
    expect_equal(c(rr$sumwgt_0, rr$sumwgt_1), c(s$`_n_0`, s$`_n_1`), tolerance = 1e-6,
                 label = paste(v, "sums of weights"))
  }
})

test_that("ps_mw_var()'s SD and PCTLDEF=1 percentiles reproduce %mw_var on its own replicates", {
  skip_without_oracle(c("mw-summary.csv", "mw-replicates.csv"))
  sas <- read_oracle("mw-summary.csv")
  reps <- read_oracle("mw-replicates.csv")
  reps <- reps[rowSums(!is.na(reps[c("y_cont", "y_bin")])) > 0L, ] # the macro's first row is empty
  for (v in c("y_cont", "y_bin")) {
    s <- sas[tolower(sas$`_outcome`) == v, ]
    expect_equal(nrow(reps), s$resample, label = paste(v, "replicates"))
    expect_equal(stats::sd(reps[[v]]), s$`_sd_trt`, tolerance = 1e-6, label = paste(v, "bootstrap SD"))
    expect_equal(.perm_percentiles(reps[[v]], type = 4), unname(unlist(s[c("p2_5", "p16", "p50", "p84", "p97_5")])),
                 tolerance = 1e-6, label = paste(v, "percentiles"))
  }
})

test_that("ps_mw_var()'s bootstrap SD is within Monte Carlo error of %mw_var's (spec section 9, item 3)", {
  skip_without_oracle("mw-summary.csv")
  sas <- read_oracle("mw-summary.csv")
  n_rep <- 2000L
  r <- ps_mw_var(fixture, treatment_col = "grp", outcomes = c("y_cont", "y_bin"), weight_col = "w",
                 n_rep = n_rep, seed = 20260930)$tables$mw_var
  for (v in c("y_cont", "y_bin")) {
    s <- sas[tolower(sas$`_outcome`) == v, ]
    # A bootstrap SD from B replicates has relative standard error about
    # 1 / sqrt(2 (B - 1)): about 5.0% for SAS's 200 and 1.6% for R's 2000, so
    # about 5.3% for their ratio. Four of those, 21%, keeps a false failure
    # below 1 in 10,000 while still catching an SD of the wrong kind.
    se_ratio <- sqrt(1 / (2 * (s$resample - 1)) + 1 / (2 * (n_rep - 1)))
    expect_lt(abs(log(r$sd[r$outcome == v] / s$`_sd_trt`)), 4 * se_ratio, label = paste(v, "bootstrap SD"))
  }
})
