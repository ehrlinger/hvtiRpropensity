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

skip_without_oracle <- function() {
  missing <- oracle_files[!file.exists(file.path(oracle_dir, oracle_files))]
  if (length(missing) > 0L) {
    skip(paste("SAS oracle output not yet returned:", paste(missing, collapse = ", ")))
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
  skip_without_oracle()
  sas <- sas_by_variable(read_oracle("sd-unweighted.csv"), "stddiff")
  expect_equal(r_stddiff()[unname(variables)], sas, tolerance = 5e-5, ignore_attr = TRUE)
})

test_that("ps_stddiff() matches %stddiff to 4 decimals, weighted", {
  skip_without_oracle()
  sas <- sas_by_variable(read_oracle("sd-weighted.csv"), "stddiff")
  expect_equal(r_stddiff(weight = "w")[unname(variables)], sas, tolerance = 5e-5, ignore_attr = TRUE)
})

# ---- %stddiffci -------------------------------------------------------------
# %stddiffci takes its permutations from the input (fgrp_1..fgrp_20), so both
# sides score the same permuted groups and the percentiles can be compared
# exactly. ps_stddiff_perm() draws its own permutations, so its percentile
# step (.perm_percentiles(type = 2)) is checked on the shared ones.

check_permutation_reference <- function(file, weighted) {
  sas <- read_oracle(file)
  observed <- r_stddiff(weight = if (weighted) "w")
  perms <- vapply(seq_len(n_perm), function(k) {
    r_stddiff(group = paste0("fgrp_", k), weight = if (weighted) paste0("w_", k))[unname(variables)]
  }, numeric(length(variables)))
  pct <- t(apply(perms, 1L, .perm_percentiles, type = 2))
  colnames(pct) <- c("p_2_5", "p_16", "p_50", "p_84", "p_97_5")
  for (v in unname(variables)) {
    row <- sas[tolower(sas$varname) == v, ]
    expect_equal(unname(observed[[v]]), row$stddiff, tolerance = 5e-5, label = paste(v, "observed"))
    expect_equal(unname(pct[v, ]), unname(unlist(row[colnames(pct)])), tolerance = 5e-5,
                 label = paste(v, "percentiles"))
  }
}

test_that("the permutation reference matches %stddiffci, unweighted", {
  skip_without_oracle()
  check_permutation_reference("ci-unweighted.csv", weighted = FALSE)
})

test_that("the permutation reference matches %stddiffci, weighted", {
  skip_without_oracle()
  check_permutation_reference("ci-weighted.csv", weighted = TRUE)
})

# ---- %mw_var ----------------------------------------------------------------
# SURVEYSELECT and R draw different bootstrap samples, so the replicates cannot
# be compared. The observed difference and group summaries can, and so can
# ps_mw_var()'s SD and percentile definitions, applied to SAS's own replicates.

test_that("ps_mw_var() matches %mw_var on the observed difference and group summaries", {
  skip_without_oracle()
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
  skip_without_oracle()
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
