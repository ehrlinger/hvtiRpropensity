## make-input.R -- writes input.csv for the %stddiff / %stddiffci / %mw_var oracle.
##
## Synthetic and deterministic: no study data. Run from this directory with
##   Rscript make-input.R
## and commit input.csv. The SAS program reads only input.csv, and the R tests
## read the same file, so both sides see identical rows, including the
## permuted group columns %stddiffci expects the caller to supply.

matching_weight <- function(p, g) round(pmin(p, 1 - p) / (p * g + (1 - p) * (1 - g)), 6)

withr::with_seed(20260930, {
  n0 <- 70L
  n1 <- 50L # unequal on purpose: a size-weighted pooled SD would differ
  n <- n0 + n1
  grp <- rep(c(0L, 1L), c(n0, n1))

  x_gauss <- round(stats::rnorm(n, mean = 60 + 4 * grp, sd = 9 + 3 * grp), 2)
  x_ord <- pmin(5L, pmax(1L, as.integer(round(stats::rnorm(n, mean = 2.6 + 0.5 * grp, sd = 1.1)))))
  x_bin <- stats::rbinom(n, 1L, 0.35 + 0.15 * grp)
  x_cat <- vapply(seq_len(n), function(i) {
    p <- if (grp[i] == 1L) c(0.20, 0.30, 0.35, 0.15) else c(0.35, 0.30, 0.20, 0.15)
    sample(c("A", "B", "C", "D"), 1L, prob = p)
  }, character(1))

  # A few missing values, so the oracle also shows how each side drops them.
  x_gauss[c(5L, 77L, 110L)] <- NA
  x_ord[c(12L, 90L)] <- NA
  x_cat[c(33L, 101L)] <- NA

  # Matching weight from a score depending on the covariates. The score is
  # kept as a column so the weights can be re-derived for a permuted group.
  lp <- -0.35 + 0.04 * (ifelse(is.na(x_gauss), 60, x_gauss) - 60) + 0.5 * x_bin
  p_score <- stats::plogis(lp)
  w <- matching_weight(p_score, grp)

  y_cont <- round(stats::rnorm(n, mean = 10 + 1.5 * grp + 0.5 * ifelse(is.na(x_ord), 3, x_ord), sd = 3), 3)
  y_bin <- stats::rbinom(n, 1L, 0.2 + 0.1 * grp)

  d <- data.frame(id = seq_len(n), grp = grp, x_gauss = x_gauss, x_ord = x_ord,
                  x_bin = x_bin, x_cat = x_cat, p_score = p_score, w = w, y_cont = y_cont, y_bin = y_bin)

  # %stddiffci reads permuted groups fgrp_1..fgrp_K and, when weighted,
  # weights w_1..w_K from the input. Each w_k is the matching weight re-derived
  # for fgrp_k from the same score, as the spec requires; ps_stddiff_perm()
  # does the same through reweight = function(d) matching_weight(d$p_score, d$grp).
  n_perm <- 20L
  for (k in seq_len(n_perm)) d[[paste0("fgrp_", k)]] <- sample(grp)
  for (k in seq_len(n_perm)) d[[paste0("w_", k)]] <- matching_weight(p_score, d[[paste0("fgrp_", k)]])
})

utils::write.csv(d, "input.csv", row.names = FALSE, na = "")
