## make-input.R -- writes input.csv for the %stddiff / %stddiffci / %mw_var oracle.
##
## Synthetic and deterministic: no study data. Run from this directory with
##   Rscript make-input.R
## and commit input.csv. The SAS program reads only input.csv, and the R tests
## read the same file, so both sides see identical rows, including the
## permuted group columns %stddiffci expects the caller to supply.

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

  # Matching weight from a score depending on the covariates.
  lp <- -0.35 + 0.04 * (ifelse(is.na(x_gauss), 60, x_gauss) - 60) + 0.5 * x_bin
  p <- stats::plogis(lp)
  w <- round(pmin(p, 1 - p) / (p * grp + (1 - p) * (1 - grp)), 6)

  y_cont <- round(stats::rnorm(n, mean = 10 + 1.5 * grp + 0.5 * ifelse(is.na(x_ord), 3, x_ord), sd = 3), 3)
  y_bin <- stats::rbinom(n, 1L, 0.2 + 0.1 * grp)

  d <- data.frame(id = seq_len(n), grp = grp, x_gauss = x_gauss, x_ord = x_ord,
                  x_bin = x_bin, x_cat = x_cat, w = w, y_cont = y_cont, y_bin = y_bin)

  # %stddiffci reads permuted groups fgrp_1..fgrp_K and, when weighted,
  # weights w_1..w_K from the input. The weights are held fixed here, which
  # ps_stddiff_perm() reproduces with reweight = function(d) d$w.
  n_perm <- 20L
  for (k in seq_len(n_perm)) d[[paste0("fgrp_", k)]] <- sample(grp)
  for (k in seq_len(n_perm)) d[[paste0("w_", k)]] <- w
})

utils::write.csv(d, "input.csv", row.names = FALSE, na = "")
