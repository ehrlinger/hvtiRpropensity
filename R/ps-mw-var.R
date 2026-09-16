###############################################################################
## ps-mw-var.R
##
## Bootstrap variance of a matching-weight treatment effect, ported from
## %mw_var in wt_mtch_boot_sd.sas (Jeevanantham Rajeswaran 2014, after Li and
## Greene 2013).
##
## Design: hvtiRtemplates dev/specs/2026-09-16-standardized-difference-design.md, section 6
##
## The macro resamples rows with replacement within each treatment group
## (PROC SURVEYSELECT method=urs samprate=1, strata group), recomputes each
## outcome's weighted group means with the ORIGINAL weights, and summarises
## the replicate differences with PROC STDIZE: the SD and PCTLDEF=1
## percentiles. The estimate is the difference on the original data;
## z = estimate / SD and p = 2 (1 - Phi(|z|)).
###############################################################################


#' Bootstrap variance of a matching-weight treatment effect
#'
#' Estimates the weighted difference in means, group 1 minus group 0, for
#' each outcome, with a bootstrap standard deviation, percentiles, z and p. A
#' port of the CCF `%mw_var` SAS macro (Rajeswaran 2014), for matching weights
#' as described by Li and Greene (2013).
#'
#' @details
#' **Estimate.** For each group, `PROC MEANS` weighted summaries: the mean
#' `sum(w y) / sum(w)`, the SD with an `n - 1` denominator, and the sum of
#' weights. The macro labels the sum of weights `N`; it is reported here as
#' `sumwgt_0` and `sumwgt_1` because it is not a count. For a 0/1 outcome the
#' means are weighted proportions.
#'
#' **Bootstrap.** Each replicate resamples rows with replacement within each
#' group, to that group's size, and recomputes the weighted means. The
#' treatment group is the only stratum, as in the macro. `sd` is the standard
#' deviation of the replicate differences, and the percentiles use
#' `quantile(type = 4)`, which is SAS `PCTLDEF=1` as the macro sets it.
#' Replicates where an outcome has no non-missing value in a group are left
#' out for that outcome, and `n_rep` counts those that remain.
#'
#' **The weights are held fixed.** Each resampled row keeps its original
#' weight; the propensity model and the weights are not re-estimated. This
#' reproduces the macro, and it ignores the uncertainty in estimating the
#' weights, so `sd` can be too small.
#'
#' **Missing values.** Rows with a missing group are dropped. A missing
#' outcome is dropped for that outcome only.
#'
#' **Random numbers.** With `seed`, the result is reproducible and the
#' caller's random number stream is restored afterwards.
#'
#' @param data A data frame.
#' @param treatment_col Name of the group column, coded 0/1 or logical, with
#'   both groups present. Default `"tavr"`.
#' @param outcomes Character vector of numeric outcome columns.
#' @param weight_col Name of the column of matching weights, which must be
#'   positive and finite.
#' @param n_rep Number of bootstrap replicates. Default `1000`, the macro's
#'   `RESAMPL=`.
#' @param seed Optional integer seed.
#'
#' @return An object of class `c("ps_mw_var", "ps_data")` with:
#' \describe{
#'   \item{`$data`}{The input data frame, unchanged.}
#'   \item{`$meta`}{Named list: `treatment_col`, `weight_col`, `outcomes`,
#'     `method`, `n_rep`, `seed`, `n_total`, `n_dropped`.}
#'   \item{`$tables`}{Named list with `mw_var`, one row per outcome:
#'     `outcome`, `label`, `sumwgt_0`, `mean_0`, `sd_0`, `sumwgt_1`,
#'     `mean_1`, `sd_1`, `estimate`, `n_rep`, `sd`, `p2_5`, `p16`, `p50`,
#'     `p84`, `p97_5`, `z` and `p`.}
#' }
#'
#' @references Li L, Greene T. A weighting analogue to pair matching in
#'   propensity score analysis. The International Journal of Biostatistics.
#'   2013;9(2):215-234.
#'
#' @seealso [ps_weight()], [ps_stddiff()]
#'
#' @examples
#' dta <- sample_ps_data(n = 200, seed = 42)
#' dta$mt_wt <- pmin(dta$prob_t, 1 - dta$prob_t) /
#'   ifelse(dta$tavr == 1, dta$prob_t, 1 - dta$prob_t)
#' obj <- ps_mw_var(dta, treatment_col = "tavr", outcomes = c("age", "ef"),
#'                  weight_col = "mt_wt", n_rep = 200, seed = 1)
#' obj$tables$mw_var
#'
#' @export
ps_mw_var <- function(data,
                      treatment_col = "tavr",
                      outcomes,
                      weight_col,
                      n_rep         = 1000L,
                      seed          = NULL) {
  .check_df(data)
  if (!is.numeric(n_rep) || length(n_rep) != 1L || !is.finite(n_rep) ||
      n_rep < 1 || n_rep > .Machine$integer.max || n_rep != round(n_rep)) {
    rlang::abort("`n_rep` must be a positive whole number.", call = NULL)
  }
  n_rep <- as.integer(n_rep)
  outcomes <- as.character(outcomes)
  if (length(outcomes) == 0L) rlang::abort("Name at least one column in `outcomes`.", call = NULL)
  .check_cols(data, c(treatment_col, outcomes, weight_col))
  for (v in outcomes) {
    if (!is.numeric(data[[v]])) {
      rlang::abort(sprintf("Outcome column `%s` must be numeric.", v), call = NULL)
    }
  }

  keep <- !is.na(data[[treatment_col]])
  base <- data[keep, , drop = FALSE]
  .check_binary(base, treatment_col)
  grp <- as.integer(base[[treatment_col]])
  if (!all(c(0L, 1L) %in% grp)) {
    rlang::abort(sprintf("Column `%s` must contain both groups, 0 and 1.", treatment_col), call = NULL)
  }
  w <- base[[weight_col]]
  if (!is.numeric(w) || any(!is.finite(w)) || any(w <= 0)) {
    rlang::abort(sprintf("Column `%s` must hold positive, finite, non-missing weights.", weight_col), call = NULL)
  }

  # Weighted mean of y over rows r, NA when none of them has a value.
  wmean <- function(y, r) {
    r <- r[!is.na(y[r])]
    if (length(r) == 0L) NA_real_ else sum(w[r] * y[r]) / sum(w[r])
  }
  i0 <- which(grp == 0L)
  i1 <- which(grp == 1L)

  summ <- lapply(outcomes, function(v) {
    y <- base[[v]]
    one <- function(r) {
      r <- r[!is.na(y[r])]
      m <- wmean(y, r)
      s <- if (length(r) < 2L) NA_real_ else sqrt(sum(w[r] * (y[r] - m)^2) / (length(r) - 1L))
      c(sum(w[r]), m, s)
    }
    c(one(i0), one(i1))
  })

  if (!is.null(seed)) withr::local_seed(seed)
  reps <- vapply(seq_len(n_rep), function(i) {
    r0 <- i0[sample.int(length(i0), replace = TRUE)]
    r1 <- i1[sample.int(length(i1), replace = TRUE)]
    vapply(outcomes, function(v) wmean(base[[v]], r1) - wmean(base[[v]], r0), numeric(1))
  }, numeric(length(outcomes)))
  reps <- matrix(reps, nrow = length(outcomes))

  rows <- lapply(seq_along(outcomes), function(k) {
    s <- summ[[k]]
    x <- reps[k, ]
    ok <- x[!is.na(x)]
    est <- s[5L] - s[2L]
    bsd <- if (length(ok) < 2L) NA_real_ else stats::sd(ok)
    pct <- .perm_percentiles(x, type = 4)
    z <- est / bsd
    lab <- attr(data[[outcomes[k]]], "label", exact = TRUE)
    data.frame(outcome = outcomes[k], label = if (is.null(lab)) NA_character_ else as.character(lab)[1L],
               sumwgt_0 = s[1L], mean_0 = s[2L], sd_0 = s[3L],
               sumwgt_1 = s[4L], mean_1 = s[5L], sd_1 = s[6L],
               estimate = est, n_rep = length(ok), sd = bsd,
               p2_5 = pct[1L], p16 = pct[2L], p50 = pct[3L], p84 = pct[4L], p97_5 = pct[5L],
               z = z, p = 2 * (1 - stats::pnorm(abs(z))), stringsAsFactors = FALSE)
  })
  tbl <- do.call(rbind, rows)
  rownames(tbl) <- NULL

  new_ps_data(
    data   = data,
    meta   = list(treatment_col = treatment_col, weight_col = weight_col, outcomes = outcomes,
                  method = "mw_var", n_rep = n_rep, seed = seed, n_total = nrow(data), n_dropped = sum(!keep)),
    tables = list(mw_var = tbl),
    subclass = "ps_mw_var"
  )
}
