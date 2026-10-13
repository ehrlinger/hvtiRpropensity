###############################################################################
## sa-sensitivity.R
##
## Sensitivity analysis functions for propensity score analyses.
##
## Four complementary methods:
##
##   sa_rosenbaum()  -- Rosenbaum gamma bounds for matched analyses; no extra
##                     dependency
##   sa_evalue()     -- E-values (VanderWeele & Ding 2017); no extra dependency
##   sa_overlap()    -- Overlap / positivity diagnostics on the PS distribution
##   sa_trim_sweep() -- IPTW trim-threshold sensitivity sweep
##
## Each returns a ps_data subclass of the same name, built by new_ps_data().
## Every one reads its column names and settings from a ps_data input through
## .resolve_meta_arg(): a supplied argument wins, a conflict with the
## object's metadata aborts, and an unsupplied argument takes the metadata.
##
###############################################################################


# ---------------------------------------------------------------------------
# Argument precedence shared by the four sa_*() functions
# ---------------------------------------------------------------------------

#' Resolve an argument against the metadata of a ps_data object
#'
#' One rule for every `sa_*()` function. An argument the caller supplied wins
#' over the object's `$meta`, and if the two disagree the call aborts and
#' names both values. An argument the caller did not supply takes the
#' metadata value, then the function's default.
#'
#' @param name      The argument name, for the error message.
#' @param value     The argument value (the default when not supplied).
#' @param supplied  `TRUE` when the caller supplied the argument, from
#'   `!missing()`. A supplied `NULL` counts as not supplied.
#' @param meta      The value recorded in the object's metadata, or `NULL`.
#' @param meta_name The name of that metadata field. Default `name`.
#' @param call_env  Environment for the error call.
#' @return The resolved value.
#' @keywords internal
.resolve_meta_arg <- function(name, value, supplied, meta, meta_name = name, call_env = rlang::caller_env()) {
  if (!supplied || is.null(value)) return(rlang::`%||%`(meta, value))
  if (is.null(meta)) return(value)
  if (!identical(as.character(value), as.character(meta))) {
    fmt <- function(v) {
      if (is.character(v)) paste0("\"", v, "\"", collapse = ", ") else paste(format(v), collapse = ", ")
    }
    rlang::abort(
      sprintf(
        paste0("`%s = %s` conflicts with the object's metadata, `$meta$%s = %s`. ",
               "Omit `%s` to use the object's value, or pass the data frame `x$data` instead."),
        name, fmt(value), meta_name, fmt(meta), name
      ),
      call = call_env
    )
  }
  value
}


# ---------------------------------------------------------------------------
# sa_rosenbaum -- Rosenbaum sensitivity bounds
# ---------------------------------------------------------------------------

#' Rosenbaum Sensitivity Bounds for Matched Analyses
#'
#' Assesses how sensitive a matched-pair analysis is to unmeasured confounding
#' using the Rosenbaum bounds framework.  For each value of the sensitivity
#' parameter Gamma -- the maximum odds ratio for differential treatment
#' assignment within a matched pair due to hidden bias -- the function reports
#' the upper and lower bounds on the p-value of the Wilcoxon signed-rank test.
#'
#' A Gamma of 1 corresponds to no hidden bias (equivalent to perfect
#' randomisation within pairs).  The *sensitivity value* is the largest Gamma
#' at which the upper p-value bound remains below `alpha`; a larger sensitivity
#' value indicates a more robust result.
#'
#' The statistic is the Wilcoxon signed-rank sum \eqn{T^+} of the ranks of
#' \eqn{|d_i|} over pairs with \eqn{d_i > 0}, where \eqn{d_i} is the treated
#' minus the control outcome and zero differences are dropped. Tied
#' \eqn{|d_i|} take their average rank. Under hidden bias Gamma each pair is
#' positive with probability at most \eqn{p^+ = \Gamma / (1 + \Gamma)}, so the
#' worst-case null distribution has mean \eqn{p^+ \sum r_i} and variance
#' \eqn{p^+ (1 - p^+) \sum r_i^2}, and the best case uses
#' \eqn{p^- = 1 / (1 + \Gamma)}. Using \eqn{\sum r_i^2} rather than
#' \eqn{n(n+1)(2n+1)/6} corrects the variance for ties; without ties the two
#' agree. The p-values use the normal approximation with a continuity
#' correction of 0.5, so at Gamma = 1 they equal the one-sided
#' [stats::wilcox.test()] p-value with `exact = FALSE` and `correct = TRUE`.
#'
#' Rosenbaum's method assumes the direction of the test is fixed before the
#' outcome is seen. Set `alternative` to `"greater"` (treated above control)
#' or `"less"` to state it. The default `"auto"` keeps the earlier behaviour
#' and picks the direction from the data: it tests the lower tail, with a
#' warning, when the median nonzero difference is negative. A direction
#' chosen from the data makes the bounds slightly optimistic, so prefer an
#' explicit `alternative` in a pre-specified analysis.
#'
#' Column names and the treated level are read from a [ps_match()] object's
#' metadata. A supplied `treatment_col`, `pair_id_col` or `treated_level`
#' must agree with that metadata; a conflict aborts.
#'
#' @param x             A `ps_match` object (returned by [ps_match()]) **or** a
#'   plain data frame whose rows are matched patients.
#' @param outcome_col   Name of the numeric outcome column.  Required.
#' @param gamma_max     Maximum Gamma to test.  Default `3`.
#' @param gamma_inc     Increment between Gamma values.  Default `0.25`.
#' @param alpha         Significance threshold.  Default `0.05`.
#' @param treatment_col Name of the binary treatment column (0/1).  Taken from
#'   `x$meta$treatment_col` when `x` is a `ps_match` object.  Required when
#'   `x` is a plain data frame.
#' @param pair_id_col   Name of the pair-ID column created by [ps_match()].
#'   Taken from `x$meta$pair_id_col` when recorded.  Default `"pair_id"`.
#' @param treated_level The value of `treatment_col` that is treated. Read
#'   from the metadata of a `ps_data` object when it records one. `NULL`
#'   (default) requires a 0/1 or logical column and takes 1 or `TRUE`.
#' @param alternative   Direction of the test: `"auto"` (default; chosen from
#'   the sign of the median difference, see Details), `"greater"` (treated
#'   outcome above control) or `"less"`.
#'
#' @return An object of class `c("sa_rosenbaum", "ps_data")`:
#' \describe{
#'   \item{`$data`}{The input data frame (`x$data` for a `ps_match` object).}
#'   \item{`$meta`}{Column names, `treated_level`, the settings (`gamma`, the
#'     grid tested, `gamma_max`, `gamma_inc`, `alpha`, `alternative`, the
#'     direction tested, and `direction_from_data`), and the results:
#'     `sensitivity_value`, the largest Gamma for which the upper p-value is
#'     below `alpha` (`NA` if the analysis is not significant at Gamma = 1),
#'     `n_pairs`, the complete matched pairs used, `n_dropped`, `n_nonzero`
#'     and the statistic `t_plus`.}
#'   \item{`$tables$bounds`}{Data frame with columns `gamma`, `p_upper`,
#'     `p_lower`, `reject_upper` (logical -- is the upper bound < `alpha`?).}
#'   \item{`$tables$pairs`}{One row per complete pair: `pair_id`, `treated`
#'     and `control` outcomes, and `difference`.}
#'   \item{`$models`}{An empty list.}
#' }
#'
#' @seealso [ps_match()], [sa_evalue()], [sa_overlap()], [sa_trim_sweep()]
#'
#' @references
#' Rosenbaum, P. R. (2002). *Observational Studies* (2nd ed.). Springer.
#'
#' @examples
#' dta <- sample_ps_data(n = 300, seed = 1)
#' m   <- ps_match(dta)
#' res <- sa_rosenbaum(m, outcome_col = "ef", gamma_max = 2, alternative = "greater")
#' res
#' res$tables$bounds
#' res$meta$sensitivity_value
#'
#' @export
sa_rosenbaum <- function(x,
                         outcome_col,
                         gamma_max     = 3,
                         gamma_inc     = 0.25,
                         alpha         = 0.05,
                         treatment_col = NULL,
                         pair_id_col   = "pair_id",
                         treated_level = NULL,
                         alternative   = c("auto", "greater", "less")) {

  alternative <- match.arg(alternative)

  # ---- Extract data frame --------------------------------------------------
  if (is_ps_data(x)) {
    treatment_col <- .resolve_meta_arg("treatment_col", treatment_col, !missing(treatment_col),
                                       x$meta$treatment_col)
    pair_id_col   <- .resolve_meta_arg("pair_id_col", pair_id_col, !missing(pair_id_col),
                                       x$meta$pair_id_col)
    treated_level <- .resolve_meta_arg("treated_level", treated_level, !missing(treated_level),
                                       x$meta$treated_level)
    dta <- x$data
  } else if (is.data.frame(x)) {
    if (is.null(treatment_col)) {
      rlang::abort(
        "`treatment_col` must be specified when `x` is a plain data frame.",
        call. = FALSE
      )
    }
    dta <- x
  } else {
    rlang::abort("`x` must be a ps_match object or a data frame.", call. = FALSE)
  }

  # ---- Validate parameters -------------------------------------------------
  if (!is.numeric(gamma_max) || length(gamma_max) != 1L || gamma_max < 1) {
    rlang::abort("`gamma_max` must be a single number >= 1.", call. = FALSE)
  }
  if (!is.numeric(gamma_inc) || length(gamma_inc) != 1L || gamma_inc <= 0) {
    rlang::abort("`gamma_inc` must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    rlang::abort("`alpha` must be a single number in (0, 1).", call. = FALSE)
  }

  # ---- Validate columns ----------------------------------------------------
  .check_cols(dta, c(outcome_col, treatment_col, pair_id_col))

  if (!is.numeric(dta[[outcome_col]])) {
    rlang::abort(
      sprintf(
        "`outcome_col` ('%s') must be numeric. Factors and characters produce wrong ranks.",
        outcome_col
      ),
      call. = FALSE
    )
  }

  # ---- Extract matched pairs -----------------------------------------------
  matched <- dta[!is.na(dta[[pair_id_col]]), , drop = FALSE]

  if (nrow(matched) == 0L) {
    rlang::abort(
      "No matched pairs found -- all `pair_id` values are NA.",
      call. = FALSE
    )
  }

  treatment   <- .treatment_indicator(matched, treatment_col, treated_level)
  matched_trt <- treatment$trt
  trt_rows <- matched[matched_trt == 1L, , drop = FALSE]
  ctl_rows <- matched[matched_trt == 0L, , drop = FALSE]

  pair_ids <- sort(unique(trt_rows[[pair_id_col]]))
  trt_out  <- trt_rows[[outcome_col]][match(pair_ids, trt_rows[[pair_id_col]])]
  ctl_out  <- ctl_rows[[outcome_col]][match(pair_ids, ctl_rows[[pair_id_col]])]

  # Drop pairs with missing outcome values (warn so the caller knows)
  complete  <- !is.na(trt_out) & !is.na(ctl_out)
  n_dropped <- sum(!complete)
  if (n_dropped > 0L) {
    rlang::warn(
      sprintf(
        "%d matched pair(s) dropped due to missing outcome values in '%s'.",
        n_dropped, outcome_col
      )
    )
  }
  pair_ids <- pair_ids[complete]
  trt_out <- trt_out[complete]
  ctl_out <- ctl_out[complete]
  n_pairs <- sum(complete)

  if (n_pairs < 2L) {
    rlang::abort(
      "At least 2 complete matched pairs are required.",
      call. = FALSE
    )
  }

  # ---- Rosenbaum bounds for Wilcoxon signed-rank statistic -----------------
  # Direct implementation of Rosenbaum (2002), Ch. 4.
  #
  # d_i = treated_i - control_i for each matched pair.
  # T+  = sum of ranks of |d_i| among nonzero pairs where d_i > 0.
  #
  # Under hidden bias Gamma (odds of differential assignment) each pair is
  # positive with probability at most p = Gamma/(1+Gamma), independently, so
  # with r_i the (average) ranks of |d_i|:
  #   max E[T+] = Gamma/(1+Gamma)   * sum(r)
  #   min E[T+] = 1/(1+Gamma)       * sum(r)
  #   Var[T+]   = Gamma/(1+Gamma)^2 * sum(r^2)
  # Without ties sum(r) = n(n+1)/2 and sum(r^2) = n(n+1)(2n+1)/6; with ties
  # sum(r^2) is the tie-corrected value (sum(r) is unchanged by midranks).
  #
  # p_upper (worst-case p-value): bias inflates E[T+], making result look
  #   less extreme.  p+ = P(T+ >= t_obs | max bias).
  # p_lower (best-case p-value):  bias deflates E[T+], making result look
  #   more extreme. p- = P(T+ >= t_obs | min bias).
  #
  # A continuity correction of 0.5 is applied (standard for discrete statistics).

  d    <- trt_out - ctl_out
  d_nz <- d[d != 0]
  n_nz <- length(d_nz)

  if (n_nz == 0L) {
    rlang::abort(
      "All matched-pair differences are zero; the signed-rank statistic is undefined.",
      call. = FALSE
    )
  }

  # ---- Orient the test ----------------------------------------------------
  # The Wilcoxon signed-rank statistic T+ counts positive-difference ranks.
  # "less" tests the negative differences, so that p_upper still represents
  # the worst-case (least significant) p-value. "auto" chooses the direction
  # from the sign of the median difference, which Rosenbaum's method assumes
  # is fixed in advance; it is kept as the default for earlier results.
  direction <- alternative
  if (alternative == "auto") {
    direction <- if (stats::median(d_nz) < 0) "less" else "greater"
    if (direction == "less") {
      rlang::warn(
        paste0(
          "Median matched-pair difference (treated - control) is negative. ",
          "The test has been re-oriented to the lower tail. ",
          "p_upper and p_lower reflect the worst- and best-case p-values ",
          "for a *harmful* treatment effect. Set `alternative` to fix the ",
          "direction in advance."
        )
      )
    }
  }
  if (direction == "less") d_nz <- -d_nz

  r      <- rank(abs(d_nz))
  # Wilcoxon signed-rank statistic T+ (sum of ranks of positive differences)
  t_plus <- sum(r[d_nz > 0])
  sum_r  <- sum(r)
  sum_r2 <- sum(r^2)

  # Round to avoid floating-point accumulation: seq(1, 3, 0.25) can produce
  # values like 2.9999999999 instead of 3.0.
  gammas <- round(seq(1, gamma_max, by = gamma_inc), digits = 10)

  bounds_list <- lapply(gammas, function(g) {
    mu_hi <- g / (1 + g) * sum_r
    mu_lo <- 1 / (1 + g) * sum_r
    sigma <- sqrt(g / (1 + g)^2 * sum_r2)

    p_up <- 1 - stats::pnorm((t_plus - 0.5 - mu_hi) / sigma)
    p_lo <- 1 - stats::pnorm((t_plus - 0.5 - mu_lo) / sigma)

    data.frame(
      gamma        = g,
      p_upper      = max(0, min(1, p_up)),
      p_lower      = max(0, min(1, p_lo)),
      reject_upper = p_up < alpha,
      stringsAsFactors = FALSE
    )
  })

  bounds <- do.call(rbind, bounds_list)
  rownames(bounds) <- NULL

  sig_gammas        <- bounds$gamma[bounds$reject_upper]
  sensitivity_value <- if (length(sig_gammas) > 0L) max(sig_gammas) else NA_real_

  pairs <- data.frame(pair_id = pair_ids, treated = trt_out, control = ctl_out,
                      difference = trt_out - ctl_out)

  new_ps_data(
    data = dta,
    meta = list(
      outcome_col         = outcome_col,
      treatment_col       = treatment_col,
      treated_level       = treatment$labels[[2L]],
      pair_id_col         = pair_id_col,
      gamma               = gammas,
      gamma_max           = gamma_max,
      gamma_inc           = gamma_inc,
      alpha               = alpha,
      alternative         = direction,
      direction_from_data = alternative == "auto",
      method              = "Rosenbaum bounds, Wilcoxon signed-rank",
      n_total             = nrow(dta),
      n_pairs             = n_pairs,
      n_dropped           = n_dropped,
      n_nonzero           = n_nz,
      t_plus              = t_plus,
      sensitivity_value   = sensitivity_value
    ),
    tables = list(bounds = bounds, pairs = pairs),
    subclass = "sa_rosenbaum"
  )
}

#' @export
print.sa_rosenbaum <- function(x, ...) {
  m <- x$meta
  cat("<sa_rosenbaum>\n")
  cat(sprintf("  Outcome     : %s (%s, %d pairs)\n", m$outcome_col,
              if (m$direction_from_data) paste(m$alternative, "chosen from data") else m$alternative, m$n_pairs))
  cat(sprintf("  Sensitivity : Gamma = %s at alpha = %g\n", format(m$sensitivity_value), m$alpha))
  print(x$tables$bounds, row.names = FALSE)
  invisible(x)
}


# ---------------------------------------------------------------------------
# sa_evalue -- E-value calculation
# ---------------------------------------------------------------------------

#' E-Value Sensitivity Analysis
#'
#' Computes the E-value for a point estimate and optionally for a confidence
#' interval bound, following VanderWeele & Ding (2017).  The E-value is the
#' minimum strength of association (on the risk-ratio scale) that an unmeasured
#' confounder would need with **both** the treatment and the outcome,
#' conditional on measured covariates, to fully explain away the observed
#' association.
#'
#' Odds ratios and hazard ratios are converted to the risk-ratio scale before
#' the E-value formula is applied, and the conversion depends on whether the
#' outcome is rare (VanderWeele & Ding 2017; VanderWeele 2017):
#'
#' * **Rare outcome** (prevalence at most about 15%): the OR or HR
#'   approximates the RR and is used directly.
#' * **Common outcome**: `RR = sqrt(OR)` for an odds ratio, and
#'   `RR = (1 - 0.5^sqrt(HR)) / (1 - 0.5^sqrt(1 / HR))` for a hazard ratio.
#'
#' There is no default for this choice: for `type = "OR"` or `"HR"` the
#' caller must set `rare`.  Rarity refers to the outcome's prevalence in the
#' study population (under about 15% at the end of follow-up), and it cannot
#' be inferred from the baseline risk `p0` alone: with OR = 4 and p0 = 0.10 the
#' exposed risk is about 0.31, and reversing the exposure coding would flip
#' any rule based on `p0`.  The same conversion is applied to the confidence
#' interval bounds.
#'
#' For risk differences a baseline outcome probability `p0` must be supplied,
#' and the conversion `RR = (RD + p0) / p0` treats `p0` as fixed: its own
#' sampling uncertainty is not carried into the CI E-value.
#'
#' @param estimate  Numeric.  The point estimate (RR, OR, HR, or RD on the
#'   original scale).  Must be > 0 for RR/OR/HR; any finite value for RD.
#' @param ci_lo     Optional numeric.  Lower confidence interval bound.
#' @param ci_hi     Optional numeric.  Upper confidence interval bound.
#'   Both `ci_lo` and `ci_hi` must be supplied together to compute the CI
#'   E-value, and must satisfy `ci_lo <= estimate <= ci_hi`.
#' @param type      One of `"RR"` (default), `"OR"`, `"HR"`, `"IRR"`, or
#'   `"RD"`.  `"IRR"` (incidence rate ratio) uses the same formula as `"RR"`;
#'   it is provided as a distinct value for clarity in count-outcome studies.
#' @param p0        Baseline outcome probability in the unexposed group.
#'   Required when `type = "RD"`, and used only then.  Must be in (0, 1).
#' @param rare      `TRUE` or `FALSE`; required when `type = "OR"` or `"HR"`.
#'   Whether the outcome is rare in the study population (prevalence under
#'   about 15% at the end of follow-up).  `TRUE` uses the OR or HR as the RR;
#'   `FALSE` applies the common-outcome conversion.  The default `NULL` aborts
#'   for OR and HR.  Ignored for other types.
#'
#' @return An object of class `c("sa_evalue", "ps_data")`. There is no input
#' data, so `$data` is an empty data frame and `$models` an empty list.
#' \describe{
#'   \item{`$meta$evalue_estimate`}{E-value for the point estimate.}
#'   \item{`$meta$evalue_ci`}{E-value for the confidence interval bound
#'     closest to the null.  `NA` if `ci_lo`/`ci_hi` are not supplied.  `1` if
#'     the CI already crosses the null.}
#'   \item{`$meta$rr_estimate`}{The point estimate converted to the
#'     risk-ratio scale (same as `estimate` when `type = "RR"`).}
#'   \item{`$meta$type`}{The measure type supplied.}
#'   \item{`$meta$rare`}{For `type = "OR"` or `"HR"`, the rare-outcome choice
#'     used for the conversion; `NA` for other types.}
#'   \item{`$meta$estimate`, `$meta$ci_lo`, `$meta$ci_hi`, `$meta$p0`}{The
#'     inputs, `NA` where not supplied.}
#'   \item{`$tables$evalue`}{Data frame with one row for the estimate and one
#'     for the CI bound closest to the null: `quantity`, `value` on the
#'     original scale, `rr` on the risk-ratio scale (`NA` when the CI crosses
#'     the null) and `evalue`.}
#' }
#'
#' @references
#' VanderWeele, T. J., & Ding, P. (2017).  Sensitivity analysis in
#' observational research: Introducing the E-value.  *Annals of Internal
#' Medicine*, 167(4), 268-274.  \doi{10.7326/M16-2607}
#'
#' VanderWeele, T. J. (2017).  On a square-root transformation of the odds
#' ratio for a common outcome.  *Epidemiology*, 28(6), e58-e60.
#' \doi{10.1097/EDE.0000000000000733}
#'
#' @seealso [sa_rosenbaum()], [sa_overlap()], [sa_trim_sweep()]
#'
#' @examples
#' # Risk ratio of 2.5 -- how robust is the result?
#' sa_evalue(2.5, ci_lo = 1.8, ci_hi = 3.5)
#'
#' # Odds ratio of 0.6 (protective effect) for a common outcome
#' sa_evalue(0.6, ci_lo = 0.4, ci_hi = 0.9, type = "OR", rare = FALSE)
#'
#' # Hazard ratio of 1.8 for an outcome under 15% at end of follow-up
#' sa_evalue(1.8, ci_lo = 1.2, ci_hi = 2.7, type = "HR", rare = TRUE)
#'
#' # Risk difference of 0.10 with baseline probability 0.20
#' ev <- sa_evalue(0.10, ci_lo = 0.02, ci_hi = 0.18, type = "RD", p0 = 0.20)
#' ev$meta$evalue_estimate
#' ev$tables$evalue
#'
#' @export
sa_evalue <- function(estimate,
                      ci_lo    = NULL,
                      ci_hi    = NULL,
                      type     = c("RR", "OR", "HR", "IRR", "RD"),
                      p0       = NULL,
                      rare     = NULL) {

  type <- match.arg(type)

  if (!is.numeric(estimate) || length(estimate) != 1L || is.na(estimate)) {
    rlang::abort("`estimate` must be a single, non-missing numeric value.",
                 call. = FALSE)
  }

  if (type %in% c("RR", "OR", "HR", "IRR") && estimate <= 0) {
    rlang::abort(
      sprintf("`estimate` must be positive for type = '%s'.", type),
      call. = FALSE
    )
  }

  if (type == "RD" && is.null(p0)) {
    rlang::abort(
      "`p0` (baseline outcome probability) must be supplied when type = 'RD'.",
      call. = FALSE
    )
  }

  if (!is.null(p0)) {
    if (!is.numeric(p0) || length(p0) != 1L || p0 <= 0 || p0 >= 1) {
      rlang::abort("`p0` must be a single number in (0, 1).", call. = FALSE)
    }
  }

  if (xor(is.null(ci_lo), is.null(ci_hi))) {
    rlang::abort(
      "Both `ci_lo` and `ci_hi` must be supplied together, or both omitted.",
      call. = FALSE
    )
  }

  # Validate CI bounds are positive for ratio measures
  for (bound_name in c("ci_lo", "ci_hi")) {
    bound_val <- get(bound_name)
    if (!is.null(bound_val) && type %in% c("RR", "OR", "HR", "IRR") &&
        bound_val <= 0) {
      rlang::abort(
        sprintf(
          "`%s` must be positive for type = '%s' (got %g).",
          bound_name, type, bound_val
        ),
        call. = FALSE
      )
    }
  }

  if (!is.null(ci_lo) && (ci_lo > estimate || estimate > ci_hi)) {
    rlang::abort(
      sprintf(
        "The CI must satisfy ci_lo <= estimate <= ci_hi (got ci_lo = %g, estimate = %g, ci_hi = %g).",
        ci_lo, estimate, ci_hi
      ),
      call. = FALSE
    )
  }

  # ---- Rare or common outcome for OR/HR ------------------------------------
  # VanderWeele & Ding (2017): a rare outcome (prevalence <= ~15%) uses the OR
  # or HR as the RR; a common outcome needs a conversion.  No silent default.
  if (!is.null(rare) && (!is.logical(rare) || length(rare) != 1L || is.na(rare))) {
    rlang::abort("`rare` must be TRUE, FALSE or NULL.", call. = FALSE)
  }
  if (type %in% c("OR", "HR")) {
    if (is.null(rare)) {
      rlang::abort(
        sprintf(
          paste0(
            "For type = '%s', set `rare = TRUE` or `rare = FALSE`. Rarity refers to the outcome's ",
            "prevalence in the study population (under about 15%% at the end of follow-up; ",
            "VanderWeele & Ding 2017). It cannot be inferred from `p0`, which is used only for type = 'RD'."
          ),
          type
        ),
        call. = FALSE
      )
    }
  } else {
    rare <- NA
  }

  # ---- Convert estimate to RR scale ----------------------------------------
  .to_rr <- function(val) {
    switch(type,
      RR  = val,
      IRR = val,                 # identical formula to RR
      OR  = if (rare) val else sqrt(val),
      HR  = if (rare) val else (1 - 0.5^sqrt(val)) / (1 - 0.5^sqrt(1 / val)),
      RD  = {
        p1 <- val + p0           # implied outcome probability in exposed
        if (p1 <= 0) {
          rlang::abort(
            sprintf(
              "RD (%g) + p0 (%g) = %g <= 0: implied outcome probability in the exposed group is non-positive.",
              val, p0, p1
            ),
            call. = FALSE
          )
        }
        if (p1 > 1) {
          rlang::abort(
            sprintf(
              "RD (%g) + p0 (%g) = %g > 1: implied outcome probability in the exposed group exceeds 1.",
              val, p0, p1
            ),
            call. = FALSE
          )
        }
        p1 / p0
      }
    )
  }

  # ---- E-value formula for RR >= 1 ------------------------------------------
  # For a protective effect (rr < 1) flip to rr > 1 first.
  .evalue_rr <- function(rr) {
    if (is.na(rr)) return(NA_real_)
    if (abs(rr - 1) < 1e-9) return(1)    # float-safe null check
    if (rr < 1) rr <- 1 / rr
    rr + sqrt(rr * (rr - 1))
  }

  rr_est <- .to_rr(estimate)
  ev_est <- .evalue_rr(rr_est)

  # ---- E-value for CI bound ------------------------------------------------
  ev_ci <- NA_real_
  rr_ci <- NA_real_
  ci_bound <- NA_real_
  if (!is.null(ci_lo) && !is.null(ci_hi)) {
    null_val <- if (type == "RD") 0 else 1

    # Pick the CI bound closer to the null
    ci_bound <- if (estimate > null_val) ci_lo else ci_hi

    # Float-safe null-crossing check
    null_crossed <- (estimate > null_val && ci_lo <= null_val) ||
                    (estimate < null_val && ci_hi >= null_val) ||
                    abs(estimate - null_val) < 1e-9

    if (null_crossed) {
      ev_ci <- 1
    } else {
      rr_ci <- .to_rr(ci_bound)
      ev_ci <- .evalue_rr(rr_ci)
    }
  }

  has_ci <- !is.null(ci_lo)
  evalue_tbl <- data.frame(
    quantity = c("estimate", "ci"),
    value    = c(estimate, ci_bound),
    rr       = c(rr_est, rr_ci),
    evalue   = c(ev_est, ev_ci),
    stringsAsFactors = FALSE
  )

  new_ps_data(
    data = data.frame(),
    meta = list(
      estimate        = estimate,
      ci_lo           = if (has_ci) ci_lo else NA_real_,
      ci_hi           = if (has_ci) ci_hi else NA_real_,
      type            = type,
      p0              = if (is.null(p0)) NA_real_ else p0,
      rare            = rare,
      method          = "E-value (VanderWeele & Ding 2017)",
      rr_estimate     = rr_est,
      evalue_estimate = ev_est,
      evalue_ci       = ev_ci
    ),
    tables = list(evalue = evalue_tbl),
    subclass = "sa_evalue"
  )
}

#' @export
print.sa_evalue <- function(x, ...) {
  m <- x$meta
  cat("<sa_evalue>\n")
  cat(sprintf("  Estimate    : %s = %g%s\n", m$type, m$estimate,
              if (is.na(m$ci_lo)) "" else sprintf(" (CI %g to %g)", m$ci_lo, m$ci_hi)))
  if (!is.na(m$rare)) cat(sprintf("  Rare outcome: %s\n", m$rare))
  cat(sprintf("  RR scale    : %g\n", m$rr_estimate))
  cat(sprintf("  E-value     : %g (estimate), %s (CI)\n", m$evalue_estimate, format(m$evalue_ci)))
  invisible(x)
}


# ---------------------------------------------------------------------------
# sa_overlap -- overlap / positivity diagnostics
# ---------------------------------------------------------------------------

#' Propensity Score Overlap and Positivity Diagnostics
#'
#' Examines the propensity score (PS) distribution for potential positivity
#' violations -- regions of covariate space where the probability of receiving
#' one treatment is near zero -- and quantifies the common support region.
#'
#' The function reports:
#' - PS summary statistics by treatment group.
#' - The empirical overlap region: \[max(min_treated, min_control),
#'   min(max_treated, max_control)\].
#' - The proportion of each group whose PS lies outside the overlap region.
#' - The number and proportion of patients with PS below `trim_threshold` or
#'   above `1 - trim_threshold` (near-positivity-violation flags).
#'
#' @param x             A `ps_data` object **or** a plain data frame.
#' @param score_col     Name of the propensity score column.  Taken from
#'   `x$meta$score_col` when `x` is a `ps_data` object.  Default `"prob_t"`.
#'   The score of a [ps_match()] object built with `scale = "logit"` is a
#'   logit, and is back-transformed with [stats::plogis()] before use.
#' @param treatment_col Name of the binary treatment column (0/1).  Taken from
#'   `x$meta$treatment_col` when `x` is a `ps_data` object`.  Default `"tavr"`.
#' @param trim_threshold Numeric in (0, 0.5).  PS values below this threshold
#'   or above `1 - trim_threshold` are flagged as near-positivity violations.
#'   Default `0.05`.
#'
#' @param treated_level The value of `treatment_col` that is treated. Read
#'   from the metadata of a `ps_data` object when it records one. `NULL`
#'   (default) requires a 0/1 or logical column and takes 1 or `TRUE`.
#'
#' @details For a `ps_data` input, a supplied `score_col`, `treatment_col` or
#'   `treated_level` must agree with the object's metadata; a conflict aborts.
#'   Arguments not supplied are read from the metadata.
#'
#' @return An object of class `c("sa_overlap", "ps_data")`:
#' \describe{
#'   \item{`$data`}{The input data frame (`x$data` for a `ps_data` object),
#'     unchanged.}
#'   \item{`$meta`}{`score_col`, `treatment_col`, `treated_level`, `scale`,
#'     `trim_threshold`, `n_total`, `n_score_missing`, `n_used` (patients
#'     with a score and a treatment), and `overlap_region`:
#'     a named numeric vector of the `lower` and `upper` bounds of the
#'     empirical common support region.}
#'   \item{`$tables$summary`}{Data frame: one row per treatment group with
#'     `group` (the value in `treatment_col`), a logical `treated`, `n`,
#'     `ps_min`, `ps_q25`, `ps_median`, `ps_mean`, `ps_q75`, `ps_max`.}
#'   \item{`$tables$outside_overlap`}{Data frame: `group`, `treated`,
#'     `n_outside`, `pct_outside` -- patients whose PS is outside the overlap
#'     region.}
#'   \item{`$tables$positivity_flags`}{Data frame: `group`, `treated`,
#'     `n_near_zero` (PS < `trim_threshold`), `n_near_one`
#'     (PS > `1 - trim_threshold`), `pct_near_zero`, `pct_near_one`.}
#'   \item{`$models`}{An empty list.}
#' }
#'
#' @seealso [ps_match()], [ps_weight()], [sa_rosenbaum()], [sa_evalue()],
#'   [sa_trim_sweep()]
#'
#' @examples
#' dta <- sample_ps_data(n = 400, seed = 7)
#' res <- sa_overlap(dta, score_col = "prob_t", treatment_col = "tavr")
#' res$meta$overlap_region
#' res$tables$positivity_flags
#'
#' # From a ps_match object
#' m   <- ps_match(dta)
#' res <- sa_overlap(m)
#' res$tables$outside_overlap
#'
#' @export
sa_overlap <- function(x,
                       score_col     = "prob_t",
                       treatment_col = "tavr",
                       trim_threshold = 0.05,
                       treated_level = NULL) {

  # ---- Extract data frame --------------------------------------------------
  scale <- "probability"
  if (is_ps_data(x)) {
    score_col     <- .resolve_meta_arg("score_col", score_col, !missing(score_col), x$meta$score_col)
    treatment_col <- .resolve_meta_arg("treatment_col", treatment_col, !missing(treatment_col),
                                       x$meta$treatment_col)
    treated_level <- .resolve_meta_arg("treated_level", treated_level, !missing(treated_level),
                                       x$meta$treated_level)
    scale <- rlang::`%||%`(x$meta$scale, scale)
    input <- x$data
    dta   <- input
    dta[[score_col]] <- .probability_score(x, score_col)
  } else if (is.data.frame(x)) {
    input <- x
    dta   <- x
  } else {
    rlang::abort("`x` must be a ps_data object or a data frame.", call. = FALSE)
  }

  # ---- Validate inputs -----------------------------------------------------
  .check_cols(dta, c(score_col, treatment_col))
  .check_probability(dta, score_col)

  if (!is.numeric(trim_threshold) || length(trim_threshold) != 1L ||
      trim_threshold <= 0 || trim_threshold >= 0.5) {
    rlang::abort("`trim_threshold` must be a single number in (0, 0.5).",
                 call. = FALSE)
  }

  ps  <- dta[[score_col]]
  treatment <- .treatment_indicator(dta, treatment_col, treated_level)
  group_labels <- treatment$labels
  trt <- treatment$trt

  # Drop rows where PS or treatment is NA before computing group sizes, so
  # that n0/n1 match the summary statistics. A missing treatment is reported
  # by .treatment_indicator(); a missing score is reported here.
  n_ps_na <- sum(is.na(ps))
  if (n_ps_na > 0L) {
    rlang::warn(
      sprintf(
        "%d patient(s) with missing '%s' values excluded from overlap diagnostics.",
        n_ps_na, score_col
      )
    )
  }
  keep <- !is.na(ps) & !is.na(trt)
  ps   <- ps[keep]
  trt  <- trt[keep]

  ps0 <- ps[trt == 0L]
  ps1 <- ps[trt == 1L]
  n0  <- length(ps0)
  n1  <- length(ps1)

  # Warn if either group is too small for reliable diagnostics
  if (n0 < 2L || n1 < 2L) {
    rlang::warn(
      sprintf(
        paste0(
          "One or both treatment groups have fewer than 2 patients after ",
          "removing missing PS values (n_control=%d, n_treated=%d). ",
          "Overlap diagnostics may be unreliable."
        ),
        n0, n1
      )
    )
  }

  # ---- Summary statistics --------------------------------------------------
  .ps_summary <- function(psv, grp, treated, n) {
    data.frame(
      group     = grp,
      treated   = treated,
      n         = n,
      ps_min    = round(min(psv,    na.rm = TRUE), 4),
      ps_q25    = round(stats::quantile(psv, 0.25, na.rm = TRUE), 4),
      ps_median = round(stats::median(psv,  na.rm = TRUE), 4),
      ps_mean   = round(mean(psv,   na.rm = TRUE), 4),
      ps_q75    = round(stats::quantile(psv, 0.75, na.rm = TRUE), 4),
      ps_max    = round(max(psv,    na.rm = TRUE), 4),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }
  summary_tbl <- rbind(.ps_summary(ps0, group_labels[[1L]], FALSE, n0),
                       .ps_summary(ps1, group_labels[[2L]], TRUE, n1))

  # ---- Overlap region ------------------------------------------------------
  overlap_lo <- max(min(ps0), min(ps1))
  overlap_hi <- min(max(ps0), max(ps1))

  if (overlap_lo > overlap_hi) {
    rlang::warn(
      sprintf(
        paste0(
          "No common support: PS distributions do not overlap ",
          "(control range [%.4f, %.4f]; treated range [%.4f, %.4f]). ",
          "All patients are outside the overlap region. ",
          "IPTW and matching results are not reliable for this data."
        ),
        min(ps0), max(ps0), min(ps1), max(ps1)
      )
    )
  }

  overlap_region <- c(lower = round(overlap_lo, 4), upper = round(overlap_hi, 4))

  # ---- Outside overlap -----------------------------------------------------
  .n_outside <- function(psv) sum(psv < overlap_lo | psv > overlap_hi)
  n_out0 <- .n_outside(ps0)
  n_out1 <- .n_outside(ps1)

  outside_overlap <- data.frame(
    group      = group_labels,
    treated    = c(FALSE, TRUE),
    n_outside  = c(n_out0, n_out1),
    pct_outside = round(100 * c(n_out0 / n0, n_out1 / n1), 1),
    stringsAsFactors = FALSE
  )

  # ---- Positivity flags ----------------------------------------------------
  nz0 <- sum(ps0 < trim_threshold, na.rm = TRUE)
  no0 <- sum(ps0 > 1 - trim_threshold, na.rm = TRUE)
  nz1 <- sum(ps1 < trim_threshold, na.rm = TRUE)
  no1 <- sum(ps1 > 1 - trim_threshold, na.rm = TRUE)

  positivity_flags <- data.frame(
    group        = group_labels,
    treated      = c(FALSE, TRUE),
    n_near_zero  = c(nz0, nz1),
    n_near_one   = c(no0, no1),
    pct_near_zero = round(100 * c(nz0 / n0, nz1 / n1), 1),
    pct_near_one  = round(100 * c(no0 / n0, no1 / n1), 1),
    stringsAsFactors = FALSE
  )

  new_ps_data(
    data = input,
    meta = list(
      score_col       = score_col,
      treatment_col   = treatment_col,
      treated_level   = group_labels[[2L]],
      scale           = scale,
      trim_threshold  = trim_threshold,
      method          = "overlap and positivity diagnostics",
      n_total         = nrow(input),
      n_score_missing = n_ps_na,
      n_used          = n0 + n1,
      overlap_region  = overlap_region
    ),
    tables = list(
      summary          = summary_tbl,
      outside_overlap  = outside_overlap,
      positivity_flags = positivity_flags
    ),
    subclass = "sa_overlap"
  )
}

#' @export
print.sa_overlap <- function(x, ...) {
  m <- x$meta
  cat("<sa_overlap>\n")
  cat(sprintf("  N total     : %d\n", m$n_total))
  cat(sprintf("  PS column   : %s\n", m$score_col))
  cat(sprintf("  Overlap     : [%.4f, %.4f]\n", m$overlap_region[["lower"]], m$overlap_region[["upper"]]))
  print(x$tables$positivity_flags, row.names = FALSE)
  invisible(x)
}


# ---------------------------------------------------------------------------
# sa_trim_sweep -- IPTW trim threshold sensitivity sweep
# ---------------------------------------------------------------------------

#' IPTW Weight-Trimming Sensitivity Sweep
#'
#' Re-computes IPTW weights across a range of winsorisation thresholds and
#' reports how the weight distribution and effective sample size (ESS) change
#' as trimming increases.  This helps identify a trimming threshold that
#' adequately controls extreme weights without an unacceptable loss of effective
#' sample size.
#'
#' @param x           A `ps_weight` object (returned by [ps_weight()]) **or**
#'   a `ps_match` object **or** a plain data frame that contains a propensity
#'   score column and a binary treatment column.
#' @param trim_seq    Numeric vector of trimming values to evaluate.  Each
#'   value is the lower tail quantile; weights outside the
#'   (`trim`, `1 - trim`) range are winsorised.  `0` means no trimming.
#'   Default `seq(0, 0.10, by = 0.01)`.
#' @param estimand    One of `"ATE"`, `"ATT"`, or `"ATC"` -- the causal
#'   estimand used to compute IPTW weights.  Taken from `x$meta$estimand`
#'   when `x` is a `ps_weight` object; defaults to `"ATE"` otherwise.
#' @param stabilise   Logical.  Whether to stabilise weights.  Taken from
#'   `x$meta$stabilised` when available; defaults to `TRUE`.  As in
#'   [ps_weight()], stabilisation applies to the ATE only.
#' @param score_col   Name of the propensity score column.  Taken from
#'   `x$meta$score_col` when `x` is a `ps_data` object.  Default `"prob_t"`.
#'   The score of a [ps_match()] object built with `scale = "logit"` is a
#'   logit, and is back-transformed with [stats::plogis()] before use.
#' @param treatment_col Name of the binary treatment column.  Taken from
#'   `x$meta$treatment_col` when `x` is a `ps_data` object`.  Default
#'   `"tavr"`.
#'
#' @param treated_level The value of `treatment_col` that is treated. Read
#'   from the metadata of a `ps_data` object when it records one. `NULL`
#'   (default) requires a 0/1 or logical column and takes 1 or `TRUE`.
#' @details For a `ps_data` input, a supplied `estimand`, `stabilise`,
#'   `score_col`, `treatment_col` or `treated_level` must agree with the
#'   object's metadata; a conflict aborts. Arguments not supplied are read
#'   from the metadata, then take their defaults.
#'
#' @return An object of class `c("sa_trim_sweep", "ps_data")`:
#' \describe{
#'   \item{`$data`}{The input data frame (`x$data` for a `ps_data` object),
#'     unchanged.}
#'   \item{`$meta`}{`score_col`, `treatment_col`, `treated_level`, `scale`,
#'     `estimand`, `stabilised` (whether stabilisation was applied, so `FALSE`
#'     for the ATT and ATC), `trim_seq` (sorted), `n_total` and `n_used`, the
#'     patients with a non-missing treatment and score.}
#'   \item{`$tables$sweep`}{Data frame with one row per trim value and
#'     columns `trim` (the winsorisation threshold), `ess_control` and
#'     `ess_treated` (effective sample sizes), `max_weight` (maximum weight
#'     after trimming), `sd_weight` (standard deviation of the weights after
#'     trimming) and `pct_trimmed` (percentage of patients whose weight was
#'     winsorised).}
#'   \item{`$models`}{An empty list.}
#' }
#'
#' @seealso [ps_weight()], [sa_rosenbaum()], [sa_evalue()], [sa_overlap()]
#'
#' @examples
#' dta    <- sample_ps_data(n = 400, seed = 3)
#' w_obj  <- ps_weight(dta, estimand = "ATE")
#' sweep  <- sa_trim_sweep(w_obj)
#' head(sweep$tables$sweep)
#'
#' # Directly from a data frame
#' sweep2 <- sa_trim_sweep(dta, score_col = "prob_t", treatment_col = "tavr",
#'                         estimand = "ATE")
#'
#' @export
sa_trim_sweep <- function(x,
                          trim_seq      = seq(0, 0.10, by = 0.01),
                          estimand      = c("ATE", "ATT", "ATC"),
                          stabilise     = TRUE,
                          score_col     = "prob_t",
                          treatment_col = "tavr",
                          treated_level = NULL) {

  estimand_supplied <- !missing(estimand)
  estimand <- match.arg(estimand)

  # ---- Extract data frame and metadata -------------------------------------
  scale <- "probability"
  if (is_ps_data(x)) {
    score_col     <- .resolve_meta_arg("score_col", score_col, !missing(score_col), x$meta$score_col)
    treatment_col <- .resolve_meta_arg("treatment_col", treatment_col, !missing(treatment_col),
                                       x$meta$treatment_col)
    estimand_use  <- .resolve_meta_arg("estimand", estimand, estimand_supplied, x$meta$estimand)
    stabilise     <- .resolve_meta_arg("stabilise", stabilise, !missing(stabilise), x$meta$stabilised,
                                       meta_name = "stabilised")
    treated_level <- .resolve_meta_arg("treated_level", treated_level, !missing(treated_level),
                                       x$meta$treated_level)
    scale <- rlang::`%||%`(x$meta$scale, scale)
    input <- x$data
    dta   <- input
    dta[[score_col]] <- .probability_score(x, score_col)
  } else if (is.data.frame(x)) {
    estimand_use <- estimand
    input <- x
    dta   <- x
  } else {
    rlang::abort("`x` must be a ps_data object or a data frame.", call. = FALSE)
  }

  estimand_use <- match.arg(estimand_use, c("ATE", "ATT", "ATC"))

  # ---- Validate inputs -----------------------------------------------------
  .check_cols(dta, c(score_col, treatment_col))
  .check_probability(dta, score_col)

  if (!is.numeric(trim_seq) || any(trim_seq < 0) || any(trim_seq >= 0.5)) {
    rlang::abort(
      "`trim_seq` must be a numeric vector with all values in [0, 0.5).",
      call. = FALSE
    )
  }

  ps  <- dta[[score_col]]
  treatment <- .treatment_indicator(dta, treatment_col, treated_level)
  trt <- treatment$trt
  # Patients with a missing treatment are excluded, as the warning from
  # .treatment_indicator() says. Patients with a missing score are excluded
  # too, with a warning here, so the weights, ESS and pct_trimmed are all
  # computed on the analysed rows.
  n_ps_na <- sum(is.na(ps))
  if (n_ps_na > 0L) {
    rlang::warn(
      sprintf("%d patient(s) with missing '%s' values excluded from the trim sweep.", n_ps_na, score_col)
    )
  }
  keep <- !is.na(trt) & !is.na(ps)
  ps  <- ps[keep]
  trt <- trt[keep]
  n   <- length(ps)

  # ---- Compute raw IPTW weights once ---------------------------------------
  raw_w <- switch(estimand_use,
    ATE = ifelse(trt == 1L, 1 / ps, 1 / (1 - ps)),
    ATT = ifelse(trt == 1L, 1,      ps / (1 - ps)),
    ATC = ifelse(trt == 1L, (1 - ps) / ps, 1)
  )

  # Guard against Inf weights from ps = 0 or ps = 1.
  # .check_probability() allows exact 0/1; flag those patients explicitly.
  n_inf <- sum(is.infinite(raw_w))
  if (n_inf > 0L) {
    rlang::warn(
      sprintf(
        paste0(
          "%d patient(s) have a propensity score of exactly 0 or 1, producing ",
          "infinite IPTW weights. Those weights cannot be trimmed and will ",
          "dominate all statistics. Consider excluding these patients or ",
          "recalibrating the propensity model."
        ),
        n_inf
      )
    )
  }

  stabilised <- isTRUE(stabilise) && estimand_use == "ATE"
  if (stabilised) {
    p_trt  <- mean(trt == 1L)
    raw_w  <- raw_w * ifelse(trt == 1L, p_trt, 1 - p_trt)
  }

  # ---- Sweep ---------------------------------------------------------------
  rows <- lapply(sort(trim_seq), function(tr) {
    # Use a machine-epsilon guard instead of exact == 0 comparison, since
    # seq(0, 0.10, by = 0.01) can produce values like 1e-17 instead of 0.
    if (tr < .Machine$double.eps) {
      w <- raw_w
      n_trimmed <- 0L
    } else {
      q_lo <- stats::quantile(raw_w, tr,       na.rm = TRUE)
      q_hi <- stats::quantile(raw_w, 1 - tr,   na.rm = TRUE)
      w    <- pmin(pmax(raw_w, q_lo), q_hi)
      n_trimmed <- sum(raw_w < q_lo | raw_w > q_hi, na.rm = TRUE)
    }

    data.frame(
      trim        = tr,
      ess_control = .effective_n(w[trt == 0L]),
      ess_treated = .effective_n(w[trt == 1L]),
      max_weight  = round(max(w, na.rm = TRUE), 4),
      sd_weight   = round(stats::sd(w, na.rm = TRUE), 4),
      pct_trimmed = round(100 * n_trimmed / n, 1),
      stringsAsFactors = FALSE
    )
  })

  new_ps_data(
    data = input,
    meta = list(
      score_col     = score_col,
      treatment_col = treatment_col,
      treated_level = treatment$labels[[2L]],
      scale         = scale,
      estimand      = estimand_use,
      stabilised    = stabilised,
      trim_seq      = sort(trim_seq),
      method        = paste0("IPTW-", estimand_use, " trim sweep"),
      n_total       = nrow(input),
      n_used        = n
    ),
    tables = list(sweep = do.call(rbind, rows)),
    subclass = "sa_trim_sweep"
  )
}

#' @export
print.sa_trim_sweep <- function(x, ...) {
  m <- x$meta
  cat("<sa_trim_sweep>\n")
  cat(sprintf("  N used      : %d of %d\n", m$n_used, m$n_total))
  cat(sprintf("  Weights     : IPTW-%s, %s\n", m$estimand, if (m$stabilised) "stabilised" else "unstabilised"))
  print(x$tables$sweep, row.names = FALSE)
  invisible(x)
}
