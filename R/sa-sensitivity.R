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
###############################################################################


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
#'   Default `"pair_id"`.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{`$bounds`}{Data frame with columns `gamma`, `p_upper`, `p_lower`,
#'     `reject_upper` (logical -- is the upper bound < `alpha`?).}
#'   \item{`$sensitivity_value`}{Largest Gamma for which the analysis remains
#'     significant (upper p-value < `alpha`), or `NA` if the initial analysis
#'     is not significant.}
#'   \item{`$n_pairs`}{Number of complete matched pairs used.}
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
#' res <- sa_rosenbaum(m, outcome_col = "ef", gamma_max = 2)
#' res$bounds
#' res$sensitivity_value
#'
#' @export
sa_rosenbaum <- function(x,
                         outcome_col,
                         gamma_max     = 3,
                         gamma_inc     = 0.25,
                         alpha         = 0.05,
                         treatment_col = NULL,
                         pair_id_col   = "pair_id") {

  # ---- Extract data frame --------------------------------------------------
  if (is_ps_data(x)) {
    if (is.null(treatment_col)) treatment_col <- x$meta$treatment_col
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

  trt_rows <- matched[as.integer(matched[[treatment_col]]) == 1L, , drop = FALSE]
  ctl_rows <- matched[as.integer(matched[[treatment_col]]) == 0L, , drop = FALSE]

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
  # Under hidden bias Gamma (odds of differential assignment):
  #   max E[T+] = Gamma/(1+Gamma) * n*(n+1)/2
  #   min E[T+] = 1/(1+Gamma)     * n*(n+1)/2
  #   Var[T+]   = Gamma/(1+Gamma)^2 * n*(n+1)*(2n+1)/6
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

  # ---- Orient the test to the observed direction --------------------------
  # The Wilcoxon signed-rank statistic T+ counts positive-difference ranks.
  # If the median difference is negative (treatment appears harmful), orient
  # the test on the negative differences so that p_upper still represents the
  # worst-case (least significant) p-value.
  if (stats::median(d_nz) < 0) {
    rlang::warn(
      paste0(
        "Median matched-pair difference (treated - control) is negative. ",
        "The test has been re-oriented to the lower tail. ",
        "p_upper and p_lower reflect the worst- and best-case p-values ",
        "for a *harmful* treatment effect."
      )
    )
    d_nz <- -d_nz
  }

  r      <- rank(abs(d_nz))
  # Wilcoxon signed-rank statistic T+ (sum of ranks of positive differences)
  t_plus <- sum(r[d_nz > 0])

  # Round to avoid floating-point accumulation: seq(1, 3, 0.25) can produce
  # values like 2.9999999999 instead of 3.0.
  gammas <- round(seq(1, gamma_max, by = gamma_inc), digits = 10)

  bounds_list <- lapply(gammas, function(g) {
    mu_hi <- g / (1 + g) * n_nz * (n_nz + 1) / 2
    mu_lo <- 1 / (1 + g) * n_nz * (n_nz + 1) / 2
    sigma <- sqrt(g * n_nz * (n_nz + 1) * (2 * n_nz + 1) / (6 * (1 + g)^2))

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

  list(
    bounds            = bounds,
    sensitivity_value = sensitivity_value,
    n_pairs           = n_pairs
  )
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
#' For odds ratios and hazard ratios the rare-outcome approximation
#' `RR_approx = sqrt(OR)` / `sqrt(HR)` is used to convert to the risk-ratio
#' scale before applying the E-value formula.  For risk differences a baseline
#' outcome probability `p0` must be supplied.
#'
#' @param estimate  Numeric.  The point estimate (RR, OR, HR, or RD on the
#'   original scale).  Must be > 0 for RR/OR/HR; any finite value for RD.
#' @param ci_lo     Optional numeric.  Lower confidence interval bound.
#' @param ci_hi     Optional numeric.  Upper confidence interval bound.
#'   Both `ci_lo` and `ci_hi` must be supplied together to compute the CI
#'   E-value.
#' @param type      One of `"RR"` (default), `"OR"`, `"HR"`, `"IRR"`, or
#'   `"RD"`.  `"IRR"` (incidence rate ratio) uses the same formula as `"RR"`;
#'   it is provided as a distinct value for clarity in count-outcome studies.
#' @param p0        Baseline outcome probability in the unexposed group.
#'   Required when `type = "RD"`.  Must be in (0, 1).
#'
#' @return A list with:
#' \describe{
#'   \item{`$evalue_estimate`}{E-value for the point estimate.}
#'   \item{`$evalue_ci`}{E-value for the confidence interval bound closest to
#'     the null.  `NA` if `ci_lo`/`ci_hi` are not supplied.  Returns `1` if
#'     the CI already crosses the null.}
#'   \item{`$rr_estimate`}{The point estimate converted to the risk-ratio
#'     scale (same as `estimate` when `type = "RR"`).}
#'   \item{`$type`}{The measure type supplied.}
#' }
#'
#' @references
#' VanderWeele, T. J., & Ding, P. (2017).  Sensitivity analysis in
#' observational research: Introducing the E-value.  *Annals of Internal
#' Medicine*, 167(4), 268-274.  \doi{10.7326/M16-2607}
#'
#' @seealso [sa_rosenbaum()], [sa_overlap()], [sa_trim_sweep()]
#'
#' @examples
#' # Risk ratio of 2.5 -- how robust is the result?
#' sa_evalue(2.5, ci_lo = 1.8, ci_hi = 3.5)
#'
#' # Odds ratio of 0.6 (protective effect)
#' sa_evalue(0.6, ci_lo = 0.4, ci_hi = 0.9, type = "OR")
#'
#' # Risk difference of 0.10 with baseline probability 0.20
#' sa_evalue(0.10, ci_lo = 0.02, ci_hi = 0.18, type = "RD", p0 = 0.20)
#'
#' @export
sa_evalue <- function(estimate,
                      ci_lo    = NULL,
                      ci_hi    = NULL,
                      type     = c("RR", "OR", "HR", "IRR", "RD"),
                      p0       = NULL) {

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

  # ---- Rare-outcome approximation warning for OR/HR ------------------------
  # sqrt(OR) and sqrt(HR) are only accurate when outcome prevalence is low.
  # If p0 is supplied (and > 0.15), warn that the approximation may be poor.
  if (type %in% c("OR", "HR") && !is.null(p0) && p0 > 0.15) {
    rlang::warn(
      sprintf(
        paste0(
          "The rare-outcome approximation RR ~ sqrt(%s) may be inaccurate when ",
          "baseline outcome probability (p0 = %.2f) exceeds ~0.15. ",
          "Consider converting to RR directly before calling sa_evalue()."
        ),
        type, p0
      )
    )
  }

  # ---- Convert estimate to RR scale ----------------------------------------
  .to_rr <- function(val) {
    switch(type,
      RR  = val,
      IRR = val,                 # identical formula to RR
      OR  = sqrt(val),           # rare-outcome approximation
      HR  = sqrt(val),           # rare-outcome approximation
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
      ev_ci <- .evalue_rr(.to_rr(ci_bound))
    }
  }

  list(
    evalue_estimate = ev_est,
    evalue_ci       = ev_ci,
    rr_estimate     = rr_est,
    type            = type
  )
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
#' @param treatment_col Name of the binary treatment column (0/1).  Taken from
#'   `x$meta$treatment_col` when `x` is a `ps_data` object`.  Default `"tavr"`.
#' @param trim_threshold Numeric in (0, 0.5).  PS values below this threshold
#'   or above `1 - trim_threshold` are flagged as near-positivity violations.
#'   Default `0.05`.
#'
#' @return A list with:
#' \describe{
#'   \item{`$summary`}{Data frame: one row per treatment group with `group`
#'     (the value in `treatment_col`), a logical `treated`, `n`, `ps_min`,
#'     `ps_q25`, `ps_median`, `ps_mean`, `ps_q75`, `ps_max`.}
#'   \item{`$overlap_region`}{Named numeric vector: `lower` and `upper` bounds
#'     of the empirical common support region.}
#'   \item{`$outside_overlap`}{Data frame: `group`, `treated`, `n_outside`,
#'     `pct_outside` -- patients whose PS is outside the overlap region.}
#'   \item{`$positivity_flags`}{Data frame: `group`, `treated`, `n_near_zero`
#'     (PS < `trim_threshold`), `n_near_one` (PS > `1 - trim_threshold`),
#'     `pct_near_zero`, `pct_near_one`.}
#' }
#'
#' @seealso [ps_match()], [ps_weight()], [sa_rosenbaum()], [sa_evalue()],
#'   [sa_trim_sweep()]
#'
#' @examples
#' dta <- sample_ps_data(n = 400, seed = 7)
#' res <- sa_overlap(dta, score_col = "prob_t", treatment_col = "tavr")
#' res$overlap_region
#' res$positivity_flags
#'
#' # From a ps_match object
#' m   <- ps_match(dta)
#' res <- sa_overlap(m)
#' res$outside_overlap
#'
#' @export
sa_overlap <- function(x,
                       score_col     = "prob_t",
                       treatment_col = "tavr",
                       trim_threshold = 0.05) {

  # ---- Extract data frame --------------------------------------------------
  treated_level <- NULL
  if (is_ps_data(x)) {
    score_col     <- rlang::`%||%`(x$meta$score_col,     score_col)
    treatment_col <- rlang::`%||%`(x$meta$treatment_col, treatment_col)
    treated_level <- x$meta$treated_level
    dta <- x$data
  } else if (is.data.frame(x)) {
    dta <- x
  } else {
    rlang::abort("`x` must be a ps_data object or a data frame.", call. = FALSE)
  }

  # ---- Validate inputs -----------------------------------------------------
  .check_cols(dta, c(score_col, treatment_col))
  .check_binary(dta, treatment_col)
  .check_probability(dta, score_col)

  if (!is.numeric(trim_threshold) || length(trim_threshold) != 1L ||
      trim_threshold <= 0 || trim_threshold >= 0.5) {
    rlang::abort("`trim_threshold` must be a single number in (0, 0.5).",
                 call. = FALSE)
  }

  ps  <- dta[[score_col]]
  # A scored object may declare 0 (or FALSE) as the treated value.
  group_labels <- .binary_group_labels(dta[[treatment_col]], treated_level)
  trt <- as.integer(as.character(dta[[treatment_col]]) == group_labels[[2L]])

  # Drop rows where PS is NA before computing group sizes, so that n0/n1
  # match the summary statistics (which use na.rm=TRUE but still exclude NAs).
  # Separately report how many were dropped.
  n_ps_na <- sum(is.na(ps))
  if (n_ps_na > 0L) {
    rlang::warn(
      sprintf(
        "%d patient(s) with missing '%s' values excluded from overlap diagnostics.",
        n_ps_na, score_col
      )
    )
  }
  keep <- !is.na(ps)
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

  list(
    summary          = summary_tbl,
    overlap_region   = overlap_region,
    outside_overlap  = outside_overlap,
    positivity_flags = positivity_flags
  )
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
#'   `x$meta$stabilised` when available; defaults to `TRUE`.
#' @param score_col   Name of the propensity score column.  Taken from
#'   `x$meta$score_col` when `x` is a `ps_data` object.  Default `"prob_t"`.
#' @param treatment_col Name of the binary treatment column.  Taken from
#'   `x$meta$treatment_col` when `x` is a `ps_data` object`.  Default
#'   `"tavr"`.
#'
#' @return A data frame with one row per trim value and columns:
#' \describe{
#'   \item{`trim`}{The winsorisation threshold.}
#'   \item{`ess_control`}{Effective sample size for the control group.}
#'   \item{`ess_treated`}{Effective sample size for the treated group.}
#'   \item{`max_weight`}{Maximum weight across all patients after trimming.}
#'   \item{`sd_weight`}{Standard deviation of weights after trimming.}
#'   \item{`pct_trimmed`}{Percentage of patients whose weight was winsorised.}
#' }
#'
#' @seealso [ps_weight()], [sa_rosenbaum()], [sa_evalue()], [sa_overlap()]
#'
#' @examples
#' dta    <- sample_ps_data(n = 400, seed = 3)
#' w_obj  <- ps_weight(dta, estimand = "ATE")
#' sweep  <- sa_trim_sweep(w_obj)
#' head(sweep)
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
                          treatment_col = "tavr") {

  # ---- Extract data frame and metadata -------------------------------------
  if (is_ps_data(x)) {
    score_col     <- rlang::`%||%`(x$meta$score_col,     score_col)
    treatment_col <- rlang::`%||%`(x$meta$treatment_col, treatment_col)
    estimand_use  <- rlang::`%||%`(x$meta$estimand,      match.arg(estimand))
    stabilise     <- rlang::`%||%`(x$meta$stabilised,    stabilise)
    dta <- x$data
  } else if (is.data.frame(x)) {
    estimand_use <- match.arg(estimand)
    dta <- x
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
  trt <- as.integer(dta[[treatment_col]])
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

  if (stabilise) {
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

  do.call(rbind, rows)
}
