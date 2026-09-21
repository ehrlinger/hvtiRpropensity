###############################################################################
## help.R
##
## Package-level documentation for hvtiRpropensity.
###############################################################################

#' hvtiRpropensity: Propensity Score Methods for Cardiac Surgery Studies
#'
#' @description
#' Provides functions for propensity score and balancing score analysis used
#' in cardiac surgery comparative-effectiveness research.  Translates the
#' Cleveland Clinic HVTI SAS template library into a consistent R workflow.
#'
#' **Step 1 — estimate scores** (new in this version):
#' - [ps_logistic()] — binary treatment propensity score (logistic regression;
#'   supports multiply-imputed data via a stacked data frame).
#' - [ps_ordinal()] — ordered treatment propensity score (proportional-odds
#'   model via `MASS::polr`; equivalent to SAS `PROC LOGISTIC` default for
#'   ordinal responses).
#' - [ps_nominal()] — nominal treatment propensity score (multinomial logistic
#'   via `nnet::multinom`; equivalent to SAS `PROC LOGISTIC link=glogit`).
#' - [bs_continuous()] — balancing score for a continuous exposure
#'   (linear regression; equivalent to SAS `PROC REG` + `%grp` macro).
#' - [bs_count()] — balancing score for a count exposure (negative-binomial
#'   or Poisson regression; equivalent to SAS `PROC GENMOD dist=nb`).
#'
#' **Step 2 — balance and use scores**:
#' - [ps_match()] — 1:1 nearest-neighbour matching on a pre-computed score.
#'   Returns a `pair_id` column in `$data` for use with [sa_rosenbaum()].
#' - [ps_weight()] — IPTW weights (ATE / ATT / ATC) from a pre-computed score.
#'
#' **Step 3 — sensitivity analysis** (optional):
#' - [sa_rosenbaum()] — Rosenbaum gamma bounds for matched analyses
#'   (requires `rbounds`).
#' - [sa_evalue()] — E-values (VanderWeele & Ding 2017); no extra dependency.
#' - [sa_overlap()] — Overlap and positivity diagnostics on the PS distribution.
#' - [sa_trim_sweep()] — IPTW trim-threshold sensitivity sweep.
#'
#' All functions return a `ps_data` S3 object whose `$data` slot holds the
#' original dataset with score / weight columns appended, and whose `$tables`
#' slot holds balance diagnostics compatible with [hvtiPlotR::hv_mirror_hist()]
#' and [hvtiPlotR::hv_balance()].
#'
#' @section Typical two-step workflow (binary treatment):
#' ```r
#' library(hvtiRpropensity)
#'
#' # Step 1: estimate propensity score
#' dta   <- sample_ps_data(n = 500, seed = 42)
#' dta$prob_t <- NULL
#' score <- ps_logistic(
#'   tavr ~ age + female + ef + diabetes + hypertension,
#'   data = dta
#' )
#' print(score)
#'
#' # Step 2a: match on the estimated score
#' m <- ps_match(score$data, score_col = score$meta$score_col)
#' matched <- m$data[m$data$match == 1L, ]
#'
#' # Step 2b: or weight
#' w <- ps_weight(score$data, score_col = score$meta$score_col)
#' ```
#'
#' @seealso
#' - [ps_logistic()] — binary propensity score estimation
#' - [ps_ordinal()] — ordinal propensity score estimation
#' - [ps_nominal()] — nominal propensity score estimation
#' - [bs_continuous()] — continuous-outcome balancing score
#' - [bs_count()] — count-outcome balancing score
#' - [ps_match()] — nearest-neighbour 1:1 matching
#' - [ps_weight()] — IPTW weighting (ATE / ATT / ATC)
#' - [sa_rosenbaum()] — Rosenbaum sensitivity bounds (matched analyses)
#' - [sa_evalue()] — E-value sensitivity analysis
#' - [sa_overlap()] — overlap and positivity diagnostics
#' - [sa_trim_sweep()] — IPTW trim-threshold sensitivity sweep
#' - [sample_ps_data()] — synthetic dataset for examples / tests
#' - [is_ps_data()] — predicate for `ps_data` objects
#'
#' @name hvtiRpropensity-package
#' @aliases hvtiRpropensity
#' @keywords internal
"_PACKAGE"
