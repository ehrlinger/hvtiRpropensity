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
#' **Step 1: estimate scores**:
#' - [ps_logistic()]: binary treatment propensity score (logistic regression;
#'   supports multiply-imputed data via a stacked data frame).
#' - [ps_forest()]: binary treatment propensity score from a random forest
#'   (`randomForestSRC`), scored out-of-bag.
#' - [ps_ordinal()]: ordered treatment propensity score (proportional-odds
#'   model via `MASS::polr`; equivalent to SAS `PROC LOGISTIC` default for
#'   ordinal responses).
#' - [ps_nominal()]: nominal treatment propensity score (multinomial logistic
#'   via `nnet::multinom`; equivalent to SAS `PROC LOGISTIC link=glogit`).
#' - [bs_continuous()]: balancing score for a continuous exposure
#'   (linear regression; equivalent to SAS `PROC REG` + `%grp` macro).
#' - [bs_count()]: balancing score for a count exposure (negative-binomial
#'   or Poisson regression; equivalent to SAS `PROC GENMOD dist=nb`).
#'
#' **Outcome models**:
#' - [fit_logistic()]: binary, ordinal or nominal outcome-model bundle, with
#'   Rubin pooling over stacked imputations and Barnard-Rubin degrees of
#'   freedom, matching `mice::pool()`.
#' - [fit_linear()]: linear outcome-model bundle, with Rubin pooling.
#' - [validate_logistic()]: applies a saved binary bundle to a validation
#'   cohort without refitting.
#'
#' **Step 2: balance and use scores**:
#' - [ps_match()]: 1:1 nearest-neighbour matching on a pre-computed score,
#'   on the probability or the logit scale. `scale = "logit", caliper_sd = 0.2`
#'   is the Austin (2011) caliper. Returns a `pair_id` column in `$data` for
#'   use with [sa_rosenbaum()] and [ps_rmst()].
#' - [ps_weight()]: IPTW weights (ATE / ATT / ATC) from a pre-computed score.
#' - [ps_support()]: compares common-support and trimming rules.
#' - [ps_rmst()]: restricted mean survival contrast, weighted, or on the
#'   matched pairs of a [ps_match()] object.
#'
#' **Balance statistics**:
#' - [ps_stddiff()]: standardized differences by variable type, optionally
#'   weighted (CCF `%stddiff`).
#' - [ps_stddiff_perm()]: permutation reference for those differences (CCF
#'   `%stddiffci`).
#' - [ps_mw_var()]: matching-weight difference in means with a bootstrap SD
#'   (CCF `%mw_var`).
#'
#' **Step 3: sensitivity analysis** (optional):
#' - [sa_rosenbaum()]: Rosenbaum gamma bounds for matched analyses; how much
#'   hidden bias would overturn the result.
#' - [sa_evalue()]: E-values (VanderWeele & Ding 2017) for unmeasured
#'   confounding. Odds and hazard ratios need `rare`.
#' - [sa_overlap()]: overlap and positivity diagnostics on the score
#'   distribution.
#' - [sa_trim_sweep()]: how the effective sample size and the weights change
#'   across IPTW trim thresholds.
#'
#' **Data**: [sample_ps_data()], [sample_ps_data_ordinal()],
#' [sample_ps_data_nominal()] and [sample_ps_data_count()] generate synthetic
#' datasets; [is_ps_data()] tests for a `ps_data` object.
#'
#' Every function except the four data generators (which return data frames)
#' and [is_ps_data()] returns a `ps_data` S3 object. Its `$data` slot holds
#' the input data with any score, weight or flag columns appended, `$meta`
#' records the columns and parameters used, and `$tables` holds the
#' diagnostics: SMD tables, group counts, effective sample sizes, or a
#' sensitivity function's bounds and sweeps.
#'
#' `hvtiPlotR::hv_mirror_hist()` reads the `$data` of a [ps_match()] or
#' [ps_weight()] object directly. `hvtiPlotR::hv_balance()` does not read
#' `$tables` as they are: it needs one long data frame with `variable`,
#' `group` and `std_diff` columns, the difference in percent, so the SMD
#' tables are stacked and multiplied by 100 first. The vignette shows the
#' reshape.
#'
#' The SAS programs this package ports ship with it, in
#' `system.file("sas-templates", package = "hvtiRpropensity")`.
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
#' # Step 2a: match on the logit of the score, caliper 0.2 SD
#' m <- ps_match(score$data, score_col = "logit_t", scale = "logit",
#'               caliper_sd = 0.2, seed = 42)
#' matched <- m$data[m$data$match == 1L, ]
#'
#' # Step 2b: or weight
#' w <- ps_weight(score$data, score_col = score$meta$score_col)
#' ```
#'
#' @references Austin PC (2011). Optimal caliper widths for propensity-score
#'   matching when estimating differences in means and differences in
#'   proportions in observational studies. *Pharmaceutical Statistics*,
#'   10(2), 150--161.  \doi{10.1002/pst.433}
#'
#' @seealso
#' - [ps_logistic()], [ps_forest()], [ps_ordinal()], [ps_nominal()]:
#'   propensity score estimation
#' - [bs_continuous()], [bs_count()]: balancing scores for continuous and
#'   count exposures
#' - [fit_logistic()], [fit_linear()], [validate_logistic()]: outcome models
#' - [ps_match()], [ps_weight()]: matching and IPTW weighting
#' - [ps_stddiff()], [ps_stddiff_perm()], [ps_mw_var()]: balance statistics
#' - [ps_support()], [ps_rmst()]: support rules and survival contrasts
#' - [sa_rosenbaum()], [sa_evalue()], [sa_overlap()], [sa_trim_sweep()]:
#'   sensitivity analysis
#' - [sample_ps_data()]: synthetic dataset for examples / tests
#' - [is_ps_data()]: predicate for `ps_data` objects
#'
#' @name hvtiRpropensity-package
#' @aliases hvtiRpropensity
#' @keywords internal
"_PACKAGE"
