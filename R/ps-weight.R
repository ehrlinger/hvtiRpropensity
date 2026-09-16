###############################################################################
## ps-weight.R
##
## Inverse-probability-of-treatment weighting (IPTW).
##
## Ports the weighting logic used in the SAS programs that produce the
## `mt_wt` column consumed by `hvtiPlotR::hv_mirror_hist()` (weighted mode).
##
## Workflow:
##   1. dta  <- sample_ps_data()            # or your own data frame
##   2. obj  <- ps_weight(dta, ...)         # build ps_weight object
##   3. print(obj) / summary(obj)           # inspect diagnostics
##   4. obj$data                            # data with weight column appended
###############################################################################


#' Inverse-Probability-of-Treatment Weighting (IPTW)
#'
#' Computes IPTW weights from a pre-computed propensity score column and
#' appends them to the dataset.  Supports ATE, ATT, and ATC estimands.
#' Returns a `ps_weight` object whose `$tables` slot contains SMD tables
#' (unweighted and weighted) and effective-N diagnostics compatible with
#' [hvtiPlotR::hv_mirror_hist()] (weighted mode) and
#' [hvtiPlotR::hv_balance()].
#'
#' @param data           A data frame.  Must contain `treatment_col` and
#'   `score_col`.
#' @param treatment_col  Name of the binary treatment column (0/1 or logical).
#'   Default `"tavr"`.
#' @param score_col      Name of the numeric propensity score column (values
#'   in \[0, 1\]).  Default `"prob_t"`.
#' @param estimand       Causal estimand.  One of:
#'   - `"ATE"` (default) — average treatment effect.
#'     Treated weights = `1/ps`; control weights = `1/(1-ps)`.
#'   - `"ATT"` — average treatment effect on the treated.
#'     Treated weights = `1`; control weights = `ps/(1-ps)`.
#'   - `"ATC"` — average treatment effect on the controls.
#'     Treated weights = `(1-ps)/ps`; control weights = `1`.
#' @param stabilise      Logical.  If `TRUE` (default), weights are
#'   stabilised by multiplying by the marginal treatment probability.
#' @param trim           Optional numeric in (0, 0.5).  If supplied, weights
#'   outside the (`trim`, `1 - trim`) quantile range are winsorised.
#'   Default `NULL` (no trimming).
#' @param covariates     Character vector of covariate columns for SMD
#'   diagnostics.  If `NULL` (default), all numeric columns other than
#'   `treatment_col` and `score_col` are used.
#' @param weight_col     Name of the output weight column appended to `$data`.
#'   Default `"iptw"`.  If the column already exists it is overwritten.
#'
#' @return An object of class `c("ps_weight", "ps_data")` with:
#' \describe{
#'   \item{`$data`}{The input data frame with `weight_col` appended.}
#'   \item{`$meta`}{Named list: `treatment_col`, `score_col`, `weight_col`,
#'     `estimand`, `stabilised`, `trim`, `method`, `n_total`.}
#'   \item{`$tables`}{Named list: `smd_unweighted`, `smd_weighted`,
#'     `group_counts`, `effective_n`.}
#' }
#'
#' @seealso [ps_match()], [sample_ps_data()],
#'   `hvtiPlotR::hv_mirror_hist()`, `hvtiPlotR::hv_balance()`
#'
#' @examples
#' dta <- sample_ps_data(n = 200, seed = 42)
#'
#' # 1. ATE weights (default)
#' obj <- ps_weight(dta)
#' print(obj)
#' summary(obj)
#'
#' # 2. ATT weights, stabilised
#' obj_att <- ps_weight(dta, estimand = "ATT", stabilise = TRUE)
#' summary(obj_att)
#'
#' # 3. Extract data with weights appended
#' head(obj$data[, c("id", "tavr", "prob_t", "iptw")])
#'
#' @export
ps_weight <- function(data,
                      treatment_col = "tavr",
                      score_col     = "prob_t",
                      estimand      = c("ATE", "ATT", "ATC"),
                      stabilise     = TRUE,
                      trim          = NULL,
                      covariates    = NULL,
                      weight_col    = "iptw") {

  # ---- Input validation ---------------------------------------------------
  estimand <- match.arg(estimand)
  .check_df(data)
  .check_cols(data, c(treatment_col, score_col))
  .check_binary(data, treatment_col)
  .check_probability(data, score_col)

  if (!is.null(trim)) {
    if (!is.numeric(trim) || length(trim) != 1L || trim <= 0 || trim >= 0.5) {
      rlang::abort("`trim` must be a number in (0, 0.5) or NULL.", call. = FALSE)
    }
  }

  # ---- Compute raw IPTW weights ------------------------------------------
  ps  <- data[[score_col]]
  trt <- as.integer(data[[treatment_col]])

  w <- switch(estimand,
    ATE = ifelse(trt == 1L, 1 / ps, 1 / (1 - ps)),
    ATT = ifelse(trt == 1L, 1,      ps / (1 - ps)),
    ATC = ifelse(trt == 1L, (1 - ps) / ps, 1)
  )

  # ---- Stabilise ----------------------------------------------------------
  if (stabilise) {
    p_trt <- mean(trt == 1L)
    p_ctl <- 1 - p_trt
    stab  <- ifelse(trt == 1L, p_trt, p_ctl)
    w     <- w * stab
  }

  # ---- Winsorise (trim) ---------------------------------------------------
  if (!is.null(trim)) {
    q_lo <- stats::quantile(w, trim,       na.rm = TRUE)
    q_hi <- stats::quantile(w, 1 - trim,   na.rm = TRUE)
    w    <- pmin(pmax(w, q_lo), q_hi)
  }
  if (any(is.infinite(w))) {
    rlang::abort(
      sprintf(paste0("%d weight(s) are infinite because `%s` is exactly 0 or 1 for a patient in the arm ",
                     "that divides by it. Set `trim` to winsorise them, or remove those patients."),
              sum(is.infinite(w)), score_col),
      call = NULL
    )
  }

  # ---- Build output data frame -------------------------------------------
  out              <- data
  out[[weight_col]] <- w

  # ---- SMD diagnostics ---------------------------------------------------
  if (is.null(covariates)) {
    covariates <- setdiff(
      names(data)[vapply(data, is.numeric, logical(1L))],
      c(treatment_col, score_col, weight_col)
    )
  }

  smd_unweighted <- .smd_table(data, treatment_col, covariates)

  # Weighted SMD uses the weight column via a weighted mean / variance
  smd_weighted <- .smd_table(out, treatment_col, covariates, weight_col = weight_col)

  group_counts <- data.frame(
    group = c("control", "treated"),
    n     = c(sum(trt == 0L), sum(trt == 1L))
  )

  effective_n <- data.frame(
    group       = c("control", "treated"),
    n_effective = c(
      .effective_n(w[trt == 0L]),
      .effective_n(w[trt == 1L])
    )
  )

  # ---- Assemble ps_data object -------------------------------------------
  new_ps_data(
    data = out,
    meta = list(
      treatment_col = treatment_col,
      score_col     = score_col,
      weight_col    = weight_col,
      estimand      = estimand,
      stabilised    = stabilise,
      trim          = trim,
      method        = paste0("IPTW-", estimand),
      n_total       = nrow(data)
    ),
    tables = list(
      smd_unweighted = smd_unweighted,
      smd_weighted   = smd_weighted,
      group_counts   = group_counts,
      effective_n    = effective_n
    ),
    subclass = "ps_weight"
  )
}


# ---------------------------------------------------------------------------
# Internal helpers for ps_weight
# ---------------------------------------------------------------------------

#' Effective sample size for a weight vector
#' @keywords internal
.effective_n <- function(w) {
  if (length(w) == 0L) return(NA_real_)
  round(sum(w)^2 / sum(w^2), 1)
}
