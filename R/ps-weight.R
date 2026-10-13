###############################################################################
## ps-weight.R
##
## Inverse-probability-of-treatment weighting (IPTW).
##
## Computes the IPTW column (`iptw` by default) that
## `hvtiPlotR::hv_mirror_hist()` reads in weighted mode. The Li-Greene matching
## weight `mt_wt` comes from ps_logistic() and
## ps_forest(), not from here.
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
#' appends them to the dataset.  `estimand` picks the target: the average
#' treatment effect in the whole sample (ATE, the default), in the treated
#' (ATT) or in the controls (ATC). For the matched population (ATM), use the
#' `mt_wt` column from [ps_logistic()] or [ps_forest()] instead.
#' Returns a `ps_weight` object whose `$tables` slot contains SMD tables
#' (unweighted and weighted) and effective-N diagnostics. Its `$data` feeds
#' `hvtiPlotR::hv_mirror_hist()` in weighted mode; `hvtiPlotR::hv_balance()`
#' needs the SMD tables stacked into long form (see the vignette).
#'
#' A patient with a missing treatment or a missing score cannot be weighted.
#' Those rows stay in `$data` with a missing `weight_col`, are left out of the
#' stabilising treatment proportion, the SMD tables, the group counts and the
#' effective sample sizes, and one warning states how many there are.
#' `$meta$n_excluded` keeps the count.
#'
#' @param data           A data frame.  Must contain `treatment_col` and
#'   `score_col`.
#' @param treatment_col  Name of the binary treatment column (0/1 or logical).
#'   Default `"tavr"`.
#' @param score_col      Name of the numeric propensity score column (values
#'   in \[0, 1\]).  Default `"prob_t"`.
#' @param estimand       Causal estimand.  One of:
#'   - `"ATE"` (default): average treatment effect.
#'     Treated weights = `1/ps`; control weights = `1/(1-ps)`.
#'   - `"ATT"`: average treatment effect on the treated.
#'     Treated weights = `1`; control weights = `ps/(1-ps)`.
#'   - `"ATC"`: average treatment effect on the controls.
#'     Treated weights = `(1-ps)/ps`; control weights = `1`.
#' @param stabilise      Logical.  If `TRUE` (default) and `estimand = "ATE"`,
#'   each weight is multiplied by the marginal probability of the patient's
#'   own arm (`P(treated)` or `P(control)`).  Stabilisation is defined for the
#'   ATE only: for `"ATT"` and `"ATC"` it has no effect, so the reference arm
#'   keeps weights of exactly 1.  `$meta$stabilised` records whether it was
#'   applied.
#' @param trim           Optional numeric in (0, 0.5).  If supplied, weights
#'   outside the (`trim`, `1 - trim`) quantile range are winsorised.
#'   Default `NULL` (no trimming).
#' @param covariates     Character vector of covariate columns for SMD
#'   diagnostics, used as given.  If `NULL` (default), every numeric column
#'   that varies is used, except `treatment_col`, `score_col`, `weight_col`,
#'   `id_col` and the fixed columns this package appends (`prob_t`,
#'   `logit_t`, `mt_wt`, `iptw`, `match`, `pair_id`, `quintile`, `decile`,
#'   `bs`, `cluster`).
#' @param weight_col     Name of the output weight column appended to `$data`.
#'   Default `"iptw"`.  If the column already exists it is overwritten.
#'
#' @param treated_level The value of `treatment_col` that is treated. `NULL`
#'   (default) requires a 0/1 or logical column and takes 1 or `TRUE`. Give a
#'   value, for example `"transcatheter"`, to use a column holding any two
#'   values.
#' @param id_col The patient identifier column, left out of the default
#'   covariates. Default `"id"`; it need not exist in `data`.
#' @return An object of class `c("ps_weight", "ps_data")` with:
#' \describe{
#'   \item{`$data`}{The input data frame with `weight_col` appended.}
#'   \item{`$meta`}{Named list: `treatment_col`, `treated_level` (the value of
#'     `treatment_col` taken as treated), `score_col`, `weight_col`,
#'     `estimand`, `stabilised`, `trim`, `method`, `n_total`, `n_excluded`
#'     (rows with a missing treatment or score).}
#'   \item{`$tables`}{Named list: `smd_unweighted`, `smd_weighted`
#'     (`variable`, `smd`), `group_counts` (`group`, `treated`, `n`) and
#'     `effective_n` (`group`, `treated`, `n_effective`, the Kish effective
#'     sample size). Group tables name each group by its value in
#'     `treatment_col` (`0`/`1`, `FALSE`/`TRUE`, or the two values when
#'     `treated_level` is given) and carry a logical `treated` column.}
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
#' # 2. ATT weights: treated patients keep a weight of 1
#' obj_att <- ps_weight(dta, estimand = "ATT")
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
                      weight_col    = "iptw",
                      treated_level = NULL,
                      id_col        = "id") {

  # ---- Input validation ---------------------------------------------------
  estimand <- match.arg(estimand)
  .check_df(data)
  .check_cols(data, c(treatment_col, score_col))
  treatment <- .treatment_indicator(data, treatment_col, treated_level, warn_na = FALSE)
  original_treatment <- data[[treatment_col]]
  data[[treatment_col]] <- treatment$trt
  .check_probability(data, score_col)

  if (!is.null(trim)) {
    if (!is.numeric(trim) || length(trim) != 1L || is.na(trim) || trim <= 0 || trim >= 0.5) {
      rlang::abort("`trim` must be a number in (0, 0.5) or NULL.", call. = FALSE)
    }
  }

  # ---- Compute raw IPTW weights ------------------------------------------
  ps  <- data[[score_col]]
  trt <- as.integer(data[[treatment_col]])
  excluded <- .warn_excluded(is.na(trt) | is.na(ps), treatment_col, score_col,
                             "weighting and the balance tables, and get a missing weight")
  kept <- !excluded

  w <- switch(estimand,
    ATE = ifelse(trt == 1L, 1 / ps, 1 / (1 - ps)),
    ATT = ifelse(trt == 1L, 1,      ps / (1 - ps)),
    ATC = ifelse(trt == 1L, (1 - ps) / ps, 1)
  )
  w[excluded] <- NA_real_

  # ---- Stabilise (ATE only) -----------------------------------------------
  # Stabilising by the marginal treatment probability is defined for the ATE.
  # For the ATT and ATC the reference arm's weights stay exactly 1.
  stabilised <- isTRUE(stabilise) && estimand == "ATE"
  if (stabilised) {
    p_trt <- mean(trt[kept] == 1L)
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
    covariates <- .default_covariates(data, c(treatment_col, score_col, weight_col, id_col))
  }

  smd_unweighted <- .smd_table(data[kept, , drop = FALSE], treatment_col, covariates)

  # Weighted SMD uses the weight column via a weighted mean / variance
  smd_weighted <- .smd_table(out[kept, , drop = FALSE], treatment_col, covariates, weight_col = weight_col)

  group_labels <- treatment$labels
  out[[treatment_col]] <- original_treatment
  group_counts <- data.frame(
    group   = group_labels,
    treated = c(FALSE, TRUE),
    n       = c(sum(kept & trt == 0L), sum(kept & trt == 1L))
  )

  effective_n <- data.frame(
    group       = group_labels,
    treated     = c(FALSE, TRUE),
    n_effective = c(
      .effective_n(w[kept & trt == 0L]),
      .effective_n(w[kept & trt == 1L])
    )
  )

  # ---- Assemble ps_data object -------------------------------------------
  new_ps_data(
    data = out,
    meta = list(
      treatment_col = treatment_col,
      treated_level = treatment$labels[[2L]],
      score_col     = score_col,
      weight_col    = weight_col,
      estimand      = estimand,
      stabilised    = stabilised,
      trim          = trim,
      method        = paste0("IPTW-", estimand),
      n_total       = nrow(data),
      n_excluded    = sum(excluded)
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
