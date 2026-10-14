###############################################################################
## ps-score.R
##
## Propensity score estimation via logistic regression (binary, ordinal,
## and nominal treatment).
##
## Ports the model-fitting step from four SAS templates:
##   tp.lm.logistic_propensity_score.sas        (binary treatment, with MI)
##   tp.lm.logistic_propensity_score.nomi.sas   (binary treatment, no MI)
##   tp.lm.logistic_propensity_ordinal.sas      (ordered treatment)
##   tp.lm.logistic_propensity.polytomous.sas   (nominal treatment)
##
## Each function fits the propensity model and returns a ps_data subclass
## whose $data slot is the original dataset with probability column(s)
## appended.  The scored data frame can then be passed directly to
## ps_match() or ps_weight().
##
## Multiple imputation (MI) is supported by all three functions via a
## stacked "long" data frame with an imputation-index column
## (e.g. mice::complete(mids, action = "long")).  The model is fit
## separately on each imputed dataset, and per-patient predicted probabilities
## are averaged across imputations, mirroring the SAS approach of
## BY _IMPUTATION_ followed by PROC SUMMARY mean=.
##
## Functions:
##   ps_logistic()  -- binary treatment (0/1)
##   ps_ordinal()   -- ordered treatment (2+ levels)
##   ps_nominal()   -- nominal (unordered) treatment (2+ levels)
###############################################################################


# ---------------------------------------------------------------------------
# Internal: rank-based quintile / decile assignment
# ---------------------------------------------------------------------------

#' Assign rank-based quintile and decile columns
#'
#' Ranks the scored patients and assigns `ceiling(rank * k / n)`, which gives
#' groups of equal size when `n` divides by `k`. This does not reproduce the
#' SAS templates exactly. They assign
#' `quintile = int(_n_ / (nobs/5)) + 1; if quintile > 5 then quintile = 5;`,
#' which puts one patient too few in the first group and folds the surplus
#' into the last: with n = 10, SAS gives quintile sizes 1, 2, 2, 2, 3 and this
#' gives 2, 2, 2, 2, 2. A patient can sit one stratum apart in the two.
#'
#' @param score Numeric score; a missing score gets missing strata.
#' @return A named list: `quintile`, `decile`.
#' @keywords internal
.assign_ps_strata <- function(score) {
  scored <- !is.na(score)
  n <- sum(scored)
  quintile <- decile <- rep(NA_integer_, length(score))
  if (n) {
    rnk <- rank(score[scored], ties.method = "first")
    quintile[scored] <- pmin(ceiling(rnk * 5L / n), 5L)
    decile[scored] <- pmin(ceiling(rnk * 10L / n), 10L)
  }
  list(quintile = quintile, decile = decile)
}


# ---------------------------------------------------------------------------
# Internal: Li-Greene matching weight
# ---------------------------------------------------------------------------

#' Li-Greene matching weight (ATM)
#'
#' `min(p, 1 - p) / p` for treated patients and `min(p, 1 - p) / (1 - p)` for
#' controls, the `mt_wt` of the SAS templates (Li & Greene, 2013). Shared by
#' [ps_logistic()] and [ps_forest()] so the two cannot drift apart.
#'
#' @param p Propensity score, P(treated).
#' @param trt Integer 0/1 treatment indicator, the same length as `p`.
#' @return Numeric vector of matching weights.
#' @keywords internal
.matching_weight <- function(p, trt) {
  pmin(p, 1 - p) / (p * trt + (1 - p) * (1L - trt))
}


# ---------------------------------------------------------------------------
# Internal: collision-free model prediction name
# ---------------------------------------------------------------------------

#' Choose a private prediction prefix absent from a data frame
#'
#' @param data A data frame.
#' @return A single column-name prefix.
#' @keywords internal
.temporary_prediction_prefix <- function(data) {
  prefix <- ".hvti_prediction"
  while (any(startsWith(names(data), prefix))) prefix <- paste0(prefix, "_")
  prefix
}

.check_output_columns <- function(data, columns) {
  if (anyNA(columns) || any(!nzchar(columns)) || anyDuplicated(columns)) {
    rlang::abort("Output column names must be unique non-empty strings.", call. = FALSE)
  }
  existing <- intersect(columns, names(data))
  if (length(existing)) {
    rlang::abort(
      sprintf("Output column(s) already exist in `data`: %s.",
              paste(existing, collapse = ", ")),
      call. = FALSE
    )
  }
  invisible(NULL)
}


# ---------------------------------------------------------------------------
# ps_logistic
# ---------------------------------------------------------------------------

#' Estimate propensity scores via binary logistic regression
#'
#' Fits `stats::glm(..., family = binomial())` and appends the estimated
#' propensity score, its logit, and a matching weight to the dataset.
#' Supports multiply-imputed data in a stacked "long" format.
#'
#' **Multiple imputation.** When `imputation_col` is supplied, the model is
#' fit separately on each imputed dataset, and each patient's predicted
#' probabilities are averaged across the imputations on the probability scale,
#' not the logit scale. `logit_col` is the logit of that average. Leyrat et
#' al. (2019) call this the "across" approach. It is not Rubin's rules, which
#' pool coefficients and their variances; `$tables$estimates` holds the
#' coefficients pooled that way. The averaging mirrors the SAS template
#' workflow of `PROC LOGISTIC ... BY _IMPUTATION_` followed by
#' `PROC SUMMARY mean=`.
#'
#' `$data` is built from the lowest-numbered imputation, and its covariate
#' values are the ones the `smd` balance table measures. That table reports
#' balance in one imputed dataset, not across all of them.
#'
#' In Leyrat et al.'s simulations of inverse probability of treatment
#' weighting, only the "within" approach was unbiased in every scenario they
#' studied, and they recommend it with the outcome included in the imputation
#' model. To follow it, call `ps_logistic()` once per imputed dataset, without
#' `imputation_col`; match or weight within each dataset; estimate the
#' treatment effect in each; then pool those estimates and their variances
#' with Rubin's rules.
#'
#' **The matching weight** appended as `weight_col` is
#' `min(p, 1-p) / (p * trt + (1-p) * (1-trt))`,
#' which equals the `mt_wt` column produced in the SAS templates (Li &
#' Greene, 2013). It targets the ATM estimand (average treatment effect among
#' the matched). Neither [ps_match()] nor [ps_weight()] reads it; both work
#' from `score_col`. Pass it as `weight_col` to [ps_stddiff()] or
#' [ps_mw_var()], or use it as the weight in an outcome model.
#'
#' A patient with a missing treatment is left out of the model fit, as
#' [stats::glm()] drops the row, and out of the strata, the balance tables and
#' the group counts. The model needs no treatment to predict, so that patient
#' still gets `score_col` and `logit_col`; `weight_col`, `quintile` and
#' `decile` are missing. One warning states how many such patients there are,
#' and `$meta$n_missing_treatment` keeps the count. A patient missing a
#' predictor gets no score; `$tables$fit_status` counts those rows.
#'
#' @param formula        A formula with the binary treatment on the
#'   left-hand side, e.g. `tavr ~ age + female + ef + diabetes`.
#' @param data           A data frame.  If `imputation_col` is set, `data`
#'   must be a stacked "long" MI data frame as produced by
#'   `mice::complete(mids, action = "long")`.
#' @param treatment_col  Name of the binary treatment column: 0/1 or
#'   logical, or any two values when `treated_level` or `treatment_levels`
#'   names them.  If `NULL` (default), extracted from the LHS of `formula`.
#' @param id_col         Patient identifier column.  Required when
#'   `imputation_col` is set.  Default `"id"`.
#' @param imputation_col Name of the imputation-index column in a stacked
#'   MI data frame.  `NULL` (default) means a single complete dataset is
#'   supplied.
#' @param score_col      Output column name for the propensity score
#'   (probability of treatment).  Default `"prob_t"`.
#' @param logit_col      Output column name for the logit of the propensity
#'   score, `log(p / (1 - p))`.  Default `"logit_t"`.
#' @param weight_col     Output column name for the Li-Greene matching
#'   weight (ATM). Default `"mt_wt"`.
#' @param covariates     Character vector of covariate column names for SMD
#'   balance diagnostics, used as given.  If `NULL` (default), every numeric
#'   column that varies is used, except `treatment_col`, `score_col`,
#'   `logit_col`, `weight_col`, `id_col` and the fixed columns this package
#'   appends (`quintile`, `decile`, `iptw`, `match`, `pair_id` and the like).
#'   A predictor named in `formula` is always kept.
#' @param treatment_levels Complete binary treatment levels. `NULL` takes them
#'   from the data: the two values present when `treated_level` is given,
#'   otherwise 0/1 or `FALSE`/`TRUE`.
#' @param treated_level Level whose probability is the propensity score. When
#'   `treatment_levels` is supplied and this is `NULL`, its last value is used.
#'   Given alone, the treatment may hold any two values (a factor such as
#'   `"TAVR"`/`"SAVR"`, for example) and this names the treated one.
#'
#' @return An object of class `c("ps_logistic", "ps_data")` with:
#' \describe{
#'   \item{`$data`}{`data` (the lowest-numbered imputation when
#'     `imputation_col` is set, without that column) with `score_col`,
#'     `logit_col`, `weight_col`, `quintile`, and `decile` columns appended.
#'     `score_col` is the probability of `treated_level`. `quintile` and
#'     `decile` rank `score_col` among patients with a treatment. `weight_col`
#'     holds the Li-Greene matching weight, `min(p, 1 - p) / p` for treated
#'     patients and `min(p, 1 - p) / (1 - p)` for controls, which targets the
#'     ATM estimand (average treatment effect among the matched; Li & Greene,
#'     2013, \doi{10.1515/ijb-2012-0030}).}
#'   \item{`$meta`}{Named list: `formula`, `treatment_col`, `id_col`,
#'     `imputation_col`, `score_col`, `logit_col`, `weight_col`,
#'     `treatment_levels`, `treated_level`, `bundle_version`, `model_family`
#'     (`"binary"`), `package_versions`, `method` (`"logistic"` or
#'     `"logistic-MI"`), `n_imputations`, `n_total`, `n_missing_treatment`.}
#'   \item{`$tables`}{Named list: `smd` (`variable`, `smd`), `group_counts`,
#'     and the model tables of [fit_logistic()]: `estimates`, `covariance`,
#'     `by_imputation` and `fit_status`. `group_counts` has one row per
#'     treatment level, named as declared, with a logical `treated` flag and
#'     the count `n`.}
#'   \item{`$models`}{The fitted [stats::glm()] models, one per imputation,
#'     named by imputation index.}
#' }
#'
#' @references Leyrat C, Seaman SR, White IR, Douglas I, Smeeth L, Kim J,
#'   Resche-Rigon M, Carpenter JR, Williamson EJ (2019). Propensity score
#'   analysis with partially observed covariates: how should multiple
#'   imputation be used? *Statistical Methods in Medical Research*, 28(1),
#'   3--19. \doi{10.1177/0962280217713032}
#'
#'   Li L, Greene T (2013). A weighting analogue to pair matching in
#'   propensity score analysis. *The International Journal of Biostatistics*,
#'   9(2), 215--234. \doi{10.1515/ijb-2012-0030}
#'
#' @seealso [ps_match()], [ps_weight()], [ps_ordinal()], [ps_nominal()],
#'   [fit_logistic()], [sample_ps_data()]
#'
#' @examples
#' dta <- sample_ps_data(n = 200, seed = 42)
#' dta$prob_t <- NULL
#'
#' # --- Single complete dataset (mirrors tp.lm.logistic_propensity_score.nomi.sas)
#' # Equivalent to: PROC LOGISTIC data=built descending; model tavr = ...;
#' obj <- ps_logistic(
#'   tavr ~ age + female + ef + diabetes + hypertension,
#'   data = dta
#' )
#' print(obj)
#' summary(obj)
#'
#' # The function appends:
#' #   prob_t   -- propensity score (p-hat; SAS _p_ / _propen_)
#' #   logit_t  -- log(p/(1-p))    (SAS _logit_)
#' #   mt_wt    -- matching weight  (SAS mt_wt = min(p,1-p)/(p*trt+(1-p)*(1-trt)))
#' #   quintile -- rank-based quintile, ceiling(rank * 5 / n)
#' #   decile   -- rank-based decile
#' head(obj$data[, c("id", "tavr", "prob_t", "logit_t", "mt_wt",
#'                   "quintile", "decile")])
#'
#' # Pass the scored data to ps_match() for downstream matching
#' matched <- ps_match(obj$data, score_col = obj$meta$score_col, seed = 42)
#' nrow(matched$data[matched$data$match == 1L, ])
#'
#' # --- Multiply-imputed data (mirrors tp.lm.logistic_propensity_score.sas)
#' # Equivalent to:
#' #   PROC LOGISTIC data=built descending; BY _IMPUTATION_; model tavr = ...;
#' #   PROC SUMMARY data=decile; class ccfid; var _p_; output out=... mean=_propen;
#' \donttest{
#' # Simulate a stacked MI dataset (2 imputations, column "_Imputation_")
#' dta_mi <- rbind(
#'   cbind(dta, `_Imputation_` = 1L),
#'   cbind(dta, `_Imputation_` = 2L)
#' )
#' names(dta_mi)[names(dta_mi) == "_Imputation_"] <- "imp"
#'
#' obj_mi <- ps_logistic(
#'   tavr ~ age + female + ef + diabetes + hypertension,
#'   data           = dta_mi,
#'   imputation_col = "imp",
#'   id_col         = "id"
#' )
#' print(obj_mi)
#' # Per-patient PS is the average across the two imputed-dataset predictions,
#' # matching the PROC SUMMARY mean= step in the SAS template.
#' head(obj_mi$data[, c("id", "tavr", "prob_t")])
#' }
#'
#' @export
ps_logistic <- function(formula,
                        data,
                        treatment_col  = NULL,
                        id_col         = "id",
                        imputation_col = NULL,
                        score_col      = "prob_t",
                        logit_col      = "logit_t",
                        weight_col     = "mt_wt",
                        covariates     = NULL,
                        treatment_levels = NULL,
                        treated_level = NULL) {

  # ---- Input validation ---------------------------------------------------
  if (!inherits(formula, "formula")) {
    rlang::abort("`formula` must be an R formula.", call. = FALSE)
  }
  .check_df(data)

  if (is.null(treatment_col)) {
    treatment_col <- as.character(formula[[2L]])
  }
  .check_cols(data, treatment_col)
  if (is.null(treatment_levels) && !is.null(treated_level)) {
    treatment_levels <- .treatment_indicator(data, treatment_col, treated_level, warn_na = FALSE)$labels
  }
  if (is.null(treatment_levels)) {
    .check_binary(data, treatment_col, warn_na = FALSE)
    treatment_levels <- if (is.logical(data[[treatment_col]])) {
      c(FALSE, TRUE)
    } else {
      c(0, 1)
    }
  }
  treatment_levels <- as.character(treatment_levels)
  if (is.null(treated_level)) treated_level <- utils::tail(treatment_levels, 1L)
  treated_level <- as.character(treated_level)

  .check_output_columns(
    data,
    c(score_col, logit_col, weight_col, "quintile", "decile")
  )

  private_prefix <- .temporary_prediction_prefix(data)
  model <- fit_logistic(
    formula = formula,
    data = data,
    family = "binary",
    outcome_col = treatment_col,
    id_col = id_col,
    imputation_col = imputation_col,
    outcome_levels = treatment_levels,
    event_level = treated_level,
    prediction_prefix = private_prefix
  )
  base_data <- model$data
  probs <- base_data[[private_prefix]]
  base_data[[private_prefix]] <- NULL

  # ---- Append score columns -----------------------------------------------
  trt <- as.integer(as.character(base_data[[treatment_col]]) == treated_level)
  missing_trt <- is.na(trt)
  if (any(missing_trt)) {
    rlang::warn(sprintf(paste(
      "Column `%s` is missing for %d patient(s). They are left out of the model fit, the strata,",
      "the balance tables and the group counts; they keep a predicted `%s`."
    ), treatment_col, sum(missing_trt), score_col))
  }

  base_data[[score_col]]  <- probs
  base_data[[logit_col]]  <- log(probs / (1 - probs))
  base_data[[weight_col]] <- .matching_weight(probs, trt)

  # ---- Quintile / decile strata -------------------------------------------
  strata <- .assign_ps_strata(replace(probs, missing_trt, NA))
  base_data[["quintile"]] <- strata$quintile
  base_data[["decile"]]   <- strata$decile

  # ---- SMD diagnostics ----------------------------------------------------
  reserved <- c(treatment_col, score_col, logit_col, weight_col,
                id_col, "quintile", "decile")
  if (is.null(covariates)) {
    covariates <- .default_covariates(base_data, reserved, keep = .retained_predictors(model$meta$formula))
  }
  diagnostic_data <- base_data
  diagnostic_data[[treatment_col]] <- trt
  smd_tbl <- .smd_table(diagnostic_data[!missing_trt, , drop = FALSE], treatment_col, covariates)

  group_counts <- data.frame(
    group   = c(setdiff(treatment_levels, treated_level), treated_level),
    treated = c(FALSE, TRUE),
    n       = c(sum(trt == 0L, na.rm = TRUE), sum(trt == 1L, na.rm = TRUE))
  )

  # ---- Assemble object ----------------------------------------------------
  new_ps_data(
    data     = base_data,
    meta     = list(
      formula        = formula,
      treatment_col  = treatment_col,
      id_col         = id_col,
      imputation_col = imputation_col,
      score_col      = score_col,
      logit_col      = logit_col,
      weight_col     = weight_col,
      treatment_levels = treatment_levels,
      treated_level = treated_level,
      bundle_version = model$meta$bundle_version,
      model_family = model$meta$model_family,
      package_versions = model$meta$package_versions,
      method         = if (is.null(imputation_col)) "logistic"
                       else "logistic-MI",
      n_imputations  = model$meta$n_imputations,
      n_total        = nrow(base_data),
      n_missing_treatment = sum(missing_trt)
    ),
    tables   = c(
      list(smd = smd_tbl, group_counts = group_counts),
      model$tables
    ),
    models = model$models,
    subclass = "ps_logistic"
  )
}


#' @export
print.ps_logistic <- function(x, ...) {
  cat("<ps_logistic>\n")
  if (!is.null(x$meta$n_total))
    cat(sprintf("  N total     : %d\n", x$meta$n_total))
  if (!is.null(x$meta$treatment_col))
    cat(sprintf("  Treatment   : %s\n", x$meta$treatment_col))
  if (!is.null(x$meta$score_col))
    cat(sprintf("  PS column   : %s\n", x$meta$score_col))
  if (!is.null(x$meta$weight_col))
    cat(sprintf("  Weight col  : %s\n", x$meta$weight_col))
  if (!is.null(x$meta$method)) {
    method_str <- x$meta$method
    if (!is.null(x$meta$n_imputations) && x$meta$n_imputations > 1L)
      method_str <- sprintf("%s (%d imputations)", method_str,
                            x$meta$n_imputations)
    cat(sprintf("  Method      : %s\n", method_str))
  }
  if (length(x$tables) > 0L)
    cat("  Tables      :", paste(names(x$tables), collapse = ", "), "\n")
  invisible(x)
}


# ---------------------------------------------------------------------------
# ps_ordinal
# ---------------------------------------------------------------------------

#' Estimate propensity scores via ordinal (cumulative logit) logistic
#' regression
#'
#' Fits a proportional-odds cumulative logit model via [MASS::polr()] for an
#' ordered treatment variable and appends one probability column per
#' treatment level to the dataset.  The model is equivalent to the default
#' `PROC LOGISTIC` behaviour for ordinal responses in SAS.
#'
#' When `imputation_col` is supplied, models are fit on each imputed dataset
#' separately and per-level probabilities are averaged per patient across
#' imputations, on the probability scale. This is the "across" approach of
#' Leyrat et al. (2019), not Rubin's rules; [ps_logistic()] describes it and
#' the "within" alternative. `$data` and the `smd` table use the
#' lowest-numbered imputation's covariates.
#'
#' @param formula          A formula with the ordered factor treatment on the
#'   LHS, e.g. `nyha_grp ~ age + female + ef`.  The treatment must be an
#'   ordered factor, whose level order is the order the model uses, unless
#'   `treatment_levels` gives that order.  Any other treatment is an error:
#'   alphabetical order is not a safe guess for an ordinal scale.
#' @param data             A data frame (single or stacked MI).
#' @param treatment_col    Name of the ordinal treatment column.  `NULL`
#'   (default) extracts from `formula`.
#' @param id_col           Patient identifier column.  Required when
#'   `imputation_col` is set.  Default `"id"`.
#' @param imputation_col   Imputation-index column for stacked MI data.
#'   `NULL` (default) assumes a single complete dataset.
#' @param score_col_prefix Prefix for the per-level output columns.  Columns
#'   are named `<prefix>_<level>` for each level of the treatment.
#'   Default `"prob"`.
#' @param covariates       Covariate columns for the balance table. `NULL`
#'   (default) uses every numeric column other than the treatment, score,
#'   identifier and strata columns.
#' @param smd_pairs        Which pairs of treatment levels the balance table
#'   compares. `"reference"` (default) compares each level with the first
#'   level; `"adjacent"` compares each level with the one below it; `"all"`
#'   compares every pair. Give more than one to combine them, for example
#'   `c("reference", "adjacent")`.
#' @param treatment_levels Complete treatment levels, lowest first. `NULL`
#'   (default) takes the levels of the treatment, which must then be an ordered
#'   factor.
#'
#' @return An object of class `c("ps_ordinal", "ps_data")` with:
#' \describe{
#'   \item{`$data`}{`data` (the lowest-numbered imputation when
#'     `imputation_col` is set, without that column) with one
#'     `<score_col_prefix>_<level>` column per treatment level, the marginal
#'     probability of that level, and `quintile` and `decile`, which rank the
#'     probability of the highest level. The treatment column is returned as
#'     an ordered factor.}
#'   \item{`$meta`}{Named list: `formula`, `treatment_col`, `id_col`,
#'     `imputation_col`, `score_cols`, `score_col_prefix`, `levels`,
#'     `smd_pairs`, `treatment_levels`, `bundle_version`, `model_family`
#'     (`"ordinal"`), `package_versions`, `cumulative_direction`, `method`
#'     (`"ordinal-logistic"` or `"ordinal-logistic-MI"`), `n_imputations`,
#'     `n_total`.}
#'   \item{`$tables`}{Named list: `smd`, `group_counts` (`group`, `n`), and
#'     the model tables of [fit_logistic()]: `estimates`, `covariance`,
#'     `by_imputation` and `fit_status`. `smd` has one row per covariate and
#'     pair of levels: `variable`, `level`, `versus`, `smd`.}
#'   \item{`$models`}{The fitted models, one per imputation, named by
#'     imputation index: [MASS::polr()] fits, or [stats::glm()] fits for a
#'     two-level treatment.}
#' }
#'
#' @references Leyrat C, Seaman SR, White IR, Douglas I, Smeeth L, Kim J,
#'   Resche-Rigon M, Carpenter JR, Williamson EJ (2019). Propensity score
#'   analysis with partially observed covariates: how should multiple
#'   imputation be used? *Statistical Methods in Medical Research*, 28(1),
#'   3--19. \doi{10.1177/0962280217713032}
#'
#' @section Package dependency:
#'   Requires the \pkg{MASS} package (a recommended package shipped with
#'   base R).  Install with `install.packages("MASS")` if missing.
#'
#' @seealso [ps_logistic()], [ps_nominal()], [ps_match()]
#'
#' @examples
#' # Mirrors the SAS template workflow:
#' #   PROC LOGISTIC (cumulative logit, default for ordinal response)
#' #   followed by decomposition of cumulative to marginal probabilities
#' #   and rank-based quintile / decile assignment.
#' \donttest{
#' if (requireNamespace("MASS", quietly = TRUE)) {
#'   dta <- sample_ps_data_ordinal(n = 300, seed = 42)
#'   obj <- ps_ordinal(
#'     nyha_grp ~ age + female + ef + diabetes,
#'     data = dta
#'   )
#'   print(obj)
#'
#'   # Each level gets its own probability column (marginal, not cumulative).
#'   # The SAS template computes: p1=col1; p2=col2-col1; p3=1-col2.
#'   # ps_ordinal() performs this decomposition internally.
#'   head(obj$data[, c("id", "nyha_grp", "prob_I", "prob_II", "prob_III")])
#'
#'   # Quintile and decile columns are appended, ordered by p(highest level).
#'   table(obj$data$quintile)
#' }
#' }
#'
#' @export
ps_ordinal <- function(formula,
                       data,
                       treatment_col    = NULL,
                       id_col           = "id",
                       imputation_col   = NULL,
                       score_col_prefix = "prob",
                       covariates       = NULL,
                       treatment_levels = NULL,
                       smd_pairs        = "reference") {

  if (!requireNamespace("MASS", quietly = TRUE)) {
    rlang::abort(
      paste0("`ps_ordinal()` requires the MASS package.\n",
             "Install it with: install.packages(\"MASS\")"),
      call. = FALSE
    )
  }

  # ---- Input validation ---------------------------------------------------
  if (!inherits(formula, "formula")) {
    rlang::abort("`formula` must be an R formula.", call. = FALSE)
  }
  .check_df(data)

  if (is.null(treatment_col)) {
    treatment_col <- as.character(formula[[2L]])
  }
  .check_cols(data, treatment_col)
  smd_pairs <- match.arg(smd_pairs, c("reference", "adjacent", "all"), several.ok = TRUE)
  if (is.null(treatment_levels)) {
    if (!is.ordered(data[[treatment_col]])) {
      rlang::abort(
        c(
          sprintf("Column `%s` must be an ordered factor, or the level order must be given.", treatment_col),
          i = sprintf("Convert it with `ordered(%s, levels = c(...))`, lowest level first.", treatment_col),
          i = "Or pass `treatment_levels`, lowest level first."
        ),
        call. = FALSE
      )
    }
    treatment_levels <- levels(factor(data[[treatment_col]]))
  }
  treatment_levels <- as.character(treatment_levels)

  score_cols <- paste0(score_col_prefix, "_", treatment_levels)
  .check_output_columns(data, c(score_cols, "quintile", "decile"))

  private_prefix <- .temporary_prediction_prefix(data)
  model <- fit_logistic(
    formula = formula,
    data = data,
    family = "ordinal",
    outcome_col = treatment_col,
    id_col = id_col,
    imputation_col = imputation_col,
    outcome_levels = treatment_levels,
    prediction_prefix = private_prefix
  )
  base_data <- model$data
  lvls <- treatment_levels

  # ---- Append probability columns -----------------------------------------
  for (i in seq_along(lvls)) {
    private_col <- paste0(private_prefix, "_", lvls[[i]])
    base_data[[score_cols[[i]]]] <- base_data[[private_col]]
    base_data[[private_col]] <- NULL
  }

  # ---- Quintile / decile strata -------------------------------------------
  strata <- .assign_ps_strata(base_data[[utils::tail(score_cols, 1L)]])
  base_data[["quintile"]] <- strata$quintile
  base_data[["decile"]]   <- strata$decile

  # ---- SMD diagnostics ----------------------------------------------------
  smd_tbl <- .smd_table_levels(
    base_data, treatment_col, lvls, smd_pairs, covariates,
    reserved = c(score_cols, id_col, imputation_col, "quintile", "decile"),
    keep = .retained_predictors(model$meta$formula)
  )

  # ---- Group counts -------------------------------------------------------
  trt_vals     <- base_data[[treatment_col]]
  group_counts <- as.data.frame(table(group = trt_vals),
                                stringsAsFactors = FALSE)
  names(group_counts) <- c("group", "n")

  # ---- Assemble object ----------------------------------------------------
  new_ps_data(
    data     = base_data,
    meta     = list(
      formula          = formula,
      treatment_col    = treatment_col,
      id_col           = id_col,
      imputation_col   = imputation_col,
      score_cols       = score_cols,
      score_col_prefix = score_col_prefix,
      levels           = lvls,
      smd_pairs        = smd_pairs,
      treatment_levels = lvls,
      bundle_version   = model$meta$bundle_version,
      model_family     = model$meta$model_family,
      package_versions = model$meta$package_versions,
      cumulative_direction = model$meta$cumulative_direction,
      method           = if (is.null(imputation_col)) "ordinal-logistic"
                         else "ordinal-logistic-MI",
      n_imputations    = model$meta$n_imputations,
      n_total          = nrow(base_data)
    ),
    tables   = c(list(smd = smd_tbl, group_counts = group_counts), model$tables),
    models = model$models,
    subclass = "ps_ordinal"
  )
}


#' @export
print.ps_ordinal <- function(x, ...) {
  cat("<ps_ordinal>\n")
  if (!is.null(x$meta$n_total))
    cat(sprintf("  N total     : %d\n", x$meta$n_total))
  if (!is.null(x$meta$treatment_col))
    cat(sprintf("  Treatment   : %s (%d levels: %s)\n",
                x$meta$treatment_col,
                length(x$meta$levels),
                paste(x$meta$levels, collapse = " < ")))
  if (!is.null(x$meta$score_cols))
    cat(sprintf("  Score cols  : %s\n",
                paste(x$meta$score_cols, collapse = ", ")))
  if (!is.null(x$meta$method)) {
    method_str <- x$meta$method
    if (!is.null(x$meta$n_imputations) && x$meta$n_imputations > 1L)
      method_str <- sprintf("%s (%d imputations)", method_str,
                            x$meta$n_imputations)
    cat(sprintf("  Method      : %s\n", method_str))
  }
  if (length(x$tables) > 0L)
    cat("  Tables      :", paste(names(x$tables), collapse = ", "), "\n")
  invisible(x)
}


# ---------------------------------------------------------------------------
# ps_nominal
# ---------------------------------------------------------------------------

#' Estimate propensity scores via nominal (generalised logit) logistic
#' regression
#'
#' Fits a multinomial logistic regression model via [nnet::multinom()] for a
#' nominal (unordered) treatment variable, equivalent to SAS
#' `PROC LOGISTIC ... / link=glogit`.  One probability column per treatment
#' level is appended to the dataset.
#'
#' When `imputation_col` is supplied, models are fit on each imputed dataset
#' and per-level probabilities are averaged per patient across imputations, on
#' the probability scale. This is the "across" approach of Leyrat et al.
#' (2019), not Rubin's rules; [ps_logistic()] describes it and the "within"
#' alternative. `$data` and the `smd` table use the lowest-numbered
#' imputation's covariates.
#'
#' @param formula          A formula with the nominal treatment on the LHS,
#'   e.g. `rtyp ~ age + female + ef + diabetes`.
#' @param data             A data frame (single or stacked MI).
#' @param treatment_col    Name of the nominal treatment column.  `NULL`
#'   (default) extracts from `formula`.
#' @param id_col           Patient identifier column.  Required when
#'   `imputation_col` is set.  Default `"id"`.
#' @param imputation_col   Imputation-index column.  `NULL` (default) means
#'   a single complete dataset.
#' @param ref_level        Reference level for the multinomial model.  `NULL`
#'   (default) uses the first factor level, matching `REF=first` in SAS.
#' @param score_col_prefix Prefix for per-level output columns (`prob_<level>`
#'   by default).  Default `"prob"`.
#' @param trace            Logical.  If `FALSE` (default), suppresses
#'   [nnet::multinom()] iteration messages.
#' @param covariates       Covariate columns for the balance table. `NULL`
#'   (default) uses every numeric column other than the treatment, score,
#'   identifier and strata columns.
#' @param smd_pairs        Which pairs of treatment levels the balance table
#'   compares. `"reference"` (default) compares each level with `ref_level`;
#'   `"all"` compares every pair. The levels have no order, so there is no
#'   `"adjacent"` choice.
#' @param treatment_levels Complete nominal treatment levels. `NULL` preserves
#'   the levels inferred by the historical interface.
#'
#' @return An object of class `c("ps_nominal", "ps_data")` with:
#' \describe{
#'   \item{`$data`}{`data` (the lowest-numbered imputation when
#'     `imputation_col` is set, without that column) with one
#'     `<score_col_prefix>_<level>` column per treatment level appended,
#'     `ref_level` first. The treatment column is returned as a factor with
#'     `ref_level` as its first level.}
#'   \item{`$meta`}{Named list: `formula`, `treatment_col`, `id_col`,
#'     `imputation_col`, `score_cols`, `score_col_prefix`, `levels` (in model
#'     order, `ref_level` first), `smd_pairs`, `ref_level`,
#'     `treatment_levels` (as declared), `bundle_version`, `model_family`
#'     (`"nominal"`), `package_versions`, `method` (`"nominal-logistic"` or
#'     `"nominal-logistic-MI"`), `n_imputations`, `n_total`.}
#'   \item{`$tables`}{Named list: `smd`, `group_counts` (`group`, `n`), and
#'     the model tables of [fit_logistic()]: `estimates`, `covariance`,
#'     `by_imputation` and `fit_status`. `smd` has one row per covariate and
#'     pair of levels: `variable`, `level`, `versus`, `smd`.}
#'   \item{`$models`}{The fitted [nnet::multinom()] models, one per
#'     imputation, named by imputation index.}
#' }
#'
#' @references Leyrat C, Seaman SR, White IR, Douglas I, Smeeth L, Kim J,
#'   Resche-Rigon M, Carpenter JR, Williamson EJ (2019). Propensity score
#'   analysis with partially observed covariates: how should multiple
#'   imputation be used? *Statistical Methods in Medical Research*, 28(1),
#'   3--19. \doi{10.1177/0962280217713032}
#'
#' @section Package dependency:
#'   Requires \pkg{nnet} (a recommended R package, usually pre-installed).
#'   Install with `install.packages("nnet")` if missing.
#'
#' @seealso [ps_logistic()], [ps_ordinal()], [ps_match()]
#'
#' @examples
#' # Mirrors the SAS template workflow:
#' #   PROC LOGISTIC ... / link=glogit  (generalised logit)
#' # The first factor level is used as the reference category, matching
#' # SAS REF=first.  Change ref_level to match a different reference.
#' \donttest{
#' if (requireNamespace("nnet", quietly = TRUE)) {
#'   dta <- sample_ps_data_nominal(n = 300, seed = 42)
#'   obj <- ps_nominal(
#'     rtyp ~ age + female + ef + diabetes,
#'     data = dta
#'   )
#'   print(obj)
#'
#'   # One probability column per treatment level (analogous to p_cos, p_per,
#'   # p_dev, p_ce from the PROC TRANSPOSE step in the SAS template).
#'   head(obj$data[, c("id", "rtyp", "prob_COS", "prob_PER", "prob_DEV", "prob_CE")])
#'
#'   # Explicitly set a different reference level
#'   obj_ce <- ps_nominal(
#'     rtyp ~ age + female + ef + diabetes,
#'     data      = dta,
#'     ref_level = "CE"    # matches REF=last in SAS
#'   )
#'   obj_ce$meta$ref_level
#' }
#' }
#'
#' @export
ps_nominal <- function(formula,
                       data,
                       treatment_col    = NULL,
                       id_col           = "id",
                       imputation_col   = NULL,
                       ref_level        = NULL,
                       score_col_prefix = "prob",
                       trace            = FALSE,
                       covariates       = NULL,
                       treatment_levels = NULL,
                       smd_pairs        = "reference") {

  if (!requireNamespace("nnet", quietly = TRUE)) {
    rlang::abort(
      paste0("`ps_nominal()` requires the nnet package.\n",
             "Install it with: install.packages(\"nnet\")"),
      call. = FALSE
    )
  }

  # ---- Input validation ---------------------------------------------------
  if (!inherits(formula, "formula")) {
    rlang::abort("`formula` must be an R formula.", call. = FALSE)
  }
  .check_df(data)

  if (is.null(treatment_col)) {
    treatment_col <- as.character(formula[[2L]])
  }
  .check_cols(data, treatment_col)
  smd_pairs <- match.arg(smd_pairs, c("reference", "all"), several.ok = TRUE)
  if (is.null(treatment_levels)) {
    treatment_levels <- levels(factor(data[[treatment_col]]))
  }
  treatment_levels <- as.character(treatment_levels)
  if (is.null(ref_level)) ref_level <- treatment_levels[[1L]]
  ref_level <- as.character(ref_level)
  model_levels <- c(ref_level, setdiff(treatment_levels, ref_level))

  score_cols <- paste0(score_col_prefix, "_", model_levels)
  .check_output_columns(data, score_cols)

  private_prefix <- .temporary_prediction_prefix(data)
  model <- fit_logistic(
    formula = formula,
    data = data,
    family = "nominal",
    outcome_col = treatment_col,
    id_col = id_col,
    imputation_col = imputation_col,
    outcome_levels = treatment_levels,
    reference_level = ref_level,
    prediction_prefix = private_prefix,
    trace = trace
  )
  base_data <- model$data
  lvls <- model_levels

  # ---- Append probability columns -----------------------------------------
  for (i in seq_along(lvls)) {
    private_col <- paste0(private_prefix, "_", lvls[[i]])
    base_data[[score_cols[[i]]]] <- base_data[[private_col]]
    base_data[[private_col]] <- NULL
  }

  # ---- SMD diagnostics ----------------------------------------------------
  smd_tbl <- .smd_table_levels(
    base_data, treatment_col, lvls, smd_pairs, covariates,
    reserved = c(score_cols, id_col, imputation_col, "quintile", "decile"),
    keep = .retained_predictors(model$meta$formula)
  )

  # ---- Group counts -------------------------------------------------------
  trt_vals     <- base_data[[treatment_col]]
  group_counts <- as.data.frame(table(group = trt_vals),
                                stringsAsFactors = FALSE)
  names(group_counts) <- c("group", "n")

  # ---- Assemble object ----------------------------------------------------
  new_ps_data(
    data     = base_data,
    meta     = list(
      formula          = formula,
      treatment_col    = treatment_col,
      id_col           = id_col,
      imputation_col   = imputation_col,
      score_cols       = score_cols,
      score_col_prefix = score_col_prefix,
      levels           = lvls,
      smd_pairs        = smd_pairs,
      ref_level        = lvls[1L],
      treatment_levels = treatment_levels,
      bundle_version   = model$meta$bundle_version,
      model_family     = model$meta$model_family,
      package_versions = model$meta$package_versions,
      method           = if (is.null(imputation_col)) "nominal-logistic"
                         else "nominal-logistic-MI",
      n_imputations    = model$meta$n_imputations,
      n_total          = nrow(base_data)
    ),
    tables   = c(list(smd = smd_tbl, group_counts = group_counts), model$tables),
    models = model$models,
    subclass = "ps_nominal"
  )
}


#' @export
print.ps_nominal <- function(x, ...) {
  cat("<ps_nominal>\n")
  if (!is.null(x$meta$n_total))
    cat(sprintf("  N total     : %d\n", x$meta$n_total))
  if (!is.null(x$meta$treatment_col))
    cat(sprintf("  Treatment   : %s (%d levels)\n",
                x$meta$treatment_col,
                length(x$meta$levels)))
  if (!is.null(x$meta$ref_level))
    cat(sprintf("  Reference   : %s\n", x$meta$ref_level))
  if (!is.null(x$meta$score_cols))
    cat(sprintf("  Score cols  : %s\n",
                paste(x$meta$score_cols, collapse = ", ")))
  if (!is.null(x$meta$method)) {
    method_str <- x$meta$method
    if (!is.null(x$meta$n_imputations) && x$meta$n_imputations > 1L)
      method_str <- sprintf("%s (%d imputations)", method_str,
                            x$meta$n_imputations)
    cat(sprintf("  Method      : %s\n", method_str))
  }
  if (length(x$tables) > 0L)
    cat("  Tables      :", paste(names(x$tables), collapse = ", "), "\n")
  invisible(x)
}
