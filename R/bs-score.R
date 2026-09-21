###############################################################################
## bs-score.R
##
## Balancing score estimation for continuous and count outcomes.
##
## Ports the model-fitting step from two SAS templates:
##   tp.rm.continuous.balncing_score.sas   (linear regression, PROC REG)
##   tp.pm.count.balncing_score.sas        (negative binomial, PROC GENMOD)
##
## The balancing score is the linear predictor from a "saturated" model
## containing all covariates thought to confound the exposure-outcome
## relationship.  Patients with similar balancing scores have similar
## expected exposure values, allowing direct within-stratum comparisons
## with reduced confounding — analogous to propensity score stratification
## for binary treatments.
##
## Workflow mirrors the SAS templates:
##   1. Fit saturated model on each imputed dataset (PROC REG / PROC GENMOD
##      BY _IMPUTATION_).
##   2. Keep the linear predictor (xbeta) or fitted value (predict).
##   3. Average across imputations per patient (PROC SUMMARY mean=).
##   4. Sort by balancing score; assign quintile, decile, and cluster
##      strata (the %grp macro in SAS).
##   5. Optionally create binary indicator columns stra_1 ... stra_k.
##
## Functions:
##   bs_continuous()  -- linear regression balancing score (PROC REG)
##   bs_count()       -- negative-binomial / Poisson balancing score (PROC GENMOD)
###############################################################################


# ---------------------------------------------------------------------------
# Internal: rank-based strata assignment
# ---------------------------------------------------------------------------

#' Assign quintile, decile, and k-cluster strata from a score vector
#'
#' Mirrors the SAS %grp macro:
#'   proc sort by balancing_score;
#'   quintile = int(_n_ / (nobs/5)) + 1;  if quintile > 5 then quintile = 5;
#'   cluster  = decile;
#'
#' @param score    Numeric vector of balancing scores.
#' @param n_strata Integer number of clusters (default 10 = deciles).
#' @return A named list: `quintile`, `decile`, `cluster`.
#' @keywords internal
.assign_bs_strata <- function(score, n_strata = 10L) {
  n   <- length(score)
  rnk <- rank(score, ties.method = "first")
  list(
    quintile = pmin(ceiling(rnk * 5L  / n), 5L),
    decile   = pmin(ceiling(rnk * 10L / n), 10L),
    cluster  = pmin(ceiling(rnk * as.integer(n_strata) / n),
                    as.integer(n_strata))
  )
}


# ---------------------------------------------------------------------------
# bs_continuous
# ---------------------------------------------------------------------------

#' Estimate a balancing score via linear regression
#'
#' Fits `stats::lm()` on a (possibly multiply-imputed) dataset, averages the
#' per-patient **fitted values** across imputations, and assigns rank-based
#' quintile, decile, and cluster strata.  Mirrors the workflow in
#' `tp.rm.continuous.balncing_score.sas` (PROC REG + PROC SUMMARY +
#' the `%grp` macro for binary stratum indicators).
#'
#' The balancing score is the expected value of the continuous exposure
#' (e.g. nadir haematocrit) given the full covariate set.  Patients in the
#' same stratum have similar expected exposure values, making within-stratum
#' comparisons less confounded.
#'
#' @param formula             A formula with the continuous outcome / exposure
#'   on the LHS, e.g.
#'   `hct_nadr ~ female + agee + in_bmi + hct_pr + iv_cpb`.
#' @param data                A data frame, or a stacked MI data frame when
#'   `imputation_col` is set.
#' @param outcome_col         Name of the continuous outcome column.  `NULL`
#'   (default) extracts from `formula`.
#' @param id_col              Patient identifier column.  Required when
#'   `imputation_col` is set.  Default `"id"`.
#' @param imputation_col      Imputation-index column for stacked MI data.
#'   `NULL` (default) assumes a single complete dataset.
#' @param score_col           Name of the output balancing score column
#'   (averaged fitted value).  Default `"bs"`.
#' @param n_strata            Number of strata (clusters).  Also controls how
#'   many binary indicator columns are created when
#'   `create_binary_strata = TRUE`.  Default `10L`.
#' @param strata_col          Name of the cluster/stratum column appended to
#'   `$data`.  Default `"cluster"`.
#' @param create_binary_strata Logical.  If `TRUE` (default), binary columns
#'   `stra_1` through `stra_<n_strata>` are appended, mirroring the `%grp`
#'   macro in the SAS template.
#' @param covariates          Covariate columns for diagnostics (not used in
#'   computation; reserved for future balance-checking helpers).
#'
#' @return An object of class `c("bs_continuous", "ps_data")` with:
#' \describe{
#'   \item{`$data`}{Base data frame with `score_col`, `"quintile"`,
#'     `"decile"`, `strata_col`, and (if `create_binary_strata`) `stra_1`
#'     through `stra_<n_strata>` appended.}
#'   \item{`$meta`}{Named list: `formula`, `outcome_col`, `id_col`,
#'     `imputation_col`, `score_col`, `n_strata`, `strata_col`, `method`,
#'     `n_imputations`, `n_total`.}
#'   \item{`$tables`}{Named list: `strata_counts`.}
#' }
#'
#' @seealso [bs_count()], [ps_logistic()], [sample_ps_data()]
#'
#' @examples
#' # Use ef (ejection fraction) as a stand-in for a continuous exposure
#' dta <- sample_ps_data(n = 300, seed = 42)
#' obj <- bs_continuous(
#'   ef ~ age + female + diabetes + hypertension,
#'   data    = dta,
#'   id_col  = "id"
#' )
#' print(obj)
#' table(obj$data$cluster)
#'
#' @export
bs_continuous <- function(formula,
                          data,
                          outcome_col          = NULL,
                          id_col               = "id",
                          imputation_col       = NULL,
                          score_col            = "bs",
                          n_strata             = 10L,
                          strata_col           = "cluster",
                          create_binary_strata = TRUE,
                          covariates           = NULL) {

  # ---- Input validation ---------------------------------------------------
  if (!inherits(formula, "formula")) {
    rlang::abort("`formula` must be an R formula.", call. = FALSE)
  }
  .check_df(data)
  if (is.null(outcome_col)) outcome_col <- as.character(formula[[2L]])
  .check_cols(data, outcome_col)
  if (!is.null(id_col))         .check_cols(data, id_col)
  if (!is.null(imputation_col)) .check_cols(data, imputation_col)
  if (!is.numeric(n_strata) || length(n_strata) != 1L || n_strata < 2L) {
    rlang::abort("`n_strata` must be an integer >= 2.", call. = FALSE)
  }
  n_strata <- as.integer(n_strata)

  # ---- Fit model(s) -------------------------------------------------------
  fit_fn <- function(df) {
    stats::lm(formula, data = df, na.action = stats::na.exclude)
  }
  pred_fn <- function(fit, df) {
    stats::predict(fit, newdata = df)
  }
  fitted <- .fit_imputations(
    data = data,
    response_col = outcome_col,
    id_col = id_col,
    imputation_col = imputation_col,
    fit_fn = fit_fn,
    predict_fn = pred_fn,
    predictor_cols = all.vars(stats::delete.response(stats::terms(formula))),
    require_stable_response = FALSE
  )
  preds <- fitted$predictions
  base_data <- fitted$data

  # ---- Append balancing score and strata ----------------------------------
  base_data[[score_col]] <- preds

  strata <- .assign_bs_strata(preds, n_strata = n_strata)
  base_data[["quintile"]] <- strata$quintile
  base_data[["decile"]]   <- strata$decile
  base_data[[strata_col]] <- strata$cluster

  if (create_binary_strata) {
    for (k in seq_len(n_strata)) {
      base_data[[paste0("stra_", k)]] <- as.integer(strata$cluster == k)
    }
  }

  # ---- Strata counts ------------------------------------------------------
  strata_counts <- as.data.frame(
    table(stratum = base_data[[strata_col]]),
    stringsAsFactors = FALSE
  )
  names(strata_counts) <- c("stratum", "n")

  # ---- Assemble object ----------------------------------------------------
  new_ps_data(
    data     = base_data,
    meta     = list(
      formula        = formula,
      outcome_col    = outcome_col,
      id_col         = id_col,
      imputation_col = imputation_col,
      score_col      = score_col,
      n_strata       = n_strata,
      strata_col     = strata_col,
      method         = if (is.null(imputation_col)) "balancing-linear"
                       else "balancing-linear-MI",
      n_imputations  = length(fitted$imputations),
      n_total        = nrow(base_data)
    ),
    tables   = list(
      strata_counts = strata_counts
    ),
    models = fitted$models,
    subclass = "bs_continuous"
  )
}


#' @export
print.bs_continuous <- function(x, ...) {
  cat("<bs_continuous>\n")
  if (!is.null(x$meta$n_total))
    cat(sprintf("  N total     : %d\n", x$meta$n_total))
  if (!is.null(x$meta$outcome_col))
    cat(sprintf("  Outcome     : %s\n", x$meta$outcome_col))
  if (!is.null(x$meta$score_col))
    cat(sprintf("  Score col   : %s\n", x$meta$score_col))
  if (!is.null(x$meta$n_strata))
    cat(sprintf("  Strata      : %d clusters (%s)\n",
                x$meta$n_strata, x$meta$strata_col))
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
# bs_count
# ---------------------------------------------------------------------------

#' Estimate a balancing score via negative-binomial or Poisson regression
#'
#' Fits a GLM with a count distribution ([MASS::glm.nb()] for negative
#' binomial, or [stats::glm()] with Poisson family) on a (possibly
#' multiply-imputed) dataset, averages the per-patient **linear predictors**
#' (log-scale, i.e. xbeta) across imputations, and assigns rank-based strata.
#' Mirrors `tp.pm.count.balncing_score.sas` (PROC GENMOD `dist=nb link=log`
#' + PROC SUMMARY + the `%grp` macro).
#'
#' The linear predictor (xbeta, log-scale) rather than the fitted count is
#' used as the balancing score because it is on an unbounded continuous scale
#' more suitable for rank-based stratification.
#'
#' @param formula             A formula with the count outcome on the LHS,
#'   e.g. `rbc_tot ~ female + agee + in_hct + hx_csurg + iv_cpb`.
#' @param data                A data frame, or a stacked MI data frame when
#'   `imputation_col` is set.
#' @param outcome_col         Name of the count outcome column.  `NULL`
#'   (default) extracts from `formula`.
#' @param id_col              Patient identifier column.  Required when
#'   `imputation_col` is set.  Default `"id"`.
#' @param imputation_col      Imputation-index column.  `NULL` (default)
#'   means a single complete dataset.
#' @param dist                Distribution for the GLM.  `"negbin"` (default)
#'   calls [MASS::glm.nb()] (negative binomial, log link);  `"poisson"` calls
#'   [stats::glm()] with `family = poisson(link = "log")`.
#' @param score_col           Name of the output balancing score column
#'   (averaged linear predictor on the log scale).  Default `"bs"`.
#' @param n_strata            Number of strata.  Default `10L`.
#' @param strata_col          Name of the cluster/stratum column.  Default
#'   `"cluster"`.
#' @param create_binary_strata Logical.  If `TRUE` (default), binary
#'   indicator columns `stra_1` through `stra_<n_strata>` are appended.
#' @param covariates          Covariate columns for diagnostics.
#'
#' @return An object of class `c("bs_count", "ps_data")` with:
#' \describe{
#'   \item{`$data`}{Base data frame with `score_col`, `"quintile"`,
#'     `"decile"`, `strata_col`, and (if `create_binary_strata`)
#'     `stra_1`...`stra_k` appended.}
#'   \item{`$meta`}{Named list: `formula`, `outcome_col`, `id_col`,
#'     `imputation_col`, `score_col`, `dist`, `n_strata`, `strata_col`,
#'     `bundle_version`, `model_family`, `package_versions`, `method`,
#'     `n_imputations`, `n_total`.}
#'   \item{`$tables`}{Named list: `strata_counts`.}
#' }
#'
#' @section Package dependency:
#'   `dist = "negbin"` (the default) requires the \pkg{MASS} package (a
#'   recommended R package, usually pre-installed).  Install with
#'   `install.packages("MASS")` if missing, or use `dist = "poisson"`.
#'
#' @seealso [bs_continuous()], [ps_logistic()]
#'
#' @examples
#' # Mirrors the SAS template workflow:
#' #   PROC GENMOD dist=nb link=log (negative binomial for overdispersed counts)
#' #   xbeta (log-scale linear predictor) used as the balancing score
#' #   (not the fitted count -- log scale is more suitable for stratification).
#' # The %grp macro binary indicators stra_1 ... stra_k are created when
#' # create_binary_strata = TRUE (the default).
#' \donttest{
#' dta <- sample_ps_data_count(n = 300, seed = 42)
#'
#' # Poisson fallback (no MASS dependency)
#' obj <- bs_count(
#'   rbc_tot ~ female + age + diabetes + hypertension,
#'   data = dta,
#'   dist = "poisson"
#' )
#' print(obj)
#' table(obj$data$cluster)    # 10 equal-ish strata by default
#'
#' # Binary stratum indicators mirror the SAS %grp macro output
#' head(obj$data[, c("id", "bs", "cluster", "stra_1", "stra_2", "stra_10")])
#'
#' # Negative binomial (accounts for overdispersion; requires MASS)
#' if (requireNamespace("MASS", quietly = TRUE)) {
#'   obj_nb <- bs_count(
#'     rbc_tot ~ female + age + diabetes + hypertension,
#'     data = dta,
#'     dist = "negbin"
#'   )
#'   print(obj_nb)
#' }
#' }
#'
#' @export
bs_count <- function(formula,
                     data,
                     outcome_col          = NULL,
                     id_col               = "id",
                     imputation_col       = NULL,
                     dist                 = c("negbin", "poisson"),
                     score_col            = "bs",
                     n_strata             = 10L,
                     strata_col           = "cluster",
                     create_binary_strata = TRUE,
                     covariates           = NULL) {

  dist <- match.arg(dist)
  model_engine <- if (identical(dist, "negbin")) "MASS" else "stats"

  if (dist == "negbin" && !requireNamespace("MASS", quietly = TRUE)) {
    rlang::abort(
      paste0(
        "`bs_count(dist = 'negbin')` requires the MASS package.\n",
        "Install it with: install.packages(\"MASS\"), ",
        "or use dist = 'poisson'."
      ),
      call. = FALSE
    )
  }

  # ---- Input validation ---------------------------------------------------
  if (!inherits(formula, "formula")) {
    rlang::abort("`formula` must be an R formula.", call. = FALSE)
  }
  .check_df(data)
  if (is.null(outcome_col)) outcome_col <- as.character(formula[[2L]])
  .check_cols(data, outcome_col)
  if (!is.null(id_col))         .check_cols(data, id_col)
  if (!is.null(imputation_col)) .check_cols(data, imputation_col)
  if (!is.numeric(n_strata) || length(n_strata) != 1L || n_strata < 2L) {
    rlang::abort("`n_strata` must be an integer >= 2.", call. = FALSE)
  }
  n_strata <- as.integer(n_strata)

  # ---- Fit function -------------------------------------------------------
  fit_fn <- if (dist == "negbin") {
    function(df) {
      MASS::glm.nb(formula, data = df, na.action = stats::na.exclude)
    }
  } else {
    function(df) {
      stats::glm(
        formula,
        data = df,
        family = stats::poisson(link = "log"),
        na.action = stats::na.exclude
      )
    }
  }
  pred_fn <- function(fit, df) {
    stats::predict(fit, newdata = df, type = "link")
  }
  fitted <- .fit_imputations(
    data = data,
    response_col = outcome_col,
    id_col = id_col,
    imputation_col = imputation_col,
    fit_fn = fit_fn,
    predict_fn = pred_fn,
    predictor_cols = all.vars(stats::delete.response(stats::terms(formula))),
    require_stable_response = FALSE
  )
  inference <- .model_inference(fitted$models, "count")
  xbeta <- fitted$predictions
  base_data <- fitted$data
  theta <- if (identical(dist, "negbin")) {
    vapply(fitted$models, `[[`, numeric(1L), "theta")
  } else {
    NULL
  }

  # ---- Append balancing score and strata ----------------------------------
  base_data[[score_col]] <- xbeta

  strata <- .assign_bs_strata(xbeta, n_strata = n_strata)
  base_data[["quintile"]] <- strata$quintile
  base_data[["decile"]]   <- strata$decile
  base_data[[strata_col]] <- strata$cluster

  if (create_binary_strata) {
    for (k in seq_len(n_strata)) {
      base_data[[paste0("stra_", k)]] <- as.integer(strata$cluster == k)
    }
  }

  # ---- Strata counts ------------------------------------------------------
  strata_counts <- as.data.frame(
    table(stratum = base_data[[strata_col]]),
    stringsAsFactors = FALSE
  )
  names(strata_counts) <- c("stratum", "n")

  # ---- Assemble object ----------------------------------------------------
  new_ps_data(
    data     = base_data,
    meta     = list(
      formula        = formula,
      outcome_col    = outcome_col,
      id_col         = id_col,
      imputation_col = imputation_col,
      score_col      = score_col,
      dist           = dist,
      theta          = theta,
      n_strata       = n_strata,
      strata_col     = strata_col,
      bundle_version = 1L,
      model_family   = "count",
      package_versions = stats::setNames(c(
        R = as.character(getRversion()),
        hvtiRpropensity = as.character(utils::packageVersion("hvtiRpropensity")),
        model = as.character(utils::packageVersion(model_engine))
      ), c("R", "hvtiRpropensity", model_engine)),
      method         = if (is.null(imputation_col))
                         paste0("balancing-", dist)
                       else paste0("balancing-", dist, "-MI"),
      n_imputations  = length(fitted$imputations),
      n_total        = nrow(base_data)
    ),
    tables   = list(
      strata_counts = strata_counts,
      estimates = inference$estimates,
      covariance = inference$covariance,
      by_imputation = inference$by_imputation,
      fit_status = fitted$status
    ),
    models = fitted$models,
    subclass = "bs_count"
  )
}


#' @export
print.bs_count <- function(x, ...) {
  cat("<bs_count>\n")
  if (!is.null(x$meta$n_total))
    cat(sprintf("  N total     : %d\n", x$meta$n_total))
  if (!is.null(x$meta$outcome_col))
    cat(sprintf("  Outcome     : %s\n", x$meta$outcome_col))
  if (!is.null(x$meta$score_col))
    cat(sprintf("  Score col   : %s\n", x$meta$score_col))
  if (!is.null(x$meta$dist))
    cat(sprintf("  Distribution: %s\n", x$meta$dist))
  if (!is.null(x$meta$n_strata))
    cat(sprintf("  Strata      : %d clusters (%s)\n",
                x$meta$n_strata, x$meta$strata_col))
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
