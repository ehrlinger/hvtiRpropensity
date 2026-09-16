###############################################################################
## utils.R
##
## Internal validation helpers shared across ps_*() constructors.
##
## None of these are exported; they are called as hvtiRpropensity:::fn()
## in tests that need to exercise them directly.
###############################################################################


# ---------------------------------------------------------------------------
# Column-existence checks
# ---------------------------------------------------------------------------

#' Assert that required columns exist in a data frame
#'
#' @param data      A data frame.
#' @param cols      Character vector of required column names.
#' @param call_env  Environment for the error call (use `rlang::caller_env()`).
#' @return Invisible NULL.  Signals a classed error on failure.
#' @keywords internal
.check_cols <- function(data, cols, call_env = rlang::caller_env()) {
  missing <- setdiff(cols, names(data))
  if (length(missing) > 0L) {
    rlang::abort(
      sprintf(
        "Required column(s) not found in `data`: %s",
        paste(missing, collapse = ", ")
      ),
      call = call_env
    )
  }
  invisible(NULL)
}


# ---------------------------------------------------------------------------
# Data-frame check
# ---------------------------------------------------------------------------

#' Assert that an object is a non-empty data frame
#'
#' @param data      Object to check.
#' @param call_env  Environment for the error call.
#' @return Invisible NULL.
#' @keywords internal
.check_df <- function(data, call_env = rlang::caller_env()) {
  if (!is.data.frame(data)) {
    rlang::abort("`data` must be a data frame.", call = call_env)
  }
  if (nrow(data) == 0L) {
    rlang::abort("`data` must contain at least one row.", call = call_env)
  }
  invisible(NULL)
}


# ---------------------------------------------------------------------------
# Binary-column check
# ---------------------------------------------------------------------------

#' Assert that a column contains only 0/1 (or logical) values
#'
#' Emits a warning (not an error) if NA values are present, because NA
#' treatment assignments silently drop patients from downstream matching and
#' weighting steps.
#'
#' @param data      A data frame.
#' @param col       Name of the column to check.
#' @param call_env  Environment for the error call.
#' @return Invisible NULL.
#' @keywords internal
.check_binary <- function(data, col, call_env = rlang::caller_env()) {
  x    <- data[[col]]
  vals <- x[!is.na(x)]
  ok   <- is.logical(x) ||
    (is.numeric(x) && all(vals %in% c(0, 1)))
  if (!ok) {
    rlang::abort(
      sprintf("Column `%s` must be binary (0/1 or logical).", col),
      call = call_env
    )
  }
  n_na <- sum(is.na(x))
  if (n_na > 0L) {
    rlang::warn(
      sprintf(
        "Column `%s` contains %d NA value(s). Those patients will be silently excluded from analysis.",
        col, n_na
      )
    )
  }
  invisible(NULL)
}


# ---------------------------------------------------------------------------
# Probability-range check
# ---------------------------------------------------------------------------

#' Assert that a column contains valid probabilities (0 <= x <= 1)
#'
#' @param data      A data frame.
#' @param col       Name of the column to check.
#' @param call_env  Environment for the error call.
#' @return Invisible NULL.
#' @keywords internal
.check_probability <- function(data, col, call_env = rlang::caller_env()) {
  x <- data[[col]]
  if (!is.numeric(x)) {
    rlang::abort(
      sprintf("Column `%s` must be numeric.", col),
      call = call_env
    )
  }
  if (any(!is.na(x) & (x < 0 | x > 1))) {
    rlang::abort(
      sprintf("Column `%s` must contain values in [0, 1].", col),
      call = call_env
    )
  }
  invisible(NULL)
}


# ---------------------------------------------------------------------------
# SMD table for a data frame
# ---------------------------------------------------------------------------

#' Build a tidy SMD table across covariates, through ps_stddiff()
#'
#' Every covariate is treated as Gaussian, so the denominator is
#' `sqrt((var1 + var0) / 2)`, and with `weight_col` the variances are the
#' `PROC MEANS` weighted ones, divided by `n - 1`. This is the one formula the
#' package uses; see [ps_stddiff()]. A subset with a group absent, such as the
#' matched rows of a match that found no pairs, gets `NA` for every covariate
#' rather than an error.
#'
#' @param data       A data frame.
#' @param treatment  Name of the binary treatment column.
#' @param covariates Character vector of covariate column names.  If `NULL`,
#'   all numeric columns other than `treatment` are used.
#' @param weight_col Optional name of a weight column.
#' @return A data frame with columns `variable` and `smd`, rounded to 4 places.
#' @keywords internal
.smd_table <- function(data, treatment, covariates = NULL, weight_col = NULL) {
  if (is.null(covariates)) {
    covariates <- setdiff(
      names(data)[vapply(data, is.numeric, logical(1))],
      c(treatment, weight_col)
    )
  }
  grp <- data[[treatment]]
  smds <- if (length(covariates) == 0L || !all(c(0, 1) %in% grp)) {
    rep(NA_real_, length(covariates))
  } else {
    ps_stddiff(data, treatment_col = treatment, gaussian = covariates, weight_col = weight_col)$tables$stddiff$stddiff
  }
  names(smds) <- covariates

  data.frame(
    variable = covariates,
    smd      = round(smds, 4L),
    stringsAsFactors = FALSE
  )
}
