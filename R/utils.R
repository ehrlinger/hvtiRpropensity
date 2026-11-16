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
#' @param warn_na   Whether to warn about missing values. A caller that counts
#'   missing treatment and missing score together, and states the total in one
#'   warning, sets this to `FALSE`.
#' @return Invisible NULL.
#' @keywords internal
.check_binary <- function(data, col, call_env = rlang::caller_env(), warn_na = TRUE) {
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
  if (warn_na && n_na > 0L) {
    rlang::warn(
      sprintf(
        "Column `%s` contains %d NA value(s). Those patients will be silently excluded from analysis.",
        col, n_na
      )
    )
  }
  invisible(NULL)
}


#' Warn once about rows excluded for a missing treatment or score
#'
#' @param excluded Logical vector, `TRUE` for a row left out of the
#'   computation.
#' @param treatment_col,score_col Column names, for the message.
#' @param what What the excluded rows are left out of, for the message.
#' @return `excluded`, invisibly.
#' @keywords internal
.warn_excluded <- function(excluded, treatment_col, score_col, what) {
  n <- sum(excluded)
  if (n > 0L) {
    rlang::warn(sprintf("%d row(s) have a missing `%s` or `%s` and are excluded from %s.",
                        n, treatment_col, score_col, what))
  }
  invisible(excluded)
}


# ---------------------------------------------------------------------------
# Probability-range check
# ---------------------------------------------------------------------------

#' Assert that a column contains valid probabilities (0 <= x <= 1)
#'
#' @param data      A data frame.
#' @param col       Name of the column to check.
#' @param call_env  Environment for the error call.
#' @param hint      Optional sentence appended to the range error.
#' @return Invisible NULL.
#' @keywords internal
.check_probability <- function(data, col, call_env = rlang::caller_env(), hint = NULL) {
  x <- data[[col]]
  if (!is.numeric(x)) {
    rlang::abort(
      sprintf("Column `%s` must be numeric.", col),
      call = call_env
    )
  }
  if (any(!is.na(x) & (x < 0 | x > 1))) {
    rlang::abort(
      paste(c(sprintf("Column `%s` must contain values in [0, 1].", col), hint), collapse = " "),
      call = call_env
    )
  }
  invisible(NULL)
}


#' The propensity score of a ps_data object on the probability scale
#'
#' A [ps_match()] object built with `scale = "logit"` names a logit column in
#' `$meta$score_col`. Readers that need a probability call this helper, which
#' applies [stats::plogis()] when `$meta$scale` is `"logit"` and otherwise
#' returns the column unchanged, then checks that every value lies in \[0, 1\].
#'
#' @param x         A `ps_data` object.
#' @param score_col Name of the score column. Default `x$meta$score_col`.
#' @param call_env  Environment for the error call.
#' @return A numeric vector of probabilities, one per row of `x$data`.
#' @keywords internal
.probability_score <- function(x, score_col = x$meta$score_col, call_env = rlang::caller_env()) {
  .check_cols(x$data, score_col, call_env)
  score <- x$data[[score_col]]
  if (identical(x$meta$scale, "logit")) score <- stats::plogis(score)
  .check_probability(stats::setNames(list(score), score_col), score_col, call_env)
  score
}


# ---------------------------------------------------------------------------
# SMD table for a data frame
# ---------------------------------------------------------------------------

#' Default covariates for a balance table
#'
#' Every numeric, logical, factor or character column except those in
#' `reserved`, the fixed output columns this package appends (`prob_t`,
#' `logit_t`, `mt_wt`, `iptw`, `match`, `pair_id`, `quintile`, `decile`, `bs`,
#' `cluster`), and constant columns. A column is constant when its non-missing
#' values take at most one value. A logical, factor or character column whose
#' non-missing values are all different is an identifier, not a covariate, and
#' is left out too. A function that generates columns under other names, such
#' as the `prob_<level>` scores of [ps_nominal()], passes those exact names in
#' `reserved`.
#'
#' Columns in `keep`, the predictors a model formula declares, are never
#' dropped by the fixed output names, the constant rule or the identifier rule;
#' only `reserved` can drop them.
#'
#' @param data     A data frame.
#' @param reserved Columns never used as covariates, such as the treatment,
#'   score and id columns.
#' @param keep     Columns protected from the fixed output names and the
#'   constant rule, typically `all.vars()` of a formula's right-hand side.
#' @return A character vector of column names.
#' @keywords internal
.default_covariates <- function(data, reserved = character(), keep = character()) {
  outputs <- c("prob_t", "logit_t", "mt_wt", "iptw", "match", "pair_id", "quintile", "decile", "bs", "cluster")
  usable <- vapply(data, function(x) is.numeric(x) || is.logical(x) || is.factor(x) || is.character(x),
                   logical(1L))
  cols <- setdiff(names(data)[usable], reserved)
  informative <- vapply(cols, function(col) {
    x <- stats::na.omit(data[[col]])
    n_values <- length(unique(x))
    n_values > 1L && (is.numeric(x) || n_values < length(x))
  }, logical(1L))
  cols[cols %in% keep | (informative & !cols %in% outputs)]
}


#' Build a tidy SMD table across covariates, through ps_stddiff()
#'
#' Each covariate is routed to the [ps_stddiff()] formula for its type. A
#' logical column, or a numeric column whose non-missing values are all 0 or
#' 1, takes the binary formula. A factor or character column takes the
#' categorical formula of Yang and Dalton (2012), one row per column, or the
#' binary formula when it has two levels. Every other numeric column is
#' Gaussian, with the denominator `sqrt((var1 + var0) / 2)`. With `weight_col`,
#' proportions and variances are weighted as [ps_stddiff()] describes. A
#' subset with a group absent, such as the matched rows of a match that found
#' no pairs, gets `NA` for every covariate rather than an error.
#'
#' @param data       A data frame.
#' @param treatment  Name of the binary treatment column.
#' @param covariates Character vector of covariate column names, used as
#'   given.  If `NULL`, the default set from `.default_covariates()` is used,
#'   excluding `treatment`, `weight_col` and `reserved`.
#' @param weight_col Optional name of a weight column.
#' @param reserved   Further columns never used as default covariates.
#' @param keep       Model predictors protected in the default set; see
#'   `.default_covariates()`.
#' @return A data frame with columns `variable` and `smd`, rounded to 4 places,
#'   one row per covariate.
#' @keywords internal
.smd_table <- function(data, treatment, covariates = NULL, weight_col = NULL, reserved = character(),
                       keep = character()) {
  if (is.null(covariates)) {
    covariates <- .default_covariates(data, c(treatment, weight_col, reserved), keep)
  }
  .check_cols(data, covariates)
  type <- vapply(covariates, function(v) {
    x <- data[[v]]
    if (is.logical(x) || (is.numeric(x) && all(x[!is.na(x)] %in% c(0, 1)))) {
      "binary"
    } else if (is.numeric(x)) {
      "gaussian"
    } else if (is.factor(x) || is.character(x)) {
      "categorical"
    } else {
      rlang::abort(sprintf(paste0("Covariate `%s` must be numeric, logical, factor or character to enter a ",
                                  "balance table; it is %s."), v, class(x)[1L]), call. = FALSE)
    }
  }, character(1L))
  grp <- data[[treatment]]
  smds <- if (length(covariates) == 0L || !all(c(0, 1) %in% grp)) {
    rep(NA_real_, length(covariates))
  } else {
    tbl <- ps_stddiff(data, treatment_col = treatment, gaussian = covariates[type == "gaussian"],
                      binary = covariates[type == "binary"], categorical = covariates[type == "categorical"],
                      weight_col = weight_col)$tables$stddiff
    tbl$stddiff[match(covariates, tbl$variable)]
  }

  data.frame(
    variable = covariates,
    smd      = round(unname(smds), 4L),
    stringsAsFactors = FALSE
  )
}


#' Winsorise IPTW weights within the arm or arms that are reweighted
#'
#' For the ATT, the treated arm is the reference and keeps weights of exactly
#' 1, so only the control weights are winsorised, at the (`trim`, `1 - trim`)
#' quantiles of the control weights. The ATC is the mirror image. For the ATE
#' both arms are reweighted and neither is a reference, so the weights are
#' winsorised at quantiles of all weights pooled across both arms, as Lee,
#' Lessler and Stuart (2011) do.
#'
#' @param w Numeric weight vector; `NA` marks an excluded patient.
#' @param trt Integer 0/1 treatment indicator, the same length as `w`.
#' @param estimand One of `"ATE"`, `"ATT"`, `"ATC"`.
#' @param trim A number in (0, 0.5).
#' @return `w`, winsorised.
#' @keywords internal
.trim_weights <- function(w, trt, estimand, trim) {
  arms <- switch(estimand, ATE = c(0L, 1L), ATT = 0L, ATC = 1L)
  i <- !is.na(w) & trt %in% arms
  if (any(i)) {
    q <- stats::quantile(w[i], c(trim, 1 - trim), names = FALSE)
    w[i] <- pmin(pmax(w[i], q[[1L]]), q[[2L]])
  }
  w
}


#' Level pairs for a multi-level balance table
#'
#' @param levels Treatment levels, reference first.
#' @param scheme Any of `"reference"`, `"adjacent"`, `"all"`.
#' @return A data frame of unique pairs with columns `level` and `versus`.
#' @keywords internal
.smd_level_pairs <- function(levels, scheme) {
  k <- length(levels)
  pairs <- list()
  if ("reference" %in% scheme && k > 1L) {
    pairs <- c(pairs, lapply(seq_len(k)[-1L], function(i) c(i, 1L)))
  }
  if ("adjacent" %in% scheme && k > 1L) {
    pairs <- c(pairs, lapply(seq_len(k)[-1L], function(i) c(i, i - 1L)))
  }
  if ("all" %in% scheme && k > 1L) {
    pairs <- c(pairs, utils::combn(k, 2L, function(ij) rev(ij), simplify = FALSE))
  }
  pairs <- unique(pairs)
  index <- do.call(rbind, pairs)
  index <- index[order(index[, 2L], index[, 1L]), , drop = FALSE]
  data.frame(level = levels[index[, 1L]], versus = levels[index[, 2L]], stringsAsFactors = FALSE)
}


#' Standardized differences between pairs of treatment levels
#'
#' Each pair is compared on the patients in those two levels only, with
#' `level` as the treated group and `versus` as the comparison group.
#'
#' @param data A data frame.
#' @param treatment Name of the treatment column.
#' @param levels Treatment levels, reference first.
#' @param scheme Any of `"reference"`, `"adjacent"`, `"all"`.
#' @param covariates Covariate columns. `NULL` uses `.default_covariates()`,
#'   excluding `treatment` and `reserved`.
#' @param reserved Columns never treated as covariates.
#' @param keep Model predictors protected in the default set; see
#'   `.default_covariates()`.
#' @return A data frame with columns `variable`, `level`, `versus`, `smd`.
#' @keywords internal
.smd_table_levels <- function(data, treatment, levels, scheme, covariates = NULL, reserved = character(),
                              keep = character()) {
  if (is.null(covariates)) {
    covariates <- .default_covariates(data, c(treatment, reserved), keep)
  }
  pairs <- .smd_level_pairs(levels, scheme)
  group <- as.character(data[[treatment]])
  tables <- lapply(seq_len(nrow(pairs)), function(i) {
    rows <- !is.na(group) & group %in% c(pairs$level[[i]], pairs$versus[[i]])
    pair_data <- data[rows, covariates, drop = FALSE]
    flag <- .temporary_prediction_prefix(pair_data)
    pair_data[[flag]] <- as.integer(group[rows] == pairs$level[[i]])
    smd <- .smd_table(pair_data, flag, covariates)
    data.frame(variable = smd$variable, level = rep(pairs$level[[i]], nrow(smd)),
               versus = rep(pairs$versus[[i]], nrow(smd)), smd = smd$smd, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, tables)
  rownames(out) <- NULL
  out
}


#' Group labels for a binary treatment column
#'
#' @param x A 0/1 or logical treatment vector.
#' @param treated_level The treated value. `NULL` means `1` or `TRUE`.
#' @return A character vector of length 2: the reference label, then the
#'   treated label.
#' @keywords internal
.binary_group_labels <- function(x, treated_level = NULL) {
  labels <- if (is.logical(x)) c("FALSE", "TRUE") else c("0", "1")
  if (is.null(treated_level)) treated_level <- labels[[2L]]
  treated_level <- as.character(treated_level)
  c(setdiff(labels, treated_level), treated_level)
}


#' Recode a two-valued treatment column to a 0/1 indicator
#'
#' With `treated_level = NULL` the column must be 0/1 or logical, and 1 or
#' `TRUE` is the treated value. With a `treated_level`, the column may hold any
#' two values and `treated_level` names the treated one. Both values must be
#' present in either case: a column holding one arm aborts.
#'
#' @param data A data frame.
#' @param col Name of the treatment column.
#' @param treated_level The treated value, or `NULL`.
#' @param call_env Environment for the error call.
#' @param warn_na Whether to warn about missing treatment values; see
#'   `.check_binary()`.
#' @return A list: `trt`, an integer 0/1 vector (`NA` where the column is
#'   missing), and `labels`, the reference then the treated value as character.
#' @keywords internal
.treatment_indicator <- function(data, col, treated_level = NULL, call_env = rlang::caller_env(),
                                 warn_na = TRUE) {
  x <- data[[col]]
  if (is.null(treated_level)) {
    .check_binary(data, col, call_env, warn_na)
    values <- unique(as.integer(x[!is.na(x)]))
    if (length(values) < 2L) {
      rlang::abort(
        sprintf("Column `%s` must hold both treatment values; it holds %s.", col,
                if (length(values)) sprintf("only %s", .binary_group_labels(x)[values + 1L]) else "none"),
        call = call_env
      )
    }
    return(list(trt = as.integer(x), labels = .binary_group_labels(x)))
  }
  if (length(treated_level) != 1L || is.na(treated_level)) {
    rlang::abort("`treated_level` must be a single non-missing value.", call = call_env)
  }
  treated_level <- as.character(treated_level)
  values <- unique(as.character(x[!is.na(x)]))
  if (length(values) != 2L) {
    rlang::abort(
      sprintf("Column `%s` must hold two treatment values; it holds %d.", col, length(values)),
      call = call_env
    )
  }
  if (!treated_level %in% values) {
    rlang::abort(
      sprintf("`treated_level` (\"%s\") is not a value of column `%s`.", treated_level, col),
      call = call_env
    )
  }
  n_na <- sum(is.na(x))
  if (warn_na && n_na > 0L) {
    rlang::warn(
      sprintf(
        "Column `%s` contains %d NA value(s). Those patients will be silently excluded from analysis.",
        col, n_na
      )
    )
  }
  list(trt = as.integer(as.character(x) == treated_level),
       labels = c(setdiff(values, treated_level), treated_level))
}


# ---------------------------------------------------------------------------
# Seeds
# ---------------------------------------------------------------------------

#' Validate a seed and return its absolute value
#'
#' Every function that draws random numbers passes its `seed` through here.
#' The draws are seeded with `withr::local_seed(abs(seed))` and
#' `randomForestSRC::rfsrc()` is given `-abs(seed)`, the sign it requires, so
#' `seed = 5` and `seed = -5` give the same result.
#'
#' @param seed `NULL`, or one whole number.
#' @param call_env Environment for the error call.
#' @return `NULL`, or `abs(seed)` as an integer.
#' @keywords internal
.seed_value <- function(seed, call_env = rlang::caller_env()) {
  if (is.null(seed)) return(NULL)
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) || seed != round(seed) ||
        abs(seed) > .Machine$integer.max) {
    rlang::abort("`seed` must be `NULL` or one whole number.", call = call_env)
  }
  as.integer(abs(seed))
}
