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
#' @return Invisible NULL.  Signals a classed error on failure.
#' @keywords internal
.check_cols <- function(data, cols) {
  missing <- setdiff(cols, names(data))
  if (length(missing) > 0L) {
    stop(errorCondition(
      paste0(.err_prefix(), sprintf(
        "Required column(s) not found in `data`: %s",
        paste(missing, collapse = ", ")
      )),
      class = "hvtiRpropensity_missing_column",
      call = NULL
    ))
  }
  invisible(NULL)
}


# ---------------------------------------------------------------------------
# Data-frame check
# ---------------------------------------------------------------------------

#' Assert that an object is a non-empty data frame
#'
#' @param data      Object to check.
#' @return Invisible NULL.
#' @keywords internal
.check_df <- function(data) {
  if (!is.data.frame(data)) {
    stop(
      paste0(.err_prefix(), "`data` must be a data frame."),
      call. = FALSE
    )
  }
  if (nrow(data) == 0L) {
    stop(
      paste0(.err_prefix(), "`data` must contain at least one row."),
      call. = FALSE
    )
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
#' @param warn_na   Whether to warn about missing values. A caller that counts
#'   missing treatment and missing score together, and states the total in one
#'   warning, sets this to `FALSE`.
#' @return Invisible NULL.
#' @keywords internal
.check_binary <- function(data, col, warn_na = TRUE) {
  x    <- data[[col]]
  vals <- x[!is.na(x)]
  ok   <- is.logical(x) ||
    (is.numeric(x) && all(vals %in% c(0, 1)))
  if (!ok) {
    stop(
      paste0(.err_prefix(), sprintf("Column `%s` must be binary (0/1 or logical).", col)),
      call. = FALSE
    )
  }
  n_na <- sum(is.na(x))
  if (warn_na && n_na > 0L) {
    warning(
      paste0(.err_prefix(), sprintf(
        "Column `%s` contains %d NA value(s). Those patients will be silently excluded from analysis.",
        col, n_na
      )),
      call. = FALSE
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
    warning(
      paste0(.err_prefix(), sprintf("%d row(s) have a missing `%s` or `%s` and are excluded from %s.",
                        n, treatment_col, score_col, what)),
      call. = FALSE
    )
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
#' @param hint      Optional sentence appended to the range error.
#' @return Invisible NULL.
#' @keywords internal
.check_probability <- function(data, col, hint = NULL) {
  x <- data[[col]]
  if (!is.numeric(x)) {
    stop(
      paste0(.err_prefix(), sprintf("Column `%s` must be numeric.", col)),
      call. = FALSE
    )
  }
  if (any(!is.na(x) & (x < 0 | x > 1))) {
    stop(
      paste0(.err_prefix(), paste(c(sprintf("Column `%s` must contain values in [0, 1].", col), hint), collapse = " ")),
      call. = FALSE
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
#' @return A numeric vector of probabilities, one per row of `x$data`.
#' @keywords internal
.probability_score <- function(x, score_col = x$meta$score_col) {
  .check_cols(x$data, score_col)
  score <- x$data[[score_col]]
  if (identical(x$meta$scale, "logit")) score <- stats::plogis(score)
  .check_probability(stats::setNames(list(score), score_col), score_col)
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
#' values take at most one value. A factor or character column whose
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
    n_values > 1L && (is.numeric(x) || is.logical(x) || n_values < length(x))
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
#' @param types      Named character vector from `.covariate_types()`. A caller
#'   that builds two tables from different rows, such as before and after a
#'   match, types the covariates once on the larger set and passes the result
#'   to both, so a column is never Gaussian in one table and binary in the
#'   other. `NULL` types them on `data`.
#' @return A data frame with columns `variable` and `smd`, rounded to 4 places,
#'   one row per covariate.
#' @keywords internal
.smd_table <- function(data, treatment, covariates = NULL, weight_col = NULL, reserved = character(),
                       keep = character(), types = NULL) {
  if (is.null(covariates)) {
    covariates <- .default_covariates(data, c(treatment, weight_col, reserved), keep)
  }
  .check_cols(data, covariates)
  type <- if (is.null(types)) .covariate_types(data, covariates) else unname(types[covariates])
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


#' Standardized-difference type of each covariate
#'
#' A logical column, or a numeric column whose non-missing values are all 0 or
#' 1, is `"binary"`; any other numeric column is `"gaussian"`; a factor or
#' character column is `"categorical"`. Any other type aborts.
#'
#' @param data A data frame.
#' @param covariates Covariate column names.
#' @return A character vector named by `covariates`.
#' @keywords internal
.covariate_types <- function(data, covariates) {
  vapply(covariates, function(v) {
    x <- data[[v]]
    if (is.logical(x) || (is.numeric(x) && all(x[!is.na(x)] %in% c(0, 1)))) {
      "binary"
    } else if (is.numeric(x)) {
      "gaussian"
    } else if (is.factor(x) || is.character(x)) {
      "categorical"
    } else {
      stop(
        paste0(.err_prefix(), sprintf(paste0("Covariate `%s` must be numeric, logical, factor or character to enter a ",
                                  "balance table; it is %s."), v, class(x)[1L])),
        call. = FALSE
      )
    }
  }, character(1L))
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
  types <- .covariate_types(data, covariates)
  group <- as.character(data[[treatment]])
  tables <- lapply(seq_len(nrow(pairs)), function(i) {
    rows <- !is.na(group) & group %in% c(pairs$level[[i]], pairs$versus[[i]])
    pair_data <- data[rows, covariates, drop = FALSE]
    flag <- .temporary_prediction_prefix(pair_data)
    pair_data[[flag]] <- as.integer(group[rows] == pairs$level[[i]])
    smd <- .smd_table(pair_data, flag, covariates, types = types)
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
#' @param warn_na Whether to warn about missing treatment values; see
#'   `.check_binary()`.
#' @return A list: `trt`, an integer 0/1 vector (`NA` where the column is
#'   missing), and `labels`, the reference then the treated value as character.
#' @keywords internal
.treatment_indicator <- function(data, col, treated_level = NULL,
                                 warn_na = TRUE) {
  x <- data[[col]]
  if (is.null(treated_level)) {
    .check_binary(data, col, warn_na = warn_na)
    values <- unique(as.integer(x[!is.na(x)]))
    if (length(values) < 2L) {
      stop(errorCondition(
        paste0(.err_prefix(), sprintf("Column `%s` must hold both treatment values; it holds %s.", col,
                if (length(values)) sprintf("only %s", .binary_group_labels(x)[values + 1L]) else "none")),
        class = "hvtiRpropensity_one_arm",
        call = NULL
      ))
    }
    return(list(trt = as.integer(x), labels = .binary_group_labels(x)))
  }
  if (length(treated_level) != 1L || is.na(treated_level)) {
    stop(
      paste0(.err_prefix(), "`treated_level` must be a single non-missing value."),
      call. = FALSE
    )
  }
  treated_level <- as.character(treated_level)
  values <- unique(as.character(x[!is.na(x)]))
  if (length(values) != 2L) {
    stop(
      paste0(.err_prefix(), sprintf("Column `%s` must hold two treatment values; it holds %d.", col, length(values))),
      call. = FALSE
    )
  }
  if (!treated_level %in% values) {
    stop(
      paste0(.err_prefix(), sprintf("`treated_level` (\"%s\") is not a value of column `%s`.", treated_level, col)),
      call. = FALSE
    )
  }
  n_na <- sum(is.na(x))
  if (warn_na && n_na > 0L) {
    warning(
      paste0(.err_prefix(), sprintf(
        "Column `%s` contains %d NA value(s). Those patients will be silently excluded from analysis.",
        col, n_na
      )),
      call. = FALSE
    )
  }
  list(trt = as.integer(as.character(x) == treated_level),
       labels = c(setdiff(values, treated_level), treated_level))
}


# ---------------------------------------------------------------------------
# Error prefix
# ---------------------------------------------------------------------------

#' Name the exported function a condition comes from
#'
#' Errors are raised with `stop(call. = FALSE)`, so the message carries the
#' function name itself, as `fn(): `. An internal helper does not know which
#' exported function called it, so it asks: this returns `"fn(): "` for the
#' outermost function on the call stack that this package exports, the one
#' the user called, or `""` when there is none (a helper called directly in a
#' test).
#'
#' @return A character string.
#' @keywords internal
.err_prefix <- function() {
  exports <- getNamespaceExports("hvtiRpropensity")
  for (call in sys.calls()) {
    fn <- call[[1L]]
    if (is.call(fn) && length(fn) == 3L && as.character(fn[[1L]]) %in% c("::", ":::")) fn <- fn[[3L]]
    if (is.name(fn) && as.character(fn) %in% exports) return(paste0(as.character(fn), "(): "))
  }
  ""
}


# ---------------------------------------------------------------------------
# Deprecations
# ---------------------------------------------------------------------------

# Deprecations already announced this session, keyed "fn:what".
.deprecation_seen <- new.env(parent = emptyenv())

#' Warn once per session about a deprecated usage
#'
#' The hvtiR family deprecates a renamed argument, or a removed default, for one
#' release: the old usage keeps working and warns once per session, with class
#' `hvtiRpropensity_deprecated`, then becomes an error in the next release.
#'
#' @param fn The exported function, for the message.
#' @param what A short key for the deprecated usage, such as `"n_rep"`.
#' @param message What changed and what to do instead.
#' @return Invisible `NULL`.
#' @keywords internal
.deprecate <- function(fn, what, message) {
  key <- paste(fn, what, sep = ":")
  if (!isTRUE(.deprecation_seen[[key]])) {
    assign(key, TRUE, envir = .deprecation_seen)
    warning(warningCondition(
      paste0(fn, "(): ", message, " This will be an error in the next release."),
      class = c("hvtiRpropensity_deprecated", "deprecatedWarning"),
      call = NULL
    ))
  }
  invisible(NULL)
}


#' The study default `treatment_col` used to have, with a deprecation warning
#'
#' `treatment_col` defaulted to `"tavr"`, a study's column name. A
#' study-specific value is never a default in the hvtiR family, so the argument
#' is now required; for one release an omitted `treatment_col` still means
#' `"tavr"`, and warns.
#'
#' @param fn The exported function, for the message.
#' @return `"tavr"`.
#' @keywords internal
.default_treatment_col <- function(fn) {
  .deprecate(fn, "treatment_col", paste(
    "`treatment_col` has no default now; `\"tavr\"` was one study's column name.",
    "Name the treatment column, for example `treatment_col = \"tavr\"`."
  ))
  "tavr"
}


#' Accept a renamed argument for one release
#'
#' @param fn The exported function, for the message.
#' @param old,new The old and new argument names.
#' @param value The value passed under the old name.
#' @param new_supplied Whether the new name was supplied too; both abort.
#' @return `value`.
#' @keywords internal
.renamed_arg <- function(fn, old, new, value, new_supplied) {
  if (new_supplied) {
    stop(sprintf("%s(): give `%s` only; `%s` is its old name.", fn, new, old), call. = FALSE)
  }
  .deprecate(fn, old, sprintf("`%s` is renamed `%s`.", old, new))
  value
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
#' @return `NULL`, or `abs(seed)` as an integer.
#' @keywords internal
.seed_value <- function(seed) {
  if (is.null(seed)) return(NULL)
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) || seed != round(seed) ||
        abs(seed) > .Machine$integer.max) {
    stop(
      paste0(.err_prefix(), "`seed` must be `NULL` or one whole number."),
      call. = FALSE
    )
  }
  as.integer(abs(seed))
}
