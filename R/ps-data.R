###############################################################################
## ps-data.R
##
## Base class infrastructure for hvtiRpropensity data objects.
##
## Every public function in this package returns an object of a specific
## subclass (e.g. "ps_match", "ps_weight") that also inherits from the
## common base class "ps_data".  This file defines:
##
##   new_ps_data()      -- internal constructor used by every ps_*() function
##   print.ps_data()    -- default print for any ps_data subclass
##   summary.ps_data()  -- default summary for any ps_data subclass
##   is_ps_data()       -- lightweight predicate for user code and tests
##
## Object structure (guaranteed for every ps_data subclass):
##
##   $data      <data.frame>  Original data with propensity scores / weights
##                            appended (columns vary by subclass).
##   $meta      <list>        Column names used, method parameters, formula,
##                            computed statistics.
##   $tables    <list>        Diagnostic tables: SMD before/after, group
##                            counts, effective N.  May be empty list.
##   $models    <list>        Fitted model objects. May be empty for methods
##                            that do not fit a model.
##
## Individual subclasses may add further named elements beyond these four,
## but $data, $meta, $tables, and $models are always present.
##
###############################################################################


# ---------------------------------------------------------------------------
# Internal constructor
# ---------------------------------------------------------------------------

#' Construct a validated ps_data object
#'
#' Internal entry point used by every `ps_*()` function to create its return
#' value.  Enforces the four-slot contract (`$data`, `$meta`, `$tables`,
#' `$models`) and
#' attaches the two-level S3 class vector.
#'
#' @param data     A data frame — the original data with propensity scores or
#'   weights appended.
#' @param meta     A named list of metadata (column names, formula, method
#'   parameters, computed statistics, etc.).
#' @param tables   A named list of diagnostic objects (SMD tables, group
#'   counts, effective N, etc.).  May be `list()`.
#' @param models   A named list of fitted model objects. May be `list()`.
#' @param subclass A single string naming the specific subclass
#'   (e.g. `"ps_match"`).
#'
#' @return A named list of class `c(subclass, "ps_data")`.
#' @keywords internal
new_ps_data <- function(data, meta, tables = list(), models = list(), subclass) {
  stopifnot(is.data.frame(data))
  stopifnot(is.list(meta))
  stopifnot(is.list(tables))
  stopifnot(is.list(models))
  stopifnot(is.character(subclass), length(subclass) == 1L, nzchar(subclass))

  structure(
    list(data = data, meta = meta, tables = tables, models = models),
    class = c(subclass, "ps_data")
  )
}


# ---------------------------------------------------------------------------
# Predicate
# ---------------------------------------------------------------------------

#' Test whether an object is a ps_data instance
#'
#' @param x Any R object.
#' @return `TRUE` if `x` inherits from `"ps_data"`, `FALSE` otherwise.
#'
#' @examples
#' dta <- sample_ps_data()
#' obj <- ps_match(dta)
#' is_ps_data(obj)   # TRUE
#' is_ps_data(list()) # FALSE
#'
#' @export
is_ps_data <- function(x) inherits(x, "ps_data")


# ---------------------------------------------------------------------------
# print
# ---------------------------------------------------------------------------

#' Print a ps_data object
#'
#' Displays a compact summary of the object class, group counts, and the
#' available diagnostic tables.
#'
#' @param x   A `ps_data` object.
#' @param ... Ignored.
#' @return `x`, invisibly.
#'
#' @export
print.ps_data <- function(x, ...) {
  cls <- class(x)[1L]
  cat(sprintf("<%s>\n", cls))

  if (!is.null(x$meta$n_total))
    cat(sprintf("  N total     : %d\n", x$meta$n_total))
  if (!is.null(x$meta$treatment_col))
    cat(sprintf("  Treatment   : %s\n", x$meta$treatment_col))
  if (!is.null(x$meta$score_col))
    cat(sprintf("  PS column   : %s\n",  x$meta$score_col))
  if (!is.null(x$meta$method))
    cat(sprintf("  Method      : %s\n",  x$meta$method))

  if (length(x$tables) > 0L) {
    cat("  Tables      :", paste(names(x$tables), collapse = ", "), "\n")
  }

  invisible(x)
}


# ---------------------------------------------------------------------------
# summary
# ---------------------------------------------------------------------------

#' Summarise a ps_data object
#'
#' Returns (and prints) a list of the available diagnostic tables.
#'
#' @param object A `ps_data` object.
#' @param ...    Ignored.
#' @return `object$tables`, invisibly.
#'
#' @export
summary.ps_data <- function(object, ...) {
  cat(sprintf("Summary of <%s>\n\n", class(object)[1L]))

  tbls <- object$tables
  if (length(tbls) == 0L) {
    cat("  (no diagnostic tables)\n")
  } else {
    for (nm in names(tbls)) {
      tbl <- tbls[[nm]]
      if (!is.null(tbl)) {
        # Convert underscores to spaces for a readable label
        label <- gsub("_", " ", nm, fixed = TRUE)
        label <- paste0(toupper(substring(label, 1L, 1L)),
                        substring(label, 2L))
        cat(label, ":\n", sep = "")
        print(tbl)
        cat("\n")
      }
    }
  }

  invisible(tbls)
}
