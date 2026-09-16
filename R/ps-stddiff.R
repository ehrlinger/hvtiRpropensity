###############################################################################
## ps-stddiff.R
##
## Standardized differences for every variable type, ported from stddiff.sas
## (Amanda Artis 2019-2020, modifying Dongsheng Yang's 2012 CCF macro).
##
## Design: hvtiRtemplates dev/specs/2026-09-16-standardized-difference-design.md
##
## Formulas, group 1 minus group 0:
##   gaussian     (m1 - m0) / sqrt((v1 + v0) / 2)
##   nong_ord     the same, on pooled ranks (ties averaged, ranks unweighted)
##   binary       (p1 - p0) / sqrt((p1(1 - p1) + p0(1 - p0)) / 2)
##   categorical  sqrt((t - c)' S^-1 (t - c)), Yang and Dalton (2012)
##   one level    0
##
## Weighted means are sum(w x) / sum(w) and weighted variances are PROC MEANS'
## sum(w (x - m)^2) / (n - 1). Weighted proportions are PROC FREQ's.
###############################################################################


#' Standardized differences between two groups
#'
#' Computes the standardized difference, group 1 minus group 0, for each
#' named variable, choosing the formula by the variable's type. A port of the
#' CCF `%stddiff` SAS macro (Artis 2019, from Yang 2012), and intended to
#' reproduce its `est.stddiff` output.
#'
#' @details
#' **Variable types.** The type decides the formula, so it is given, not
#' inferred: an ordinal NYHA class stored as an integer would otherwise be
#' treated as Gaussian.
#'
#' * `gaussian`: difference in means over `sqrt((var1 + var0) / 2)`. The
#'   denominator averages the two variances; it does not pool them by sample
#'   size, so it differs from a Cohen's d when the groups are unequal.
#' * `nong_ord`: the same formula on ranks of the pooled sample, ties
#'   averaged, as `PROC RANK` does by default. Ranks are unweighted; only
#'   their means and variances use `weight_col`.
#' * `binary` and `categorical` are routed by the number of levels observed,
#'   as the macro routes its `BINARY=` and `CATG=` lists: one level gives 0,
#'   two use the binary formula on the proportion at the last level in sort
#'   order (`1` for a 0/1 variable), and more than two use the Mahalanobis
#'   form of Yang and Dalton (2012). For that form the last level in sort
#'   order is the one dropped, and `S` has diagonal
#'   `(t_i (1 - t_i) + c_i (1 - c_i)) / 2` and off-diagonal
#'   `-(t_i t_j + c_i c_j) / 2`. A categorical difference is never negative.
#'   When `S` is singular, as when the two groups share no levels, the
#'   difference is `NA` with a warning naming the variable.
#'
#' **Weights.** Means are `sum(w x) / sum(w)`; variances are
#' `sum(w (x - m)^2) / (n - 1)`, which is what `PROC MEANS` reports with a
#' `WEIGHT` statement. Weights must be positive.
#'
#' **Missing values.** Rows with a missing group are dropped. Otherwise a
#' missing value is dropped for that variable only.
#'
#' **Units.** The result is a proportion of a standard deviation, not a
#' percent. The SAS output labels the column "(%)" but does not multiply by
#' 100.
#'
#' @param data A data frame.
#' @param treatment_col Name of the group column, coded 0/1 or logical, with
#'   both groups present. Default `"tavr"`.
#' @param gaussian,nong_ord,binary,categorical Character vectors of column
#'   names, by type. See Details. At least one variable is required, and a
#'   variable may appear under one type only.
#' @param weight_col Optional name of a column of positive weights, for
#'   example matching weights from [ps_weight()].
#'
#' @return An object of class `c("ps_stddiff", "ps_data")` with:
#' \describe{
#'   \item{`$data`}{The input data frame, unchanged.}
#'   \item{`$meta`}{Named list: `treatment_col`, `weight_col`, `variables`
#'     (a list by type), `method`, `n_total`, `n_dropped`.}
#'   \item{`$tables`}{Named list with `stddiff`, a data frame with one row
#'     per variable in the order given: `variable`, `label` (the column's
#'     `"label"` attribute, or `NA`), `type` (`"gaussian"`, `"nong_ord"`,
#'     `"binary"`, `"categorical"` or `"onelevel"`) and `stddiff`, which is
#'     `NA` when the pooled variance is zero.}
#' }
#'
#' @references Yang D, Dalton JE. A unified approach to measuring the effect
#'   size between two groups using SAS. SAS Global Forum 2012, paper 335-2012.
#'
#' @seealso [ps_match()], [ps_weight()]
#'
#' @examples
#' dta <- sample_ps_data(n = 200, seed = 42)
#' obj <- ps_stddiff(dta, treatment_col = "tavr",
#'                   gaussian = c("age", "ef"), binary = c("female", "diabetes"))
#' obj$tables$stddiff
#'
#' @export
ps_stddiff <- function(data,
                       treatment_col = "tavr",
                       gaussian      = NULL,
                       nong_ord      = NULL,
                       binary        = NULL,
                       categorical   = NULL,
                       weight_col    = NULL) {
  .check_df(data)
  vars <- list(gaussian = gaussian, nong_ord = nong_ord, binary = binary, categorical = categorical)
  vars <- lapply(vars, function(v) if (is.null(v)) character(0) else as.character(v))
  all_vars <- unlist(vars, use.names = FALSE)
  if (length(all_vars) == 0L) {
    rlang::abort("Name at least one variable in `gaussian`, `nong_ord`, `binary` or `categorical`.", call = NULL)
  }
  dup <- unique(all_vars[duplicated(all_vars)])
  if (length(dup) > 0L) {
    rlang::abort(sprintf("Variable(s) named under more than one type: %s", paste(dup, collapse = ", ")), call = NULL)
  }
  .check_cols(data, c(treatment_col, all_vars, weight_col))
  for (v in c(vars$gaussian, vars$nong_ord)) {
    if (!is.numeric(data[[v]])) {
      rlang::abort(sprintf("Column `%s` must be numeric to be summarised by mean or rank.", v), call = NULL)
    }
  }

  keep <- !is.na(data[[treatment_col]])
  base <- data[keep, , drop = FALSE]
  .check_binary(base, treatment_col)
  grp <- as.integer(base[[treatment_col]])
  if (!all(c(0L, 1L) %in% grp)) {
    rlang::abort(sprintf("Column `%s` must contain both groups, 0 and 1.", treatment_col), call = NULL)
  }
  w <- if (is.null(weight_col)) rep(1, nrow(base)) else base[[weight_col]]
  if (!is.null(weight_col) && (!is.numeric(w) || anyNA(w) || any(w <= 0))) {
    rlang::abort(sprintf("Column `%s` must hold positive, non-missing weights.", weight_col), call = NULL)
  }

  one <- function(v) {
    x  <- base[[v]]
    ok <- !is.na(x)
    type <- names(vars)[vapply(vars, function(s) v %in% s, logical(1))]
    if (type %in% c("gaussian", "nong_ord")) {
      if (type == "nong_ord") x[ok] <- rank(x[ok], ties.method = "average")
      est <- .stddiff_mean(x[ok], grp[ok], w[ok])
    } else {
      lv <- sort(unique(x[ok]))
      if (length(lv) == 1L) {
        type <- "onelevel"
        est  <- 0
      } else {
        prop <- function(g) {
          i <- ok & grp == g
          vapply(lv, function(l) sum(w[i][x[i] == l]) / sum(w[i]), numeric(1))
        }
        p1 <- prop(1L)
        p0 <- prop(0L)
        if (length(lv) == 2L) {
          type <- "binary"
          est  <- .stddiff_ratio(p1[2L] - p0[2L], (p1[2L] * (1 - p1[2L]) + p0[2L] * (1 - p0[2L])) / 2)
        } else {
          type <- "categorical"
          est  <- .stddiff_categorical(p1[-length(lv)], p0[-length(lv)])
          if (is.na(est)) {
            rlang::warn(sprintf(
              "`%s`: the groups share too few levels for a categorical standardized difference; returning NA.", v))
          }
        }
      }
    }
    lab <- attr(data[[v]], "label", exact = TRUE)
    data.frame(variable = v, label = if (is.null(lab)) NA_character_ else as.character(lab)[1L],
               type = type, stddiff = est, stringsAsFactors = FALSE)
  }
  tbl <- do.call(rbind, lapply(all_vars, one))
  rownames(tbl) <- NULL

  new_ps_data(
    data   = data,
    meta   = list(treatment_col = treatment_col, weight_col = weight_col, variables = vars,
                  method = "stddiff", n_total = nrow(data), n_dropped = sum(!keep)),
    tables = list(stddiff = tbl),
    subclass = "ps_stddiff"
  )
}


#' Mean difference over sqrt of the averaged weighted variances
#' @keywords internal
#' @noRd
.stddiff_mean <- function(x, grp, w) {
  stat <- function(g) {
    xi <- x[grp == g]
    wi <- w[grp == g]
    m  <- sum(wi * xi) / sum(wi)
    c(m, if (length(xi) < 2L) NA_real_ else sum(wi * (xi - m)^2) / (length(xi) - 1L))
  }
  s1 <- stat(1L)
  s0 <- stat(0L)
  .stddiff_ratio(s1[1L] - s0[1L], (s1[2L] + s0[2L]) / 2)
}


#' Difference over the square root of a variance, NA when that variance is 0
#' @keywords internal
#' @noRd
.stddiff_ratio <- function(diff, var) {
  if (is.na(var) || var <= 0) return(NA_real_)
  diff / sqrt(var)
}


#' Yang and Dalton's Mahalanobis standardized difference
#'
#' `t` and `c` are the group 1 and group 0 proportions with the last level
#' already removed.
#' @keywords internal
#' @noRd
.stddiff_categorical <- function(t, c) {
  s_mat <- -0.5 * (outer(t, t) + outer(c, c))
  diag(s_mat) <- 0.5 * (t * (1 - t) + c * (1 - c))
  d <- t - c
  # S is singular when the groups share no levels, e.g. group 1 all "a"/"b"
  # and group 0 all "c": t sums to 1 and c is all zero. No finite difference
  # is defined then, so return NA and let the caller name the variable, rather
  # than let solve() stop the whole table. Decided 2026-09-16.
  if (anyNA(s_mat) || qr(s_mat)$rank < length(d)) return(NA_real_)
  sqrt(sum(d * solve(s_mat, d)))
}
