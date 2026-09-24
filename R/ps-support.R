#' Compare support (overlap) rules
#'
#' Flags patients with weak support under several rules and measures how well
#' the rules agree. Built-in rules use the propensity score of a `ps_data`
#' object; any other rule (for example an isolation-forest tail flag) can be
#' supplied as a logical vector.
#'
#' @param x A `ps_data` object with a score and treatment column.
#' @param flags Optional named list of extra logical vectors, one element per
#'   row of `x$data`, `TRUE` meaning the rule would drop that patient.
#' @param trim Length-2 numeric, or `NULL`. Adds the rule "score outside
#'   `[trim[1], trim[2]]`". Default `c(0.1, 0.9)`.
#' @param common Logical; add the rule "score outside the range occupied by
#'   both groups" (common support).
#'
#' @return A `ps_support` / `ps_data` object. `$data` is `x$data` plus one
#'   logical column `weak_<rule>` per rule. `$tables$counts` gives the number
#'   flagged per rule; `$tables$agreement` gives, for every pair of rules,
#'   `n_both`, `n_neither`, `pct_agree`, Cohen `kappa` and `jaccard` (overlap
#'   of the two dropped sets). `$meta$common_support` is the interval used.
#'
#' @examples
#' dta <- sample_ps_data(n = 200, seed = 3)[, c("id", "tavr", "age", "ef")]
#' sup <- ps_support(ps_logistic(tavr ~ age + ef, dta))
#' sup$tables$agreement
#' @export
ps_support <- function(x, flags = list(), trim = c(0.1, 0.9), common = TRUE) {
  if (!is_ps_data(x)) rlang::abort("`x` must be a ps_data object.", call. = FALSE)
  score <- x$data[[x$meta$score_col]]
  trt <- as.integer(as.character(x$data[[x$meta$treatment_col]]) == as.character(x$meta$treated_level))
  n <- nrow(x$data)

  rules <- list()
  bounds <- NULL
  if (isTRUE(common)) {
    bounds <- c(max(min(score[trt == 1L]), min(score[trt == 0L])),
                min(max(score[trt == 1L]), max(score[trt == 0L])))
    rules$common <- score < bounds[1L] | score > bounds[2L]
  }
  if (!is.null(trim)) {
    if (!is.numeric(trim) || length(trim) != 2L || trim[1L] >= trim[2L]) {
      rlang::abort("`trim` must be two increasing numbers.", call. = FALSE)
    }
    rules$trim <- score < trim[1L] | score > trim[2L]
  }
  if (length(flags)) {
    if (is.null(names(flags)) || any(!nzchar(names(flags)))) {
      rlang::abort("`flags` must be a named list.", call. = FALSE)
    }
    if (any(lengths(flags) != n)) {
      rlang::abort("Every element of `flags` must have one value per row of `x$data`.", call. = FALSE)
    }
    rules <- c(rules, lapply(flags, as.logical))
  }
  if (length(rules) < 1L) rlang::abort("No support rules requested.", call. = FALSE)

  out <- x$data
  for (nm in names(rules)) out[[paste0("weak_", nm)]] <- rules[[nm]]

  pairs <- utils::combn(names(rules), 2L, simplify = FALSE)
  agreement <- if (length(pairs)) {
    do.call(rbind, lapply(pairs, function(p) {
      cbind(rule_a = p[1L], rule_b = p[2L], .support_agreement(rules[[p[1L]]], rules[[p[2L]]]),
            stringsAsFactors = FALSE)
    }))
  } else {
    data.frame()
  }
  counts <- data.frame(rule = names(rules), n_weak = vapply(rules, sum, numeric(1L), na.rm = TRUE),
                       n = n, row.names = NULL)

  new_ps_data(
    data = out,
    meta = c(x$meta[c("treatment_col", "score_col", "treated_level")],
             list(rules = names(rules), common_support = bounds, trim = trim, n_total = n)),
    tables = list(counts = counts, agreement = agreement),
    subclass = "ps_support"
  )
}

# Cohen kappa + Jaccard for two logical "weak support" flags.
.support_agreement <- function(flag_a, flag_b) {
  ok <- !is.na(flag_a) & !is.na(flag_b)
  a <- flag_a[ok]
  b <- flag_b[ok]
  n <- length(a)
  n11 <- sum(a & b)
  n00 <- sum(!a & !b)
  po <- (n11 + n00) / n
  pe <- mean(a) * mean(b) + mean(!a) * mean(!b)
  union_n <- sum(a | b)
  data.frame(
    n = n, n_a = sum(a), n_b = sum(b), n_both = n11, n_neither = n00,
    pct_agree = 100 * po,
    kappa = if (isTRUE(all.equal(pe, 1))) NA_real_ else (po - pe) / (1 - pe),
    jaccard = if (union_n == 0L) NA_real_ else n11 / union_n
  )
}

#' @export
print.ps_support <- function(x, ...) {
  cat("<ps_support>\n")
  cat(sprintf("  N total : %d\n", x$meta$n_total))
  cat("  Rules   :", paste(x$meta$rules, collapse = ", "), "\n")
  print(x$tables$counts, row.names = FALSE)
  invisible(x)
}
