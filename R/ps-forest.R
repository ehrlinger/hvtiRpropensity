#' Random-forest propensity score (out-of-bag)
#'
#' Estimates P(treated) with a random-forest classifier
#' ([randomForestSRC::rfsrc()]) and returns the **out-of-bag** class
#' probability as the score, so no patient is scored by trees that saw them.
#' The result has the same contract as [ps_logistic()], so it feeds
#' [ps_weight()], [ps_match()], [sa_overlap()] and
#' `hvtiPlotR::hv_balance()` / `hvtiPlotR::hv_mirror_hist()` unchanged.
#'
#' @param formula Two-sided formula, `treatment ~ x1 + x2 + ...`.
#' @param data Data frame with no missing values in the formula variables.
#' @param treatment_col Treatment column name; defaults to the formula
#'   response. Binary, 0/1 or logical.
#' @param id_col Optional patient identifier column, carried through.
#' @param score_col,logit_col,weight_col Names of the appended score, logit and
#'   overlap-weight columns.
#' @param treated_level Value of `treatment_col` that is "treated" (the score
#'   is P(treated)). Defaults to `1` / `TRUE`.
#' @param ntree Number of trees.
#' @param seed Optional integer seed passed to [randomForestSRC::rfsrc()].
#' @param clip Scores are clipped to `[clip, 1 - clip]` before the logit and
#'   weight are computed. Forest probabilities can be exactly 0 or 1.
#' @param ... Further arguments for [randomForestSRC::rfsrc()].
#'
#' @return A `ps_forest` / `ps_data` object. `$data` is `data` plus
#'   `score_col`, `logit_col`, `weight_col` (overlap weights), `quintile` and
#'   `decile`; `$meta` follows [ps_logistic()] with `method = "forest-oob"`;
#'   `$tables` holds `smd` and `group_counts`; `$models$forest` is the fitted
#'   forest.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("randomForestSRC", quietly = TRUE)) {
#'   dta <- sample_ps_data(n = 150, seed = 1)[, c("id", "tavr", "age", "ef")]
#'   obj <- ps_forest(tavr ~ age + ef, dta, ntree = 100, seed = 1)
#'   summary(obj$data$prob_t)
#' }
#' }
#' @export
ps_forest <- function(formula,
                      data,
                      treatment_col = NULL,
                      id_col        = "id",
                      score_col     = "prob_t",
                      logit_col     = "logit_t",
                      weight_col    = "mt_wt",
                      treated_level = NULL,
                      ntree         = 500L,
                      seed          = NULL,
                      clip          = 1e-3,
                      ...) {
  rlang::check_installed("randomForestSRC", reason = "to fit a forest propensity score.")
  if (!inherits(formula, "formula")) {
    rlang::abort("`formula` must be an R formula.", call. = FALSE)
  }
  .check_df(data)
  if (is.null(treatment_col)) treatment_col <- as.character(formula[[2L]])
  .check_cols(data, treatment_col)
  .check_binary(data, treatment_col)
  .check_output_columns(data, c(score_col, logit_col, weight_col, "quintile", "decile"))

  vars <- all.vars(formula)
  .check_cols(data, vars)
  if (anyNA(data[, vars, drop = FALSE])) {
    rlang::abort("`data` has missing values in the formula variables; impute or drop them first.",
                 call. = FALSE)
  }

  if (is.null(treated_level)) treated_level <- if (is.logical(data[[treatment_col]])) TRUE else 1
  treated_level <- as.character(treated_level)
  trt <- as.integer(as.character(data[[treatment_col]]) == treated_level)

  fit_data <- data[, vars, drop = FALSE]
  fit_data[[treatment_col]] <- factor(trt, levels = c(0L, 1L))
  if (!is.null(seed)) seed <- -abs(as.integer(seed))   # rfsrc wants a negative seed
  args <- c(list(formula = formula, data = fit_data, ntree = ntree), if (!is.null(seed)) list(seed = seed),
            list(...))
  fit <- do.call(randomForestSRC::rfsrc, args)

  probs <- pmin(pmax(fit$predicted.oob[, "1"], clip), 1 - clip)
  out <- data
  out[[score_col]]  <- probs
  out[[logit_col]]  <- log(probs / (1 - probs))
  out[[weight_col]] <- ifelse(trt == 1L, 1 - probs, probs)
  strata <- .assign_ps_strata(probs)
  out[["quintile"]] <- strata$quintile
  out[["decile"]]   <- strata$decile

  covariates <- setdiff(vars, treatment_col)
  covariates <- covariates[vapply(data[covariates], is.numeric, logical(1L))]
  diag <- out
  diag[[treatment_col]] <- trt
  smd_tbl <- .smd_table(diag, treatment_col, covariates)

  new_ps_data(
    data = out,
    meta = list(
      formula = formula, treatment_col = treatment_col, id_col = id_col,
      score_col = score_col, logit_col = logit_col, weight_col = weight_col,
      treated_level = treated_level, treatment_levels = c("0", "1"),
      method = "forest-oob", ntree = ntree, n_total = nrow(out)
    ),
    tables = list(
      smd = smd_tbl,
      group_counts = data.frame(group = c("control", "treated"),
                                n = c(sum(trt == 0L), sum(trt == 1L)))
    ),
    models = list(forest = fit),
    subclass = "ps_forest"
  )
}

#' @export
print.ps_forest <- function(x, ...) {
  cat("<ps_forest>\n")
  cat(sprintf("  N total     : %d\n", x$meta$n_total))
  cat(sprintf("  Treatment   : %s\n", x$meta$treatment_col))
  cat(sprintf("  PS column   : %s (out-of-bag, %d trees)\n", x$meta$score_col, x$meta$ntree))
  cat("  Tables      :", paste(names(x$tables), collapse = ", "), "\n")
  invisible(x)
}
