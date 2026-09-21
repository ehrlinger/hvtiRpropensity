###############################################################################
## lm-validate.R
##
## Validation of saved binary logistic model bundles without refitting.
###############################################################################


#' Validate a saved binary logistic model
#'
#' Scores a declared validation cohort with every fitted model retained in a
#' binary [fit_logistic()] bundle. It does not refit or modify the source
#' bundle.
#'
#' @param model A binary `lm_fit` bundle with bundle version 1.
#' @param data Validation data containing the stored predictors and outcome.
#' @param outcome_col Outcome column. Defaults to the bundle's outcome column.
#' @param prediction_col Output column for averaged predicted probabilities.
#' @param groups Number of equal-frequency calibration groups.
#'
#' @return An object of class `c("lm_validation", "ps_data")`.
#' @export
validate_logistic <- function(model, data, outcome_col = NULL,
                              prediction_col = "predicted", groups = 10L) {
  if (!inherits(model, "lm_fit")) {
    rlang::abort("`model` must be an `lm_fit` object.", call. = FALSE)
  }
  if (!identical(model$meta$bundle_version, 1L)) {
    rlang::abort("The `model` bundle version must be 1.", call. = FALSE)
  }
  if (!identical(model$meta$model_family, "binary")) {
    rlang::abort("`validate_logistic()` accepts only a binary model bundle.",
                 call. = FALSE)
  }
  .check_df(data)
  if (!is.numeric(groups) || length(groups) != 1L || is.na(groups) ||
        groups < 2 || groups != as.integer(groups)) {
    rlang::abort("`groups` must be one integer greater than or equal to 2.",
                 call. = FALSE)
  }
  groups <- as.integer(groups)
  if (!is.character(prediction_col) || length(prediction_col) != 1L ||
        is.na(prediction_col) || !nzchar(prediction_col)) {
    rlang::abort("`prediction_col` must be one non-empty string.",
                 call. = FALSE)
  }
  if (prediction_col %in% names(data)) {
    rlang::abort(sprintf("Column `%s` already exists in `data`.", prediction_col),
                 call. = FALSE)
  }
  if (is.null(outcome_col)) outcome_col <- model$meta$outcome_col
  if (!is.character(outcome_col) || length(outcome_col) != 1L ||
        is.na(outcome_col) || !nzchar(outcome_col)) {
    rlang::abort("`outcome_col` must be one non-empty column name.",
                 call. = FALSE)
  }
  .check_cols(data, outcome_col)
  missing_predictors <- setdiff(model$meta$predictors, names(data))
  if (length(missing_predictors)) {
    rlang::abort(
      sprintf("Validation data is missing required predictor(s): %s.",
              paste(missing_predictors, collapse = ", ")),
      call. = FALSE
    )
  }
  .check_outcome_levels(
    data, outcome_col, model$meta$outcome_levels, imputation_col = NULL
  )
  if (!is.list(model$models) || !length(model$models)) {
    rlang::abort("The model bundle does not contain fitted models.",
                 call. = FALSE)
  }

  predictions <- lapply(seq_along(model$models), function(i) {
    out <- tryCatch(
      stats::predict(model$models[[i]], newdata = data, type = "response"),
      error = function(error) {
        rlang::abort(
          sprintf("Stored model %s could not score the validation data: %s",
                  names(model$models)[[i]], conditionMessage(error)),
          call. = FALSE
        )
      }
    )
    if (!is.numeric(out) || length(out) != nrow(data)) {
      rlang::abort(
        sprintf("Stored model %s did not return one numeric prediction per row.",
                names(model$models)[[i]]),
        call. = FALSE
      )
    }
    out
  })
  missing_pattern <- is.na(predictions[[1L]])
  if (any(vapply(predictions[-1L], function(x) {
    !identical(is.na(x), missing_pattern)
  }, logical(1L)))) {
    rlang::abort("Stored models disagree on which validation rows can be scored.",
                 call. = FALSE)
  }
  predicted <- rowMeans(do.call(cbind, predictions))
  if (any(!is.na(predicted) & !is.finite(predicted))) {
    rlang::abort("The stored models returned non-finite predictions.",
                 call. = FALSE)
  }

  event <- model$meta$event_level
  observed <- as.integer(as.character(data[[outcome_col]]) == event)
  analyzed <- !is.na(observed) & !is.na(predicted)
  n_analyzed <- sum(analyzed)
  if (groups > n_analyzed) {
    rlang::abort("`groups` cannot exceed the number of analyzed observations.",
                 call. = FALSE)
  }
  y <- observed[analyzed]
  p <- predicted[analyzed]
  if (length(unique(y)) != 2L) {
    rlang::abort("The analyzed validation cohort must contain both outcome levels.",
                 call. = FALSE)
  }

  analyzed_rows <- which(analyzed)
  ranked_rows <- analyzed_rows[order(predicted[analyzed_rows], analyzed_rows)]
  calibration_group <- rep.int(NA_integer_, nrow(data))
  calibration_group[ranked_rows] <- pmin(
    ceiling(seq_along(ranked_rows) * groups / length(ranked_rows)),
    groups
  )
  calibration <- do.call(rbind, lapply(seq_len(groups), function(group) {
    rows <- calibration_group == group & !is.na(calibration_group)
    data.frame(
      group = group,
      n = sum(rows),
      observed = sum(observed[rows]),
      expected = sum(predicted[rows])
    )
  }))

  n_event <- sum(y == 1L)
  n_nonevent <- sum(y == 0L)
  ranks <- rank(p, ties.method = "average")
  auc <- (sum(ranks[y == 1L]) - n_event * (n_event + 1) / 2) /
    (n_event * n_nonevent)
  expected_events <- sum(p)
  performance <- data.frame(
    n = n_analyzed,
    observed = n_event,
    expected = expected_events,
    oe_ratio = n_event / expected_events,
    auc = auc,
    brier = mean((y - p)^2)
  )

  scored <- data
  scored[[prediction_col]] <- unname(predicted)
  meta <- list(
    bundle_version = model$meta$bundle_version,
    model_family = model$meta$model_family,
    formula = model$meta$formula,
    outcome_col = outcome_col,
    outcome_levels = model$meta$outcome_levels,
    event_level = event,
    predictors = model$meta$predictors,
    prediction_col = prediction_col,
    groups = groups,
    n_input = nrow(data),
    n_analyzed = n_analyzed,
    n_excluded = nrow(data) - n_analyzed,
    package_versions = model$meta$package_versions
  )
  new_ps_data(
    scored,
    meta,
    tables = list(calibration = calibration, performance = performance),
    models = model$models,
    subclass = "lm_validation"
  )
}
