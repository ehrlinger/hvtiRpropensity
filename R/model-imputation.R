###############################################################################
## model-imputation.R
##
## Strict model fitting over one dataset or stacked imputations.
###############################################################################


#' Report whether a fitted model converged
#'
#' @param fit A fitted model object.
#' @return A single logical value.
#' @keywords internal
.model_converged <- function(fit) {
  if (!is.null(fit$converged)) return(isTRUE(fit$converged))
  if (!is.null(fit$convergence)) return(isTRUE(fit$convergence == 0L))
  TRUE
}


#' Count rows used by a fitted model
#'
#' @param fit A fitted model object.
#' @param fallback Count to return when the model has no `nobs()` method.
#' @return A single integer.
#' @keywords internal
.model_nobs <- function(fit, fallback) {
  out <- tryCatch(stats::nobs(fit), error = function(e) fallback)
  as.integer(out)
}


#' Validate predictions against their input partition
#'
#' @param predictions Numeric vector or matrix.
#' @param n Expected number of rows.
#' @param imputation Imputation label used in errors.
#' @param columns Expected matrix column names, or `NULL` for the first fit.
#' @return Matrix column names, or `NULL` for vector predictions.
#' @keywords internal
.check_model_predictions <- function(predictions, n, imputation, columns = NULL) {
  actual <- if (is.matrix(predictions)) nrow(predictions) else length(predictions)
  if (actual != n) {
    rlang::abort(
      sprintf("Imputation %s returned %d predictions for %d rows.",
              imputation, actual, n),
      call. = FALSE
    )
  }
  if (!is.numeric(predictions)) {
    rlang::abort(
      sprintf("Imputation %s predictions must be numeric.", imputation),
      call. = FALSE
    )
  }
  if (!is.matrix(predictions)) return(NULL)

  current <- colnames(predictions)
  if (is.null(current) || anyNA(current) || any(!nzchar(current))) {
    rlang::abort(
      sprintf("Imputation %s prediction columns must be named.", imputation),
      call. = FALSE
    )
  }
  if (!is.null(columns) && !identical(current, columns)) {
    rlang::abort(
      sprintf("Imputation %s prediction columns differ from the first imputation.",
              imputation),
      call. = FALSE
    )
  }
  current
}


#' Compare response values, treating paired missing values as equal
#'
#' @param x Current response vector.
#' @param reference Reference response vector.
#' @return Logical vector.
#' @keywords internal
.same_response <- function(x, reference) {
  same <- x == reference
  paired_na <- is.na(x) & is.na(reference)
  same[is.na(same)] <- paired_na[is.na(same)]
  same
}


#' Fit one model per imputation with strict patient alignment
#'
#' @param data A complete data frame or stacked long-form imputations.
#' @param response_col Response column whose patient values must remain stable.
#' @param id_col Patient identifier column.
#' @param imputation_col Imputation-index column, or `NULL` for one dataset.
#' @param fit_fn Function accepting one data-frame partition.
#' @param predict_fn Function accepting a fit and its aligned partition.
#' @param predictor_cols Predictor columns whose categorical levels must agree
#'   across imputations.
#' @param require_stable_response Whether each patient's response must agree
#'   across imputations.
#' @return A list containing base data, averaged predictions, fitted models,
#'   fit status and imputation labels.
#' @keywords internal
.fit_imputations <- function(data, response_col, id_col, imputation_col,
                             fit_fn, predict_fn, predictor_cols = NULL,
                             require_stable_response = TRUE) {
  .check_df(data)
  .check_cols(data, response_col)
  if (!is.null(id_col)) .check_cols(data, id_col)
  if (length(predictor_cols)) .check_cols(data, predictor_cols)

  if (is.null(imputation_col)) {
    fit <- fit_fn(data)
    predictions <- predict_fn(fit, data)
    .check_model_predictions(predictions, nrow(data), 1L)
    n_analyzed <- .model_nobs(fit, nrow(data))
    return(list(
      data = data,
      predictions = predictions,
      models = list("1" = fit),
      status = data.frame(
        imputation = 1L,
        converged = .model_converged(fit),
        n_input = nrow(data),
        n_analyzed = n_analyzed,
        n_excluded = nrow(data) - n_analyzed
      ),
      imputations = 1L
    ))
  }

  if (is.null(id_col)) {
    rlang::abort("`id_col` is required when `imputation_col` is supplied.",
                 call. = FALSE)
  }
  .check_cols(data, imputation_col)
  if (anyNA(data[[imputation_col]])) {
    rlang::abort("`imputation_col` must not contain missing values.",
                 call. = FALSE)
  }
  imputations <- sort(unique(data[[imputation_col]]))
  if (length(imputations) < 2L) {
    rlang::abort(
      sprintf("Column `%s` must contain at least 2 distinct imputation indices.",
              imputation_col),
      call. = FALSE
    )
  }

  partitions <- lapply(imputations, function(imp) {
    data[data[[imputation_col]] == imp, , drop = FALSE]
  })
  for (i in seq_along(partitions)) {
    if (anyNA(partitions[[i]][[id_col]])) {
      rlang::abort(
        sprintf("Imputation %s contains a missing patient key.", imputations[[i]]),
        call. = FALSE
      )
    }
    duplicated_ids <- duplicated(partitions[[i]][[id_col]])
    if (any(duplicated_ids)) {
      rlang::abort(
        sprintf("Imputation %s contains %d duplicate patient key(s).",
                imputations[[i]], sum(duplicated_ids)),
        call. = FALSE
      )
    }
  }

  reference_ids <- partitions[[1L]][[id_col]]
  reference_response <- partitions[[1L]][[response_col]]
  for (i in seq_along(partitions)) {
    ids <- partitions[[i]][[id_col]]
    missing_ids <- setdiff(reference_ids, ids)
    extra_ids <- setdiff(ids, reference_ids)
    if (length(missing_ids) || length(extra_ids)) {
      detail <- c(
        if (length(missing_ids)) sprintf("missing %d patient(s)", length(missing_ids)),
        if (length(extra_ids)) sprintf("has %d extra patient(s)", length(extra_ids))
      )
      rlang::abort(
        sprintf("Imputation %s %s relative to imputation %s.",
                imputations[[i]], paste(detail, collapse = " and "),
                imputations[[1L]]),
        call. = FALSE
      )
    }
    partitions[[i]] <- partitions[[i]][match(reference_ids, ids), , drop = FALSE]
    same <- .same_response(partitions[[i]][[response_col]], reference_response)
    if (require_stable_response && any(!same)) {
      rlang::abort(
        sprintf("Imputation %s response differs for %d patient(s) from imputation %s.",
                imputations[[i]], sum(!same), imputations[[1L]]),
        call. = FALSE
      )
    }
  }

  categorical_predictors <- predictor_cols[vapply(
    data[predictor_cols],
    function(x) is.factor(x) || is.character(x) || is.logical(x),
    logical(1L)
  )]
  for (predictor in categorical_predictors) {
    reference_levels <- sort(unique(as.character(
      partitions[[1L]][[predictor]][!is.na(partitions[[1L]][[predictor]])]
    )))
    for (i in seq_along(partitions)[-1L]) {
      current_levels <- sort(unique(as.character(
        partitions[[i]][[predictor]][!is.na(partitions[[i]][[predictor]])]
      )))
      if (!identical(current_levels, reference_levels)) {
        rlang::abort(
          sprintf("Predictor `%s` has different observed levels in imputation %s.",
                  predictor, imputations[[i]]),
          call. = FALSE
        )
      }
    }
  }

  models <- vector("list", length(imputations))
  predictions <- vector("list", length(imputations))
  status <- vector("list", length(imputations))
  prediction_columns <- NULL
  for (i in seq_along(imputations)) {
    partition <- partitions[[i]]
    fit <- fit_fn(partition)
    predicted <- predict_fn(fit, partition)
    current_columns <- .check_model_predictions(
      predicted, nrow(partition), imputations[[i]], prediction_columns
    )
    if (i == 1L) prediction_columns <- current_columns
    n_analyzed <- .model_nobs(fit, nrow(partition))
    models[[i]] <- fit
    predictions[[i]] <- predicted
    status[[i]] <- data.frame(
      imputation = imputations[[i]],
      converged = .model_converged(fit),
      n_input = nrow(partition),
      n_analyzed = n_analyzed,
      n_excluded = nrow(partition) - n_analyzed
    )
  }
  names(models) <- as.character(imputations)

  averaged <- if (is.matrix(predictions[[1L]])) {
    Reduce(`+`, predictions) / length(predictions)
  } else {
    rowMeans(do.call(cbind, predictions))
  }
  base_data <- partitions[[1L]]
  base_data[[imputation_col]] <- NULL
  rownames(base_data) <- NULL

  list(
    data = base_data,
    predictions = averaged,
    models = models,
    status = do.call(rbind, status),
    imputations = imputations
  )
}
