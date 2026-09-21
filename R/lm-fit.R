###############################################################################
## lm-fit.R
##
## General logistic model bundles.
###############################################################################


#' Validate a declared outcome-level contract
#'
#' @param data Model data.
#' @param outcome_col Outcome column.
#' @param outcome_levels Complete declared levels as character strings.
#' @param imputation_col Imputation column, or `NULL`.
#' @return Invisible `NULL`.
#' @keywords internal
.check_outcome_levels <- function(data, outcome_col, outcome_levels,
                                  imputation_col = NULL) {
  check_one <- function(values, label = NULL) {
    observed <- unique(as.character(values[!is.na(values)]))
    undeclared <- setdiff(observed, outcome_levels)
    absent <- setdiff(outcome_levels, observed)
    where <- if (is.null(label)) "" else sprintf(" in imputation %s", label)
    if (length(undeclared)) {
      rlang::abort(
        sprintf("The outcome contains undeclared level(s)%s: %s.",
                where, paste(undeclared, collapse = ", ")),
        call. = FALSE
      )
    }
    if (length(absent)) {
      rlang::abort(
        sprintf("The declared outcome level(s) %s are absent%s.",
                paste(absent, collapse = ", "), where),
        call. = FALSE
      )
    }
  }

  if (is.null(imputation_col)) {
    check_one(data[[outcome_col]])
  } else {
    for (imputation in sort(unique(data[[imputation_col]]))) {
      rows <- data[[imputation_col]] == imputation
      check_one(data[[outcome_col]][rows], imputation)
    }
  }
  invisible(NULL)
}


#' Fit a general logistic model bundle
#'
#' Fits binary, proportional-odds ordinal, or generalized-logit nominal models
#' under an explicit outcome-level contract. Stacked long-form imputations are
#' fitted separately and combined with Rubin's rules.
#'
#' @param formula Model formula.
#' @param data A data frame containing the outcome and predictors.
#' @param family Model family.
#' @param outcome_col Outcome column. By default it is read from `formula`.
#' @param id_col Patient identifier column.
#' @param imputation_col Imputation-index column, or `NULL` for one dataset.
#' @param outcome_levels Complete outcome levels in their declared order.
#' @param event_level Event level for a binary model.
#' @param reference_level Reference level for a nominal model.
#' @param prediction_prefix Name for a binary probability column, or prefix for
#'   ordinal and nominal level-specific probability columns.
#' @param trace Whether model engines should print fitting progress.
#'
#' @return An object of class `c("lm_fit", "ps_data")`.
#' @export
fit_logistic <- function(formula, data,
                         family = c("binary", "ordinal", "nominal"),
                         outcome_col = NULL, id_col = "id",
                         imputation_col = NULL, outcome_levels,
                         event_level = NULL, reference_level = NULL,
                         prediction_prefix = "prob", trace = FALSE) {
  family <- match.arg(family)
  .check_df(data)
  if (!inherits(formula, "formula")) {
    rlang::abort("`formula` must be a model formula.", call. = FALSE)
  }
  response_vars <- all.vars(formula[[2L]])
  if (length(response_vars) != 1L) {
    rlang::abort("`formula` must contain one untransformed outcome column.",
                 call. = FALSE)
  }
  formula_outcome <- response_vars[[1L]]
  if (is.null(outcome_col)) outcome_col <- formula_outcome
  if (!is.character(outcome_col) || length(outcome_col) != 1L ||
        !nzchar(outcome_col) || !identical(outcome_col, formula_outcome)) {
    rlang::abort("`outcome_col` must name the outcome on the left side of `formula`.",
                 call. = FALSE)
  }
  .check_cols(data, all.vars(formula))
  if (!is.null(id_col)) .check_cols(data, id_col)
  if (!is.null(imputation_col)) .check_cols(data, imputation_col)
  if (!is.character(prediction_prefix) || length(prediction_prefix) != 1L ||
        is.na(prediction_prefix) || !nzchar(prediction_prefix)) {
    rlang::abort("`prediction_prefix` must be one non-empty string.",
                 call. = FALSE)
  }
  declared_levels <- as.character(outcome_levels)
  minimum_levels <- 2L
  valid_length <- if (identical(family, "binary")) {
    length(declared_levels) == 2L
  } else {
    length(declared_levels) >= minimum_levels
  }
  if (!valid_length || anyNA(outcome_levels) ||
        any(!nzchar(declared_levels)) || anyDuplicated(declared_levels)) {
    requirement <- if (identical(family, "binary")) {
      "exactly two"
    } else {
      "at least two"
    }
    rlang::abort(
      sprintf("%s `outcome_levels` must contain %s unique non-missing values.",
              family, requirement),
      call. = FALSE
    )
  }
  declared_event <- NULL
  declared_reference <- NULL
  if (identical(family, "binary")) {
    declared_event <- as.character(event_level)
    if (length(declared_event) != 1L || is.na(event_level) ||
          !declared_event %in% declared_levels) {
      rlang::abort("`event_level` must be one of the declared `outcome_levels`.",
                   call. = FALSE)
    }
    model_levels <- c(setdiff(declared_levels, declared_event), declared_event)
  } else if (identical(family, "ordinal")) {
    model_levels <- declared_levels
  } else {
    declared_reference <- as.character(reference_level)
    if (length(declared_reference) != 1L || is.na(reference_level) ||
          !declared_reference %in% declared_levels) {
      rlang::abort(
        "`reference_level` must be one of the declared `outcome_levels`.",
        call. = FALSE
      )
    }
    model_levels <- c(
      declared_reference,
      setdiff(declared_levels, declared_reference)
    )
  }
  two_level <- length(model_levels) == 2L
  prediction_names <- if (identical(family, "binary")) {
    prediction_prefix
  } else {
    paste0(prediction_prefix, "_", model_levels)
  }
  existing_predictions <- intersect(prediction_names, names(data))
  if (length(existing_predictions)) {
    label <- if (length(existing_predictions) == 1L) "Column" else "Columns"
    rlang::abort(
      sprintf("%s %s already exist%s in `data`.",
              label,
              paste0("`", existing_predictions, "`", collapse = ", "),
              if (length(existing_predictions) == 1L) "s" else ""),
      call. = FALSE
    )
  }
  .check_outcome_levels(data, outcome_col, declared_levels, imputation_col)

  fit_fn <- function(partition) {
    model_data <- partition
    ordered_outcome <- identical(family, "ordinal") && !two_level
    model_data[[outcome_col]] <- factor(
      as.character(model_data[[outcome_col]]),
      levels = model_levels,
      ordered = ordered_outcome
    )
    model <- switch(
      family,
      binary = stats::glm(
        formula,
        data = model_data,
        family = stats::binomial(),
        na.action = stats::na.exclude
      ),
      ordinal = if (two_level) {
        stats::glm(
          formula,
          data = model_data,
          family = stats::binomial(),
          na.action = stats::na.exclude
        )
      } else {
        MASS::polr(
          formula,
          data = model_data,
          Hess = TRUE,
          na.action = stats::na.exclude
        )
      },
      nominal = nnet::multinom(
        formula,
        data = model_data,
        Hess = TRUE,
        trace = trace,
        na.action = stats::na.exclude
      )
    )
    if (identical(family, "ordinal") && two_level) {
      attr(model, "hvti_outcome_levels") <- model_levels
    }
    model
  }
  predict_fn <- function(fit, partition) {
    prediction_type <- if (identical(family, "binary") ||
                             (identical(family, "ordinal") && two_level)) {
      "response"
    } else {
      "probs"
    }
    predictions <- stats::predict(fit, newdata = partition, type = prediction_type)
    if (!identical(family, "binary") && two_level && is.null(dim(predictions))) {
      predictions <- cbind(1 - predictions, predictions)
      colnames(predictions) <- model_levels
    }
    predictions
  }
  fitted <- .fit_imputations(
    data = data,
    response_col = outcome_col,
    id_col = id_col,
    imputation_col = imputation_col,
    fit_fn = fit_fn,
    predict_fn = predict_fn,
    predictor_cols = all.vars(stats::delete.response(stats::terms(formula)))
  )
  inference <- .model_inference(fitted$models, family)
  scored <- fitted$data
  if (identical(family, "binary")) {
    scored[[prediction_names]] <- unname(fitted$predictions)
  } else {
    for (i in seq_along(model_levels)) {
      scored[[prediction_names[[i]]]] <- unname(fitted$predictions[, model_levels[[i]]])
    }
    scored[[outcome_col]] <- factor(
      as.character(scored[[outcome_col]]),
      levels = model_levels,
      ordered = identical(family, "ordinal")
    )
  }
  status <- fitted$status
  model_engine <- switch(
    family,
    binary = "stats",
    ordinal = if (two_level) "stats" else "MASS",
    nominal = "nnet"
  )

  meta <- list(
    bundle_version = 1L,
    formula = formula,
    outcome_col = outcome_col,
    model_family = family,
    outcome_levels = declared_levels,
    event_level = declared_event,
    reference_level = declared_reference,
    cumulative_direction = if (identical(family, "ordinal")) {
      "P(Y <= level) = logistic(threshold - linear predictor)"
    } else {
      NULL
    },
    predictors = all.vars(stats::delete.response(stats::terms(formula))),
    id_col = id_col,
    imputation_col = imputation_col,
    n_imputations = as.integer(length(fitted$imputations)),
    n_input = as.integer(sum(status$n_input)),
    n_analyzed = as.integer(sum(status$n_analyzed)),
    n_excluded = as.integer(sum(status$n_excluded)),
    prediction_prefix = prediction_prefix,
    package_versions = stats::setNames(c(
      R = as.character(getRversion()),
      hvtiRpropensity = as.character(utils::packageVersion("hvtiRpropensity")),
      model = as.character(utils::packageVersion(model_engine))
    ), c("R", "hvtiRpropensity", model_engine))
  )
  tables <- list(
    estimates = inference$estimates,
    covariance = inference$covariance,
    by_imputation = inference$by_imputation,
    fit_status = status
  )
  new_ps_data(scored, meta, tables, fitted$models, subclass = "lm_fit")
}
