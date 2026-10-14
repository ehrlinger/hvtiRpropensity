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
#' fitted separately, one model per imputation. The coefficients are combined
#' with Rubin's rules; the predicted probabilities are averaged per patient.
#'
#' @details
#' **Ordinal standard errors.** An ordinal model is fitted with [MASS::polr()]
#' at a tight `optim()` tolerance (`reltol = 1e-14`), and its covariance in
#' `$tables` is the inverse of the expected (Fisher) information, the matrix
#' SAS `PROC LOGISTIC` reports. `polr()` itself, and `vcov()` on the model in
#' `$models`, invert the observed information instead. For the cumulative logit
#' the two differ at the estimate, here by a fraction of a percent in the
#' standard errors, so the bundle's ordinal standard errors and p-values match
#' SAS rather than `summary()` of the stored model. Binary and nominal models
#' use canonical links, where the two matrices are the same.
#'
#' **Multiple imputation.** With `imputation_col` set, one model is fitted to
#' each imputed dataset. The two outputs are then combined differently.
#' `$tables$estimates` and `$tables$covariance` pool the coefficients and their
#' variances with Rubin's rules, with Barnard and Rubin's degrees of freedom
#' (see "Pooled degrees of freedom" below). The probability columns in `$data` are the
#' per-patient mean, on the probability scale, of the predictions from each
#' imputation's model. Rubin's rules play no part there: averaging the
#' predictions is what Leyrat et al. (2019) call the "across" approach. For a
#' propensity score they favour the "within" approach instead; see
#' [ps_logistic()] for the recipe.
#'
#' `$data` is built from the lowest-numbered imputation: its rows, ordered as
#' they appear there, its covariate values, and no `imputation_col`.
#'
#' @section Pooled degrees of freedom:
#' With `imputation_col` set, each pooled coefficient is tested against a t
#' distribution, and the `df` column of `$tables$estimates` gives its degrees
#' of freedom. We use Barnard and Rubin's (1999) small-sample degrees of
#' freedom, computed the way `mice::pool()` computes them, so a pooled table
#' here and one from `mice` agree.
#'
#' Rubin's (1987) original degrees of freedom are
#' \eqn{\nu_{old} = (m - 1) / \lambda^2}, where \eqn{m} is the number of
#' imputations and \eqn{\lambda = (1 + 1/m) B / T} is the proportion of a
#' coefficient's total variance \eqn{T} due to missing data, through the
#' between-imputation variance \eqn{B}. (This is `lambda` in `mice`, not its
#' finite-sample fraction of missing information, `fmi`.) When little is
#' missing, \eqn{\lambda} is small and \eqn{\nu_{old}} can run into the
#' millions, larger than the sample could ever support, and the intervals come
#' out too narrow. Barnard and Rubin cap it with the complete-data degrees of
#' freedom \eqn{\nu_{com}}, the value the analysis would have had with no
#' missing data:
#' \deqn{\nu_{obs} = \frac{\nu_{com} + 1}{\nu_{com} + 3} \nu_{com} (1 - \lambda), \qquad
#'   \nu = \left(\frac{1}{\nu_{old}} + \frac{1}{\nu_{obs}}\right)^{-1}.}{
#'   nu_obs = (nu_com + 1) / (nu_com + 3) * nu_com * (1 - lambda),
#'   nu = 1 / (1 / nu_old + 1 / nu_obs).}
#' So \eqn{\nu} never exceeds \eqn{\nu_{com}}, and it approaches
#' \eqn{\nu_{old}} as \eqn{\nu_{com}} grows.
#'
#' \eqn{\nu_{com}} is \eqn{n - k}, the residual degrees of freedom of the
#' first imputation's model: \eqn{n} analysed patients minus \eqn{k} estimated
#' parameters. That holds for every family, linear, logistic, ordinal, nominal
#' and count, as it does in `mice`. [nnet::multinom()] reports no residual
#' degrees of freedom, so for a nominal model \eqn{k} is its effective degrees
#' of freedom, again as `mice` does.
#' Pooling stops with an error when \eqn{\nu_{com}} is zero or negative, as it
#' is for a saturated model, because no small-sample degrees of freedom exist
#' then.
#'
#' Two consequences are worth knowing. First, a single logistic, ordinal,
#' nominal or count fit is tested against the normal distribution (`df` is
#' `Inf`), but the same model pooled over imputations reports a finite `df`
#' and a t test. A single linear fit uses its residual degrees of freedom.
#' Second, SAS `PROC MIANALYZE` reports \eqn{\nu_{old}} unless you give it the
#' complete-data degrees of freedom with `EDF=`. To reproduce these `df`, and
#' the p-values and intervals that depend on them, run `PROC MIANALYZE` with
#' `EDF=` set to \eqn{n - k}.
#'
#' Barnard J, Rubin DB (1999). Small-sample degrees of freedom with multiple
#' imputation. *Biometrika*, 86(4), 948--955. \doi{10.1093/biomet/86.4.948}
#'
#' Rubin DB (1987). *Multiple Imputation for Nonresponse in Surveys*. Wiley.
#' \doi{10.1002/9780470316696}
#'
#' @references Leyrat C, Seaman SR, White IR, Douglas I, Smeeth L, Kim J,
#'   Resche-Rigon M, Carpenter JR, Williamson EJ (2019). Propensity score
#'   analysis with partially observed covariates: how should multiple
#'   imputation be used? *Statistical Methods in Medical Research*, 28(1),
#'   3--19. \doi{10.1177/0962280217713032}
#'
#' @param formula Model formula with the outcome column, untransformed, on the
#'   left. A `.` on the right-hand side expands to the data's other columns,
#'   excluding the outcome, `id_col` and `imputation_col`.
#' @param data A data frame containing the outcome and predictors.
#' @param family Model family: `"binary"` (default, [stats::glm()] with a logit
#'   link), `"ordinal"` (proportional-odds cumulative logit, [MASS::polr()]; a
#'   two-level outcome is fitted with [stats::glm()]) or `"nominal"`
#'   (generalized logit, [nnet::multinom()]).
#' @param outcome_col Outcome column. By default it is read from `formula`.
#' @param id_col Patient identifier column. Default `"id"`. Required, and
#'   checked, only when `imputation_col` is set, where it aligns each patient's
#'   rows across imputations.
#' @param imputation_col Imputation-index column, or `NULL` for one dataset.
#' @param outcome_levels Complete outcome levels in their declared order.
#' @param event_level Event level for a binary model.
#' @param reference_level Reference level for a nominal model.
#' @param prediction_prefix Name for a binary probability column, or prefix for
#'   ordinal and nominal level-specific probability columns.
#' @param trace Whether model engines should print fitting progress.
#'
#' @return An object of class `c("lm_fit", "ps_data")` with:
#' \describe{
#'   \item{`$data`}{`data` (the lowest-numbered imputation when
#'     `imputation_col` is set, without that column) with the predicted
#'     probabilities appended: one column named `prediction_prefix` for a
#'     binary model (`prob` by default), the probability of `event_level`; or
#'     one `<prediction_prefix>_<level>` column per outcome level for an
#'     ordinal or nominal model, whose outcome column is returned as a factor
#'     (ordered for an ordinal model). A row missing a predictor has missing
#'     probabilities.}
#'   \item{`$meta`}{Named list: `bundle_version` (1), `formula`, `outcome_col`,
#'     `model_family`, `outcome_levels`, `event_level` (binary only, else
#'     `NULL`), `reference_level` (nominal only, else `NULL`),
#'     `cumulative_direction` (ordinal only, else `NULL`), `predictors`,
#'     `id_col`, `imputation_col`, `n_imputations`, `n_input`, `n_analyzed`,
#'     `n_excluded` (summed over imputations), `prediction_prefix` and
#'     `package_versions` (R, this package and the model engine).}
#'   \item{`$tables`}{Named list: `estimates` (one row per term: `term`,
#'     `estimate`, `std.error`, `statistic`, `df`, `p.value`, `conf.low`,
#'     `conf.high`, `odds_ratio` (`NA` for an ordinal threshold), and
#'     `pooled`, `TRUE` when Rubin's rules
#'     combined imputations; a pooled `df` is Barnard and Rubin's, see
#'     "Pooled degrees of freedom"), `covariance` (the coefficient covariance matrix,
#'     pooled when imputed), `by_imputation` (`imputation`, `term`,
#'     `estimate`, `std.error`) and `fit_status` (one row per imputation:
#'     `imputation`, `converged`, `n_input`, `n_analyzed`, `n_excluded`).}
#'   \item{`$models`}{The fitted models, one per imputation, named by
#'     imputation index (`"1"` for a single dataset). For binary bundles,
#'     [validate_logistic()] scores new data with them; it rejects ordinal and
#'     nominal bundles.}
#' }
#'
#' @examples
#' dta <- sample_ps_data(n = 150, seed = 42)
#' fit <- fit_logistic(
#'   tavr ~ age + female + ef,
#'   data           = dta,
#'   family         = "binary",
#'   outcome_levels = c(0, 1),
#'   event_level    = 1
#' )
#' fit$tables$estimates[, c("term", "estimate", "odds_ratio")]
#' head(fit$data[, c("id", "tavr", "prob")])
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
  if (!is.name(formula[[2L]])) {
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
  if (!is.null(imputation_col)) .check_cols(data, c(id_col, imputation_col))
  formula <- .expand_dot_formula(formula, data, c(id_col, imputation_col))
  .check_cols(data, all.vars(formula))
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
  if (identical(family, "ordinal") && !two_level) {
    rlang::check_installed("MASS", reason = "to fit an ordinal logistic model.")
  }
  if (identical(family, "nominal")) {
    rlang::check_installed("nnet", reason = "to fit a nominal logistic model.")
  }

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
        # polr()'s default optim() tolerance stops with a gradient near 1e-4;
        # a tighter one reaches the maximum SAS PROC LOGISTIC reports.
        MASS::polr(
          formula,
          data = model_data,
          Hess = TRUE,
          na.action = stats::na.exclude,
          control = list(reltol = 1e-14, maxit = 1000L)
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
    predictor_cols = .retained_predictors(formula)
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
    predictors = .retained_predictors(formula),
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


#' Stop when a fitted linear model has an aliased term
#'
#' `stats::lm()` returns `NA` for a term it cannot estimate, which would
#' otherwise surface later as an anonymous non-finite coefficient.
#'
#' @param models Named list of fitted models.
#' @return Invisible `NULL`.
#' @keywords internal
.check_aliased <- function(models) {
  for (i in seq_along(models)) {
    aliased <- names(stats::coef(models[[i]]))[is.na(stats::coef(models[[i]]))]
    if (length(aliased)) {
      label <- if (is.null(names(models))) i else names(models)[[i]]
      rlang::abort(
        sprintf("Term(s) aliased in imputation %s: %s. Remove them from the formula.",
                label, paste(aliased, collapse = ", ")),
        call. = FALSE
      )
    }
  }
  invisible(NULL)
}


#' Resolve `outcome_col` from the formula's response
#'
#' The response must read exactly one column, whether or not `outcome_col` is
#' supplied: `I(count + other)` or `cbind(count, other)` would fit one quantity
#' while the stable-response check, and the metadata, used another. A
#' transformed response such as `log(ef)` reads one column, `ef`. A supplied
#' `outcome_col` must name that column.
#'
#' @param formula Two-sided model formula.
#' @param outcome_col Supplied outcome column name, or `NULL`.
#' @return The outcome column name.
#' @keywords internal
.check_outcome_response <- function(formula, outcome_col = NULL) {
  response <- all.vars(formula[[2L]])
  if (length(response) != 1L) {
    rlang::abort("The formula's response must name exactly one column.", call. = FALSE)
  }
  if (is.null(outcome_col)) return(response)
  if (!is.character(outcome_col) || length(outcome_col) != 1L ||
        !identical(outcome_col, response)) {
    rlang::abort(
      sprintf("`outcome_col` (\"%s\") is not the formula's response (%s).",
              paste(outcome_col, collapse = ", "), paste(response, collapse = ", ")),
      call. = FALSE
    )
  }
  response
}


#' Variables the model's retained terms read
#'
#' `all.vars()` on the right-hand side also returns variables a formula
#' subtracts, such as `dropme` in `y ~ . - dropme`. These are the variables of
#' the terms the model keeps, plus any `offset()`.
#'
#' @param formula Model formula.
#' @return Character vector of column names, possibly empty.
#' @keywords internal
.retained_predictors <- function(formula) {
  rhs <- stats::delete.response(stats::terms(formula))
  kept <- lapply(attr(rhs, "term.labels"), str2lang)
  offsets <- as.list(attr(rhs, "variables"))[-1L][attr(rhs, "offset")]
  unique(as.character(unlist(lapply(c(kept, offsets), all.vars))))
}


#' Expand a dot formula over the model columns only
#'
#' `stats::terms()` expands `.` to every other column, which would fit the
#' patient identifier and imputation index as covariates. The expansion is done
#' once, against `data` without those columns and without the response
#' variables, and the explicit formula is returned.
#'
#' @param formula Two-sided model formula.
#' @param data Model data.
#' @param exclude Column names that `.` must not expand to, such as the
#'   identifier and imputation columns. `NULL` entries are ignored.
#' @return `formula` unchanged when its right-hand side has no `.`, otherwise
#'   the expanded formula, in the original formula's environment.
#' @keywords internal
.expand_dot_formula <- function(formula, data, exclude = NULL) {
  if (!"." %in% all.vars(formula[[3L]])) return(formula)
  # Name the columns the formula spells out first: terms() only warns about one
  # that is missing, and a subtracted typo would otherwise be dropped silently.
  .check_cols(data, setdiff(all.vars(formula), "."))
  keep <- setdiff(names(data), c(exclude, all.vars(formula[[2L]])))
  if (!length(keep)) {
    rlang::abort("`.` in `formula` leaves no columns once the response, `id_col` and `imputation_col` are set aside.",
                 call. = FALSE)
  }
  expanded <- stats::formula(stats::terms(formula, data = data[keep]))
  environment(expanded) <- environment(formula)
  expanded
}


#' Fit an ordinary least-squares model bundle
#'
#' Fits `stats::lm()` to one dataset, or once per imputation of stacked
#' long-form data, and combines imputations with Rubin's rules. A single fit is
#' inferred on the t distribution with its residual degrees of freedom, as
#' `PROC REG` reports it. Pooled coefficients use Barnard and Rubin's
#' degrees of freedom, described below.
#'
#' @inheritSection fit_logistic Pooled degrees of freedom
#'
#' @param formula Model formula. The response may be transformed, e.g.
#'   `log(los) ~ age`. A `.` on the right-hand side expands to the data's other
#'   columns, excluding the response, `id_col` and `imputation_col`.
#' @param data A data frame, or stacked imputations when `imputation_col` is set.
#' @param outcome_col Outcome column. By default the one variable in the
#'   formula's response; if supplied, it must name that variable. The
#'   response must read exactly one column.
#' @param id_col Patient identifier column. Required only when `imputation_col` is set.
#' @param imputation_col Imputation-index column, or `NULL` for one dataset.
#' @param prediction_col Name of the fitted-value column added to `$data`,
#'   averaged across imputations.
#'
#' @return An object of class `c("rm_fit", "ps_data")`. `$tables$fit_quality`
#'   gives each fit's residual degrees of freedom, residual standard error and
#'   R-squared; R-squared is `NA` for imputed data, where Rubin's rules define
#'   none.
#'
#' @seealso [fit_logistic()], [bs_continuous()]
#'
#' @examples
#' dta <- sample_ps_data(n = 150, seed = 42)
#' fit <- fit_linear(ef ~ age + female, data = dta)
#' fit$tables$estimates[, c("term", "estimate", "std.error", "p.value")]
#' @export
fit_linear <- function(formula, data, outcome_col = NULL, id_col = "id",
                       imputation_col = NULL, prediction_col = "fitted") {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    rlang::abort("`formula` must be a two-sided R formula.", call. = FALSE)
  }
  .check_df(data)
  outcome_col <- .check_outcome_response(formula, outcome_col)
  .check_cols(data, outcome_col)
  if (!is.numeric(data[[outcome_col]])) {
    rlang::abort(sprintf("The outcome `%s` must be numeric for a linear model.", outcome_col),
                 call. = FALSE)
  }

  if (!is.character(prediction_col) || length(prediction_col) != 1L ||
        is.na(prediction_col) || !nzchar(prediction_col)) {
    rlang::abort("`prediction_col` must be one non-empty string.", call. = FALSE)
  }
  if (prediction_col %in% names(data)) {
    rlang::abort(sprintf("Column `%s` already exists in `data`.", prediction_col), call. = FALSE)
  }
  formula <- .expand_dot_formula(formula, data, c(id_col, imputation_col))
  .check_cols(data, all.vars(formula))

  fitted <- .fit_imputations(
    data = data,
    response_col = outcome_col,
    id_col = id_col,
    imputation_col = imputation_col,
    fit_fn = function(df) stats::lm(formula, data = df, na.action = stats::na.exclude),
    predict_fn = function(fit, df) stats::predict(fit, newdata = df),
    predictor_cols = .retained_predictors(formula)
  )
  .check_aliased(fitted$models)
  inference <- .model_inference(fitted$models, "linear")
  imputed <- length(fitted$models) > 1L

  scored <- fitted$data
  scored[[prediction_col]] <- unname(fitted$predictions)
  fit_quality <- data.frame(
    imputation = names(fitted$models),
    df_residual = vapply(fitted$models, function(m) m$df.residual, numeric(1L)),
    residual_se = vapply(fitted$models, stats::sigma, numeric(1L)),
    # Rubin's rules define no pooled R-squared, so none is reported for imputed data.
    r_squared = if (imputed) NA_real_ else summary(fitted$models[[1L]])$r.squared,
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  meta <- list(
    formula = formula,
    outcome_col = outcome_col,
    id_col = id_col,
    imputation_col = imputation_col,
    prediction_col = prediction_col,
    bundle_version = 1L,
    model_family = "linear",
    package_versions = c(
      R = as.character(getRversion()),
      hvtiRpropensity = as.character(utils::packageVersion("hvtiRpropensity")),
      stats = as.character(utils::packageVersion("stats"))
    ),
    method = if (imputed) "ols-MI" else "ols",
    n_imputations = length(fitted$imputations),
    n_total = nrow(scored)
  )
  tables <- list(
    estimates = inference$estimates,
    covariance = inference$covariance,
    by_imputation = inference$by_imputation,
    fit_status = fitted$status,
    fit_quality = fit_quality
  )
  new_ps_data(scored, meta, tables, fitted$models, subclass = "rm_fit")
}
