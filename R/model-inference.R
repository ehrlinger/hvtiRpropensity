###############################################################################
## model-inference.R
##
## Parameter normalization and Rubin inference for fitted model bundles.
###############################################################################


#' Validate an inference confidence level
#'
#' @param conf_level Confidence level.
#' @return Invisible `NULL`.
#' @keywords internal
.check_conf_level <- function(conf_level) {
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
        is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    rlang::abort("`conf_level` must be one number strictly between 0 and 1.",
                 call. = FALSE)
  }
  invisible(NULL)
}


#' Validate coefficient and covariance inputs for one imputation
#'
#' @param estimate Named numeric coefficient vector.
#' @param covariance Named covariance matrix.
#' @param terms Reference term order.
#' @param imputation Imputation label used in errors.
#' @return Invisible `NULL`.
#' @keywords internal
.check_inference_input <- function(estimate, covariance, terms, imputation) {
  if (!is.numeric(estimate) || is.null(names(estimate)) ||
        !identical(names(estimate), terms)) {
    rlang::abort(
      sprintf("The coefficient names in imputation %s differ from imputation 1.",
              imputation),
      call. = FALSE
    )
  }
  if (any(!is.finite(estimate))) {
    rlang::abort(
      sprintf("A non-finite coefficient was found in imputation %s.", imputation),
      call. = FALSE
    )
  }
  if (!is.matrix(covariance) || nrow(covariance) != ncol(covariance)) {
    rlang::abort(
      sprintf("The covariance matrix must be square in imputation %s.", imputation),
      call. = FALSE
    )
  }
  if (!identical(rownames(covariance), terms) ||
        !identical(colnames(covariance), terms)) {
    rlang::abort(
      sprintf("The covariance names in imputation %s differ from its coefficients.",
              imputation),
      call. = FALSE
    )
  }
  if (any(!is.finite(covariance))) {
    rlang::abort(
      sprintf("A non-finite covariance was found in imputation %s.", imputation),
      call. = FALSE
    )
  }
  invisible(NULL)
}


#' Build a final inference table
#'
#' @param terms Parameter names.
#' @param estimate Parameter estimates.
#' @param covariance Final covariance matrix.
#' @param df Per-parameter degrees of freedom.
#' @param conf_level Confidence level.
#' @param pooled Whether estimates were pooled.
#' @return A data frame.
#' @keywords internal
.inference_table <- function(terms, estimate, covariance, df, conf_level, pooled) {
  standard_error <- sqrt(diag(covariance))
  statistic <- estimate / standard_error
  alpha <- 1 - conf_level
  critical <- vapply(df, function(value) {
    if (is.infinite(value)) stats::qnorm(1 - alpha / 2) else stats::qt(1 - alpha / 2, value)
  }, numeric(1L))
  p_value <- vapply(seq_along(statistic), function(i) {
    if (is.infinite(df[[i]])) {
      2 * stats::pnorm(abs(statistic[[i]]), lower.tail = FALSE)
    } else {
      2 * stats::pt(abs(statistic[[i]]), df = df[[i]], lower.tail = FALSE)
    }
  }, numeric(1L))
  data.frame(
    term = terms,
    estimate = unname(estimate),
    std.error = unname(standard_error),
    statistic = unname(statistic),
    df = unname(df),
    p.value = unname(p_value),
    conf.low = unname(estimate - critical * standard_error),
    conf.high = unname(estimate + critical * standard_error),
    odds_ratio = unname(exp(estimate)),
    pooled = rep.int(pooled, length(terms)),
    stringsAsFactors = FALSE
  )
}


#' Combine coefficient vectors and covariance matrices with Rubin's rules
#'
#' The degrees of freedom are Barnard and Rubin's (1999) small-sample
#' adjustment, computed as `mice::pool()` computes them. With `df_com = Inf`
#' they reduce to Rubin's (1987) large-sample degrees of freedom.
#'
#' @param estimates List of identically named coefficient vectors.
#' @param covariances List of identically named covariance matrices.
#' @param conf_level Confidence level.
#' @param df_com Complete-data degrees of freedom, from [.complete_df()]. Must be
#'   positive; `Inf` gives Rubin's (1987) degrees of freedom.
#' @return A list with final estimates, covariance and per-imputation values.
#' @keywords internal
.rubin_pool <- function(estimates, covariances, conf_level = 0.95, df_com = Inf) {
  .check_conf_level(conf_level)
  if (!is.numeric(df_com) || length(df_com) != 1L || is.na(df_com) || df_com <= 0) {
    rlang::abort(
      sprintf(paste("Pooling over imputations needs positive complete-data degrees of freedom",
                    "(n - k); the model has %s. A saturated model leaves none."),
              format(df_com)),
      call. = FALSE
    )
  }
  if (!is.list(estimates) || !is.list(covariances) ||
        length(estimates) < 2L || length(estimates) != length(covariances)) {
    rlang::abort(
      "Rubin pooling requires equal lists of estimates and covariances from at least 2 imputations.",
      call. = FALSE
    )
  }
  terms <- names(estimates[[1L]])
  if (is.null(terms) || anyNA(terms) || any(!nzchar(terms))) {
    rlang::abort("Coefficient names are required for Rubin pooling.",
                 call. = FALSE)
  }
  imputation_labels <- names(estimates)
  if (is.null(imputation_labels) || any(!nzchar(imputation_labels))) {
    imputation_labels <- seq_along(estimates)
  }
  for (i in seq_along(estimates)) {
    .check_inference_input(
      estimates[[i]], covariances[[i]], terms, imputation_labels[[i]]
    )
  }

  m <- length(estimates)
  q <- do.call(rbind, estimates)
  qbar <- colMeans(q)
  ubar <- Reduce(`+`, covariances) / m
  between <- stats::cov(q)
  if (length(terms) == 1L) {
    between <- matrix(between, 1L, 1L, dimnames = list(terms, terms))
  } else {
    dimnames(between) <- list(terms, terms)
  }
  total <- ubar + (1 + 1 / m) * between
  dimnames(total) <- list(terms, terms)

  # lambda is the fraction of the total variance due to missing data.
  lambda <- (1 + 1 / m) * diag(between) / diag(total)
  df <- if (is.infinite(df_com)) {
    (m - 1) / lambda^2
  } else {
    observed <- (1 - lambda) * (1 + df_com) * df_com
    (m - 1) * observed / ((df_com + 3) * (m - 1) + lambda^2 * observed)
  }
  by_imputation <- do.call(rbind, lapply(seq_len(m), function(i) {
    data.frame(
      imputation = imputation_labels[[i]],
      term = terms,
      estimate = unname(estimates[[i]]),
      std.error = sqrt(diag(covariances[[i]])),
      stringsAsFactors = FALSE
    )
  }))

  list(
    estimates = .inference_table(terms, qbar, total, df, conf_level, TRUE),
    covariance = total,
    by_imputation = by_imputation
  )
}


#' Covariance of a proportional-odds fit from the expected information
#'
#' [MASS::polr()] inverts the observed information (the numerical Hessian of
#' the log likelihood). SAS `PROC LOGISTIC` reports the inverse of the expected
#' (Fisher) information, and for the cumulative logit, which is not a
#' canonical link, the two differ at the estimate. This computes the expected
#' information analytically, so ordinal standard errors match SAS:
#' \eqn{I = \sum_i w_i \sum_j g_{ij} g_{ij}^T / p_{ij}}, where \eqn{p_{ij}} is
#' the fitted probability of category \eqn{j} and \eqn{g_{ij}} its gradient in
#' the coefficients and thresholds.
#'
#' The category probabilities are computed without cancellation, so extreme
#' linear predictors give finite contributions, and an `offset()` enters through
#' the fitted linear predictor.
#'
#' @param fit A [MASS::polr()] fit with the logistic link and its model frame.
#' @return The covariance matrix, named and ordered as `vcov(fit)`.
#' @keywords internal
.polr_expected_vcov <- function(fit) {
  if (!identical(fit$method, "logistic") || is.null(fit$model)) {
    rlang::abort("Expected information needs a logistic polr() fit that kept its model frame.",
                 call. = FALSE)
  }
  frame <- fit$model
  x <- stats::model.matrix(stats::delete.response(fit$terms), frame)
  x <- x[, colnames(x) != "(Intercept)", drop = FALSE]
  beta <- stats::coef(fit)
  # An intercept-only model has no columns and no coefficients.
  if (ncol(x) != length(beta) || !all(colnames(x) == names(beta))) {
    rlang::abort("The polr() model matrix does not match its coefficients.", call. = FALSE)
  }
  weights <- frame[["(weights)"]]
  if (is.null(weights)) weights <- rep.int(1, nrow(x))
  # polr()'s linear predictor includes any offset().
  eta <- fit$lp
  cuts <- c(-Inf, fit$zeta, Inf)
  n_cuts <- length(fit$zeta)
  info <- matrix(0, length(beta) + n_cuts, length(beta) + n_cuts)
  for (j in seq_len(n_cuts + 1L)) {
    upper <- cuts[[j + 1L]] - eta
    lower <- cuts[[j]] - eta
    # F(u) - F(l) = F(u) S(l) (1 - exp(l - u)) has no cancellation: the naive
    # difference rounds to 0 when both tail probabilities round to 1.
    prob <- stats::plogis(upper) * stats::plogis(-lower) * -expm1(lower - upper)
    density_upper <- stats::dlogis(upper)
    density_lower <- stats::dlogis(lower)
    gradient_cuts <- matrix(0, nrow(x), n_cuts)
    if (j <= n_cuts) gradient_cuts[, j] <- density_upper
    if (j > 1L) gradient_cuts[, j - 1L] <- -density_lower
    gradient <- cbind(-x * (density_upper - density_lower), gradient_cuts)
    scale <- ifelse(prob > 0, sqrt(weights / prob), 0)
    info <- info + crossprod(gradient * scale)
  }
  covariance <- solve(info)
  dimnames(covariance) <- list(c(names(beta), names(fit$zeta)), c(names(beta), names(fit$zeta)))
  covariance
}


#' Extract normalized coefficients and covariance from a fitted model
#'
#' @param fit Fitted model object.
#' @param family Model family.
#' @return A list with `estimate` and `covariance`.
#' @keywords internal
.model_components <- function(fit, family) {
  if (family %in% c("binary", "count", "linear")) {
    estimate <- stats::coef(fit)
    covariance <- stats::vcov(fit)
    terms <- names(estimate)
  } else if (identical(family, "ordinal") && inherits(fit, "glm")) {
    raw_estimate <- stats::coef(fit)
    raw_covariance <- stats::vcov(fit)
    predictor_terms <- setdiff(names(raw_estimate), "(Intercept)")
    threshold_term <- paste0(
      "threshold:",
      paste(attr(fit, "hvti_outcome_levels"), collapse = "|")
    )
    order <- c(match(predictor_terms, names(raw_estimate)), match("(Intercept)", names(raw_estimate)))
    signs <- c(rep.int(1, length(predictor_terms)), -1)
    estimate <- raw_estimate[order] * signs
    terms <- c(predictor_terms, threshold_term)
    names(estimate) <- terms
    covariance <- raw_covariance[order, order, drop = FALSE] * outer(signs, signs)
  } else if (identical(family, "ordinal")) {
    estimate <- c(stats::coef(fit), fit$zeta)
    terms <- c(names(stats::coef(fit)), paste0("threshold:", names(fit$zeta)))
    names(estimate) <- terms
    covariance <- .polr_expected_vcov(fit)
  } else if (identical(family, "nominal")) {
    coefficients <- stats::coef(fit)
    if (is.null(dim(coefficients))) {
      coefficients <- matrix(coefficients, nrow = 1L,
                             dimnames = list(fit$lev[[2L]], names(coefficients)))
    }
    terms <- unlist(lapply(rownames(coefficients), function(level) {
      paste0(level, ":", colnames(coefficients))
    }), use.names = FALSE)
    estimate <- as.vector(t(coefficients))
    names(estimate) <- terms
    covariance <- stats::vcov(fit)
  } else {
    rlang::abort(sprintf("Unsupported model family: %s", family), call. = FALSE)
  }

  if (!is.matrix(covariance) || nrow(covariance) != length(estimate) ||
        ncol(covariance) != length(estimate)) {
    rlang::abort(
      sprintf("The %s model returned a covariance matrix incompatible with its coefficients.",
              family),
      call. = FALSE
    )
  }
  dimnames(covariance) <- list(terms, terms)
  list(estimate = estimate, covariance = covariance)
}


#' Complete-data degrees of freedom of a fitted model
#'
#' The residual degrees of freedom, n minus the number of estimated parameters,
#' which is the complete-data value `mice::pool()` uses. `nnet::multinom()`
#' defines no `df.residual()` method, so n minus its effective degrees of
#' freedom is used, as `mice` does.
#'
#' @param model A fitted model.
#' @return One number; `Inf` when the model reports no residual degrees of freedom.
#' @keywords internal
.complete_df <- function(model) {
  df <- stats::df.residual(model)
  if (is.null(df) && inherits(model, "multinom")) df <- nrow(model$residuals) - model$edf
  if (is.null(df)) Inf else as.numeric(df)
}


#' Compute normalized inference for a model bundle
#'
#' @param models Non-empty list of fitted models.
#' @param family One of `"binary"`, `"ordinal"`, `"nominal"`, `"count"`, or `"linear"`.
#' @param conf_level Confidence level.
#' @return A list with final estimates, covariance and per-imputation values.
#' @keywords internal
.model_inference <- function(models, family, conf_level = 0.95) {
  .check_conf_level(conf_level)
  family <- match.arg(family, c("binary", "ordinal", "nominal", "count", "linear"))
  if (!is.list(models) || !length(models)) {
    rlang::abort("`models` must be a non-empty list of fitted models.",
                 call. = FALSE)
  }
  convergence <- vapply(models, .model_converged, logical(1L))
  if (any(!convergence)) {
    failed <- names(models)[!convergence]
    if (is.null(failed) || any(!nzchar(failed))) failed <- which(!convergence)
    rlang::abort(
      sprintf("Model fitting failed to converge in imputation(s): %s.",
              paste(failed, collapse = ", ")),
      call. = FALSE
    )
  }

  components <- lapply(models, .model_components, family = family)
  estimates <- lapply(components, `[[`, "estimate")
  covariances <- lapply(components, `[[`, "covariance")
  if (length(models) > 1L) {
    result <- .rubin_pool(estimates, covariances, conf_level, .complete_df(models[[1L]]))
  } else {
    terms <- names(estimates[[1L]])
    .check_inference_input(estimates[[1L]], covariances[[1L]], terms, 1L)
    result <- list(
      estimates = .inference_table(
        terms, estimates[[1L]], covariances[[1L]],
        # An OLS fit is inferred on the t distribution with its residual degrees
        # of freedom, as PROC REG reports it; the other families are large-sample.
        if (identical(family, "linear")) {
          rep.int(models[[1L]]$df.residual, length(terms))
        } else {
          rep.int(Inf, length(terms))
        },
        conf_level, FALSE
      ),
      covariance = covariances[[1L]],
      by_imputation = data.frame(
        imputation = 1L,
        term = terms,
        estimate = unname(estimates[[1L]]),
        std.error = sqrt(diag(covariances[[1L]])),
        stringsAsFactors = FALSE
      )
    )
  }
  if (identical(family, "ordinal")) {
    threshold <- startsWith(result$estimates$term, "threshold:")
    result$estimates$odds_ratio[threshold] <- NA_real_
  }
  if (identical(family, "count")) {
    # exp() of a log-link coefficient is a rate ratio, not an odds ratio
    names(result$estimates)[names(result$estimates) == "odds_ratio"] <- "rate_ratio"
  }
  if (identical(family, "linear")) {
    # An identity-link coefficient has no ratio scale; exp() of it means nothing.
    result$estimates$odds_ratio <- NULL
  }
  imputation_labels <- names(models)
  if (is.null(imputation_labels) || any(!nzchar(imputation_labels))) {
    imputation_labels <- seq_along(models)
  }
  result$by_imputation$imputation <- rep(
    imputation_labels,
    each = nrow(result$estimates)
  )
  result
}
