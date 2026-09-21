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
#' @param estimates List of identically named coefficient vectors.
#' @param covariances List of identically named covariance matrices.
#' @param conf_level Confidence level.
#' @return A list with final estimates, covariance and per-imputation values.
#' @keywords internal
.rubin_pool <- function(estimates, covariances, conf_level = 0.95) {
  .check_conf_level(conf_level)
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

  between_diag <- diag(between)
  within_diag <- diag(ubar)
  df <- ifelse(
    between_diag == 0,
    Inf,
    (m - 1) * (1 + within_diag / ((1 + 1 / m) * between_diag))^2
  )
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


#' Extract normalized coefficients and covariance from a fitted model
#'
#' @param fit Fitted model object.
#' @param family Model family.
#' @return A list with `estimate` and `covariance`.
#' @keywords internal
.model_components <- function(fit, family) {
  if (family %in% c("binary", "count")) {
    estimate <- stats::coef(fit)
    covariance <- stats::vcov(fit)
    terms <- names(estimate)
  } else if (identical(family, "ordinal")) {
    estimate <- c(stats::coef(fit), fit$zeta)
    terms <- c(names(stats::coef(fit)), paste0("threshold:", names(fit$zeta)))
    names(estimate) <- terms
    covariance <- stats::vcov(fit)
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


#' Compute normalized inference for a model bundle
#'
#' @param models Non-empty list of fitted models.
#' @param family One of `"binary"`, `"ordinal"`, `"nominal"`, or `"count"`.
#' @param conf_level Confidence level.
#' @return A list with final estimates, covariance and per-imputation values.
#' @keywords internal
.model_inference <- function(models, family, conf_level = 0.95) {
  .check_conf_level(conf_level)
  family <- match.arg(family, c("binary", "ordinal", "nominal", "count"))
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
    result <- .rubin_pool(estimates, covariances, conf_level)
  } else {
    terms <- names(estimates[[1L]])
    .check_inference_input(estimates[[1L]], covariances[[1L]], terms, 1L)
    result <- list(
      estimates = .inference_table(
        terms, estimates[[1L]], covariances[[1L]],
        rep.int(Inf, length(terms)), conf_level, FALSE
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
