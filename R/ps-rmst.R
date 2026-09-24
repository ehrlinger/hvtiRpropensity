#' Observed-outcome restricted mean survival contrast, weighted
#'
#' Design-based comparison of restricted mean survival time (RMST) between
#' treated and control patients using only the observed follow-up: weighted
#' Kaplan-Meier curves (weights from the propensity score) integrated to
#' `tau`, differenced, with a stratified bootstrap interval.
#'
#' @param x A `ps_data` object with a score and treatment column.
#' @param time_col,event_col Follow-up time and event indicator (0/1) columns
#'   of `x$data`.
#' @param tau Horizon, in the units of `time_col`.
#' @param weights Weighting schemes to run: any of `"unweighted"`, `"ato"`
#'   (overlap), `"att"`, `"ate"`.
#' @param subsets Optional named list of logical vectors (one value per row);
#'   each is run as a restricted analysis (for example a common-support
#'   subset). The full data is always included as `"all"`.
#' @param n_boot Bootstrap draws. Resampling is stratified by arm.
#' @param refit If `FALSE` (default) the score is held fixed in the bootstrap,
#'   which understates uncertainty. If `TRUE` the score model is refit on each
#'   draw (supported for `ps_logistic` and `ps_forest` inputs).
#' @param seed Bootstrap seed.
#' @param time_unit_days Days per unit of `time_col`, for the `diff_days`
#'   column. Default `365.2425` (time in years).
#' @param clip Scores are clipped to `[clip, 1 - clip]` before weighting.
#'
#' @return A `ps_rmst` / `ps_data` object. `$tables$estimates` has one row
#'   per subset-and-weighting estimator: `estimator`, `subset`, `weighting`,
#'   `n`, `ess_treated`, `ess_control`, `rmst_treated`, `rmst_control`, `diff`
#'   (treated minus control, positive favours treated), `diff_days`,
#'   `lo_days`, `hi_days` (2.5 and 97.5 percent bootstrap) and `n_failed`.
#'   `$tables$curves` has the weighted Kaplan-Meier steps (`estimator`, `arm`,
#'   `time`, `surv`) for plotting; `$data` is `x$data` plus one weight column
#'   `w_<weighting>` per scheme.
#'
#' @examples
#' dta <- sample_ps_data(n = 200, seed = 5)[, c("id", "tavr", "age", "ef")]
#' dta$t <- stats::rexp(nrow(dta), 0.2)
#' dta$e <- as.integer(dta$t < 5)
#' dta$t <- pmin(dta$t, 5)
#' res <- ps_rmst(ps_logistic(tavr ~ age + ef, dta), "t", "e", tau = 4, n_boot = 20)
#' res$tables$estimates
#' @export
ps_rmst <- function(x, time_col, event_col, tau,
                    weights = c("unweighted", "ato"),
                    subsets = list(),
                    n_boot = 200L, refit = FALSE, seed = 1024L,
                    time_unit_days = 365.2425, clip = 1e-3) {
  rlang::check_installed("survival", reason = "to fit Kaplan-Meier curves.")
  if (!is_ps_data(x)) rlang::abort("`x` must be a ps_data object.", call. = FALSE)
  weights <- match.arg(weights, c("unweighted", "ato", "att", "ate"), several.ok = TRUE)
  .check_cols(x$data, c(time_col, event_col))
  if (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) || tau <= 0) {
    rlang::abort("`tau` must be one positive number.", call. = FALSE)
  }
  if (length(subsets) && (is.null(names(subsets)) || any(!nzchar(names(subsets))) ||
                          any(lengths(subsets) != nrow(x$data)))) {
    rlang::abort("`subsets` must be a named list of logical vectors, one value per row.", call. = FALSE)
  }
  if (refit && !inherits(x, c("ps_logistic", "ps_forest"))) {
    rlang::abort("`refit = TRUE` needs a ps_logistic or ps_forest object.", call. = FALSE)
  }

  d <- x$data
  arm <- as.character(d[[x$meta$treatment_col]]) == as.character(x$meta$treated_level)
  time <- d[[time_col]]
  event <- as.numeric(d[[event_col]])
  score <- d[[x$meta$score_col]]

  wts <- function(sc, a, type) {
    p <- pmin(pmax(sc, clip), 1 - clip)
    switch(type,
      unweighted = rep(1, length(a)),
      ato = ifelse(a, 1 - p, p),
      att = ifelse(a, 1, p / (1 - p)),
      ate = ifelse(a, 1 / p, 1 / (1 - p))
    )
  }
  for (type in weights) d[[paste0("w_", type)]] <- wts(score, arm, type)

  subsets <- c(list(all = rep(TRUE, nrow(d))), lapply(subsets, function(s) as.logical(s) & !is.na(s)))
  grid <- expand.grid(subset = names(subsets), weighting = weights, stringsAsFactors = FALSE)

  est_one <- function(idx, type, sc = score) {
    a <- arm[idx]
    w <- wts(sc[idx], a, type)
    ra <- .weighted_km(time[idx][a], event[idx][a], w[a], tau)
    rb <- .weighted_km(time[idx][!a], event[idx][!a], w[!a], tau)
    list(ra = ra, rb = rb)
  }
  refit_scores <- function(idx) {
    dd <- x$data[idx, , drop = FALSE]
    dd[c(x$meta$score_col, x$meta$logit_col, x$meta$weight_col, "quintile", "decile")] <- NULL
    obj <- if (inherits(x, "ps_forest")) {
      ps_forest(x$meta$formula, dd, treatment_col = x$meta$treatment_col, id_col = x$meta$id_col,
                score_col = x$meta$score_col, treated_level = x$meta$treated_level,
                ntree = x$meta$ntree, clip = clip)
    } else {
      ps_logistic(x$meta$formula, dd, treatment_col = x$meta$treatment_col, id_col = x$meta$id_col,
                  score_col = x$meta$score_col, treated_level = x$meta$treated_level,
                  treatment_levels = x$meta$treatment_levels)
    }
    obj$data[[x$meta$score_col]]
  }

  set.seed(seed)
  rows <- vector("list", nrow(grid))
  curves <- list()
  for (i in seq_len(nrow(grid))) {
    idx <- which(subsets[[grid$subset[i]]])
    type <- grid$weighting[i]
    label <- paste(grid$subset[i], type, sep = "/")
    pt <- est_one(idx, type)
    a <- arm[idx]
    w <- wts(score[idx], a, type)
    ess <- function(g) if (any(g)) sum(w[g])^2 / sum(w[g]^2) else 0
    diff_pt <- (pt$ra$rmst - pt$rb$rmst)

    idx_a <- idx[arm[idx]]
    idx_b <- idx[!arm[idx]]
    draws <- rep(NA_real_, n_boot)
    for (b in seq_len(n_boot)) {
      bi <- c(sample(idx_a, length(idx_a), replace = TRUE), sample(idx_b, length(idx_b), replace = TRUE))
      draws[b] <- tryCatch({
        sc <- if (refit) {
          tmp <- score
          tmp[bi] <- refit_scores(bi)
          tmp
        } else {
          score
        }
        r <- est_one(bi, type, sc)
        (r$ra$rmst - r$rb$rmst) * time_unit_days
      }, error = function(e) NA_real_)
    }
    ci <- stats::quantile(draws, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
    rows[[i]] <- data.frame(
      estimator = label, subset = grid$subset[i], weighting = type, n = length(idx),
      ess_treated = ess(a), ess_control = ess(!a),
      rmst_treated = pt$ra$rmst, rmst_control = pt$rb$rmst,
      diff = diff_pt, diff_days = diff_pt * time_unit_days,
      lo_days = ci[1L], hi_days = ci[2L], n_failed = sum(is.na(draws))
    )
    curves[[i]] <- rbind(
      data.frame(estimator = label, arm = "treated", pt$ra$curve),
      data.frame(estimator = label, arm = "control", pt$rb$curve)
    )
  }

  new_ps_data(
    data = d,
    meta = c(x$meta[c("treatment_col", "score_col", "treated_level")],
             list(time_col = time_col, event_col = event_col, tau = tau, n_boot = n_boot, refit = refit,
                  time_unit_days = time_unit_days, n_total = nrow(d))),
    tables = list(estimates = do.call(rbind, rows), curves = do.call(rbind, curves)),
    subclass = "ps_rmst"
  )
}

# Weighted Kaplan-Meier -> list(rmst, curve). Zero-weight rows are dropped.
.weighted_km <- function(time, event, weights, tau) {
  keep <- !is.na(time) & !is.na(event) & !is.na(weights) & weights > 0
  dd <- data.frame(t = time[keep], e = event[keep], w = weights[keep])
  km <- survival::survfit(survival::Surv(t, e) ~ 1, data = dd, weights = dd$w)
  curve <- data.frame(time = c(0, km$time), surv = c(1, km$surv))
  list(rmst = .rmst_step(curve$time, curve$surv, tau), curve = curve)
}

# Area under a right-continuous step curve on [0, tau].
.rmst_step <- function(time, surv, tau) {
  ord <- order(time)
  time <- time[ord]
  surv <- surv[ord]
  if (max(time) < tau) {
    time <- c(time, tau)
    surv <- c(surv, surv[length(surv)])
  } else if (max(time) > tau) {
    keep <- time <= tau
    last <- surv[keep][sum(keep)]
    time <- c(time[keep], tau)
    surv <- c(surv[keep], last)
  }
  sum(diff(time) * surv[-length(surv)])
}

#' @export
print.ps_rmst <- function(x, ...) {
  cat(sprintf("<ps_rmst>  tau = %g, %d bootstrap draws%s\n", x$meta$tau, x$meta$n_boot,
              if (x$meta$refit) " (score refit)" else " (score fixed)"))
  e <- x$tables$estimates
  print(data.frame(estimator = e$estimator, n = e$n, diff_days = round(e$diff_days, 1),
                   lo = round(e$lo_days, 1), hi = round(e$hi_days, 1)), row.names = FALSE)
  invisible(x)
}
