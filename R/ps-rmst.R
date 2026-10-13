#' Observed-outcome restricted mean survival contrast, weighted
#'
#' Design-based comparison of restricted mean survival time (RMST) between
#' treated and control patients using only the observed follow-up: weighted
#' Kaplan-Meier curves (weights from the propensity score) integrated to
#' `tau`, differenced, with a stratified bootstrap interval.
#'
#' @param x A `ps_data` object with a score and treatment column. A
#'   [ps_match()] object is analysed on its matched pairs only, and a
#'   [ps_weight()] object contributes its stored weight column; see Details.
#'   A score on the logit scale (`scale = "logit"`, as a [ps_match()] object
#'   can carry) is back-transformed with [stats::plogis()] before any weights
#'   are built.
#' @param time_col,event_col Follow-up time and event indicator (0/1) columns
#'   of `x$data`.
#' @param tau Horizon, in the units of `time_col`. It must not exceed the
#'   largest observed time among the rows that enter either arm's curve, for
#'   any analysed subset and weighting.
#' @param weights Weighting schemes to run: any of `"unweighted"`, `"ato"`
#'   (overlap), `"att"`, `"ate"`. For a [ps_match()] input the default is
#'   `"unweighted"` alone, and no other scheme is accepted.
#' @param subsets Optional named list of logical vectors (one value per row);
#'   each is run as a restricted analysis (for example a common-support
#'   subset). The full data is always included as `"all"`.
#' @param n_boot Bootstrap draws. Resampling is by matched pair for a
#'   [ps_match()] input and by row, stratified by arm, otherwise.
#' @param refit If `FALSE` (default) the score is held fixed in the bootstrap,
#'   so the intervals ignore the uncertainty in estimating the propensity
#'   score and understate the total uncertainty. If `TRUE` the score model is
#'   refit on each draw (supported for `ps_logistic` and `ps_forest` inputs).
#' @param seed Bootstrap seed: `NULL` (default) or one whole number. With a seed, the
#'   random number stream is set to `abs(seed)` immediately before the draws,
#'   the result is reproducible whatever ran earlier in the session, the
#'   caller's stream is restored afterwards, and the seed is kept in
#'   `$meta$seed`. With `NULL` the draws come from the caller's stream and
#'   advance it, and `$meta$seed` is `NA`.
#' @param time_unit_days Days per unit of `time_col`, for the `diff_days`
#'   column. Default `365.2425` (time in years).
#' @param clip Scores are clipped to `[clip, 1 - clip]` before weighting.
#'
#' @details
#' **Weights.** `"unweighted"` gives every row weight 1. `"ato"` weights
#' treated rows by `1 - p` and controls by `p`, targeting the overlap
#' population; `"att"` weights controls by `p / (1 - p)`, targeting the
#' treated; `"ate"` weights by `1 / p` and `1 / (1 - p)`, targeting the whole
#' cohort. Here `p` is the score clipped to `[clip, 1 - clip]`.
#'
#' **Matched input.** For a [ps_match()] object only the rows with
#' `match_col == 1` are analysed, so `n` is twice the number of matched pairs,
#' and the estimand is the one the match targets. Weighting a matched set by
#' the score adjusts twice for the same confounders, so only `"unweighted"`
#' is accepted. The bootstrap resamples whole pairs: it draws `pair_id`
#' values with replacement, as many as there are pairs, and takes both members
#' of each drawn pair, so a pair drawn twice enters twice. This keeps the
#' within-pair correlation that a row bootstrap ignores. A subset that splits
#' a pair keeps the member it contains, which is then resampled alone.
#'
#' **Weighted input.** For a [ps_weight()] object, the scheme that matches its
#' estimand (`"ate"` for `estimand = "ATE"`, `"att"` for `"ATT"`) uses the
#' stored column `x$meta$weight_col` instead of recomputing weights from the
#' score, so any trimming is kept. Stabilisation multiplies each arm by a
#' constant, which leaves each arm's Kaplan-Meier curve and effective sample
#' size unchanged. The other schemes are computed from the score, and an
#' `estimand = "ATC"` object has no matching scheme. In the bootstrap each row
#' keeps its stored weight.
#'
#' **Horizon.** The weighted Kaplan-Meier curve is not defined past the last
#' observed time in an arm, counting only the rows that enter the curve (time
#' and event present, weight above zero), so a `tau` beyond it stops with an error rather
#' than carrying the last survival value forward. A bootstrap draw can still
#' end short of `tau`, and its curve is then carried flat to `tau`.
#'
#' **Failures.** A bootstrap draw that errors, or that leaves an arm empty,
#' gives `NA`; the interval is taken over the other draws and the number of
#' `NA` draws is reported as `n_failed`.
#'
#' @return A `ps_rmst` / `ps_data` object. `$tables$estimates` has one row
#'   per subset-and-weighting estimator: `estimator`, `subset`, `weighting`,
#'   `n` (rows analysed; for a matched input, twice the matched pairs),
#'   `ess_treated`, `ess_control`, `rmst_treated`, `rmst_control`, `diff`
#'   (treated minus control, positive favours treated), `diff_days`,
#'   `lo_days`, `hi_days` (2.5 and 97.5 percent bootstrap) and `n_failed`
#'   (bootstrap draws that gave no estimate).
#'   `$tables$curves` has the weighted Kaplan-Meier steps (`estimator`, `arm`,
#'   `time`, `surv`) for plotting; `$data` is `x$data` plus one weight column
#'   `w_<weighting>` per scheme. `$meta$seed` is the seed used (`NA` when none
#'   was given). `$meta$scale` is the scale of `score_col`, carried from `x`
#'   (`"probability"` when `x` records none). `$meta$resampling` is `"pair"`
#'   for a matched input and `"row"` (stratified by arm) otherwise, and
#'   `$meta$stored_weighting` names the scheme that used a [ps_weight()]
#'   object's stored weights (`NA` when none did).
#'
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   dta <- sample_ps_data(n = 200, seed = 5)[, c("id", "tavr", "age", "ef")]
#'   dta$t <- stats::rexp(nrow(dta), 0.2)
#'   dta$e <- as.integer(dta$t < 5)
#'   dta$t <- pmin(dta$t, 5)
#'   res <- ps_rmst(ps_logistic(tavr ~ age + ef, dta), "t", "e", tau = 4, n_boot = 20,
#'                   seed = 1024)
#'   res$tables$estimates
#' }
#' @export
ps_rmst <- function(x, time_col, event_col, tau,
                    weights = c("unweighted", "ato"),
                    subsets = list(),
                    n_boot = 200L, refit = FALSE, seed = NULL,
                    time_unit_days = 365.2425, clip = 1e-3) {
  rlang::check_installed("survival", reason = "to fit Kaplan-Meier curves.")
  if (!is_ps_data(x)) rlang::abort("`x` must be a ps_data object.", call. = FALSE)
  matched <- inherits(x, "ps_match")
  if (matched && missing(weights)) weights <- "unweighted"
  weights <- match.arg(weights, c("unweighted", "ato", "att", "ate"), several.ok = TRUE)
  if (matched && any(weights != "unweighted")) {
    rlang::abort(paste0("A ps_match object is analysed unweighted on its matched pairs; ",
                        "weighting the matched set by the score adjusts twice. Use `weights = \"unweighted\"`."),
                 call. = FALSE)
  }
  .check_cols(x$data, c(time_col, event_col))
  seed <- .seed_value(seed)
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
  score <- .probability_score(x)

  stored_type <- NA_character_
  stored_w <- NULL
  if (inherits(x, "ps_weight")) {
    stored_w <- d[[x$meta$weight_col]]  # read before any w_<type> column can overwrite it
    stored_type <- unname(c(ATE = "ate", ATT = "att")[x$meta$estimand])
    if (is.na(stored_type) || !stored_type %in% weights) stored_type <- NA_character_
  }
  in_match <- if (matched) d[[x$meta$match_col]] %in% 1L else rep(TRUE, nrow(d))
  pair <- if (matched) d[[x$meta$pair_id_col]] else NULL

  wts <- function(sc, a, type) {
    p <- pmin(pmax(sc, clip), 1 - clip)
    switch(type,
      unweighted = rep(1, length(a)),
      ato = ifelse(a, 1 - p, p),
      att = ifelse(a, 1, p / (1 - p)),
      ate = ifelse(a, 1 / p, 1 / (1 - p))
    )
  }
  # Weights for rows `idx`: the stored ps_weight column for its own scheme, else from the score.
  wts_rows <- function(idx, type, sc = score) {
    if (identical(type, stored_type)) stored_w[idx] else wts(sc[idx], arm[idx], type)
  }
  for (type in weights) d[[paste0("w_", type)]] <- wts_rows(seq_len(nrow(d)), type)

  subsets <- c(list(all = rep(TRUE, nrow(d))), lapply(subsets, function(s) as.logical(s) & !is.na(s)))
  subsets <- lapply(subsets, function(s) s & in_match)
  grid <- expand.grid(subset = names(subsets), weighting = weights, stringsAsFactors = FALSE)
  # Horizon guard over the rows that enter each weighted curve, as .weighted_km() keeps them.
  for (i in seq_len(nrow(grid))) {
    idx <- which(subsets[[grid$subset[i]]])
    w <- wts_rows(idx, grid$weighting[i])
    used <- !is.na(time[idx]) & !is.na(event[idx]) & !is.na(w) & w > 0
    for (g in c(TRUE, FALSE)) {
      tg <- time[idx][used & arm[idx] == g]
      if (length(tg) && tau > max(tg)) {
        rlang::abort(sprintf(paste0("`tau` = %g is past the last observed time (%g) in the %s arm of ",
                                    "estimator \"%s/%s\". Choose a `tau` within follow-up."),
                             tau, max(tg), if (g) "treated" else "control", grid$subset[i], grid$weighting[i]),
                     call. = FALSE)
      }
    }
  }

  est_one <- function(idx, type, sc = score) {
    a <- arm[idx]
    w <- wts_rows(idx, type, sc)
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

  if (!is.null(seed)) withr::local_seed(seed)
  # One bootstrap draw of row indices from `idx`: whole pairs for a matched input, else rows within each arm.
  resample <- function(idx) {
    draw <- function(v) v[sample.int(length(v), length(v), replace = TRUE)]
    if (matched) {
      members <- split(idx, pair[idx])
      unlist(members[draw(names(members))], use.names = FALSE)
    } else {
      c(draw(idx[arm[idx]]), draw(idx[!arm[idx]]))
    }
  }
  rows <- vector("list", nrow(grid))
  curves <- list()
  for (i in seq_len(nrow(grid))) {
    idx <- which(subsets[[grid$subset[i]]])
    type <- grid$weighting[i]
    label <- paste(grid$subset[i], type, sep = "/")
    pt <- est_one(idx, type)
    a <- arm[idx]
    w <- wts_rows(idx, type)
    has_w <- !is.na(w)
    ess <- function(g) {
      g <- g & has_w
      if (any(g)) sum(w[g])^2 / sum(w[g]^2) else 0
    }
    diff_pt <- (pt$ra$rmst - pt$rb$rmst)

    draws <- rep(NA_real_, n_boot)
    for (b in seq_len(n_boot)) {
      bi <- resample(idx)
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
    ci <- if (all(is.na(draws))) {
      c(NA_real_, NA_real_)
    } else {
      stats::quantile(draws, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
    }
    rows[[i]] <- data.frame(
      estimator = label, subset = grid$subset[i], weighting = type, n = sum(has_w),
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
             list(scale = rlang::`%||%`(x$meta$scale, "probability"),
                  time_col = time_col, event_col = event_col, tau = tau, n_boot = n_boot, refit = refit,
                  seed = if (is.null(seed)) NA_integer_ else seed, time_unit_days = time_unit_days, n_total = nrow(d),
                  resampling = if (matched) "pair" else "row", stored_weighting = stored_type)),
    tables = list(estimates = do.call(rbind, rows), curves = do.call(rbind, curves)),
    subclass = "ps_rmst"
  )
}

# Weighted Kaplan-Meier -> list(rmst, curve). Zero-weight rows are dropped.
.weighted_km <- function(time, event, weights, tau) {
  keep <- !is.na(time) & !is.na(event) & !is.na(weights) & weights > 0
  dd <- data.frame(t = time[keep], e = event[keep], w = weights[keep])
  if (!nrow(dd)) return(list(rmst = NA_real_, curve = data.frame(time = 0, surv = 1)))
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
