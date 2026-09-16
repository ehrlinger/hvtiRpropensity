###############################################################################
## ps-stddiff-perm.R
##
## Permutation reference for standardized differences, ported from
## stddiffci.sas (Amanda Artis 2020, corrected 2021).
##
## Design: hvtiRtemplates dev/specs/2026-09-16-standardized-difference-design.md, section 5
##
## The macro calls %stddiff once on the observed groups and NPERM times on
## permuted groups, then takes PROC UNIVARIATE percentiles 2.5, 16, 50, 84 and
## 97.5 of each variable's permuted values. PROC UNIVARIATE's default
## percentile definition is PCTLDEF=5, which is R's quantile(type = 2).
###############################################################################


#' Permutation reference for standardized differences
#'
#' Computes [ps_stddiff()] on the observed groups and on `n_perm` random
#' permutations of the group labels, and reports the observed value beside
#' percentiles of its permutation distribution. A port of the CCF
#' `%stddiffci` SAS macro (Artis 2020).
#'
#' @details
#' **What the percentiles are.** They describe the standardized difference
#' you would see if group membership carried no information about the
#' variable. They are not a confidence interval for the observed value,
#' whatever the SAS macro's name suggests: an observed difference outside the
#' 2.5 to 97.5 range is larger than label shuffling alone tends to produce.
#' The 16 and 84 percentiles bound the central 68%, the range the macro's 2021
#' correction fixed.
#'
#' **Percentile definition.** `quantile(type = 2)`, which is SAS
#' `PCTLDEF=5`, the `PROC UNIVARIATE` default the macro relies on.
#' Permutations where a variable's difference is `NA` are left out of that
#' variable's percentiles.
#'
#' **Weights.** A matching or IPTW weight depends on the group, so a weight
#' carried over unchanged onto permuted groups describes neither the observed
#' design nor the null. `reweight` must therefore be given whenever
#' `weight_col` is: a function taking a data frame and returning one weight
#' per row. It is called once on the unpermuted data, for the observed row,
#' and once on each permuted data set, and its result is written to
#' `weight_col` before [ps_stddiff()] runs. The macro instead expects the
#' permuted weights to exist as columns before it is called.
#'
#' **Missing groups.** Rows with a missing group keep it and are not
#' permuted into either group.
#'
#' **Random numbers.** With `seed`, the permutations are reproducible and the
#' caller's random number stream is restored afterwards. Without it, the
#' current stream is used and advanced.
#'
#' @inheritParams ps_stddiff
#' @param weight_col Optional name of the column [ps_stddiff()] weights by.
#'   Requires `reweight`, which supplies its values.
#' @param reweight A function of a data frame returning one positive, finite
#'   weight per row, for example one that refits a propensity model and
#'   returns matching weights. Required when `weight_col` is given, and
#'   rejected without it.
#' @param n_perm Number of permutations. Default `1000`, the macro's.
#' @param seed Optional integer seed.
#'
#' @return An object of class `c("ps_stddiff_perm", "ps_data")` with:
#' \describe{
#'   \item{`$data`}{The input data frame, unchanged.}
#'   \item{`$meta`}{Named list: `treatment_col`, `weight_col`, `variables`,
#'     `method`, `n_perm`, `seed`, `n_total`.}
#'   \item{`$tables`}{Named list with `stddiff_perm`, one row per variable:
#'     `variable`, `label`, `type`, `observed` (as [ps_stddiff()] reports it)
#'     and the permutation percentiles `p2_5`, `p16`, `p50`, `p84`,
#'     `p97_5`.}
#' }
#'
#' @seealso [ps_stddiff()]
#'
#' @examples
#' dta <- sample_ps_data(n = 200, seed = 42)
#' obj <- ps_stddiff_perm(dta, treatment_col = "tavr", gaussian = c("age", "ef"),
#'                        binary = "female", n_perm = 200, seed = 1)
#' obj$tables$stddiff_perm
#'
#' @export
ps_stddiff_perm <- function(data,
                            treatment_col = "tavr",
                            gaussian      = NULL,
                            nong_ord      = NULL,
                            binary        = NULL,
                            categorical   = NULL,
                            weight_col    = NULL,
                            reweight      = NULL,
                            n_perm        = 1000L,
                            seed          = NULL) {
  .check_df(data)
  if (!is.numeric(n_perm) || length(n_perm) != 1L || !is.finite(n_perm) ||
      n_perm < 1 || n_perm > .Machine$integer.max || n_perm != round(n_perm)) {
    rlang::abort("`n_perm` must be a positive whole number.", call = NULL)
  }
  n_perm <- as.integer(n_perm)
  if (!is.null(weight_col) && !is.function(reweight)) {
    rlang::abort(paste0("`weight_col` needs `reweight`: weights that depend on the group must be recomputed ",
                        "for each permutation."), call = NULL)
  }
  if (is.null(weight_col) && !is.null(reweight)) {
    rlang::abort("`reweight` was given without `weight_col` to write its weights to.", call = NULL)
  }
  .check_cols(data, treatment_col)

  if (!is.null(seed)) {
    # Seed before anything that may draw, reweight on the observed data
    # included, and restore the caller's stream on exit.
    rng <- ".Random.seed"
    had_seed <- exists(rng, envir = globalenv(), inherits = FALSE)
    old_seed <- if (had_seed) get(rng, envir = globalenv(), inherits = FALSE)
    on.exit(
      if (had_seed) assign(rng, old_seed, envir = globalenv()) else rm(list = rng, envir = globalenv()),
      add = TRUE
    )
    set.seed(seed)
  }

  run <- function(d, quiet = FALSE) {
    if (!is.null(weight_col)) {
      w <- reweight(d)
      if (!is.numeric(w) || length(w) != nrow(d)) {
        rlang::abort("`reweight` must return one weight per row of the data it is given.", call = NULL)
      }
      d[[weight_col]] <- w
    }
    withCallingHandlers(
      ps_stddiff(d, treatment_col = treatment_col, gaussian = gaussian, nong_ord = nong_ord,
                 binary = binary, categorical = categorical, weight_col = weight_col)$tables$stddiff,
      # A permutation can leave a categorical variable with no shared levels.
      # Its NA is dropped from the percentiles, so that one warning would only
      # repeat; every other warning, reweight's included, still surfaces.
      warning = function(w) {
        if (quiet && grepl("share too few levels", conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning")
      }
    )
  }
  obs <- run(data)

  idx  <- which(!is.na(data[[treatment_col]]))
  sims <- vapply(seq_len(n_perm), function(i) {
    d <- data
    d[[treatment_col]][idx] <- data[[treatment_col]][idx][sample.int(length(idx))]
    run(d, quiet = TRUE)$stddiff
  }, numeric(nrow(obs)))
  sims <- matrix(sims, nrow = nrow(obs))

  pct <- t(apply(sims, 1L, .perm_percentiles))
  colnames(pct) <- c("p2_5", "p16", "p50", "p84", "p97_5")
  tbl <- data.frame(obs[c("variable", "label", "type")], observed = obs$stddiff, pct,
                    stringsAsFactors = FALSE, row.names = NULL)

  new_ps_data(
    data   = data,
    meta   = list(treatment_col = treatment_col, weight_col = weight_col,
                  variables = list(gaussian = gaussian, nong_ord = nong_ord,
                                   binary = binary, categorical = categorical),
                  method = "stddiff_perm", n_perm = n_perm, seed = seed, n_total = nrow(data)),
    tables = list(stddiff_perm = tbl),
    subclass = "ps_stddiff_perm"
  )
}


#' PROC UNIVARIATE percentiles 2.5, 16, 50, 84, 97.5 (PCTLDEF=5, R type 2)
#' @keywords internal
#' @noRd
.perm_percentiles <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(rep(NA_real_, 5L))
  stats::quantile(x, c(0.025, 0.16, 0.5, 0.84, 0.975), type = 2, names = FALSE)
}
