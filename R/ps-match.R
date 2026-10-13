###############################################################################
## ps-match.R
##
## Propensity score nearest-neighbour matching (binary treatment).
##
## Ports the matching logic used in the SAS programs that produce the
## `match` column consumed by `hvtiPlotR::hv_mirror_hist()`.
##
## Workflow:
##   1. dta  <- sample_ps_data()            # or your own data frame
##   2. obj  <- ps_match(dta, ...)          # build ps_match object
##   3. print(obj) / summary(obj)           # inspect diagnostics
##   4. obj$data                            # matched data frame (match == 1)
###############################################################################


#' Propensity Score Nearest-Neighbour Matching
#'
#' Performs 1:1 nearest-neighbour matching without replacement on a
#' pre-computed propensity score column.  Returns a `ps_match` object whose
#' `$data` slot is the **full** original dataset with a `match` column set to
#' `1` for matched pairs and `0` for unmatched patients.
#'
#' The `$tables` slot contains standardised mean difference (SMD) tables
#' before and after matching, and group counts. These are the diagnostics
#' that `hvtiPlotR::hv_mirror_hist()` and `hvtiPlotR::hv_balance()` expect.
#'
#' @param data           A data frame.  Must contain `treatment_col`,
#'   `score_col`, and `id_col`.
#' @param treatment_col  Name of the binary treatment column (0/1 or
#'   logical).  Default `"tavr"`.
#' @param score_col      Name of the numeric propensity score column, on the
#'   scale named by `scale`: probabilities in \[0, 1\] for
#'   `scale = "probability"`, or logits (any finite value, such as the
#'   `logit_t` column [ps_logistic()] appends) for `scale = "logit"`.
#'   Default `"prob_t"`.
#' @param id_col         Name of the patient identifier column.
#'   Default `"id"`.
#' @param caliper        Maximum allowable difference in `score_col`, on its
#'   own scale, for a match to be accepted.  `NULL` (default) imposes no
#'   caliper unless `caliper_sd` is given.  Give at most one of `caliper` and
#'   `caliper_sd`.
#' @param covariates     Character vector of covariate columns for SMD
#'   diagnostics.  If `NULL` (default), all numeric columns other than
#'   `treatment_col` and `score_col` are used.
#' @param match_col      Name of the output match-indicator column appended
#'   to `$data`.  Default `"match"`.  If the column already exists it is
#'   overwritten.
#' @param ratio          Matching ratio (controls : treated).  Currently
#'   only `1` (1:1) is supported.  Default `1`.
#' @param seed           Seed for the random shuffling of treated-patient
#'   order that reduces order-dependence. `NULL` (default) or one whole number. With a
#'   seed, the random number stream is set to `abs(seed)` immediately before the draws,
#'   the result is reproducible whatever ran earlier in the session, the
#'   caller's stream is restored afterwards, and the seed is kept in
#'   `$meta$seed`. With `NULL` the draws come from the caller's stream and
#'   advance it, and `$meta$seed` is `NA`.
#'
#' @param treated_level The value of `treatment_col` that is treated. `NULL`
#'   (default) requires a 0/1 or logical column and takes 1 or `TRUE`. Give a
#'   value, for example `"transcatheter"`, to use a column holding any two
#'   values.
#' @param scale          The scale `score_col` is on, and so the scale distances
#'   and the caliper are measured on: `"probability"` (default) or `"logit"`.
#'   The recommended caliper, 0.2 standard deviations of the logit of the
#'   propensity score (Austin 2011), needs `scale = "logit"`.
#' @param caliper_sd     Caliper in pooled standard-deviation units of
#'   `score_col`: the caliper applied is `caliper_sd * sqrt((var_t + var_c) / 2)`,
#'   where `var_t` and `var_c` are the variances of `score_col` among treated
#'   and control patients (missing values dropped).  `NULL` (default) uses
#'   `caliper` instead.  `scale = "logit", caliper_sd = 0.2`
#'   is the Austin (2011) recommendation.
#' @references Austin PC (2011). Optimal caliper widths for propensity-score
#'   matching when estimating differences in means and differences in
#'   proportions in observational studies. *Pharmaceutical Statistics*,
#'   10(2), 150--161.  \doi{10.1002/pst.433}
#' @return An object of class `c("ps_match", "ps_data")` with:
#' \describe{
#'   \item{`$data`}{The full input data frame with `match_col` and `pair_id`
#'     appended.  `match_col` is `1` for matched patients and `0` otherwise.
#'     `pair_id` is an integer identifying which matched pair each patient
#'     belongs to (`NA` for unmatched patients).  Filter to `match == 1` to
#'     obtain the matched subset.}
#'   \item{`$meta`}{Named list: `treatment_col`, `score_col`, `id_col`,
#'     `match_col`, `pair_id_col`, `scale`, `caliper` (the absolute caliper
#'     applied, on `scale`), `caliper_sd`, `method`, `seed`, `n_total`,
#'     `n_matched`, `n_unmatched`.}
#'   \item{`$tables`}{Named list: `smd_before`, `smd_after`,
#'     `group_counts_before`, `group_counts_after`. Group tables name each
#'     group by its value in `treatment_col` (`0`/`1` or `FALSE`/`TRUE`) and
#'     carry a logical `treated` column.}
#' }
#'
#' @seealso [ps_weight()], [sample_ps_data()],
#'   `hvtiPlotR::hv_mirror_hist()`, `hvtiPlotR::hv_balance()`
#'
#' @examples
#' dta <- sample_ps_data(n = 200, seed = 42)
#'
#' # 1. Build the matched object
#' obj <- ps_match(dta)
#'
#' # 2. Inspect diagnostics
#' print(obj)
#' summary(obj)
#'
#' # 3. Extract the full data (with match indicator)
#' head(obj$data)
#' table(obj$data$match)
#'
#' # 4. Matched subset only
#' matched <- obj$data[obj$data$match == 1L, ]
#' nrow(matched)
#'
#' # 5. Match on the logit with a 0.2 SD caliper (Austin 2011)
#' dta$logit_t <- stats::qlogis(dta$prob_t)
#' obj_logit <- ps_match(dta, score_col = "logit_t", scale = "logit", caliper_sd = 0.2)
#' obj_logit$meta$caliper
#'
#' @export
ps_match <- function(data,
                     treatment_col = "tavr",
                     score_col     = "prob_t",
                     id_col        = "id",
                     caliper       = NULL,
                     covariates    = NULL,
                     match_col     = "match",
                     ratio         = 1L,
                     seed          = NULL,
                     treated_level = NULL,
                     scale         = c("probability", "logit"),
                     caliper_sd    = NULL) {

  # ---- Input validation ---------------------------------------------------
  .check_df(data)
  seed <- .seed_value(seed)
  .check_cols(data, c(treatment_col, score_col, id_col))
  treatment <- .treatment_indicator(data, treatment_col, treated_level)
  original_treatment <- data[[treatment_col]]
  data[[treatment_col]] <- treatment$trt
  scale <- match.arg(scale)
  if (scale == "probability") {
    .check_probability(data, score_col,
                       hint = "For a score on the logit scale, such as `logit_t`, set `scale = \"logit\"`.")
  } else if (!is.numeric(data[[score_col]]) || any(is.infinite(data[[score_col]]))) {
    rlang::abort(sprintf("Column `%s` must contain finite numeric logits.", score_col), call. = FALSE)
  } else if (all(data[[score_col]] >= 0 & data[[score_col]] <= 1, na.rm = TRUE)) {
    rlang::warn(sprintf(paste(
      "Every value of `%s` lies in [0, 1], so it looks like a probability, but `scale = \"logit\"`",
      "expects logits, such as `score_col = \"logit_t\"`."
    ), score_col))
  }

  if (!is.null(caliper) &&
      (!is.numeric(caliper) || length(caliper) != 1L || caliper <= 0)) {
    rlang::abort("`caliper` must be a positive number or NULL.", call. = FALSE)
  }
  if (!is.null(caliper_sd)) {
    if (!is.numeric(caliper_sd) || length(caliper_sd) != 1L || is.na(caliper_sd) || caliper_sd <= 0) {
      rlang::abort("`caliper_sd` must be a positive number or NULL.", call. = FALSE)
    }
    if (!is.null(caliper)) {
      rlang::abort("Give at most one of `caliper` and `caliper_sd`.", call. = FALSE)
    }
    is_t <- as.integer(data[[treatment_col]]) == 1L
    pooled_sd <- sqrt((stats::var(data[[score_col]][is_t], na.rm = TRUE) +
                         stats::var(data[[score_col]][!is_t], na.rm = TRUE)) / 2)
    caliper <- caliper_sd * pooled_sd
  }
  if (!identical(as.integer(ratio), 1L)) {
    rlang::abort("Only 1:1 matching (`ratio = 1`) is currently supported.",
                 call. = FALSE)
  }

  # ---- Warn if pair_id column already exists in input data ---------------
  if ("pair_id" %in% names(data)) {
    rlang::warn(
      "Input data already contains a column named 'pair_id'. It will be overwritten."
    )
  }

  # ---- Split treated / control -------------------------------------------
  trt <- as.integer(data[[treatment_col]])
  ps  <- data[[score_col]]

  idx_t <- which(trt == 1L)   # treated rows
  idx_c <- which(trt == 0L)   # control rows

  ps_t <- ps[idx_t]
  ps_c <- ps[idx_c]

  # ---- Nearest-neighbour matching (greedy, without replacement) ----------
  # Pre-allocate to maximum possible size to avoid O(n^2) vector growth.
  n_t       <- length(idx_t)
  matched_c <- logical(length(idx_c))    # TRUE = this control already used
  pair_t    <- integer(n_t)
  pair_c    <- integer(n_t)
  k         <- 0L                        # matched pairs found so far

  # Randomise treated order to reduce order-dependence.
  # Use seed if supplied for reproducibility.
  if (!is.null(seed)) withr::local_seed(seed)
  order_t <- sample(seq_len(n_t))

  for (i in order_t) {
    remaining_idx <- which(!matched_c)
    if (length(remaining_idx) == 0L) break

    dists     <- abs(ps_t[i] - ps_c[remaining_idx])
    best_pos  <- which.min(dists)           # position within remaining_idx
    best      <- remaining_idx[best_pos]    # index into idx_c
    best_dist <- dists[best_pos]

    if (!is.null(caliper) && best_dist > caliper) next

    k              <- k + 1L
    matched_c[best] <- TRUE
    pair_t[k]      <- idx_t[i]
    pair_c[k]      <- idx_c[best]
  }

  # Trim to actual matched count
  if (k > 0L) {
    pair_t <- pair_t[seq_len(k)]
    pair_c <- pair_c[seq_len(k)]
  } else {
    pair_t <- integer(0)
    pair_c <- integer(0)
  }

  # ---- Build output data frame -------------------------------------------
  out <- data
  out[[match_col]] <- 0L
  out[[match_col]][c(pair_t, pair_c)] <- 1L

  # Pair identifier: integer 1..n_pairs for matched rows; NA for unmatched
  out[["pair_id"]] <- NA_integer_
  if (length(pair_t) > 0L) {
    pair_ids <- seq_along(pair_t)
    out[["pair_id"]][pair_t] <- pair_ids
    out[["pair_id"]][pair_c] <- pair_ids
  }

  # ---- SMD diagnostics ---------------------------------------------------
  if (is.null(covariates)) {
    covariates <- setdiff(
      names(data)[vapply(data, is.numeric, logical(1L))],
      c(treatment_col, score_col, id_col, match_col)
    )
  }

  smd_before <- .smd_table(data, treatment_col, covariates)

  matched_rows <- out[[match_col]] == 1L
  smd_after    <- .smd_table(out[matched_rows, ], treatment_col, covariates)

  n_matched_t <- length(pair_t)
  n_matched_c <- length(pair_c)

  group_labels <- treatment$labels
  out[[treatment_col]] <- original_treatment
  group_counts_before <- data.frame(
    group   = group_labels,
    treated = c(FALSE, TRUE),
    n       = c(length(idx_c), length(idx_t))
  )
  group_counts_after <- data.frame(
    group   = group_labels,
    treated = c(FALSE, TRUE),
    n       = c(n_matched_c, n_matched_t)
  )

  # ---- Assemble ps_data object -------------------------------------------
  new_ps_data(
    data = out,
    meta = list(
      treatment_col = treatment_col,
      treated_level = treatment$labels[[2L]],
      score_col     = score_col,
      id_col        = id_col,
      match_col     = match_col,
      pair_id_col   = "pair_id",
      scale         = scale,
      caliper       = caliper,
      caliper_sd    = caliper_sd,
      method        = "nearest-neighbour 1:1",
      seed          = if (is.null(seed)) NA_integer_ else seed,
      n_total       = nrow(data),
      n_matched     = n_matched_t,   # matched pairs (each arm)
      n_unmatched   = length(idx_t) - n_matched_t
    ),
    tables = list(
      smd_before           = smd_before,
      smd_after            = smd_after,
      group_counts_before  = group_counts_before,
      group_counts_after   = group_counts_after
    ),
    subclass = "ps_match"
  )
}
