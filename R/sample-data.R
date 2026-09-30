###############################################################################
## sample-data.R
##
## Reproducible synthetic dataset for examples and tests.
##
## `sample_ps_data()` generates a data frame that mimics the structure of
## a typical cardiac-surgery SAS export used for propensity score analyses:
##   - binary treatment indicator (`tavr`)
##   - true propensity score (`prob_t`)  computed from a logistic model
##   - several continuous and binary covariates
##   - a match indicator column (`match`) populated with 0s (filled in by
##     `ps_match()`)
###############################################################################


#' Generate sample data for propensity score examples
#'
#' Produces a synthetic data frame that resembles a typical SAS export from a
#' cardiac-surgery propensity score analysis.  Propensity scores are generated
#' from a logistic model driven by the supplied covariates, so the resulting
#' `prob_t` column has realistic overlap and the SMD values are non-trivial.
#'
#' @param n         Number of patients per treatment arm.  Total rows =
#'   `2 * n`.  Must be a positive integer.  Default `500`.
#' @param seed      Random seed for reproducibility.  Default `42L`. The fixed default gives the
#'   same demo dataset on every call; `NULL` gives a fresh one.
#' @param separation Numeric scalar controlling how well the covariates
#'   separate the two groups.  Higher values give less overlap.  Default `1.0`.
#'
#' @return A data frame with `2 * n` rows and the following columns:
#'   \describe{
#'     \item{`id`}{Integer patient identifier.}
#'     \item{`tavr`}{Binary treatment indicator: `0` = SAVR, `1` = TF-TAVR.}
#'     \item{`age`}{Continuous: age in years.}
#'     \item{`female`}{Binary: `1` = female.}
#'     \item{`ef`}{Continuous: ejection fraction (%).}
#'     \item{`diabetes`}{Binary: `1` = diabetic.}
#'     \item{`hypertension`}{Binary: `1` = hypertensive.}
#'     \item{`prob_t`}{Estimated propensity score (numeric, 0 to 1).}
#'     \item{`match`}{Match indicator, initialised to `0`.  Populated by
#'       [ps_match()].}
#'   }
#'
#' @seealso [ps_match()], [ps_weight()]
#'
#' @examples
#' dta <- sample_ps_data(n = 200, seed = 42)
#' head(dta)
#' table(dta$tavr)
#' hist(dta$prob_t, breaks = 20, main = "Propensity score distribution")
#'
#' @export
sample_ps_data <- function(n = 500L, seed = 42L, separation = 1.0) {
  if (!is.numeric(n) || length(n) != 1L || n < 1L)
    stop("`n` must be a positive integer.", call. = FALSE)
  if (!is.numeric(separation) || length(separation) != 1L || separation <= 0)
    stop("`separation` must be a positive number.", call. = FALSE)

  withr::local_seed(seed)
  n <- as.integer(n)

  # Control arm (tavr = 0): older, lower EF
  age0      <- stats::rnorm(n, mean = 75 + separation * 3, sd = 8)
  female0   <- stats::rbinom(n, 1, 0.40)
  ef0       <- stats::rnorm(n, mean = 55 - separation * 3, sd = 10)
  diabetes0 <- stats::rbinom(n, 1, 0.30)
  hyp0      <- stats::rbinom(n, 1, 0.65)

  # Treatment arm (tavr = 1): younger, higher EF
  age1      <- stats::rnorm(n, mean = 75 - separation * 3, sd = 8)
  female1   <- stats::rbinom(n, 1, 0.45)
  ef1       <- stats::rnorm(n, mean = 55 + separation * 3, sd = 10)
  diabetes1 <- stats::rbinom(n, 1, 0.25)
  hyp1      <- stats::rbinom(n, 1, 0.60)

  dta <- data.frame(
    id           = seq_len(2L * n),
    tavr         = c(rep(0L, n), rep(1L, n)),
    age          = c(age0, age1),
    female       = c(female0, female1),
    ef           = c(ef0, ef1),
    diabetes     = c(diabetes0, diabetes1),
    hypertension = c(hyp0, hyp1)
  )

  # Fit a quick logistic model to produce realistic propensity scores
  fit <- stats::glm(
    tavr ~ age + female + ef + diabetes + hypertension,
    data   = dta,
    family = stats::binomial()
  )
  dta$prob_t <- as.numeric(stats::predict(fit, type = "response"))
  dta$match  <- 0L

  dta
}


# ---------------------------------------------------------------------------
# sample_ps_data_ordinal
# ---------------------------------------------------------------------------

#' Generate sample data with an ordinal treatment for ps_ordinal() examples
#'
#' Produces a synthetic data frame with a three-level ordered treatment
#' variable (`nyha_grp`: NYHA functional classes I, II, III), mimicking a
#' typical cardiac-surgery dataset used with
#' `tp.lm.logistic_propensity_ordinal.sas`.
#'
#' @param n    Number of patients per treatment level.  Total rows = `3 * n`.
#'   Default `300L`.
#' @param seed Random seed.  Default `42L`. The fixed default gives the
#'   same demo dataset on every call; `NULL` gives a fresh one.
#'
#' @return A data frame with `3 * n` rows and columns: `id`, `nyha_grp`
#'   (ordered factor: I < II < III), `age`, `female`, `ef`, `diabetes`,
#'   `hypertension`.
#'
#' @seealso [ps_ordinal()]
#'
#' @examples
#' dta <- sample_ps_data_ordinal(n = 100, seed = 7)
#' table(dta$nyha_grp)
#' str(dta)
#'
#' @export
sample_ps_data_ordinal <- function(n = 300L, seed = 42L) {
  if (!is.numeric(n) || length(n) != 1L || n < 1L)
    stop("`n` must be a positive integer.", call. = FALSE)

  withr::local_seed(seed)
  n <- as.integer(n)

  # NYHA I: younger, higher EF, fewer comorbidities
  age0      <- stats::rnorm(n, mean = 62, sd = 9)
  female0   <- stats::rbinom(n, 1, 0.38)
  ef0       <- stats::rnorm(n, mean = 63, sd = 8)
  diabetes0 <- stats::rbinom(n, 1, 0.18)
  hyp0      <- stats::rbinom(n, 1, 0.50)

  # NYHA II: intermediate
  age1      <- stats::rnorm(n, mean = 68, sd = 9)
  female1   <- stats::rbinom(n, 1, 0.42)
  ef1       <- stats::rnorm(n, mean = 57, sd = 9)
  diabetes1 <- stats::rbinom(n, 1, 0.26)
  hyp1      <- stats::rbinom(n, 1, 0.60)

  # NYHA III: older, lower EF, more comorbidities
  age2      <- stats::rnorm(n, mean = 74, sd = 9)
  female2   <- stats::rbinom(n, 1, 0.46)
  ef2       <- stats::rnorm(n, mean = 49, sd = 10)
  diabetes2 <- stats::rbinom(n, 1, 0.34)
  hyp2      <- stats::rbinom(n, 1, 0.70)

  data.frame(
    id           = seq_len(3L * n),
    nyha_grp     = factor(
      c(rep("I", n), rep("II", n), rep("III", n)),
      levels  = c("I", "II", "III"),
      ordered = TRUE
    ),
    age          = c(age0, age1, age2),
    female       = c(female0, female1, female2),
    ef           = c(ef0,  ef1,  ef2),
    diabetes     = c(diabetes0, diabetes1, diabetes2),
    hypertension = c(hyp0, hyp1, hyp2)
  )
}


# ---------------------------------------------------------------------------
# sample_ps_data_nominal
# ---------------------------------------------------------------------------

#' Generate sample data with a nominal treatment for ps_nominal() examples
#'
#' Produces a synthetic data frame with a four-level nominal treatment
#' variable (`rtyp`: tricuspid repair types COS, PER, DEV, CE), mimicking
#' the dataset used in `tp.lm.logistic_propensity.polytomous.sas`.
#'
#' @param n    Number of patients per treatment level.  Total rows = `4 * n`.
#'   Default `200L`.
#' @param seed Random seed.  Default `42L`. The fixed default gives the
#'   same demo dataset on every call; `NULL` gives a fresh one.
#'
#' @return A data frame with `4 * n` rows and columns: `id`, `rtyp`
#'   (unordered factor), `age`, `female`, `ef`, `diabetes`, `hypertension`.
#'
#' @seealso [ps_nominal()]
#'
#' @examples
#' dta <- sample_ps_data_nominal(n = 100, seed = 7)
#' table(dta$rtyp)
#' str(dta)
#'
#' @export
sample_ps_data_nominal <- function(n = 200L, seed = 42L) {
  if (!is.numeric(n) || length(n) != 1L || n < 1L)
    stop("`n` must be a positive integer.", call. = FALSE)

  withr::local_seed(seed)
  n <- as.integer(n)

  # COS (commissurotomy): younger, lower EF
  age0      <- stats::rnorm(n, mean = 58, sd = 10)
  female0   <- stats::rbinom(n, 1, 0.55)
  ef0       <- stats::rnorm(n, mean = 56, sd = 9)
  diabetes0 <- stats::rbinom(n, 1, 0.20)
  hyp0      <- stats::rbinom(n, 1, 0.48)

  # PER (pericardial): intermediate
  age1      <- stats::rnorm(n, mean = 63, sd = 10)
  female1   <- stats::rbinom(n, 1, 0.48)
  ef1       <- stats::rnorm(n, mean = 53, sd = 9)
  diabetes1 <- stats::rbinom(n, 1, 0.25)
  hyp1      <- stats::rbinom(n, 1, 0.55)

  # DEV (DeVega): older, lower EF
  age2      <- stats::rnorm(n, mean = 68, sd = 10)
  female2   <- stats::rbinom(n, 1, 0.44)
  ef2       <- stats::rnorm(n, mean = 50, sd = 10)
  diabetes2 <- stats::rbinom(n, 1, 0.30)
  hyp2      <- stats::rbinom(n, 1, 0.62)

  # CE (Carpentier-Edwards ring): oldest, more comorbidities
  age3      <- stats::rnorm(n, mean = 72, sd = 10)
  female3   <- stats::rbinom(n, 1, 0.41)
  ef3       <- stats::rnorm(n, mean = 48, sd = 11)
  diabetes3 <- stats::rbinom(n, 1, 0.35)
  hyp3      <- stats::rbinom(n, 1, 0.68)

  data.frame(
    id           = seq_len(4L * n),
    rtyp         = factor(
      c(rep("COS", n), rep("PER", n), rep("DEV", n), rep("CE", n)),
      levels = c("COS", "PER", "DEV", "CE")
    ),
    age          = c(age0, age1, age2, age3),
    female       = c(female0, female1, female2, female3),
    ef           = c(ef0,  ef1,  ef2,  ef3),
    diabetes     = c(diabetes0, diabetes1, diabetes2, diabetes3),
    hypertension = c(hyp0, hyp1, hyp2, hyp3)
  )
}


# ---------------------------------------------------------------------------
# sample_ps_data_count
# ---------------------------------------------------------------------------

#' Generate sample data with a count outcome for bs_count() examples
#'
#' Produces a synthetic data frame with a Poisson-distributed count outcome
#' (`rbc_tot`, total intra-operative RBC units transfused), along with
#' baseline covariates and an optional multiply-imputed stacked format.
#' Mirrors the dataset used in `tp.pm.count.balncing_score.sas`.
#'
#' @param n             Number of patients.  Default `500L`.
#' @param seed          Random seed.  Default `42L`. The fixed default gives the
#'   same demo dataset on every call; `NULL` gives a fresh one.
#' @param n_imputations If `> 1`, returns a stacked multiply-imputed data
#'   frame with an `_IMPUTATION_` index column containing `n_imputations`
#'   replicates (each with slightly jittered covariate values to mimic
#'   imputed datasets).  Default `1L` (single complete dataset).
#'
#' @return A data frame with columns: `id`, `rbc_tot` (count), `age`,
#'   `female`, `diabetes`, `hypertension`, `hct_pr` (pre-op haematocrit),
#'   `cpb_time` (cardiopulmonary bypass time, minutes), and optionally
#'   `_IMPUTATION_` when `n_imputations > 1`.
#'
#' @seealso [bs_count()], [bs_continuous()]
#'
#' @examples
#' dta <- sample_ps_data_count(n = 200, seed = 7)
#' hist(dta$rbc_tot, main = "RBC units distribution")
#'
#' # Stacked multiply-imputed version (3 imputations)
#' dta_mi <- sample_ps_data_count(n = 200, seed = 7, n_imputations = 3)
#' table(dta_mi[["_IMPUTATION_"]])
#'
#' @export
sample_ps_data_count <- function(n = 500L, seed = 42L, n_imputations = 1L) {
  if (!is.numeric(n) || length(n) != 1L || n < 1L)
    stop("`n` must be a positive integer.", call. = FALSE)
  if (!is.numeric(n_imputations) || length(n_imputations) != 1L ||
      n_imputations < 1L)
    stop("`n_imputations` must be a positive integer.", call. = FALSE)

  withr::local_seed(seed)
  n             <- as.integer(n)
  n_imputations <- as.integer(n_imputations)

  .one_dataset <- function(seed_offset = 0L) {
    if (!is.null(seed)) withr::local_seed(seed + seed_offset)
    age          <- stats::rnorm(n, mean = 68, sd = 10)
    female       <- stats::rbinom(n, 1, 0.40)
    diabetes     <- stats::rbinom(n, 1, 0.28)
    hypertension <- stats::rbinom(n, 1, 0.60)
    hct_pr       <- stats::rnorm(n, mean = 39, sd = 5)    # pre-op haematocrit
    cpb_time     <- stats::rnorm(n, mean = 90, sd = 30)   # bypass minutes

    # Log-linear mean for RBC: more units for anaemic, long bypass, older
    log_mu <- 0.5 +
      (-0.03) * (hct_pr - 39) +
      0.01  * (cpb_time - 90) +
      0.02  * (age - 68) +
      0.15  * diabetes
    rbc_tot <- stats::rpois(n, lambda = exp(log_mu))

    data.frame(
      id           = seq_len(n),
      rbc_tot      = rbc_tot,
      age          = age,
      female       = female,
      diabetes     = diabetes,
      hypertension = hypertension,
      hct_pr       = hct_pr,
      cpb_time     = cpb_time
    )
  }

  if (n_imputations == 1L) {
    return(.one_dataset(0L))
  }

  # Stacked MI format: same IDs, slightly jittered covariates per imputation
  dsets <- lapply(seq_len(n_imputations), function(m) {
    d                    <- .one_dataset(seed_offset = m * 100L)
    d[["_IMPUTATION_"]] <- m
    d
  })
  do.call(rbind, dsets)
}
