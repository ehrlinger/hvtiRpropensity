# Deterministic synthetic input for the rm SAS oracle. Contains no study data.
# Run from the package root: Rscript tests/testthat/fixtures/rm-sas/generate-input.R
set.seed(20261013)
n <- 200L
base <- data.frame(
  id = seq_len(n),
  age = round(stats::rnorm(n, 62, 11), 1),
  female = stats::rbinom(n, 1L, 0.4)
)
base$y <- round(30 + 0.25 * base$age - 3 * base$female + stats::rnorm(n, 0, 6), 3)
imputations <- lapply(1:5, function(m) {
  d <- base
  # Imputation varies a covariate only; the response is the same in every imputation.
  d$age <- round(d$age + stats::rnorm(n, 0, 1.5), 1)
  cbind(`_Imputation_` = m, d)
})
utils::write.csv(do.call(rbind, imputations), "tests/testthat/fixtures/rm-sas/input.csv", row.names = FALSE)
