###############################################################################
## helper_barnard_rubin.R
##
## The SAS oracles ran PROC MIANALYZE without EDF=, so their DF column is
## Rubin's (1987) large-sample df. The package reports Barnard and Rubin's
## (1999) df, as mice::pool() does. Both are closed-form in the fraction of
## missing information, so the Barnard-Rubin df, interval and p-value are
## rebuilt here from SAS's own Estimate, StdErr and DF.
###############################################################################

barnard_rubin_from_sas <- function(sas, df_com, conf_level = 0.95) {
  m <- sas$NImpute
  lambda <- sqrt((m - 1) / sas$DF)
  observed <- (1 - lambda) * (1 + df_com) * df_com
  df <- (m - 1) * observed / ((df_com + 3) * (m - 1) + lambda^2 * observed)
  critical <- stats::qt(1 - (1 - conf_level) / 2, df)
  list(
    df = df,
    conf.low = sas$Estimate - critical * sas$StdErr,
    conf.high = sas$Estimate + critical * sas$StdErr,
    p.value = 2 * stats::pt(abs(sas$Estimate / sas$StdErr), df, lower.tail = FALSE)
  )
}
