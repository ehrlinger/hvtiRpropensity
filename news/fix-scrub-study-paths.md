* The SAS templates in `inst/sas-templates/` no longer carry the study they
  were taken from: each study path is now `/studies/<study>`, and the study
  title and cohort line read `<study description>` and
  `CCF, 19xx to 20xx, n=xxxx`, as the `%let STUDYDESC` and `%let STUDYPOP`
  lines already did. The analysis code is unchanged.
