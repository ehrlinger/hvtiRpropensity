# Synthetic SAS OLS oracles

`input.csv` is deterministic (`generate-input.R`), contains no study data, and
is the sole input to `generate-fixtures.sas`. The SAS output is the independent
oracle for `fit_linear()` coefficients, covariance, fitted values and Rubin
pooling. Do not substitute R-generated values.

## Controlled run

1. Copy a folder holding `program/generate-fixtures.sas`, `input/input.csv` and
   an empty `out/` to a path visible to SAS.
2. In the program, set only `%let root=` to that folder's absolute server path,
   without a trailing slash.
3. Submit it. Return the complete `out/` directory without opening and resaving
   any CSV.

The last log marker must be:

```text
NOTE: RM_ORACLE_COMPLETE expected_outputs=8.
```

The eight CSVs are `single-estimates`, `single-covariance`, `single-anova`
(whose `Error` row carries the residual degrees of freedom), `single-scores`,
`pooled-estimates`, `pooled-covariance`, `mi-scores` and `sas-environment`.
Before acceptance, confirm that the log has no `ERROR:` and no `WARNING:`, and
that the ODS tables `ParameterEstimates`, `CovB`, `ANOVA` and `TCov` were all
created. Values are exported unformatted, so a CSV holding `<.0001` or a
five-decimal estimate means the format strip did not run. A
missing ODS table is a NOTE, not an error, so read for it. The raw log stays
outside Git because its header carries licensing metadata;
`sas-run-summary.txt` records its SHA-256, its diagnostic counts and the
completion marker.

## Degrees of freedom

`PROC MIANALYZE` ran without `EDF=`, so `pooled-estimates.csv` holds Rubin's
(1987) large-sample `DF`, and the `LCLMean`, `UCLMean` and `Probt` built on it.
The package reports Barnard and Rubin's (1999) degrees of freedom, as
`mice::pool()` does, with complete-data degrees of freedom n - k. The tests
rebuild those values from SAS's own `Estimate`, `StdErr` and `DF`
(`helper_barnard_rubin.R`). A regenerated fixture should pass `EDF=` set to
n - k, after which the tests can compare against `DF` directly.
