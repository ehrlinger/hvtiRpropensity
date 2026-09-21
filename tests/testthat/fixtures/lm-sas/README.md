# Synthetic SAS model oracles

`input.csv` is deterministic, contains no study data, and is the sole input to
`generate-fixtures.sas`. The SAS output is the independent oracle for response
coding, parameter names, covariance ordering, predicted values, and Rubin
pooling. Do not substitute R-generated values.

## Controlled run

1. Copy the complete `lm-oracle-20260921` folder to a path visible to SAS.
2. In `program/generate-fixtures.sas`, set only `%let root=` to that folder's
   absolute server path. Do not add a trailing slash.
3. Confirm that `input/input.csv` and the empty `out/` directory exist.
4. Submit the SAS program. Return the complete `out/` directory, including the
   log, without editing or opening and resaving any CSV.

The last log marker must be:

```text
NOTE: LM_ORACLE_COMPLETE expected_outputs=19.
```

The run must produce these 19 CSVs:

- binary, ordinal, nominal, and count: `*-estimates.csv`,
  `*-covariance.csv`, and `*-scores.csv` (12 files);
- multiply-imputed binary: `mi-estimates.csv`, `mi-covariance.csv`,
  `mi-outest.csv`, and `mi-scores.csv` (4 files);
- pooled binary: `pooled-estimates.csv` and `pooled-covariance.csv` (2 files);
- `sas-environment.csv` (1 file).

It also writes `generate-fixtures.log`. The raw log is retained with the
controlled run outside Git because its SAS header contains institutional
licensing metadata. `sas-run-summary.txt` records its SHA-256 checksum, clean
diagnostic counts and completion marker. The log is not included in the
19-CSV count. Before acceptance, confirm that it has no `ERROR:` and no warning
about convergence, separation, Hessian, covariance, or `PROC MIANALYZE`. The
environment CSV records the exact SAS product/version, platform, and run time.

## Why raw tables are returned

The program exports raw ODS and scored tables rather than renaming parameters
inside SAS. Once the controlled result is returned, the R tests will map the
observed SAS schema explicitly and compare estimates, covariance, scores, and
pooled inference. This prevents either implementation from defining the other
one's answer.

Status on 2026-09-21: the controlled run completed under SAS
9.04.01M8P02222023 on Linux. The 19 raw CSVs and the non-sensitive run summary
are committed; the raw log remains in the controlled-run folder.
`test_lm_sas_oracle.R` checks coefficients, covariance, predictions, and pooled
inference against them. The ordinal covariance check uses an absolute tolerance
because SAS and `MASS::polr()` use different Hessian calculations; their
coefficient, standard-error, and probability results agree at the documented
cross-engine tolerances.
