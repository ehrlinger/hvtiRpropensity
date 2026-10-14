# Synthetic SAS model oracles

`input.csv` is deterministic, contains no study data, and is the sole input to
`generate-fixtures.sas`. The SAS output is the independent oracle for response
coding, parameter names, covariance ordering, predicted values, and Rubin
pooling. Do not substitute R-generated values.

## Controlled run

1. Copy a folder holding `program/generate-fixtures.sas`, `input/input.csv` and
   an empty `out/` to a path visible to SAS.
2. In `program/generate-fixtures.sas`, set only `%let root=` to that folder's
   absolute server path. Do not add a trailing slash.
3. Confirm that `input/input.csv` and the empty `out/` directory exist.
4. Submit the SAS program. Return the complete `out/` directory, including the
   log, without editing or opening and resaving any CSV.

The last log marker must be:

```text
NOTE: LM_ORACLE_COMPLETE expected_outputs=21.
```

The run must produce these 21 CSVs:

- binary, ordinal, nominal, and count: `*-estimates.csv`,
  `*-covariance.csv`, and `*-scores.csv` (12 files);
- multiply-imputed binary: `mi-estimates.csv`, `mi-covariance.csv`,
  `mi-outest.csv`, and `mi-scores.csv` (4 files);
- pooled binary: `pooled-estimates.csv` and `pooled-covariance.csv` (2 files);
- pooled binary with `EDF=`: `pooled-estimates-edf.csv` and `complete-df.csv`
  (2 files);
- `sas-environment.csv` (1 file).

It also writes `generate-fixtures.log`. The raw log is retained with the
controlled run outside Git because its SAS header contains institutional
licensing metadata. `sas-run-summary.txt` records its SHA-256 checksum, clean
diagnostic counts and completion marker. The log is not included in the
21-CSV count. Before acceptance, confirm that it has no `ERROR:` and no warning
about convergence, separation, Hessian, covariance, or `PROC MIANALYZE`. The
environment CSV records the exact SAS product/version, platform, and run time.

## Why raw tables are returned

The program exports raw ODS and scored tables rather than renaming parameters
inside SAS. Once the controlled result is returned, the R tests will map the
observed SAS schema explicitly and compare estimates, covariance, scores, and
pooled inference. This prevents either implementation from defining the other
one's answer.

Status on 2026-10-14: the controlled run (`lm-oracle-20261014c`) completed
under SAS 9.04.01M8P02222023 on Linux. The 21 raw CSVs and the non-sensitive
run summary are committed; the raw log remains in the controlled-run folder.
The program strips display formats before export, as the rm oracle does, so
every value carries full precision, and every `PROC LOGISTIC` fit runs at
`GCONV=1e-12`; the log reports that criterion satisfied for all six fits.

Two earlier runs that day isolated the causes of the remaining gaps. The run
with the default `GCONV=` (1e-8, a relative gradient criterion) left binary
coefficients 7.8e-5 from `glm()`; at 1e-12 they agree to about 1e-9, and the
binary, nominal and pooled checks are tested near that. The `PROC GENMOD`
Poisson fit has no `GCONV=`; it already agrees to 1e-6 or better and is tested
there too. Two ordinal gaps are
not SAS convergence and stay at looser tolerances. `polr()`'s own default
tolerance moves its coefficients by about 2e-5. The ordinal covariance check
uses an absolute tolerance because SAS and `MASS::polr()` use different
information matrices: SAS's covariance is the inverse expected Fisher
information, `polr()` inverts the observed Hessian. Tightening either fit's
convergence leaves that gap (about 3e-3) unchanged.

## Degrees of freedom

`pooled-estimates.csv` is `PROC MIANALYZE` without `EDF=`: Rubin's (1987)
large-sample `DF`, and the `LCLMean`, `UCLMean` and `Probt` built on it. The
package reports Barnard and Rubin's (1999) degrees of freedom, as
`mice::pool()` does, with complete-data degrees of freedom n - k.
`generate-fixtures.sas` therefore runs `PROC MIANALYZE` a second time with
`EDF=` set to n - k for imputation 1 and writes `pooled-estimates-edf.csv`,
with the n, k and EDF it used in `complete-df.csv`.

The controlled run on 2026-10-14 produced both files. The tests compare `df`,
`p.value`, `conf.low` and `conf.high` against `pooled-estimates-edf.csv`
directly, check `complete-df.csv` against `df.residual()` of the first
imputation's model, and check that `helper_barnard_rubin.R`, applied to the
default run's `DF`, reproduces the `EDF=` degrees of freedom.
