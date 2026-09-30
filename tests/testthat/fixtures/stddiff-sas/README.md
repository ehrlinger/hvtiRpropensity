# Synthetic SAS oracle for the standardized-difference ports

`input.csv` is synthetic and deterministic (`make-input.R`), holds no study
data, and is the only input to `generate-fixtures.sas`. That program runs the
production CCF macros `%stddiff`, `%stddiffci` and `%mw_var` from `!MACROS`, and
the CSVs it writes are the oracle for `test_stddiff_sas_oracle.R`. Do not
substitute R-generated values. Tracked in issue #34.

## What the input exercises
- Unequal groups: 70 in group 0, 50 in group 1.
- One variable of each `%stddiff` type: `x_gauss` (Gaussian), `x_ord`
  (non-Gaussian or ordinal, ranked), `x_bin` (binary) and `x_cat` (categorical,
  four levels).
- Missing values in `x_gauss`, `x_ord` and `x_cat`, so the oracle also shows how
  each side drops them.
- A propensity score, `p_score`, and the matching weight derived from it, `w`.
- Two outcomes for `%mw_var`: `y_cont` and `y_bin`.
- Twenty permuted groups, `fgrp_1` to `fgrp_20`, with weights `w_1` to `w_20`.
  `%stddiffci` reads its permutations from the input rather than drawing them,
  so R scores the same permuted groups and its percentiles can be compared
  exactly. Each `w_k` is the matching weight re-derived from `p_score` for the
  permuted group `fgrp_k`, as the spec requires: fixed weights would answer a
  different question. The R test re-derives them the same way, through the
  `reweight` step `ps_stddiff_perm()` uses, and checks they equal the columns
  SAS reads.

## Controlled run
1. Stage this folder as `/studies/general/_development/stddiff-oracle-20260930/`
   with `program/generate-fixtures.sas`, `input/input.csv` and an empty `out/`.
2. In the program, check `%let root=` matches that path. No trailing slash.
3. Submit the program. The last log line must be:

   ```text
   NOTE: STDDIFF_ORACLE_COMPLETE expected_outputs=10.
   ```

4. Copy these from `out/` into this directory, without opening or resaving any
   CSV: `sd-unweighted.csv`, `sd-weighted.csv`, `ci-unweighted.csv`,
   `ci-weighted.csv`, `mw-summary.csv`, `mw-replicates.csv`,
   `sas-environment.csv`.

The run also writes three `macro-*.sas` files, the macro sources as they ran,
plus `mw-var.rtf` and `generate-fixtures.log`. The log stays with the controlled
run, outside Git, because its SAS header carries licensing metadata. Before
accepting the run, confirm it has no `ERROR:`.

## What is compared
| SAS output | R | tolerance |
|---|---|---|
| `%stddiff`, unweighted and weighted | `ps_stddiff()` | 4 decimals |
| `%stddiffci` observed value and 2.5/16/50/84/97.5 percentiles | `ps_stddiff()` on each permuted group, then `.perm_percentiles(type = 2)` | 4 decimals |
| `%mw_var` difference, group means, SDs and sums of weights | `ps_mw_var()` | 1e-6 |
| `%mw_var` bootstrap SD and percentiles | `sd()` and `.perm_percentiles(type = 4)` on SAS's own replicates | 1e-6 |
| `%mw_var` bootstrap SD | `ps_mw_var()` with 2,000 of its own replicates | 4 Monte Carlo standard errors of the SD ratio (about 21%) |

`%mw_var` draws its bootstrap with `PROC SURVEYSELECT`, which R cannot
reproduce, so the replicates themselves are not compared.

The tests skip, naming the missing files, if the output CSVs are absent.

## Status
Run 2026-09-30 under SAS 9.04.01M8P02222023 on Linux, twice. The first run
held the permutation weights fixed; `input.csv` now re-derives them (`w_1` to
`w_20` changed, `p_score` was added, every other column is unchanged), and the
second run, 16:57, is the one committed here. Its log ends
`NOTE: STDDIFF_ORACLE_COMPLETE expected_outputs=10` and its SHA-256 is
`5f7e0520f241ad320f24b931df60e29b06d6d97548fcfe3844884392ce0b11ab`. It has no `ERROR:` except the `%mw_var` RTF report, which cannot
render in a batch session without fonts; its datasets are exported
afterwards. The three macro sources that ran match `~/Documents/macro.library`
line for line, apart from line endings.

Between the two runs, every output that does not depend on the new columns or
on random draws came out byte-identical: both `%stddiff` runs, the unweighted
`%stddiffci`, and `%mw_var`'s group means, SDs and sums of weights. The
`%mw_var` bootstrap differs, as it should with `SEED=-1`.

Every comparison in `test_stddiff_sas_oracle.R` passes. The largest
difference in a standardized difference is 3.5e-7, for the categorical
variable, whose value `%stddiff` passes through a macro variable and so rounds
to six significant digits. The `%mw_var` difference agrees to 3e-10.

Two things about the macros this run showed:
- `%mw_var` passes the same `SEED=` to `PROC SURVEYSELECT` on every replicate,
  so a positive seed draws one sample repeatedly and the bootstrap SD is 0. Only
  the default `SEED=-1` gives a bootstrap. `ps_mw_var()` seeds once.
- `%mw_var`'s report step prints `_LABEL_`, so unlabelled outcomes stop it.
