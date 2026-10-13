# Assign rank-based quintile and decile columns

Ranks the scored patients and assigns `ceiling(rank * k / n)`, which
gives groups of equal size when `n` divides by `k`. This does not
reproduce the SAS templates exactly. They assign
`quintile = int(_n_ / (nobs/5)) + 1; if quintile > 5 then quintile = 5;`,
which puts one patient too few in the first group and folds the surplus
into the last: with n = 10, SAS gives quintile sizes 1, 2, 2, 2, 3 and
this gives 2, 2, 2, 2, 2. A patient can sit one stratum apart in the
two.

## Usage

``` r
.assign_ps_strata(score)
```

## Arguments

- score:

  Numeric score; a missing score gets missing strata.

## Value

A named list: `quintile`, `decile`.
