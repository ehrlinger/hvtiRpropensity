# Assign quintile, decile, and k-cluster strata from a score vector

Ranks the scored patients and assigns `ceiling(rank * k / n)`, which
gives groups of equal size when `n` divides by `k`. This does not
reproduce the SAS %grp macro exactly. The macro assigns
`quintile = int(_n_ / (nobs/5)) + 1; if quintile > 5 then quintile = 5;`,
which puts one patient too few in the first group and folds the surplus
into the last: with n = 10, SAS gives quintile sizes 1, 2, 2, 2, 3 and
this gives 2, 2, 2, 2, 2. A patient can sit one stratum apart in the
two.

## Usage

``` r
.assign_bs_strata(score, n_strata = 10L)
```

## Arguments

- score:

  Numeric vector of balancing scores.

- n_strata:

  Integer number of clusters (default 10 = deciles).

## Value

A named list: `quintile`, `decile`, `cluster`. A missing score gets
missing strata, and the strata divide the scored patients only.
