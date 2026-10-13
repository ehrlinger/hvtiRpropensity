# Assign quintile, decile, and k-cluster strata from a score vector

Mirrors the SAS %grp macro: proc sort by balancing_score; quintile =
int(*n* / (nobs/5)) + 1; if quintile \> 5 then quintile = 5; cluster =
decile;

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
