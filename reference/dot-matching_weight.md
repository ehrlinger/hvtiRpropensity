# Li-Greene matching weight (ATM)

`min(p, 1 - p) / p` for treated patients and `min(p, 1 - p) / (1 - p)`
for controls, the `mt_wt` of the SAS templates (Li & Greene, 2013).
Shared by
[`ps_logistic()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_logistic.md)
and
[`ps_forest()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_forest.md)
so the two cannot drift apart.

## Usage

``` r
.matching_weight(p, trt)
```

## Arguments

- p:

  Propensity score, P(treated).

- trt:

  Integer 0/1 treatment indicator, the same length as `p`.

## Value

Numeric vector of matching weights.
