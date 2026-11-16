# Winsorise IPTW weights within the arm or arms that are reweighted

For the ATT, the treated arm is the reference and keeps weights of
exactly 1, so only the control weights are winsorised, at the (`trim`,
`1 - trim`) quantiles of the control weights. The ATC is the mirror
image. For the ATE both arms are reweighted and neither is a reference,
so the weights are winsorised at quantiles of all weights pooled across
both arms, as Lee, Lessler and Stuart (2011) do.

## Usage

``` r
.trim_weights(w, trt, estimand, trim)
```

## Arguments

- w:

  Numeric weight vector; `NA` marks an excluded patient.

- trt:

  Integer 0/1 treatment indicator, the same length as `w`.

- estimand:

  One of `"ATE"`, `"ATT"`, `"ATC"`.

- trim:

  A number in (0, 0.5).

## Value

`w`, winsorised.
