# Complete-data degrees of freedom of a fitted model

The residual degrees of freedom, n minus the number of estimated
parameters, which is the complete-data value `mice::pool()` uses.
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) defines
no [`df.residual()`](https://rdrr.io/r/stats/df.residual.html) method,
so n minus its effective degrees of freedom is used, as `mice` does.

## Usage

``` r
.complete_df(model)
```

## Arguments

- model:

  A fitted model.

## Value

One number; `Inf` when the model reports no residual degrees of freedom.
