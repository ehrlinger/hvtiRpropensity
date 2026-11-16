# Stop when a fitted linear or logistic model has an aliased term

[`stats::lm()`](https://rdrr.io/r/stats/lm.html) and
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) return `NA` for a
term they cannot estimate, which would otherwise surface later as an
anonymous non-finite coefficient. Models from other engines are skipped:
[`MASS::polr()`](https://rdrr.io/pkg/MASS/man/polr.html) drops an
aliased column itself, and a
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html)
coefficient matrix carries no `NA`.

## Usage

``` r
.check_aliased(models)
```

## Arguments

- models:

  Named list of fitted models.

## Value

Invisible `NULL`.
