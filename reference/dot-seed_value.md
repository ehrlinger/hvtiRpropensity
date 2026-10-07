# Validate a seed and return its absolute value

Every function that draws random numbers passes its `seed` through here.
The draws are seeded with `withr::local_seed(abs(seed))` and
[`randomForestSRC::rfsrc()`](https://www.randomforestsrc.org//reference/rfsrc.html)
is given `-abs(seed)`, the sign it requires, so `seed = 5` and
`seed = -5` give the same result.

## Usage

``` r
.seed_value(seed, call_env = rlang::caller_env())
```

## Arguments

- seed:

  `NULL`, or one whole number.

- call_env:

  Environment for the error call.

## Value

`NULL`, or `abs(seed)` as an integer.
