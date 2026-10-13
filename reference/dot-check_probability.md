# Assert that a column contains valid probabilities (0 \<= x \<= 1)

Assert that a column contains valid probabilities (0 \<= x \<= 1)

## Usage

``` r
.check_probability(data, col, call_env = rlang::caller_env(), hint = NULL)
```

## Arguments

- data:

  A data frame.

- col:

  Name of the column to check.

- call_env:

  Environment for the error call.

- hint:

  Optional sentence appended to the range error.

## Value

Invisible NULL.
