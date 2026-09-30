# Group labels for a binary treatment column

Group labels for a binary treatment column

## Usage

``` r
.binary_group_labels(x, treated_level = NULL)
```

## Arguments

- x:

  A 0/1 or logical treatment vector.

- treated_level:

  The treated value. `NULL` means `1` or `TRUE`.

## Value

A character vector of length 2: the reference label, then the treated
label.
