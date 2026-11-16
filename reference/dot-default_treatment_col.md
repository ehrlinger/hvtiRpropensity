# The study default `treatment_col` used to have, with a deprecation warning

`treatment_col` defaulted to `"tavr"`, a study's column name. A
study-specific value is never a default in the hvtiR family, so the
argument is now required; for one release an omitted `treatment_col`
still means `"tavr"`, and warns.

## Usage

``` r
.default_treatment_col(fn)
```

## Arguments

- fn:

  The exported function, for the message.

## Value

`"tavr"`.
