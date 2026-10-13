# Resolve an argument against the metadata of a ps_data object

One rule for every `sa_*()` function. An argument the caller supplied
wins over the object's `$meta`, and if the two disagree the call aborts
and names both values. An argument the caller did not supply takes the
metadata value, then the function's default.

## Usage

``` r
.resolve_meta_arg(
  name,
  value,
  supplied,
  meta,
  meta_name = name,
  call_env = rlang::caller_env()
)
```

## Arguments

- name:

  The argument name, for the error message.

- value:

  The argument value (the default when not supplied).

- supplied:

  `TRUE` when the caller supplied the argument, from `!missing()`. A
  supplied `NULL` counts as not supplied.

- meta:

  The value recorded in the object's metadata, or `NULL`.

- meta_name:

  The name of that metadata field. Default `name`.

- call_env:

  Environment for the error call.

## Value

The resolved value.
