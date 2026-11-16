# Accept a renamed argument for one release

Accept a renamed argument for one release

## Usage

``` r
.renamed_arg(fn, old, new, value, new_supplied)
```

## Arguments

- fn:

  The exported function, for the message.

- old, new:

  The old and new argument names.

- value:

  The value passed under the old name.

- new_supplied:

  Whether the new name was supplied too; both abort.

## Value

`value`.
