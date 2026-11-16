# Name the exported function a condition comes from

Errors are raised with `stop(call. = FALSE)`, so the message carries the
function name itself, as `fn(): `. An internal helper does not know
which exported function called it, so it asks: this returns `"fn(): "`
for the outermost function on the call stack that this package exports,
the one the user called, or `""` when there is none (a helper called
directly in a test).

## Usage

``` r
.err_prefix()
```

## Value

A character string.
