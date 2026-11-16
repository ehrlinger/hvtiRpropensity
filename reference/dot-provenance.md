# The hvtiR family provenance record

A named list with the seven fields every hvtiR package records, in this
order: `package` and `version` (character), `sha` (the installed commit,
`RemoteSha`, else `NA`), `method` (what was run), `seed` (integer, `NA`
when none was given), `rng_kind`
([`RNGkind()`](https://rdrr.io/r/base/Random.html), three strings) and
`time` (POSIXct). A field that does not apply is `NA` of its type, so a
provenance table shows a gap rather than a plausible value. The family's
spec, not a shared helper: each package keeps its own copy and a
conformance test.

## Usage

``` r
.provenance(method, seed = NULL)
```

## Arguments

- method:

  What was run, for example `"logistic-binary"`.

- seed:

  The seed used, or `NULL`.

## Value

A named list.
