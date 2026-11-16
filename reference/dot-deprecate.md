# Warn once per session about a deprecated usage

The hvtiR family deprecates a renamed argument, or a removed default,
for one release: the old usage keeps working and warns once per session,
with class `hvtiRpropensity_deprecated`, then becomes an error in the
next release.

## Usage

``` r
.deprecate(fn, what, message)
```

## Arguments

- fn:

  The exported function, for the message.

- what:

  A short key for the deprecated usage, such as `"n_rep"`.

- message:

  What changed and what to do instead.

## Value

Invisible `NULL`.
