# Compute normalized inference for a model bundle

Compute normalized inference for a model bundle

## Usage

``` r
.model_inference(models, family, conf_level = 0.95)
```

## Arguments

- models:

  Non-empty list of fitted models.

- family:

  One of `"binary"`, `"ordinal"`, `"nominal"`, or `"count"`.

- conf_level:

  Confidence level.

## Value

A list with final estimates, covariance and per-imputation values.
