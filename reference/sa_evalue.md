# E-Value Sensitivity Analysis

Computes the E-value for a point estimate and optionally for a
confidence interval bound, following VanderWeele & Ding (2017). The
E-value is the minimum strength of association (on the risk-ratio scale)
that an unmeasured confounder would need with **both** the treatment and
the outcome, conditional on measured covariates, to fully explain away
the observed association.

## Usage

``` r
sa_evalue(
  estimate,
  ci_lo = NULL,
  ci_hi = NULL,
  type = c("RR", "OR", "HR", "IRR", "RD"),
  p0 = NULL,
  rare = NULL
)
```

## Arguments

- estimate:

  Numeric. The point estimate (RR, OR, HR, or RD on the original scale).
  Must be \> 0 for RR/OR/HR; any finite value for RD.

- ci_lo:

  Optional numeric. Lower confidence interval bound.

- ci_hi:

  Optional numeric. Upper confidence interval bound. Both `ci_lo` and
  `ci_hi` must be supplied together to compute the CI E-value, and must
  satisfy `ci_lo <= estimate <= ci_hi`.

- type:

  One of `"RR"` (default), `"OR"`, `"HR"`, `"IRR"`, or `"RD"`. `"IRR"`
  (incidence rate ratio) uses the same formula as `"RR"`; it is provided
  as a distinct value for clarity in count-outcome studies.

- p0:

  Baseline outcome probability in the unexposed group. Required when
  `type = "RD"`, and used only then. Must be in (0, 1).

- rare:

  `TRUE` or `FALSE`; required when `type = "OR"` or `"HR"`. Whether the
  outcome is rare in the study population (prevalence under about 15% at
  the end of follow-up). `TRUE` uses the OR or HR as the RR; `FALSE`
  applies the common-outcome conversion. The default `NULL` aborts for
  OR and HR. Ignored for other types.

## Value

An object of class `c("sa_evalue", "ps_data")`. There is no input data,
so `$data` is an empty data frame and `$models` an empty list.

- `$meta$evalue_estimate`:

  E-value for the point estimate.

- `$meta$evalue_ci`:

  E-value for the confidence interval bound closest to the null. `NA` if
  `ci_lo`/`ci_hi` are not supplied. `1` if the CI already crosses the
  null.

- `$meta$rr_estimate`:

  The point estimate converted to the risk-ratio scale (same as
  `estimate` when `type = "RR"`).

- `$meta$type`:

  The measure type supplied.

- `$meta$rare`:

  For `type = "OR"` or `"HR"`, the rare-outcome choice used for the
  conversion; `NA` for other types.

- `$meta$estimate`, `$meta$ci_lo`, `$meta$ci_hi`, `$meta$p0`:

  The inputs, `NA` where not supplied.

- `$tables$evalue`:

  Data frame with one row for the estimate and one for the CI bound
  closest to the null: `quantity`, `value` on the original scale, `rr`
  on the risk-ratio scale (`NA` when the CI crosses the null) and
  `evalue`.

## Details

Odds ratios and hazard ratios are converted to the risk-ratio scale
before the E-value formula is applied, and the conversion depends on
whether the outcome is rare (VanderWeele & Ding 2017; VanderWeele 2017):

- **Rare outcome** (prevalence at most about 15%): the OR or HR
  approximates the RR and is used directly.

- **Common outcome**: `RR = sqrt(OR)` for an odds ratio, and
  `RR = (1 - 0.5^sqrt(HR)) / (1 - 0.5^sqrt(1 / HR))` for a hazard ratio.

There is no default for this choice: for `type = "OR"` or `"HR"` the
caller must set `rare`. Rarity refers to the outcome's prevalence in the
study population (under about 15% at the end of follow-up), and it
cannot be inferred from the baseline risk `p0` alone: with OR = 4 and p0
= 0.10 the exposed risk is about 0.31, and reversing the exposure coding
would flip any rule based on `p0`. The same conversion is applied to the
confidence interval bounds.

For risk differences a baseline outcome probability `p0` must be
supplied, and the conversion `RR = (RD + p0) / p0` treats `p0` as fixed:
its own sampling uncertainty is not carried into the CI E-value.

## References

VanderWeele, T. J., & Ding, P. (2017). Sensitivity analysis in
observational research: Introducing the E-value. *Annals of Internal
Medicine*, 167(4), 268-274.
[doi:10.7326/M16-2607](https://doi.org/10.7326/M16-2607)

VanderWeele, T. J. (2017). On a square-root transformation of the odds
ratio for a common outcome. *Epidemiology*, 28(6), e58-e60.
[doi:10.1097/EDE.0000000000000733](https://doi.org/10.1097/EDE.0000000000000733)

## See also

[`sa_rosenbaum()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_rosenbaum.md),
[`sa_overlap()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_overlap.md),
[`sa_trim_sweep()`](https://ehrlinger.github.io/hvtiRpropensity/reference/sa_trim_sweep.md)

## Examples

``` r
# Risk ratio of 2.5 -- how robust is the result?
sa_evalue(2.5, ci_lo = 1.8, ci_hi = 3.5)
#> <sa_evalue>
#>   Estimate    : RR = 2.5 (CI 1.8 to 3.5)
#>   RR scale    : 2.5
#>   E-value     : 4.43649 (estimate), 3 (CI)

# Odds ratio of 0.6 (protective effect) for a common outcome
sa_evalue(0.6, ci_lo = 0.4, ci_hi = 0.9, type = "OR", rare = FALSE)
#> <sa_evalue>
#>   Estimate    : OR = 0.6 (CI 0.4 to 0.9)
#>   Rare outcome: FALSE
#>   RR scale    : 0.774597
#>   E-value     : 1.90392 (estimate), 1.292878 (CI)

# Hazard ratio of 1.8 for an outcome under 15% at end of follow-up
sa_evalue(1.8, ci_lo = 1.2, ci_hi = 2.7, type = "HR", rare = TRUE)
#> <sa_evalue>
#>   Estimate    : HR = 1.8 (CI 1.2 to 2.7)
#>   Rare outcome: TRUE
#>   RR scale    : 1.8
#>   E-value     : 3 (estimate), 1.689898 (CI)

# Risk difference of 0.10 with baseline probability 0.20
ev <- sa_evalue(0.10, ci_lo = 0.02, ci_hi = 0.18, type = "RD", p0 = 0.20)
ev$meta$evalue_estimate
#> [1] 2.366025
ev$tables$evalue
#>   quantity value  rr   evalue
#> 1 estimate  0.10 1.5 2.366025
#> 2       ci  0.02 1.1 1.431662
```
