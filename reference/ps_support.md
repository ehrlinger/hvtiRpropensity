# Compare support (overlap) rules

Flags patients with weak support under several rules and measures how
well the rules agree. Built-in rules use the propensity score of a
`ps_data` object; any other rule (for example an isolation-forest tail
flag) can be supplied as a logical vector.

## Usage

``` r
ps_support(x, flags = list(), trim = c(0.1, 0.9), common = TRUE)
```

## Arguments

- x:

  A `ps_data` object with a score and treatment column. The score of a
  [`ps_match()`](https://ehrlinger.github.io/hvtiRpropensity/reference/ps_match.md)
  object built with `scale = "logit"` is a logit, and is
  back-transformed with
  [`stats::plogis()`](https://rdrr.io/r/stats/Logistic.html), so `trim`
  and the common-support interval are on the probability scale.

- flags:

  Optional named list of extra logical vectors, one element per row of
  `x$data`, `TRUE` meaning the rule would drop that patient. The names
  `"trim"` and `"common"` belong to the built-in rules and are refused.

- trim:

  Length-2 numeric, or `NULL`. Adds the rule "score outside
  `[trim[1], trim[2]]`". Both values must lie in \[0, 1\], increasing.
  Default `c(0.1, 0.9)`.

- common:

  Logical; add the rule "score outside the range occupied by both
  groups" (common support).

## Value

A `ps_support` / `ps_data` object. `$data` is `x$data` plus one logical
column `weak_<rule>` per rule. `$tables$counts` gives the number flagged
per rule; `$tables$agreement` gives, for every pair of rules, `n_both`,
`n_neither`, `pct_agree`, Cohen `kappa` and `jaccard` (overlap of the
two dropped sets). `$meta$common_support` is the interval used, on the
probability scale. `$meta$scale` is the scale of `score_col`, carried
from `x` (`"probability"` when `x` records none).

## Examples

``` r
dta <- sample_ps_data(n = 200, seed = 3)[, c("id", "tavr", "age", "ef")]
sup <- ps_support(ps_logistic(tavr ~ age + ef, dta))
sup$tables$agreement
#>   rule_a rule_b   n n_a n_b n_both n_neither pct_agree     kappa   jaccard
#> 1 common   trim 400  17   7      7       383      97.5 0.5727409 0.4117647
```
