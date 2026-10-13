* `ps_match()` gains `scale = c("probability", "logit")` and `caliper_sd`, so
  pairs can be matched on the logit of the propensity score with a caliper in
  standard-deviation units. `scale = "logit", caliper_sd = 0.2` is the
  Austin (2011) recommendation. The default stays `"probability"`, which keeps
  the \[0, 1\] check on `score_col`. `$meta` now records `scale` and
  `caliper_sd`, and `$meta$caliper` is the absolute caliper applied.
