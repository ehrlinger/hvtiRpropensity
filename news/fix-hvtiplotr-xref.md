* The help pages for `ps_match()`, `ps_weight()` and the package now name
  `hvtiPlotR::hv_mirror_hist()` and `hvtiPlotR::hv_balance()` as code rather
  than linking to them. 'hvtiPlotR' is not a declared dependency, so
  `R CMD check` could not resolve the links on a machine without it and
  reported "Unknown package 'hvtiPlotR' in Rd xrefs".
