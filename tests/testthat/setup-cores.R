# Cap parallel workers at two for the whole test run. randomForestSRC uses
# every core by default; CRAN allows at most two. A setup file rather than a
# helper, because teardown_env() exists only once the test run has started.
withr::local_options(rf.cores = 2L, mc.cores = 2L, .local_envir = testthat::teardown_env())
