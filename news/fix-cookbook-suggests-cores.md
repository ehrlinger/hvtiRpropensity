* `fit_logistic()` checks that 'MASS' (ordinal) or 'nnet' (nominal) is
  installed before it fits, and the tests that use either package skip when
  it is absent (#89).
* The test suite and the `ps_forest()` example cap 'randomForestSRC' at two
  cores (#89).
