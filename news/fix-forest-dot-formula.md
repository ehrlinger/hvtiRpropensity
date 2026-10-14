* `ps_forest()` accepts a dot formula such as `tavr ~ .`. It stopped with
  "Required column(s) not found in `data`: ." before. The dot expands to every
  column except the treatment and `id_col`, as it does in `ps_logistic()`.
  A column the formula subtracts, as in `. - x3`, is no longer required, fitted,
  checked for missing values or reported in the SMD table.
* `ps_logistic()`, `ps_ordinal()` and `ps_nominal()` leave a constant column
  the formula subtracts out of the SMD table.
