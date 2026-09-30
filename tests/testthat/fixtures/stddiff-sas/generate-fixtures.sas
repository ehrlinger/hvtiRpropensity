/*
 * Synthetic oracle for hvtiRpropensity ps_stddiff(), ps_stddiff_perm() and
 * ps_mw_var(), ported from the CCF macros %stddiff, %stddiffci and %mw_var.
 *
 * This program reads only input/input.csv, which is synthetic and holds no
 * study data. It runs the production macros from !MACROS and exports their
 * raw output datasets, so no R code decides the oracle's values. It also
 * copies the macro sources it ran into out/, to record which versions they
 * were.
 *
 * Tracked in https://github.com/ehrlinger/hvtiRpropensity/issues/34
 */

/* EDIT: absolute server path to this folder. No trailing slash. */
%let root=/studies/general/_development/stddiff-oracle-20260930;
%let input_dir=&root./input;
%let output_dir=&root./out;

filename orcllog "&output_dir./generate-fixtures.log";
proc printto log=orcllog new;
run;

%put NOTE: STDDIFF oracle root is &root.;
%put NOTE: SAS version is &sysvlong4.;
%put NOTE: Platform is &sysscpl.;

/* ---- Record the macro sources that run ---------------------------------- */
%macro copy_source(name);
  filename msrc "!MACROS/&name..sas";
  filename mdst "&output_dir./macro-&name..sas";
  data _null_;
    infile msrc lrecl=32767;
    file mdst lrecl=32767;
    input;
    put _infile_;
  run;
  filename msrc clear;
  filename mdst clear;
%mend copy_source;

%copy_source(stddiff);
%copy_source(stddiffci);
%copy_source(wt_mtch_boot_sd);

/* stddiffci.sas includes stddiff.sas itself. */
%include "!MACROS/stddiffci.sas";
%include "!MACROS/wt_mtch_boot_sd.sas";

/* ---- Input ---------------------------------------------------------------- */
proc import datafile="&input_dir./input.csv"
    out=fixture dbms=csv replace;
  guessingrows=max;
  getnames=yes;
run;

/* ---- %stddiff: observed standardized differences -------------------------- */
%stddiff(DSN=fixture,
         GROUP=grp,
         GAUSSIAN=x_gauss,
         NONG_ORD=x_ord,
         BINARY=x_bin,
         CATG=x_cat,
         BLACKBOX=1,
         OUT=sd_unweighted);

%stddiff(DSN=fixture,
         GROUP=grp,
         GAUSSIAN=x_gauss,
         NONG_ORD=x_ord,
         BINARY=x_bin,
         CATG=x_cat,
         WEIGHT=w,
         BLACKBOX=1,
         OUT=sd_weighted);

/* ---- %stddiffci: permutation reference ------------------------------------ */
/* The macro reads permuted groups fgrp_1..fgrp_20 and weights w_1..w_20 from   */
/* the input; input.csv supplies them, so R can use the same permutations.     */
%stddiffci(DSN=fixture,
           GROUP=grp,
           GAUSSIAN=x_gauss,
           NONG_ORD=x_ord,
           BINARY=x_bin,
           CATG=x_cat,
           NPERM=20,
           BLACKBOX=1,
           OUT=ci_unweighted);

%stddiffci(DSN=fixture,
           GROUP=grp,
           GAUSSIAN=x_gauss,
           NONG_ORD=x_ord,
           BINARY=x_bin,
           CATG=x_cat,
           WEIGHT=w,
           NPERM=20,
           BLACKBOX=1,
           OUT=ci_weighted);

/* ---- Export the %stddiff and %stddiffci results -------------------------- */
/* Before %mw_var, whose report step can fail in batch (see below).             */
%macro export(ds, file);
  proc export data=&ds.
      outfile="&output_dir./&file."
      dbms=csv replace;
  run;
%mend export;

%export(sd_unweighted,   sd-unweighted.csv);
%export(sd_weighted,     sd-weighted.csv);
%export(ci_unweighted,   ci-unweighted.csv);
%export(ci_weighted,     ci-weighted.csv);

/* ---- %mw_var: matching-weight difference, bootstrap SD -------------------- */
/* %mw_var sorts and rewrites its input, so it gets a copy. Its report step     */
/* prints _LABEL_, which exists only when the outcomes carry labels.            */
data mwin;
  set fixture(keep=id grp w y_cont y_bin);
  label y_cont='Continuous outcome'
        y_bin='Binary outcome';
run;

/* SEED is left at the macro default (-1). The macro passes the same SEED to   */
/* PROC SURVEYSELECT on every replicate, so a positive seed draws the same     */
/* sample 200 times and the bootstrap SD is 0. The tests compare SAS's own     */
/* replicates, so they need not be reproducible.                               */
%mw_var(indat=mwin,
        group=grp,
        varlist=y_cont y_bin,
        mt_weight=w,
        RESAMPL=200,
        out=&output_dir./mw-var.rtf);

/* The macro leaves NONOTES set, and its RTF report can fail in a batch session */
/* without fonts, which puts SAS into syntax-check mode. Its output datasets    */
/* exist by then; recover so they are exported.                                 */
options notes obs=max nosyntaxcheck;
%let syscc=0;

/* ---- Environment ---------------------------------------------------------- */
data sas_environment;
  length item $32 value $200;
  item='sysvlong4'; value="&sysvlong4."; output;
  item='sysscpl';   value="&sysscpl.";   output;
  item='run_datetime'; value=put(datetime(), e8601dt19.); output;
run;

/* ---- Export the %mw_var datasets and the environment ----------------------- */
%export(_myout_01,       mw-summary.csv);
%export(_mwtrt,          mw-replicates.csv);
%export(sas_environment, sas-environment.csv);

/* ---- Completion check ----------------------------------------------------- */
%macro check_outputs;
  %let expected=sd-unweighted.csv sd-weighted.csv ci-unweighted.csv ci-weighted.csv
                mw-summary.csv mw-replicates.csv sas-environment.csv
                macro-stddiff.sas macro-stddiffci.sas macro-wt_mtch_boot_sd.sas;
  %let n_expected=%sysfunc(countw(&expected., %str( )));
  %let n_found=0;
  %do i=1 %to &n_expected.;
    %let f=%scan(&expected., &i., %str( ));
    %if %sysfunc(fileexist(&output_dir./&f.)) %then %let n_found=%eval(&n_found. + 1);
    %else %put ERROR: STDDIFF oracle output missing: &f.;
  %end;
  %if &n_found. = &n_expected. %then
    %put NOTE: STDDIFF_ORACLE_COMPLETE expected_outputs=&n_expected.;
  %else
    %put ERROR: STDDIFF oracle incomplete, &n_found. of &n_expected. outputs;
%mend check_outputs;

%check_outputs;

proc printto;
run;
