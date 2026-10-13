/* rm SAS oracle: PROC REG on one dataset, and per imputation pooled by PROC MIANALYZE.
   Set only ROOT. Input is input/input.csv; output goes to out/. */
%let root=SET_ME;

proc import datafile="&root/input/input.csv" out=work.input dbms=csv replace;
  guessingrows=max;
run;

proc sort data=work.input; by _Imputation_ id; run;

/* One dataset: imputation 1 only. */
ods output ParameterEstimates=work.single_est CovB=work.single_cov;
proc reg data=work.input(where=(_Imputation_=1));
  model y = age female / covb;
  output out=work.single_scores p=fitted;
run; quit;

/* Per imputation, then Rubin's rules. */
proc reg data=work.input outest=work.mi_outest covout noprint;
  by _Imputation_;
  model y = age female;
  output out=work.mi_pred p=fitted;
run; quit;

ods output ParameterEstimates=work.pooled_est TCov=work.pooled_cov;
proc mianalyze data=work.mi_outest tcov;
  modeleffects Intercept age female;
run;

/* The balancing score is the per-patient fitted value averaged over imputations. */
proc summary data=work.mi_pred nway;
  class id;
  var fitted;
  output out=work.mi_scores(drop=_type_ _freq_) mean=fitted;
run;

data work.sas_environment;
  length product $40 version $40 platform $40 run_at $25;
  product = "SAS"; version = "&sysvlong"; platform = "&sysscpl";
  run_at = put(datetime(), e8601dt19.);
run;

%macro export(ds, file);
  proc export data=work.&ds outfile="&root/out/&file" dbms=csv replace; run;
%mend;
%export(single_est, single-estimates.csv)
%export(single_cov, single-covariance.csv)
%export(single_scores, single-scores.csv)
%export(pooled_est, pooled-estimates.csv)
%export(pooled_cov, pooled-covariance.csv)
%export(mi_scores, mi-scores.csv)
%export(sas_environment, sas-environment.csv)

%put NOTE: RM_ORACLE_COMPLETE expected_outputs=7.;
