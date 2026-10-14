/* rm SAS oracle: PROC REG on one dataset, and per imputation pooled by PROC MIANALYZE.
   Set only ROOT. Input is input/input.csv; output goes to out/. */
%let root=SET_ME;

proc import datafile="&root/input/input.csv" out=work.input dbms=csv replace;
  guessingrows=max;
run;

proc sort data=work.input; by _Imputation_ id; run;

/* One dataset: imputation 1 only. */
ods output ParameterEstimates=work.single_est CovB=work.single_cov ANOVA=work.single_anova;
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

/* Barnard-Rubin (1999) degrees of freedom. EDF= is the complete-data df,
   n - k: the Error DF of imputation 1's fit, the value R takes from
   df.residual() of the first imputation's model. */
proc sql noprint;
  select DF into :edf trimmed from work.single_anova where Source = "Error";
quit;
%if %length(&edf) = 0 %then %do;
  %put ERROR: The complete-data degrees of freedom were not computed.;
  %abort cancel;
%end;
%put NOTE: Complete-data degrees of freedom EDF=&edf..;

data work.complete_df;
  edf = &edf;
run;

ods output ParameterEstimates=work.pooled_est_edf;
proc mianalyze data=work.mi_outest edf=&edf;
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

/* PROC EXPORT writes formatted values, so strip the display formats first:
   an oracle needs every digit, and a p-value printed as <.0001 is not a number. */
%macro export(ds, file);
  data work._export;
    set work.&ds;
    format _numeric_ best32.;
  run;
  proc export data=work._export outfile="&root/out/&file" dbms=csv replace; run;
%mend;
%export(single_est, single-estimates.csv)
%export(single_cov, single-covariance.csv)
%export(single_anova, single-anova.csv)
%export(single_scores, single-scores.csv)
%export(pooled_est, pooled-estimates.csv)
%export(pooled_cov, pooled-covariance.csv)
%export(pooled_est_edf, pooled-estimates-edf.csv)
%export(complete_df, complete-df.csv)
%export(mi_scores, mi-scores.csv)
%export(sas_environment, sas-environment.csv)

%put NOTE: RM_ORACLE_COMPLETE expected_outputs=10.;
