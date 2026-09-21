/*
 * Synthetic model oracles for hvtiRpropensity.
 *
 * This program uses only the committed deterministic input.csv. It exports
 * raw ODS and scored tables so no R implementation decides the oracle's
 * parameter names, ordering, covariance layout, or predicted probabilities.
 */

/* EDIT: absolute server path to the staged lm-oracle folder. */
%let root=;
%let input_dir=&root./input;
%let output_dir=&root./out;

filename orcllog "&output_dir./generate-fixtures.log";
proc printto log=orcllog new;
run;

%put NOTE: LM oracle root is &root.;
%put NOTE: SAS version is &sysvlong4.;
%put NOTE: Platform is &sysscpl.;

proc import datafile="&input_dir./input.csv"
    out=fixture dbms=csv replace;
  guessingrows=max;
  getnames=yes;
run;

proc format;
  value ordinal_order_format
    1='low'
    2='mid'
    3='high';
run;

data fixture;
  set fixture;
  length ordinal_order 8;
  select (ordinal);
    when ('low')  ordinal_order=1;
    when ('mid')  ordinal_order=2;
    when ('high') ordinal_order=3;
    otherwise call missing(ordinal_order);
  end;
  format ordinal_order ordinal_order_format.;
run;

data one;
  set fixture;
  where imputation=1;
run;

proc sort data=fixture;
  by imputation;
run;

ods exclude all;

/* Binary logistic outcome model. */
ods output ParameterEstimates=binary_estimates CovB=binary_covariance;
proc logistic data=one order=internal;
  model binary(event='1') = x1_imp x2 / covb;
  output out=binary_scores p=probability xbeta=linear_predictor;
run;

/* Proportional-odds cumulative logistic outcome model. */
ods output ParameterEstimates=ordinal_estimates CovB=ordinal_covariance;
proc logistic data=one order=internal;
  model ordinal_order = x1_imp x2 / covb;
  output out=ordinal_scores predprobs=individual;
run;

/* Nominal generalized-logit outcome model, reference level a. */
ods output ParameterEstimates=nominal_estimates CovB=nominal_covariance;
proc logistic data=one order=internal;
  class nominal(ref='a') / param=ref;
  model nominal(ref='a') = x1_imp x2 / link=glogit covb;
  output out=nominal_scores predprobs=individual;
run;

/* Poisson count model with log-link balancing score. */
ods output ParameterEstimates=count_estimates CovB=count_covariance;
proc genmod data=one;
  model count = x1_imp x2 / dist=poisson link=log covb;
  output out=count_scores pred=predicted_mean xbeta=linear_predictor;
run;

/* Binary logistic fits by imputation plus legacy-compatible OUTEST input to
 * PROC MIANALYZE. */
ods output ParameterEstimates=mi_estimates CovB=mi_covariance;
proc logistic data=fixture order=internal covout outest=mi_outest;
  by imputation;
  model binary(event='1') = x1_imp x2 / covb;
  output out=mi_scores p=probability xbeta=linear_predictor;
run;

data mi_estimates_mianalyze;
  set mi_estimates;
  rename imputation=_Imputation_;
run;

data mi_covariance_mianalyze;
  set mi_covariance;
  rename imputation=_Imputation_;
run;

ods output ParameterEstimates=pooled_estimates TCov=pooled_covariance;
proc mianalyze parms=mi_estimates_mianalyze
    covb(effectvar=stacking)=mi_covariance_mianalyze tcov mult;
  modeleffects Intercept x1_imp x2;
run;

ods exclude none;

data sas_environment;
  length product_version platform $200 run_datetime $32;
  product_version="&sysvlong4";
  platform="&sysscpl";
  run_datetime=put(datetime(), e8601dt.);
run;

%macro export_raw(data=, file=);
  proc export data=&data
      outfile="&output_dir./&file..csv" dbms=csv replace;
  run;
%mend;

%export_raw(data=binary_estimates, file=binary-estimates);
%export_raw(data=binary_covariance, file=binary-covariance);
%export_raw(data=binary_scores, file=binary-scores);

%export_raw(data=ordinal_estimates, file=ordinal-estimates);
%export_raw(data=ordinal_covariance, file=ordinal-covariance);
%export_raw(data=ordinal_scores, file=ordinal-scores);

%export_raw(data=nominal_estimates, file=nominal-estimates);
%export_raw(data=nominal_covariance, file=nominal-covariance);
%export_raw(data=nominal_scores, file=nominal-scores);

%export_raw(data=count_estimates, file=count-estimates);
%export_raw(data=count_covariance, file=count-covariance);
%export_raw(data=count_scores, file=count-scores);

%export_raw(data=mi_estimates, file=mi-estimates);
%export_raw(data=mi_covariance, file=mi-covariance);
%export_raw(data=mi_outest, file=mi-outest);
%export_raw(data=mi_scores, file=mi-scores);

%export_raw(data=pooled_estimates, file=pooled-estimates);
%export_raw(data=pooled_covariance, file=pooled-covariance);
%export_raw(data=sas_environment, file=sas-environment);

%put NOTE: LM_ORACLE_COMPLETE expected_outputs=19.;

proc printto;
run;
