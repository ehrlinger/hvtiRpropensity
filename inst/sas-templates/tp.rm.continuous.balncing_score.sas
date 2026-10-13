%macro skkip;
// JOB rm.hct_nadr.balncing_score.sas
sas9.2 rm.hct_nadr.balncing_score.sas
//cp rm.hct_nadr.balncing_score.l*\
//  /studies/<study>/analyses/.
//spool cont printer 'lptxt -l110 -s6 -f LetterGothic-Bold'
spool cont to email
splfile rm.hct_nadr.balncing_score.l*
// FILE rm.hct_nadr.balncing_score.sas
%mend;

*______________________________________________________________________________;
*                                                                              ;
* /analyses/rm.hct_nadr.balncing_score.sas
*______________________________________________________________________________;
*                                                                              ;
  %let STUDY=/studies/<study>;
  %let STUDYDESC = xxxxxxxxxxxxxxxxxxxxxxxx;
  %let STUDYPOP  = CCF, 19xx to 20xx, n=xxxx;   
*______________________________________________________________________________;
*
* <study description>
* (CCF, 19xx to 20xx, n=xxxx) 
*
* Multivariable analysis of nadir HCT
* Linear regression on  nadir HCT
** Estimating Balacing score based on Nadir HCT;
*______________________________________________________________________________;
*                                                                              ;
  options pagesize=107 linesize=132;
  libname library "&STUDY/datasets";
  libname est     "&STUDY/estimates";
  
   
  data built; set library.built_multimpt; run;
  
  title1 "&STUDYDESC";
  title2 "&STUDYPOP";
  title3 "Multivariable Analysis of LOS operative";
  title4 "Estimating Balacing score based on Nadir HCT";
  
*******************************************************************************;
* Data transformations                                                         ;
/*
  filename vars "&STUDY/datasets/vars.sas";
  %inc vars; %vars;
*/ 
  
*______________________________________________________________________________;
*                                                                              ;
*                M U L T I V A R I A B L E   A N A L Y S I S                   ;
*______________________________________________________________________________;
*                                                                              ;
%macro model; 
 proc reg data=built (where=(_imputation_=2));
       model hct_nadr=
       
 /* exposure groups */
           group_1 group_2 group_3 group_4 
       
  /* Demographic */
       female  male  age  ln_age  in_age  age2  agee  
        ht  ln_ht  wt  wt2  in_wt  ln_wt 
        bmi  bmi2  in_bmi  ln_bmi  race_bl  race_wh  

 
  /* Symptoms */
       nyha_pr  ln_nyha  in_nyha   emgsrg 

  /* Ventricular function */
        hx_mi  hdef  ln_hdef  in_hdef  hdef2 
      

  /* Cardiac comorbidity */
       afib_pr  chb_pr  varr_pr  hx_csurg  hx_endo  hx_chf 
       iabp_pr  CarShock  

  /* Non-cardiac comorbidity */
       hx_smoke  hx_pad   hx_cardz   hx_copd  hx_htn 
       hx_iddm  hx_niddm  hx_dmtrt  hx_rnldz  hx_cva 
       chol_pr  ln_chol  chol2  in_chol  
       hct_pr  ln_hct  in_hct  hct2 
       bun_pr  ln_bun  in_bun  bun2 
       gfr_pr  ln_gfr  in_gfr  gfr2 
       blrbn_pr  ln_blrbn  in_blrbn  blrbn2 
      
  /* Coronary anatomy */
       vd0 vd3 
       lmt70  lmt50  lmtany 
       lad50  lad70  ladany 
       lcx50  lcx70  lcxany 
       rca50  rca70  rcaany 

  /* Experience */
       iv_opyrs  ln_opyrs  in_opyrs  ivopyrs2 
       
   /* Procedure */
     ita_none  ita_sin  ita_bil  ita_any  ita_num  
     iso_cabg  iso_valv  iso_vlcg  iso_oth 

  /* Support */
      iv_aocc  ln_aocc  in_aocc  aocc2  
      iv_cpb  ln_cpb  in_cpb  cpb2 

       /selection=stepwise sle=0.07 sls=0.05 details;
       ; run;
 
  %mend;
*______________________________________________________________________________;
*                                                                              ;
*              F I N A L   M U L T I V A R I A B L  E   M O D E L              ;
*              Model checking  With one dataset
*______________________________________________________________________________;
*                                                                              ;
  %macro check;
   title5 " Parsimonious model";
   
   proc reg data=built(where=(_imputation_=2)) covout outest=outest edf;
       model hct_nadr =afib_pr hx_csurg lmt70 lcx70 in_blrbn 
                       hct_pr iv_cpb female agee in_bmi
                       iso_valv;
       output out=regout residual=resid predicted=predict;
    run;
  
data _outest1;
    set outest; keep _TYPE_ _EDF_ _DEPVAR_ _RMSE_; run; 
    
data _outest1; set _outest1;  noob=_n_; run;    
 
*______________________________________________________________________________;
*                                                                              ;
*                       G O O D N E S S   O F   F I T						   ;
*______________________________________________________________________________;
*                                                                              ;
  title5 "Goodness of Fit";
  proc univariate data=regout plot normal;
       var resid; run;
  proc plot data=regout;
  plot predict*hct_nadr resid*hct_nadr resid*predict;
  run;
  %mend check;
*******************************************************************************;
%check;

*______________________________________________________________________________;
*                                                                              ;
*               F I N A L   M U L T I V A R I A B L  E   M O D E L             ;
*                     P A R S I M O N I O U S    M O D E L                     ;
*                U S I N G   M U L T I P L E   I M P U T A T I O N             ;
*______________________________________________________________________________;
*                                                                              ;
  %macro parsimonious;
  title5 "Regression on 5 completed datasets";
  proc sort data=built; by _imputation_;run;
  proc reg  data=built covout outest=outest noprint;
       model  hct_nadr =afib_pr hx_csurg lmt70 lcx70 in_blrbn 
                       hct_pr iv_cpb female agee in_bmi
                       iso_valv;
       by _Imputation_;
    run;
 %mend  parsimonious;
 %parsimonious; 
*******************************************************************************;
*******************************************************************************;
* Combine the estimates using PROC MIANALYZE;
* since there are lots of variables in the model;
* we use the following macro to extract the variable list;
*NOTE: input the response variable;

%macro combine(_resp=hct_nadr);
 
  %local _modvar;
   data _outest_; set outest;
   drop  &_resp _MODEL_ _NAME_ _RMSE_ _Imputation_ _DEPVAR_;
   
   
 proc SQL noprint;
    select name into: _modvar separated by ' '
        from sashelp.vcolumn
                where   libname="WORK" AND 
                        memname=%upcase("_outest_") AND
                        upcase(memtype)='DATA' AND
                        (substr(name,1,1) ~='_' 
                            OR substr(name,length(name),1)~='_');                        
  quit;
 
%put &_modvar;

title5 "Combined estimates";
 ods output  ParameterEstimates=parms;
 ods output TCov=tcov;
 PROC MIANALYZE data=outest edf=9107  MULT tcov ;
      modeleffects &_modvar ;                                  
 run;
 
 %mend combine;

%combine;

*******************************************************************************;
%macro save_est;
*******************************************************************************;
* Rearrange parameter estimates and variance-covariance matrix suitable for    ;
* Linear regression plot macro;
*******************************************************************************;
* Rearrange parameter estimates;
*******************************************************************************;
 * proc print data=parms;
  data parms1; set parms;
     keep parm estimate;
  proc transpose data= parms1 out=parms1;
     id parm;   
  run;
  
  data parms1; set parms1;
  _TYPE_="PARMS";
  drop _NAME_; 
 proc print data=parms1;run;
*******************************************************************************; 
* Rearrange the variance-covariance matrix of parameter estimates              ;
*******************************************************************************;
*proc print data=tcov;run;
data covout1; set tcov;
  _TYPE_="Cov";
  _NAME_=Variable;
  drop variable;
proc print data=covout1;run;
*******************************************************************************;
* Merge all the estimates;
*******************************************************************************;
data outesti1; set parms1 covout1;
noob=_n_;run;

proc sort data=outesti1; by noob;
proc sort data=_outest1; by noob;

data outesti; merge outesti1 _outest1; by noob;run;

data est.hctnr_par; set outesti; drop noob; run; 
proc print data=est.hctnr_par;run;
  
%mend save_est;
*******************************************************************************;
*******************************************************************************;
%save_est; 

*______________________________________________________________________________;
*                                                                              ;
*               F I N A L   M U L T I V A R I A B L  E   M O D E L             ;
*                        S A T U R A T E D   M O D E L                         ;
*                    B A L A N C I N G    SC O R E    M O D E L                ;
*______________________________________________________________________________;
*                                                                              ;
%macro saturated;

  proc sort data=built; by _IMPUTATION_; run;
  proc reg  data=built ;
  
        by _IMPUTATION_;
       model hct_nadr= 
      
      /* Demographic */
       female  agee  in_bmi race_wh 

 
  /* Symptoms */
       nyha_pr emgsrg 

  /* Ventricular function */
        hx_mi  hdef2       

  /* Cardiac comorbidity */
       afib_pr  chb_pr  varr_pr  hx_csurg  hx_endo  hx_chf 
       iabp_pr  CarShock  

  /* Non-cardiac comorbidity */
       hx_smoke  hx_pad   hx_cardz   hx_copd  hx_htn 
       hx_iddm  hx_niddm  hx_dmtrt  hx_rnldz  hx_cva 
       ln_chol hct_pr bun_pr BUN2 ln_gfr  in_blrbn 
     
  /* Coronary anatomy */
       lmt50 lad50  lcx50 rca50 

  /* Experience */
       surgy_05  surgy_06 surgy_07  surgy_08  surgy_09
       
   /* Procedure */
     ita_none  ita_sin  
     iso_cabg  iso_valv  iso_vlcg 

  /* Support */
       iv_cpb;
       
       output out=_hctpred predicted=predict;       
       
      
      run;
      
%mend saturated;
******************************************************************************;
%saturated;
*_____________________________________________________________________________;
*                                                                             ;
*                    B A L A N C I N G   S C O R E                            ;
*_____________________________________________________________________________;              

%MACRO BLC_SCR;     
******************************************************************************;
* Keep  estimated CONTINUOUS RESPONSE and ccfid                               ;
******************************************************************************;
                                             
data _hctpred; set _hctpred;
keep ccfid predict _IMPUTATION_; run;                                             

proc sort data=_hctpred; by _IMPUTATION_ ccfid; run;
 
*****************************************************************************;     
* Combine  the balancing score estimates (average of m complete datasets);
*****************************************************************************;

proc sort data=_hctpred; by ccfid;

proc summary data=_hctpred; class ccfid; var predict;
       output out=_hct_balnc mean=prd_hctn; run;

data _hct_balnc; set _hct_balnc;
 if ccfid=. then delete;  run;       
*_________________________________________________________________________________;
*                                                                                 ;
**      Create Strata using the the Balancing Score                               ;
*_________________________________________________________________________________;
** making strata using the balancing score: 10 strata                             ;
 
  proc sort data=_hct_balnc; by prd_hctn;
  data _hct_balnc; set _hct_balnc nobs=num;
  _nobs_=num;
  quint=_nobs_/5;
  quintile=int(_n_/quint) +1;
    if quintile>5 then quintile=5;
  desci=_nobs_/10;
  decile=int(_n_/desci) +1;
    if decile>10 then decile=10;  

data _hct_balnc; set _hct_balnc;
cluster=decile;

data _hct_balnc; set  _hct_balnc;
keep ccfid cluster prd_hctn decile quintile; run;


*********************************************************************************;
* Create binary variables for each strata: for future use                        ;
*********************************************************************************;

%macro grp(in=_hct_balnc, out=_hct_balnc, _strata=cluster, ngrps=10);
data &out; set &in;

  %do i=1 %to &ngrps;
      stra_&i=0;
      if &_strata=&i then stra_&i=1;
  %end;
%mend;

%grp;   

%mend blc_scr;

%blc_scr; 

*_________________________________________________________________________________;
*                                                                                 ;
**      Check the goodness of the Balancing Score                                 ;
*_________________________________________________________________________________;
** Checking to see if the covariates are balanced with respect to the             ;
*  exposure variable  (nadr hct, in this case)                                    ;
**********************************************************************************;
%macro blnc_checking;
**********************************************************************************;
* joint the strata variable with patient data;
data built; set library.built_multimpt(where=(_imputation_=3)); run;     
proc sql;
   create table hctpred as
   select * from built left join _hct_balnc
   on built.ccfid =_hct_balnc.ccfid;

**********************************************************************************;
* keep only the variables in my saturated model;
data built1; set hctpred; keep hct_nadr cluster
       female  agee  in_bmi race_wh 
       nyha_pr emgsrg hx_mi  hdef2       
       afib_pr  chb_pr  varr_pr  hx_csurg  hx_endo  hx_chf 
       iabp_pr  CarShock  
       hx_smoke  hx_pad   hx_cardz   hx_copd  hx_htn 
       hx_iddm  hx_niddm  hx_dmtrt  hx_rnldz  hx_cva 
       ln_chol hct_pr bun_pr BUN2 ln_gfr  in_blrbn 
       lmt50 lad50  lcx50 rca50 
       iv_opyrs
       /*surgy_05  surgy_06 surgy_07  surgy_08  surgy_09 */
       ita_none iso_cabg  iso_valv  iso_vlcg iv_cpb; run;

**********************************************************************************;
**Initialize an output datasets: for estimates and p-value                        ;
**********************************************************************************; 
data built2; set built1; keep cluster;
proc sort data=built2 nodupkey;
          by cluster;run; 
          
**Initialize estimate dataset;
data _outest; set built2; keep cluster;

**Initialize Pvalue dataset;
data _outpval;set built2; keep cluster; run;  

********************************************************************************;
** A macro for balance checking;
********************************************************************************;
* For each covariate, This macro use a univariate linear regression (PROC REG)  ;
* INPUT: dataset with response variable and the covariate                       ;
*        Name of the response variable                                          ;                            
********************************************************************************;
%macro blc_ckmod(_indat=_indat, _resp=_resp);
 
 data _indat_; set &_indat;
   drop  &_resp cluster;
   
 %local _modvar;

 proc SQL noprint;
    select name into: _modvar separated by ' '
        from sashelp.vcolumn
                where   libname="WORK" AND 
                        memname=%upcase("_indat_") AND
                        upcase(memtype)='DATA' AND
                        (substr(name,1,1) ~='_' 
                            OR substr(name,length(name),1)~='_');                        
  quit;
 
%put &_modvar;



%let i=1;
%let item=%scan(&_modvar, &i);

%do %until (&item.=);

/**********************************/
   ods output "Parameter Estimates"=parest;
   proc sort data=&_indat; by cluster;
   proc reg data=&_indat; by cluster;
            model  &_resp = &item ;
            run;

 data parest; set parest; keep cluster variable Probt Estimate;

 data parest; set parest; where variable="&item"; run;  
 data parestP; set parest; keep cluster Probt; rename Probt=&item;
 data parestE; set parest; keep cluster Estimate; rename Estimate=&item;

/***********************************/

data _outest; merge _outest parestE; by cluster;
data _outpval; merge _outpval parestP; by cluster; 
 
            
 %let i=%eval(&i+1);
 %let item=%scan(&_modvar, &i);
%end;

%mend;

*********************************************************************************;
%blc_ckmod(_indat=built1, _resp=hct_nadr);


title6 "Covariate balance in HCT nadir model: Parameter Estimates";
proc print data=_outest; proc means data=_outest; run;
title6 "Covariate balance in HCT nadir model: P-values";
proc print data=_outpval; run; proc means data=_outpval; run;

 %mend blnc_checking;
********************************************************************************;
%blnc_checking;

**********************************************************************************;
**Save the balancing score the strata variables to a final dataset;
**********************************************************************************;
%macro final;
data est.hct_blnc_fin; set _hct_balnc; run;
%mend final;
%final;