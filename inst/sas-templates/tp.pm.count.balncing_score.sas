%macro skipp;
//JOB pm.rbctot.balncing_score.sas
sas9.2 pm.rbctot.balncing_score.sas
//cp  splfile pm.rbctot.balncing_score.l*\
//  /studies/<study>/analyses/.
//spool cont printer 'lptxt -l115 -s6 -f LetterGothic-Bold'
spool cont to email
splfile pm.rbctot.balncing_score.l*
// FILE pm.rbctot.balncing_score.sas
%mend;
*______________________________________________________________________________;
*                                                                              ;
* /analyses/gm.rbctot.balcing_score.sas
*______________________________________________________________________________;
*                                                                              ;
  %let STUDY=/studies/<study>;
  %let STUDYDESC = xxxxxxxxxxxxxxxxxxxxxxxx;
  %let STUDYPOP  = CCF, 19xx to 20xx, n=xxxx;    
*______________________________________________________________________________;
*                                                                              ;
* <study description>
* (CCF, 19xx to 20xx, n=xxxx)
* Poisson regression analysis on total number of intrasop RBC units
* Predictors of number of units of intra-op RBC trasfusion
* Estimating Balacing score based on number of units of post-op RBC;
*******************************************************************************;  
*   P O I S S O N   R E G R E S S I O N;
*______________________________________________________________________________;
*                                                                              ;
*                                                                              ;
  options pagesize=100 linesize=132 ;
  libname library "&STUDY/datasets";
  libname est     "&STUDY/estimates";
  data built; set library.built_multimpt;run; 
  title1 "&STUDYDESC";
  title2 "&STUDYPOP";
  title3 "Predictors of number of units of post-op RBC trasfusion using Poisson regression";
  title4 "Estimating Balacing score based on number of units of post-op RBC";
  
*******************************************************************************;
* Data transformations                                                         ;
 filename vars "&STUDY/datasets/vars.sas";
  %inc vars; %vars(transf=0);
   filename vars "&STUDY/datasets/vars.sas";
   data built; set built;
   nyha1=(nyha_pr<1.5);
   nyha2=(1.5<=nyha_pr<2.5);
   nyha3=(2.5<=nyha_pr<3.5);
   nyha4=(nyha_pr>=3.5);
   
   proc freq; tables nyha_pr nyha1 nyha2 nyha3 nyha4; run;
*______________________________________________________________________________;
*                                                                              ;
* Because PROC GENMOD does not include a system for variable selection, use    ;
* PROC REG for poisson response to do variable screening of correlates.        ;  
*______________________________________________________________________________;
*                                                                              ;
*______________________________________________________________________________;
*                                                                              ;
*              F I N A L   M U L T I V A R I A B L  E   M O D E L              ;
*                        Model checking  With one dataset                      ;
*______________________________________________________________________________;

*******************************************************************************;
* Since the poisson fit has over dispersion we use a negative binomial;
* distribution to fit the data;
*******************************************************************************;
%macro modelcheck;
title5 "FINAL Model: RBC";
 ods trace on; 
  proc genmod data=built(where=(_imputation_=3));
       model rbc_tot = agee  female in_bmi
                       in_hct cpb2 hx_csurg emgsrg
                       iso_valv iso_oth afib_pr hx_cva 
                       surgy_05 surgy_06 surgy_07 surgy_08 surgy_09
      /dist = nb link=log type3; 
       output out=predrbc p=pred   RESDEV=resdev l=cll u=clu xbeta=xbeta
                 stdxbeta=std hesswgt=hesswgt reschi=reschi 
                 stdresdev=stresdev stdreschi=streschi ;
       
  run;
 
*******************************************************************************;
*  D I A G N O S T I C   P L O T S;
*******************************************************************************;
* Print estimates                                                              ;
  title4 "Individual Estimates";
  data predrbc1; set predrbc; 
     keep ccfid dt_surg rbc_tot pred resdev cll clu xbeta ;
   run;
  proc print data=predrbc1 (firstobs=1 obs=20);
  
  title4 "Predicted value Vs Actual Value";
  proc plot data=predrbc;
       plot pred*rbc_tot /overlay;
  
  title4 "Deviance Residual Vs Actual Value";    
       plot resdev*rbc_tot;
  
  title4 "Deviance Residual Vs Predicted  Value";    
       plot resdev*pred;     
  run;
*******************************************************************************;
data residuals; set predrbc; 
data residuals;     set residuals;
  h=hesswgt*std*2; /* diagonal element of hession matrix */
  cookd=(streschi**2)*h/((1-h)*22); /* Cook s D statistics */
  adjpred=2*sqrt(pred); /* Adjusted predicted Value  */
  adjlinp=xbeta+(rbc_tot - pred)/pred; /* Adjusted linear predictor */
  absres=abs(stresdev); /* Absolute value of standardized deviance residuals */
run;

*******************************************************************************;
* suggested diagnostic plots;
*******************************************************************************;
title5 "Standardized residuals Vs. The predicted Values";
proc plot data=residuals;
       plot stresdev*pred;

title5 "Standardized residuals Vs. Adjusted predicted Values";       
       plot stresdev*adjpred;
       
title5 "Absolute of Standardized residuals Vs. Adjusted predicted Values";       
       plot absres*adjpred;

title5 "Adjusted linear predictor Vs. Estimated linear predictor";       
       plot adjlinp*xbeta;       
  
title5 "Adjusted linear predictor Vs. Estimated linear predictor";       
       plot adjlinp*xbeta;       

title5 "Cook's D statistics Vs. Estimated predicted values";       
       plot cookd*pred; run;
       
  proc print data=residuals; 
       var  ccfid rbc_tot pred cookd; 
       where cookd>1.00 or pred>30;     run;
*******************************************************************************;
%mend modelcheck;

%modelcheck;

*______________________________________________________________________________;
*                                                                              ;
*               F I N A L   M U L T I V A R I A B L  E   M O D E L             ;
*                     P A R S I M O N I O U S    M O D E L                     ;
*                U S I N G   M U L T I P L E   I M P U T A T I O N             ;
*______________________________________________________________________________;
*                                                                              ;

%macro PARSIMONIOUS;
title5 "Parsimonious Model: RBC";
proc sort data=built; by _imputation_;run;

 ods output parameterestimates=peML_mul;
 ods output CovB=pecov_mul; 
  proc genmod data=built; by  _imputation_;
       model rbc_tot = agee  female in_bmi
                       in_hct cpb2 hx_csurg emgsrg
                       iso_valv iso_oth afib_pr hx_cva 
                       surgy_05 surgy_06 surgy_07 surgy_08 surgy_09
      /dist = nb link=log type3 covb; 
       
       
  run;
%mend PARSIMONIOUS; 

%PARSIMONIOUS;
*******************************************************************************;
* Combine the estimates using PROC MIANALYZE;

%macro combine;
 
title5 "Combined estimates";
 ods output  ParameterEstimates=parms;
 PROC MIANALYZE parms=peML_mul ;
      modeleffects agee  female in_bmi
                       in_hct cpb2 hx_csurg emgsrg
                       iso_valv iso_oth afib_pr hx_cva 
                       surgy_05 surgy_06 surgy_07 surgy_08 surgy_09;                                  
 run;

data est.rbcparms; set parms; run;
 
%mend combine;

%combine;
*______________________________________________________________________________;
*                                                                              ;
*               F I N A L   M U L T I V A R I A B L  E   M O D E L             ;
*                        S A T U R A T E D   M O D E L                         ;
*                    B A L A N C I N G    SC O R E    M O D E L                ;
*______________________________________________________________________________;
*                                                                              ;

%macro saturated;
title5 "Balancing Score Model: intraop RBC";
proc sort data=built; by _imputation_;
 ods trace on; 
  proc genmod data=built; by _imputation_;
       model rbc_tot = 
         /* Demographic */
       female agee in_bmi race_wh 

 
  /* Symptoms */
       nyha2 nyha3 nyha4  emgsrg 

  /* Ventricular function */
        hx_mi  hdef   

  /* Cardiac comorbidity */
       afib_pr  chb_pr  varr_pr  hx_csurg  hx_endo  hx_chf 
       iabp_pr  CarShock  

  /* Non-cardiac comorbidity */
       hx_smoke  hx_pad   hx_cardz   hx_copd  hx_htn 
       hx_iddm   hx_dmtrt  hx_rnldz  hx_cva 
       ln_chol in_hct bun_pr gfr_pr in_blrbn 
     
  /* Coronary anatomy */
      lmt50 lad50 lcx50 rca50 

  /* Experience */
        surgy_05 surgy_06 surgy_07 surgy_08 surgy_09
       
   /* Procedure */
     ita_none  ita_sin
     iso_cabg  iso_valv  iso_vlcg 

  /* Support */
      iv_cpb           
      /dist = nb link=log;       
       output out=predrbc p=pred  xbeta=xbeta;
       
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
* Keep  estimated XBETA and ccfid                               ;
******************************************************************************;
                                             
data _predrbc; set predrbc;
keep ccfid xbeta _IMPUTATION_; run;                                             

proc sort data=_predrbc; by _IMPUTATION_ ccfid; run;
 
*****************************************************************************;     
* Combine  the balancing score estimates (linear predictors) (average of m complete datasets);
*****************************************************************************;

proc sort data=_predrbc; by ccfid;

proc summary data=_predrbc; class ccfid; var xbeta;
       output out=_rbc_balnc mean=prd_xbeta; run; 

data _rbc_balnc; set _rbc_balnc;
 if ccfid=. then delete; 
*_________________________________________________________________________________;
*                                                                                 ;
**      Create Strata using the the Balancing Score                               ;
*_________________________________________________________________________________;
** making strata using the balancing score: 10 strata                             ;
 
  proc sort data=_rbc_balnc; by prd_xbeta;
  data _rbc_balnc; set _rbc_balnc nobs=num;
  _nobs_=num;
  quint=_nobs_/5;
  quintile=int(_n_/quint) +1;
    if quintile>5 then quintile=5;
  desci=_nobs_/10;
  decile=int(_n_/desci) +1;
    if decile>10 then decile=10;  

data _rbc_balnc; set _rbc_balnc;
cluster=decile;

data _rbc_balnc; set  _rbc_balnc;
keep ccfid cluster prd_xbeta decile quintile; run;


*********************************************************************************;
* Create binary variables for each strata: for future use                        ;
*********************************************************************************;

%macro grp(in=_rbc_balnc, out=_rbc_balnc, _strata=cluster, ngrps=10);
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
   create table rbcpred as
   select * from built left join _rbc_balnc
   on built.ccfid =_rbc_balnc.ccfid;

**********************************************************************************;
* keep only the variables in my saturated model;
data built1; set rbcpred; keep rbc_tot cluster
        female agee in_bmi race_wh 
       nyha_pr  emgsrg 
        hx_mi  hdef   
       afib_pr  chb_pr  varr_pr  hx_csurg  hx_endo  hx_chf 
       iabp_pr  CarShock  
       hx_smoke  hx_pad   hx_cardz   hx_copd  hx_htn 
       hx_iddm   hx_dmtrt  hx_rnldz  hx_cva 
       ln_chol in_hct bun_pr gfr_pr in_blrbn 
      lmt50 lad50 lcx50 rca50 
      iv_opyrs ita_none iso_cabg  iso_valv  iso_vlcg iv_cpb; run;

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
   ods output  ParameterEstimates=parest;
   proc sort data=&_indat; by cluster;
    proc genmod data=&_indat; by cluster;
            model  &_resp = &item /dist = nb link=log;   
            run;

 data parest; set parest; keep cluster Parameter ProbChiSq Estimate;

 data parest; set parest; where Parameter="&item"; run;  
 data parestP; set parest; keep cluster ProbChiSq; rename ProbChiSq=&item;
 data parestE; set parest; keep cluster Estimate; rename Estimate=&item;

/***********************************/

data _outest; merge _outest parestE; by cluster;
data _outpval; merge _outpval parestP; by cluster; 
 
            
 %let i=%eval(&i+1);
 %let item=%scan(&_modvar, &i);
%end;

%mend blc_ckmod;

%blc_ckmod(_indat=built1, _resp=rbc_tot);



title6 "Covariate balance in RBC total: Parameter Estimates";
proc print data=_outest; proc means data=_outest; run;
title6 "Covariate balance in RBC total: P-values";
proc print data=_outpval; run; proc means data=_outpval; run;


 %mend blnc_checking;
********************************************************************************;
%blnc_checking;

**********************************************************************************;
**Save the balancing score the strata variables to a final dataset;
**********************************************************************************;
%macro final;
data est.rbc_blnc_fin; set _rbc_balnc; run;
%mend final;
%final;