%macro skip;
// JOB lm.mip.sas
sas9.1 lm.mip.sas
//cp lm.mip.l*\
//   /studies/<study>/analysis/.
//spool cont printer 'lptxt -l110 -s6 -f LetterGothic-Bold'
spool cont to email
splfile lm.mip.l*
// FILE lm.mip.sas 
%mend;
*______________________________________________________________________________;
*                                                                              ;
* /analyses/lm.mip.sas
*______________________________________________________________________________;
*                                                                              ;
  %let STUDY=/studies/<study>;
  %let STUDYDESC = xxxxxxxxxxxxxxxxxxxxxxxx;
  %let STUDYPOP  = CCF, 19xx to 20xx, n=xxxx;   
*______________________________________________________________________________;
*                                                                              ;
* <study description>
* (CCF, 19xx to 20xx, n=xxxx)
*
* Multivariable logistic regression analysis of mitral valve mini vs. full surgery
* Parsimoneous and propensity models on mitralprocedure alone patients
*______________________________________________________________________________;
*                                                                              ;
  options pagesize=107 linesize=132 pageno=1;
  libname library "&STUDY/datasets";
  libname est v8  "&STUDY/estimates";
  filename deciles "!MACROS/deciles.new"; %inc deciles;
  
  title1 "&STUDYDESC";
  title2 "&STUDYPOP";
  title3 "Multivariable Analysis of <study description>";
*******************************************************************************;
*******************************************************************************;
* read in the imputed complete datasets;
data built; set library.built_mip_imput; run;
*******************************************************************************;
* Data Transformation                                                          ;                                                      
  filename vars "&STUDY/datasets/vars.sas";
  %inc vars; %vars; 
  

*______________________________________________________________________________;
*                                                                              ;
*       M U L T I V A R I A B L E   L O G I S T I C   R E G R E S S I O N      ;
*______________________________________________________________________________;
*                                                                              ;
* Use one complete data set for initial stepwise modeling;
%macro select;
  title4 "Selection Model";
  proc logistic data=built(where=(_IMPUTATION_=1)) covout outest=outest descending;
       output out=decile p=_y_;
   model mip=
  
  /* Demographic */
  /* Symptoms */
  /* Ventricular function */
  /* Cardiac comorbidity */
  /* Non-cardiac comorbidity */
  /* Coronary anatomy */
  /* Experience */
     
      /selection=stepwise details maxstep=15; run;   

%mend; 
*______________________________________________________________________________;
*                                                                              ;
*                G O O D N E S S    O F    F I T                               ; 
*     Check goodness of fit of model to the data (Hosmer-Lemeshow statistics)  ;
*______________________________________________________________________________;
* Use one complete data set to check the goodness of fit;
%macro model_fit;
 proc logistic data=built(where=(_IMPUTATION_=1)) descending ;
        output out=decile p=phat;
       model mip = ln_opyrs ave_dege ave_rheu sgn_jose sgn_gost sgn_mcca lvfc2 tvrgsev
                  bun_pr wt2 ;run;

  %deciles(in=decile, _event_=mip,_p_=phat,
           vars=ln_opyrs ave_dege ave_rheu sgn_jose sgn_gost lvfc2 tvrgsev
                  bun_pr wt2);
%mend;
%model_fit;
*______________________________________________________________________________;
*                                                                              ;
*               F I N A L   M U L T I V A R I A B L  E   M O D E L             ;
*                     P A R S I M O N I O U S    M O D E L                     ;
*______________________________________________________________________________;
*                                                                              ;
%macro pars;
  proc sort data=built; by _IMPUTATION_; run;
  proc logistic data=built descending covout outest=outest noprint;
        by _IMPUTATION_;
       model mip =ln_opyrs ave_dege ave_rheu sgn_jose sgn_mcca sgn_gost lvfc2 tvrgsev
                  bun_pr wt2;run;
  
/* C-statistic average of 5 models: 0.73*/       
*******************************************************************************;
* Combine the estimates using PROC MIANALYZE;       
 ods output ParameterEstimates=parms;
 ods output TCov=tcov;
 PROC MIANALYZE data=outest tcov mult;
      modeleffects Intercept ln_opyrs ave_dege ave_rheu sgn_jose sgn_mcca sgn_gost lvfc2 tvrgsev
                  bun_pr wt2;               
 run;

*******************************************************************************;
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
data outesti; set parms1 covout1;

proc print; run;

data est.mip; set outesti; run;
********************************************************************************;
%mend;

%pars;

*______________________________________________________________________________;
*                                                                              ;
*               F I N A L   M U L T I V A R I A B L  E   M O D E L             ;
*                         P R O P E N S I T Y    M O D E L                     ;
*______________________________________________________________________________;
*                                                                              ;
%macro propen;

  proc sort data=built; by _IMPUTATION_; run;
  proc logistic data=built descending covout outest=outest noprint;
        by _IMPUTATION_;
        output out=decile p=_p_;
       model mip = 

  /* Demographic */
  /* Symptoms */
  /* Ventricular function */
  /* Pathology */
  /* Cardiac comorbidity */
  /* Non-cardiac comorbidity */
  /* Coronary anatomy */
  /* Experience */
  /* Interaction terms:  Patient */
  /* Missing value flags */
     
  ;
      
      run;
/* C-STATITICS= average of the c-stat from the 5 models=0.770 */        
*******************************************************************************;
* Combine  the propensity score estimates (average of m complete datasets);

proc sort data=decile; by ccfid;

proc summary data=decile; class ccfid; var _p_;
       output out=_decile1 mean=_propen; run;

data _decile1; set _decile1;
 if ccfid=. then delete; 
 keep ccfid _propen; run;   

data decile; set decile(where=(_IMPUTATION_=1)); drop _p_; run;     
proc sql;
   create table decile as
   select * from decile left join _decile1
   on decile.ccfid = _decile1.ccfid;


*______________________________________________________________________________;
*                                                                              ;
*                   P R O P E N S I T Y   S T A T I S T I C S                  ;
*              C H E C K I N G   P R O P E N S I T Y  M O D E L                ;
*______________________________________________________________________________;
*                                                                              ;
* Form quintiles and deciles                                                   ;
  proc sort data=decile; by _propen;
  data decile; set decile nobs=num;
  _nobs_=num;
  quint=_nobs_/5;
  quintile=int(_n_/quint) +1;
    if quintile>5 then quintile=5;
  dec=_nobs_/10;
  decile=int(_n_/dec) +1;
    if decile>10 then decile=10;
  run;
*******************************************************************************;
* Contingency table statistics                                                 ;
  title5 "Quintile Statistics: Categorical";
  proc freq; tables mip*quintile;
  proc freq data=decile; by quintile; tables mip*(
    /* Demographic */
    /* Symptoms */
    /* Ventricular function */
    /* Pathology */
    /* Cardiac comorbidity */
    /* Non-cardiac comorbidity */
    /* Coronary anatomy */
    /* Experience */
    /* Interaction terms:  Patient */
    /* Missing value flags */
       )/chisq;
  run;
*******************************************************************************;
* T-testing                                                                    ;
  title5 "Quintile Statistics: Continuous";
  proc ttest data=decile; class mip; by quintile; var
     /* Demographic */
     /* Symptoms */
     /* Ventricular function */
     /* Pathology */
     /* Cardiac comorbidity */
     /* Non-cardiac comorbidity */
     /* Coronary anatomy */
     /* Experience */
     /* Interaction terms:  Patient */
     /* Missing value flags */   
  ;
  run;
*______________________________________________________________________________;
*                                                                              ;
*             E S T I M A T E    M A T C H I N G  W E I G H T S                ;
*______________________________________________________________________________;

data decile; set decile;
 mt_wt=min(_p_, (1-_p_))/(_p_*mip + (1-_p_)*(1-mip) );                                    
run;  
*______________________________________________________________________________;
*                                                                              ;
*     S T O R E   P R O P E N S I T Y   S C O R E   F O R   E A C H   O B S    ;
*______________________________________________________________________________;
*                                                                              ;
* Put individual propensity scores to file                                     ;
  data library.propen_mip; set decile;
  _p_=_propen;
  _logit_=log(_p_/(1-_p_));
  prob=_p_;
  
  keep ccfid _p_ _propen  _logit_ mt_wt mip prob quintile decile;
  proc means n mean std min max sum;
  run;

  %mend;
 %propen;