//JOB lm.rtyp.propensity.sas
sas8.2 lm.rtyp.propensity.sas
//spool cont to email
spool cont printer 'lptxt -l110 -s6 -f LetterGothic-Bold'
cp  lm.rtyp.propensity.l*\
   /studies/<study>/analyses/.
splfile lm.rtyp.propensity.l*
// FILE lm.rtyp.propensity.sas
*______________________________________________________________________________;
*                                                                              ;
* /analysis/lm.rtyp.propensity.sas
*______________________________________________________________________________;
*                                                                              ;
  %let STUDY=/studies/<study>;
  %let STUDYDESC = xxxxxxxxxxxxxxxxxxxxxxxx;
  %let STUDYPOP  = CCF, 19xx to 20xx, n=xxxx;    
*______________________________________________________________________________;
*                                                                              ;
* <study description>                                                        ;
* CCF, 19xx to 20xx, n=xxxx                                                    ;
*                                                                              ;
* Multivariable analysis of tricuspid valve reoperation for regurgitation      ;
* Propensity analysis of repair types                                              ;
*______________________________________________________________________________;
*                                                                              ;
  options pagesize=107 linesize=132;
  libname library "&STUDY/datasets";
  libname est     "&STUDY/estimates";
  data built; set library.built; run;
       
  title1 "&STUDYDESC";
  title2 "&STUDYPOP";
  title3 "Propensity analysis of repair types ";
*______________________________________________________________________________;
*                                                                              ;
* Invoke vars to get newly organized variables                                 ;
*______________________________________________________________________________;
*                                                                              ;
  filename vars "&STUDY/datasets/vars.sas";
  %inc vars; %vars(in=built, out=built,  missing=1);   

*NOTE: PROC LOGISTIC is fitting the generalized logit model;
* The logits modeled contrast each response level against the reference;
*Use the response variable option REF= first or last in PROC option;
* if you want to change the reference level.
  
*______________________________________________________________________________;

*                                                                              ;
* Logistic Model: Stepwise Selection: With generalized logit link with         ;
* CE as the reference category                                                 ;
*______________________________________________________________________________;
*                                                                              ;
%macro select;
 title5 'Model for Predictors of Repair types';
 
  proc logistic data=built ref=first;
            model rtyp1= 
  /* Demographic */
       female male age ln_age in_age age2 agee ht ln_ht wt wt2 in_wt
       ln_wt wtht bsa ln_bsa bmi bmi2 in_bmi ln_bmi

  /* Symptoms */
       nyhc ln_nyhc in_nyhc 

  /* Ventricular functon */
       lvfcath ln_lvfc lvf2 lvf34 lvf3 lvf32 lvfnorml ekginf hxmi
       ejfcath ln_ejfc in_ejfc

  /* Pathology and etiology */
       mvrgrg mvrgsev mvsten  ln_mr
       tr_prtte ln_tr rvsp_pre ln_rvspr in_rvspr rvspre2
   

  /* Cardiac comorbidity */
       afib chb ventarr  prev_cs  endo hx_ppm

  /* Non-cardiac comorbidity */
       eversmok  pvd 
       carotid popdz copd htn iddm niddm diabetes rnldz  creat
       ln_creat in_creat creat2 blrbn ln_blrbn in_blrbn blrbn2
       hct ln_hct in_hct hct2 bun bun2 ln_bun in_bun

  /* Coronary anatomy */
       sys_dis sysdis2 vd0 vd1 vd0_1 vd2 vd3 vd2_3
       lmt lmt50 lmtany
       lad lad50 lad70 ladany
       lcx lcx50 lcx70 lcxany
       rca rca50 rca70 rcaany 
 
  /* Experience 
       opyrs ln_opyrs in_opyrs in2opyrs opyrs2 */

 
  /* Interaction terms:  Patient */
       age_diab lagediab iagediab i2agdiab age2diab
      
       
        
 
             /link=glogit  selection=stepwise sle=0.1  sls=0.07  details  ;
                  
                  run;
                  
/*prev_cs in_creat mvsten hct ejfcath tr_prtte lad agee in_nyhc ln_ejfc  */
/* ln_mr age_diab niddm ekginf htn  ladany hx_ppm chb female */

%mend; 
*_____________________________________________________________________________;
*
*   PARSIMONIOUS MODEL                                                        ;
*_____________________________________________________________________________;
%macro parsimonious;

title5 'PARSIMONIOUS MODEL: Model for Predictors of Repair types';
 
  proc logistic data=built outest=outest covout;
            model rtyp1= prev_cs in_creat ln_mr  hct2 tr_prtte 
                         ln_age ln_nyhc  ln_ejfc 
                         ekginf htn lad female hx_ppm
                          / link=glogit ;   run; 
 
data est.lprtyp;
       set outest;
  run;

%mend;  
*_____________________________________________________________________________;
*
*   PROPENSITY   MODEL                                                        ;
*_____________________________________________________________________________;
title5 'PROPENSITY MODEL: ';
 
  proc logistic data=built outest=outest covout;
            output out=probs p=_p_;
            model rtyp1= prev_cs in_creat hct2  
                         ln_age   ln_ejfc 
                         ekginf htn lad female  hx_ppm 
                         
                         in_bmi nyhc_1 nyhc_2 nyhc_3      
                         lvf_1 lvf_2 lvf_3   
                         hxmi
                         mvrg_0 mvrg_1 mvrg_2 mvrg_3 m_mvrgs
                         tvrg_2 tvrg_3 tvrg_4   ln_rvspr
                         afib chb ventarr endo 
                         eversmok pvd carotid popdz copd 
                         iddm niddm diabetes  rnldz ln_blrbn in_bun
                         vd0 vd1 vd2 lmtany 
                          m_lvfca  m_mvrgs m_htn m_sys_d   m_rvsp_
                                    / link=glogit ;   run;    
                                             
                                             
******************************************************************************;
* Keep  estimated probabilities and ccfid                                     ;
******************************************************************************;
                                             
data probs; set probs;
keep ccfid _p_; run;                                             

proc sort data=probs; by ccfid; run;

proc transpose data=probs out=propen; by ccfid; var _p_; run;

******************************************************************************;
* Rename the estimated probabilities ;
******************************************************************************;
data propen; set propen;
rename
   col1=p_cos
   col2=p_per
   col3=p_dev 
   col4=p_ce;
   keep ccfid col1 col2 col3  col4; run;

*******************************************************************************;
* Save in a permanent data set;
*******************************************************************************;   
data library.propen; set propen; run;   

*******************************************************************************;
   



                                                                                       
                                                         