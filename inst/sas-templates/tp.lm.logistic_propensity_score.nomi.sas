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
**  OLD TEMPLATE (WITHOUT MULTIPLE IMPUTATION);
*______________________________________________________________________________;
*                                                                              ;
* /analyses/lm.mip.sas
*______________________________________________________________________________;
*                                                                              ;
  %let STUDY=/studies/<study>;
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
  
  title1 "Aortic valve alone surgery";
  title2 "(CCF, 19xx to 20xx, n=xxxx)";
  title3 "Multivariable Analysis of <study description>";
*******************************************************************************;
*******************************************************************************;
* read in the imputed complete datasets;
data built; set library.built; run;
*******************************************************************************;
* Data Transformation                                                          ;                                                      
  filename vars "&STUDY//datasets/vars.sas";
  %inc vars; %vars; 
  

*______________________________________________________________________________;
*                                                                              ;
*       M U L T I V A R I A B L E   L O G I S T I C   R E G R E S S I O N      ;
*______________________________________________________________________________;
*                                                                              ;
%macro select;
  title4 "Selection Model";
  proc logistic data=built covout outest=outest descending;
       output out=decile p=_y_;
   model mip=
  
  /* Demographic */
  /* Symptoms */
  /* Ventricular function */
  /* Cardiac comorbidity */
  /* Non-cardiac comorbidity */
  /* Coronary anatomy */
  /* Experience */
     
      /selection=stepwise details maxstep=15; 
   *selection=forward details sle=0.07 sls=0.05 start=13 stop=14;
   run;   

%mend; 
*______________________________________________________________________________;
*                                                                              ;
*                G O O D N E S S    O F    F I T                               ; 
*     Check goodness of fit of model to the data (Hosmer-Lemeshow statistics)  ;
*______________________________________________________________________________;
* Use complete data set to check the goodness of fit;
%macro model_fit;
 proc logistic data=built descending ;
        output out=decile p=phat;
       model mip = ln_opyrs ave_dege ave_rheu sgn_jose sgn_gost sgn_mcca lvfc2 tvrgsev
                  bun_pr wt2 ; run;

  %deciles(in=decile, _event_=mip,_p_=phat,
           vars=ln_opyrs ave_dege ave_rheu sgn_jose sgn_gost lvfc2 tvrgsev
                bun_pr wt2);
%mend;
*%model_fit;
*______________________________________________________________________________;
*                                                                              ;
*               F I N A L   M U L T I V A R I A B L  E   M O D E L             ;
*                     P A R S I M O N I O U S    M O D E L                     ;
*______________________________________________________________________________;
*                                                                              ;
%macro pars;
ods rtf body= "&STUDY/analyses/lm.mip.rtf";
  title4 "Parsimoneous Model";
  proc logistic data=built covout outest=outest descending;
      output out=decile p=_y_;
      model mip = agee iv_opyrs  mp_cal
                   /*mp_c_ru  mp_c_el hx_cva app_mini */  mg_la_pr
                   sgn_cosg  sgn_gost sgn_tayl sgn_stew
                    ;
       run;

       run;
  ods rtf close; 
%mend;
%pars;

*******************************************************************************;
* Put estimates to file                                                        ;
  data est.lmmip; set outest; run;
*******************************************************************************;
* Check goodness of fit of model to the data (Hosmer-Lemeshow statistics)      ;
  %deciles(in=decile, _event_=mip, _p_=_y_,
           vars=age iv_opyrs female nyha_pr mp_an_dl sgn_cosg mg_laopr mg_lappr);
*______________________________________________________________________________;
*                                                                              ;
*               F I N A L   M U L T I V A R I A B L  E   M O D E L             ;
*                         P R O P E N S I T Y    M O D E L                     ;
*______________________________________________________________________________;
*                                                                              ;
%macro propen;
  title4 "Propensity Model";
  proc logistic data=built covout outest=outest descending;
       output out=decile p=_y_;
      model mip = agee iv_opyrs hx_cva  mp_cal  mg_la_pr
                    sgn_cosg sgn_gost  sgn_tayl sgn_stew chb_pr

       /* Other variables */
  /* Demographic */
       female /*ln_ht ln_wt*/ ln_bsa
  /* Ventricular function */
       lvfc2   lad50 
  /* Symptoms */   nyha_pr
  /* Pathology */  tvrgsev
  /* Cardiac comorbidity */
  /* Non-cardiac comorbidity */
  /* Coronary anatomy */
  /* Experience */
  /* Interaction terms:  Patient */
  /* Missing value flags */

       ;
      run;


*******************************************************************************;
* Put estimates to file                                                        ;
  data est.mip; set outest;
      run;
*______________________________________________________________________________;
*                                                                              ;
*                   P R O P E N S I T Y   S T A T I S T I C S                  ;
*______________________________________________________________________________;
*                                                                              ;
* Form quintiles and deciles                                                   ;
  proc sort data=decile; by _y_;
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
  proc freq data=built; table mip*
     (nyha_pr mp_endo mp_an_dl mp_cal mp_cporu
      mg_la_pr mg_lappr  mp_c_el mp_ca_el mp_cp_el mp_c_ru mp_di_lv hx_cva
      mp_c_ru  sgn_cosg sgn_tayl sgn_stew sgn_gost vd0 cad_sys vd123);
  run;

  title5 "Quintile Statistics: Categorical";
  proc freq data=decile; by quintile; tables mip*(
  /* Demographic */
       female
  /* Symptoms */
       nyha_pr  cac_pr emgsrg
  /* Ventricular function */
       lvfcath  hx_mi
  /* Cardiac comorbidity */
       hx_fcad afib_pr chb_pr varr_pr hx_csurg surg_num hx_endo hx_chf
  /* Non-cardiac comorbidity */
       hx_smoke hx_pvd hx_cardz hx_popdz hx_copd hx_htn
       hx_iddm hx_niddm hx_dmtrt hx_rnldz hx_cva
  /* Coronary anatomy */
       cad_sys lmt50 lad50 lcx50 rca50
  /* Specific mitral valve pathology & etiology*/
       mg_la_pr mg_lappr  mp_c_el  mp_c_ru
  /* Procedure */
       ita_num cabg tvrpr tvrpl
  /* Support */
       onpump circarr
  /* Postoperative management */
       iabp_po
  /* Postoperative complications */
       mb_bld mb_cva mb_rnl mb_rsp mb_sep
       mb_afib mb_stok mb_tamp hdeath
       dead reop1

       )/chisq;
  run;
*******************************************************************************;
* T-testing                                                                    ;
  title5 "Quintile Statistics: Continuous";
  proc ttest data=decile; class mip; by quintile; var
  /* Demographic */
       age ht wt bsa bmi
  /* Ventricular functon */
       lvefcath lvefecho
  /* Non-cardiac comorbidity */
       chol_pr trig_pr hdl_pr ldl_pr bun_pr creat_pr blrbn_pr hct_pr
  /* Experience */
       iv_opyrs
  ;
  run;
*______________________________________________________________________________;
*                                                                              ;
*     S T O R E   P R O P E N S I T Y   S C O R E   F O R   E A C H   O B S    ;
*______________________________________________________________________________;
*                                                                              ;
* Put individual propensity scores to file                                     ;
  data library.propen; set decile;
  _p_=_propen;
  _logit_=log(_p_/(1-_p_));
  prob=_p_;
  keep ccfid _propen prob _logit_ mip quintile decile;
  label _p_   ="Propensity for MIP"
        prob  ="Propensity for MIP"
       _logit_="Logit for MIP"
       ;
  proc means n mean std min max sum;
  run;
*******************************************************************************;
* Check goodness of fit of model to the data (Hosmer-Lemeshow statistics)      ;
  title5 "Goodness of Fit";
  %deciles(in=decile, _event_=mip, _p_=_y_,
           vars=age iv_opyrs female nyha_pr sgn_cosg mp_cal age mg_la_pr);
  run;
*******************************************************************************;

  %mend;
 %propen;