%macro skip;
// JOB lm.nyha_grp.sas
sas8.2 lm.nyha_grp.sas
//cp lm.nyha_grp.l*\
//   /studies/<study>/analyses/.
//spool cont printer 'lptxt -l110 -s6 -f LetterGothic-Bold'
spool cont to email
splfile lm.nyha_grp.l*
// FILE lm.nyha_grp.sas
%mend;
*______________________________________________________________________________;
*                                                                              ;
* /analyses/lm.nyha_grp.sas
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
* Multivariable logistic regression analysis of NYHA I vs. II vs. III/IV
* Parsimoneous and propensity models for ordinal response;
* when a category variable have more than 2 categories;
* PROC LOGISTIC by default will fit a cumulative logit model;
*______________________________________________________________________________;
*                                                                              ;
  options pagesize=90 linesize=120 pageno=1;
  libname library "&STUDY/datasets";
  libname est     "&STUDY/estimates";
  filename deciles "!MACROS/deciles.new"; %inc deciles;
  data built; set library.built;
  title1 "&STUDYDESC";
  title2 "&STUDYPOP";
  title3 "Multivariable Analysis of NYHA I vs. II vs. III/IV";
*******************************************************************************;
* Data transformations                                                         ;
  filename vars "&STUDY/datasets/vars.sas";
  %inc vars; %vars(missing=1);
*______________________________________________________________________________;
*                                                                              ;
*       M U L T I V A R I A B L E   L O G I S T I C   R E G R E S S I O N      ;
*    Stepwise Selection:
*          With cululative logit link   ;
*______________________________________________________________________________;
*                                                                              ;
  %macro model;
  proc logistic data=built ;

       model nyha_grp =
                       agee ln_opyrs iv_opyrs  hx_chf hdl_pr
                       afib_pr female hx_copd hx_htn hx_endo
                       sgn_cosg  mp_fib ln_bmi

  /* Demographic */
       female male age ln_age in_age age2 agee
       ht ln_ht wt wt2 in_wt ln_wt wtht
       bsa ln_bsa bmi bmi2 in_bmi ln_bmi
  /* Symptoms */
       /*nyha_pr ln_nyha in_nyha cac_pr emgsrg
       nyha1 nyha2 nyha3 nyha4 nyha34 nyha_grp*/
  /* Ventricular function */
       lvfcath ln_lvfc in_lvfc lvfc2 in2lvfc lvf34 lvf3 lvfnorml
       ekginf hx_mi lvefcath lvfecho /*lvefecho lvaneur lvglobal*/
  /* General valve pathology */
       avrgrg avrgsev avsten mvrgrg mvrgsev mvsten tvrgsev tvsten
  /* Specific mitral valve pathology & etiology*/
       mp_ca_ru mp_cp_ru  mp_ca_el  mp_cp_el  mp_cu_el
       mp_la_pr mp_lp_pr  mp_p_el  /*mp_sv_fu*/  mp_lf_pf mp_lf_cl
      /*mp_an_fu mp_lf_th  mp_lf_ca  mp_an_ca mp_al_ca mp_lf_pr*/
       mp_an_dl mp_endo   mp_deg  mp_cal  mp_fib
       mp_regur mp_di_lv  mp_lappr  mp_laopr mp_lpopr
       mp_capru  mp_caoru  mp_cporu  mp_c_el  mp_c_ru
       mg_la_pr mg_lp_pr mg_lappr  mg_laopr mg_lpopr
  /* Cardiac comorbidity */
       hx_fcad afib_pr chb_pr varr_pr  hx_endo hx_chf
  /* Non-cardiac comorbidity */
       hx_smoke hx_pvd  hx_cardz hx_popdz hx_copd hx_htn
       hx_iddm hx_niddm hx_dmtrt hx_rnldz hx_cva
       chol_pr ln_cholp chol2 in_cholp trig_pr hdl_pr ldl_pr
       bun_pr creat_pr  blrbn_pr hct_pr
       ln_bun ln_blrbn
  /* Coronary anatomy */
       cad_sys cadsys2 vd0 vd1 vd0_1 vd2 vd3 vd2_3 vd123
       lmt lmt70 lmt50 lmtany lad lad50 lad70 ladany
       lcx lcx50 lcx70 lcxany rca rca50 rca70 rcaany
  /* Surgeon */
       sgn_effl sgn_grov sgn_loop sgn_tayl sgn_cosg sgn_lytl
       sgn_stew sgn_mcca sgn_smed sgn_sabk sgn_gilv sgn_banb
       sgn_gost sgn_othe
  /* Experience */
       iv_opyrs ln_opyrs in_opyrs in2opyrs
  /* Interaction terms:  Patient */
       age_diab lagediab iagediab i2agdiab age2diab
       /*mplppr_r*/
  /* Missing value flags */
       /*ms_size ms_nyhap ms_cacp ms_lvfc ms_lvfe ms_emerg ms_lvefc
       ms_lvefe ms_hxmi ms_fcad ms_endo ms_smoke ms_cholp ms_trigp
       ms_hdlp ms_ldlp ms_bunp ms_creap ms_blrbp ms_hctp ms_varrp
       ms_dmtx ms_lmt ms_lad ms_lcx ms_rca ms_cadsy ms_ivcpb ms_ivaoc*/

     /* pladia plvidd plvisd ppwt pivswt pef_echo pavpkgra
           pavmngra pav_area pmv_sten pmvpkgra pmvmngra ptv_rvel plasarea
           pladarea prvsp plvedvol plvesvol pav_reg ptv_reg pmv_reg */
           
/ selection=forward details sle=0.07 sls=0.05 start=13 stop=14;


     run;

*******************************************************************************;
* Check goodness of fit of model to the data (Hosmer-Lemeshow statistics)      ;
  %deciles(in=decile, _event_=nyha_grp,_p_=phat,
           vars=age iv_opyrs female nyha_pr mp_an_dl sgn_cosg mg_laopr
          mg_lappr mg_lpopr);


  proc freq data=built;
     table (mp_regur mp_endo mp_fib hx_endo hx_chf female mg_lpopr)*nyha_grp;
  run;

  %mend;
*______________________________________________________________________________;
*                                                                              ;
*                   P A R S I M O N E O U S   M O D E L                        ;
*______________________________________________________________________________;
*                                                                              ;
 
 *ods rtf body="&STUDY/xxx.rtf"; 
  title4 "Parsimoneous Model";
  proc logistic data=built covout outest=outest ;
      model nyha_grp = agee  female bmi2 bsa hx_mi 
               afib_pr hdl_pr  hx_copd hx_htn hx_endo mp_fib
               ln_opyrs iv_opyrs  sgn_cosg
                ;
       run;
 *ods rtf close;
 
 
 
 
 
 
 
 

*******************************************************************************;
* Put estimates to file                                                        ;
  data est.tlmnygrp;
    set outest;
  run;
*______________________________________________________________________________;
*                                                                              ;
*                    P R O P E N S I T Y   M O D E L                           ;
*______________________________________________________________________________;
*                                                                              ;
  title4 "Propensity Model";
  proc logistic data=built covout outest=outest ;
        output out=probs p=_p_;
        model nyha_grp = agee ln_opyrs iv_opyrs
               hx_chf /*keep? removed from parsimoneous model*/
               hdl_pr afib_pr female hx_copd hx_htn hx_endo
               sgn_cosg  mp_fib hx_mi

              /* Other variables */
                  ln_bmi
                  cad_sys
                  emgsrg  mp_regur chb_pr
                  hx_pvd  tvrgsev
                  ;

          run;
******************************************************************************;
* Keep  estimated cumulative probabilities and ccfid                         ;
******************************************************************************;
*****_p_ is estimate of the cumulative probability;

 data probs1; set probs;
 keep ccfid _p_; run;

 proc sort data=probs1; by ccfid; run;

 proc transpose data=probs1 out=propen; by ccfid; var _p_; run;

******************************************************************************; 
* calculate the individual probabilities from cumulative probabilities;
* NOTE: NOTE: If decending option is used in the PROC  LOGISTIC statement;
* cumulatice probabilities are given from the last category; 
******************************************************************************;
data propen; set propen;
p1=col1;
p2=col2-col1;
p3=1-col2; run;

******************************************************************************;
* Rename the estimated probabilities                                          ;
******************************************************************************;
 data propen1; set propen;
 rename
   p1=p_ny1
   p2=p_ny2
   p3=p_ny34
   ;
   keep ccfid p1 p2 p3; run;

*******************************************************************************;
* Form quintiles and deciles                                                   ;
  proc sort data=propen1; by p_ny34;

  data propen2; set propen1 nobs=num;
  _nobs_=num;
  quint=_nobs_/5;
  quintile=int(_n_/quint) +1;
    if quintile>5 then quintile=5;
  dec=_nobs_/10;
  decile=int(_n_/dec) +1;
    if decile>10 then decile=10;
  logtny34=log(p_ny34/(1-p_ny34));
  label quintile="Propensity quintile based on p_ny34"
        decile= "Propensity decile based on p_ny34"
     ;
  rename quintile=quinty34 decile=declny34;
  keep ccfid p_ny1 p_ny2 p_ny34 logtny34 quintile decile;
  run;
  proc sort data=propen2; by ccfid; run;


*******************************************************************************;
* Save in a permanent data set;
*******************************************************************************;
 data library.propen_nygrpt; set propen2; run;

  proc contents; run;

*******************************************************************************;
