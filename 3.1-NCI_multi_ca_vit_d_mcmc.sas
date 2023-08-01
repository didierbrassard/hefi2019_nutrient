 /*************************************************************************/
 /*                                                                       */
 /*                     CCHS 2015 - Nutrition (PUMF)                      */
 /*                                                                       */
 /*      Estimate usual intakes based on the NCI multivariate method      */
 /*                 Code for outcome: CALCIUM + VITAMIN D                 */
 /*                                                                       */
 /*                        Author: Didier Brassard                        */
 /*                                                                       */
 /*                               Version 1                               */
 /*                                OCT2022                                */
 /*                                                                       */
 /*************************************************************************/

 /*************************************************************************/
 /*                    SET MACRO VARIABLES AND LIBRARIES                  */
 /*************************************************************************/

/* indicate file location - warning: case sensitive */
%global path suffix;
	%let path = /home/DIBRA22/hefi2019_nutrients/ ;
	%let pathCCHS = &path./Data/Raw/CCHS_Nutrition_2015_PUMF/;
	%let suffix = _ca_vit_d ;

/* indicate library folder */
	options dlcreatedir;
	libname Fmtdata "&path./Data/Processed/";	/* input data set preNCI */
	libname NCI "&path./NCI/";			/* create NCI folder, if need be */
	
/* NCI macros */
	/* boxcox svy */ %include "&path./Macros/boxcox_survey.macro.v1.2.sas";
	/* std boxcox */ %include "&path./Macros/std_cov_boxcox24hr_conday_minamt_macro_v2.0.sas";
	/* multi dist */ %include "&path./Macros/multivar_distrib_macro_v2.1.sas";
	/* percentiles*/ %include "&path./Macros/percentiles_survey.macro.v1.1.sas";

	/* Available at: https://prevention.cancer.gov/research-groups/biometry/measurement-error-impact/software-measurement-error/several-regularly-consumed-or-0 */

/* include HEFI-2019 scoring macro */
	%include "&path./Macros/hefi2019.scoring.macro.sas" ;

/* auto assign proper libraries */
	%macro mcmclib(suffix);
	options dlcreatedir;
	libname MCMC "&path./NCI/MCMC&suffix./";		  		/* create MCMC folder, if need be */
	libname chain "&path./NCI/MCMC&suffix./Markovchains/";	/* folder for trace plots */
	libname reslib "&path./NCI/MCMC&suffix./Results/";		/* results */
	libname baselib "&path./NCI/MCMC&suffix./Model/";		/* base estimates */
	libname bootlib "&path./NCI/MCMC&suffix./Bootlib/";		/* bootstrap replicates */
	libname reg "&path./NCI/MCMC&suffix./Reg/";			  	 		/* base estimate for linear reg */
	libname regb "&path./NCI/MCMC&suffix./Reg/Bootstraps/";  		/* bootstrap data for linear reg */
	
	/* multi MCMC */ %include "&path./Macros/multivar_mcmc_macro_v2.1.sas";
	%mend;
	
	%mcmclib(&suffix.);

 /*************************************************************************/
 /*                                                                       */
 /*         Prepare an input data set for the multivariate method         */
 /*                                                                       */
 /*************************************************************************/


/* 0) Import bootstrap weights */
	data work.BSW(drop=FWGT);
		%let datafid = "&pathCCHS./Bootstrap/Data_Donnee/b5.txt";
		%include "&pathCCHS./Bootstrap/SAS_SPSS/b5_i.sas";
	run;

	proc sort;
		by adm_rno;
	run;

/* 1) Combine dietary intakes data with sociodemo data */
	data intake_and_sociodeom(drop=nonzero_energy);
		merge fmtdata.intake_per24hr(in=a)
			  fmtdata.hs_nci(in=b keep=adm_rno suppid wts_p r24_weekend drig age sex );
		by adm_rno suppid;
	* remove 24-h recall with 0 energy intake ;
		if energy=0 then delete ;
	* keep adm_rno present in both data + single 24-h recall only ;
		if (a and b) then output;
	run;

	proc sort;
		by adm_rno;
	run;
	
	/* note: sample size of respondents 2y+ for both 24-h recall = 27,529 */
	
/* 2) Re-code covariates as dummy variables, select final sample and output data*/
	data preNCI;
		set intake_and_sociodeom (in=intake rename=(r24_weekend=_r24_weekend)) ;
	* Dummy variable for sequence of 24-h recalls ;
	if suppid=2 then seq2=1;
		else seq2=0;
	
	* Dummy variable for weekend 24-h recalls ;
	if _r24_weekend = 1 then r24_weekend=1;
		else r24_weekend =0;
	drop _r24_weekend;

	* make age categories according to drig groupings;
	length agec $ 6;
	
	else if drig in (12 13) then agec="65to70";
	else if drig in (14 15) then agec="71plus";
	
	* Make a dummy variable for 71 or older ;
	if agec="71plus" then agec_71plus=1 ;
	else agec_71plus=0;

	* Change variable name for sampling weights to be consistent with bootstrap weights ;
	rename wts_p = bsw0 ;
		
	* Filter to keep only respondents 65 years or older;
	if age <65 then delete ;
	run;
	
/* 3) Confirm age categories recoding */
	
	proc freq data=preNCI;
	table drig drig * agec: sex:;
	run;

/* 4) Look at intakes on a given day */

	/* 4.0) define dietary constituents */
	%let daily_list  = vf otherfoods pfab water mufa pufa sfa freesugars sodium energy calcium vit_d ;
	%let episo_list  = wg pfpb otherbevs milk rg;

	/* 4.1) Define format for zeros */
		proc format;
		value zerofmt
			0-<0.001 = "Zero"
			0.001-high = "Non-zero"
			;
		run;

	/* 4.2) Descriptive statistics */
	ods select none;
		proc freq data=preNCI ;
		format &daily_list. &episo_list. zerofmt. ; 
		table &daily_list. &episo_list. ;
		where suppid=1 ; * first 24-h recall only, non-zero energy intake;
		ods output OneWayFreqs=pZero(where=(CumPercent < 99.9));
		run;
	ods select all;
	
		proc sort data=pZero;
			by descending percent;
		run;
		
		proc print data=pZero;
		title1 "Proportion of zero intake among HEFI-2019 dietary constituents for first 24-h recall, CCHS 2015 - Nutrition";
		id table;
		var Frequency Percent ;
		run;
		title1;
	
		proc means data=preNCI n mean min p25 p50 p75 max stackods maxdec=1;
		class suppid; 
		var &daily_list. &episo_list. ;
		weight bsw0 ;
		run;

	
 /*************************************************************************/
 /*                                                                       */
 /*            NCI Multivariate method based on MCMC algorithm            */
 /*                                                                       */
 /*************************************************************************/

/*************************************************************************/
/* Stratum selection based on STRATA=sex (1/2)                           */
/*************************************************************************/
 
	 * Output data of participants (adm_rno) in stratum 1 only ;
	data _tempstratum1 ;
	set preNCI ;
	if (sex = 1) then output;
	run;
	 
	proc sort ;
	by adm_rno ;
	run;
 
/*************************************************************************/
/* Box-Cox transformations (loopboxcox2, std_cov_boxcox24hr...)          */
/*************************************************************************/

	* Output data for the first 24-h dietary recall completed;
	proc sort data = _tempstratum1 nodupkey out=inboxcox;
	by adm_rno ;
	run;
	 
	* Make sure the macro variable <best_lambda> is available outside <boxcox_survey>;
	 %global best_lambda ;
 
 /********************************************************************/
 /*   Dietary constituents #1 - wg: min. amount and lambda   */
 /********************************************************************/
 
* 1) Get minimum non-zero consumption amount ;
proc means data=_tempstratum1 min noprint;
var wg ;
where wg > 0;
output out=_min1(keep=minamount) min=minamount ;
run;
 
data _min1;
set _min1;
minamount=minamount/2;
tran_paramindex=1;
length varname $ 32;
varname="wg";
run;
 
data _temp1;
set inboxcox(keep=adm_rno bsw0 wg seq2 r24_weekend agec_71plus );
if (wg > 0) then output;
run;
 
* 2) Call the <boxcox_survey> macro to find best normal transformation ;
 %boxcox_survey(
  data    = _temp1, 
  subject = adm_rno, 
  var     = wg,  
  covars  = seq2 r24_weekend  agec_71plus ,  
  weight  = bsw0,  
  print   = N, 
  plot    = N, 
  ntitle  = 4 ); 
 

 
*3) Save <best_lambda> macro variable from <boxcox_survey> for current dietary constituent ;
data _x1;
length varname $ 32;
tran_paramindex=1;
varname="wg";
 tran_lambda=&best_lambda ;
 * The value of best_lambda was ... ;
 *
tran_lambda=0.21;
run;
 
 
 /********************************************************************/
 /*   Dietary constituents #2 - pfpb: min. amount and lambda         */
 /********************************************************************/
 
* 1) Get minimum non-zero consumption amount ;
proc means data=_tempstratum1 min noprint;
var pfpb ;
where pfpb > 0;
output out=_min2(keep=minamount) min=minamount ;
run;
 
data _min2;
set _min2;
minamount=minamount/2;
tran_paramindex=2;
length varname $ 32;
varname="pfpb";
run;
 
data _temp2;
set inboxcox(keep=adm_rno bsw0 pfpb seq2 r24_weekend agec_71plus );
if (pfpb > 0) then output;
run;
 
* 2) Call the <boxcox_survey> macro to find best normal transformation ;
 %boxcox_survey(
  data    = _temp2, 
  subject = adm_rno, 
  var     = pfpb,  
  covars  = seq2 r24_weekend  agec_71plus ,  
  weight  = bsw0,  
  print   = N, 
  plot    = N, 
  ntitle  = 4 ); 
 

 
*3) Save <best_lambda> macro variable from <boxcox_survey> for current dietary constituent ;
data _x2;
length varname $ 32;
tran_paramindex=2;
varname="pfpb";
 tran_lambda=&best_lambda ;
 * The value of best_lambda was ... ;
 *
tran_lambda=0.18;
run;
 
 
 /********************************************************************/
 /*   Dietary constituents #3 - otherbevs: min. amount and lambda   */
 /********************************************************************/
 
* 1) Get minimum non-zero consumption amount ;
proc means data=_tempstratum1 min noprint;
var otherbevs ;
where otherbevs > 0;
output out=_min3(keep=minamount) min=minamount ;
run;
 
data _min3;
set _min3;
minamount=minamount/2;
tran_paramindex=3;
length varname $ 32;
varname="otherbevs";
run;
 
data _temp3;
set inboxcox(keep=adm_rno bsw0 otherbevs seq2 r24_weekend agec_71plus );
if (otherbevs > 0) then output;
run;
 
* 2) Call the <boxcox_survey> macro to find best normal transformation ;
 %boxcox_survey(
  data    = _temp3, 
  subject = adm_rno, 
  var     = otherbevs,  
  covars  = seq2 r24_weekend  agec_71plus ,  
  weight  = bsw0,  
  print   = N, 
  plot    = N, 
  ntitle  = 4 ); 
 

 
*3) Save <best_lambda> macro variable from <boxcox_survey> for current dietary constituent ;
data _x3;
length varname $ 32;
tran_paramindex=3;
varname="otherbevs";
 tran_lambda=&best_lambda ;
 * The value of best_lambda was ... ;
 *
tran_lambda=0.28;
run;
 
 
 /********************************************************************/
 /*   Dietary constituents #4 - milk: min. amount and lambda   */
 /********************************************************************/
 
* 1) Get minimum non-zero consumption amount ;
proc means data=_tempstratum1 min noprint;
var milk ;
where milk > 0;
output out=_min4(keep=minamount) min=minamount ;
run;
 
data _min4;
set _min4;
minamount=minamount/2;
tran_paramindex=4;
length varname $ 32;
varname="milk";
run;
 
data _temp4;
set inboxcox(keep=adm_rno bsw0 milk seq2 r24_weekend agec_71plus );
if (milk > 0) then output;
run;
 
* 2) Call the <boxcox_survey> macro to find best normal transformation ;
 %boxcox_survey(
  data    = _temp4, 
  subject = adm_rno, 
  var     = milk,  
  covars  = seq2 r24_weekend  agec_71plus ,  
  weight  = bsw0,  
  print   = N, 
  plot    = N, 
  ntitle  = 4 ); 
 

 
*3) Save <best_lambda> macro variable from <boxcox_survey> for current dietary constituent ;
data _x4;
length varname $ 32;
tran_paramindex=4;
varname="milk";
 tran_lambda=&best_lambda ;
 * The value of best_lambda was ... ;
 *
tran_lambda=0.3;
run;
 
 
 /********************************************************************/
 /*   Dietary constituents #5 - rg: min. amount and lambda   */
 /********************************************************************/
 
* 1) Get minimum non-zero consumption amount ;
proc means data=_tempstratum1 min noprint;
var rg ;
where rg > 0;
output out=_min5(keep=minamount) min=minamount ;
run;
 
data _min5;
set _min5;
minamount=minamount/2;
tran_paramindex=5;
length varname $ 32;
varname="rg";
run;
 
data _temp5;
set inboxcox(keep=adm_rno bsw0 rg seq2 r24_weekend agec_71plus );
if (rg > 0) then output;
run;
 
* 2) Call the <boxcox_survey> macro to find best normal transformation ;
 %boxcox_survey(
  data    = _temp5, 
  subject = adm_rno, 
  var     = rg,  
  covars  = seq2 r24_weekend  agec_71plus ,  
  weight  = bsw0,  
  print   = N, 
  plot    = N, 
  ntitle  = 4 ); 
 

 
*3) Save <best_lambda> macro variable from <boxcox_survey> for current dietary constituent ;
data _x5;
length varname $ 32;
tran_paramindex=5;
varname="rg";
 tran_lambda=&best_lambda ;
 * The value of best_lambda was ... ;
 *
tran_lambda=0.37;
run;
 
 
 /********************************************************************/
 /*   Dietary constituents #6 - vf: min. amount and lambda   */
 /********************************************************************/
 
* 1) Get minimum non-zero consumption amount ;
proc means data=_tempstratum1 min noprint;
var vf ;
where vf > 0;
output out=_min6(keep=minamount) min=minamount ;
run;
 
data _min6;
set _min6;
minamount=minamount/2;
tran_paramindex=6;
length varname $ 32;
varname="vf";
run;
 
data _temp6;
set inboxcox(keep=adm_rno bsw0 vf seq2 r24_weekend agec_71plus );
if (vf > 0) then output;
run;
 
* 2) Call the <boxcox_survey> macro to find best normal transformation ;
 %boxcox_survey(
  data    = _temp6, 
  subject = adm_rno, 
  var     = vf,  
  covars  = seq2 r24_weekend  agec_71plus ,  
  weight  = bsw0,  
  print   = N, 
  plot    = N, 
  ntitle  = 4 ); 
 

 
*3) Save <best_lambda> macro variable from <boxcox_survey> for current dietary constituent ;
data _x6;
length varname $ 32;
tran_paramindex=6;
varname="vf";
 tran_lambda=&best_lambda ;
 * The value of best_lambda was ... ;
 *
tran_lambda=0.19;
run;
 
 
 /********************************************************************/
 /*   Dietary constituents #7 - otherfoods: min. amount and lambda   */
 /********************************************************************/
 
* 1) Get minimum non-zero consumption amount ;
proc means data=_tempstratum1 min noprint;
var otherfoods ;
where otherfoods > 0;
output out=_min7(keep=minamount) min=minamount ;
run;
 
data _min7;
set _min7;
minamount=minamount/2;
tran_paramindex=7;
length varname $ 32;
varname="otherfoods";
run;
 
data _temp7;
set inboxcox(keep=adm_rno bsw0 otherfoods seq2 r24_weekend agec_71plus );
if (otherfoods > 0) then output;
run;
 
* 2) Call the <boxcox_survey> macro to find best normal transformation ;
 %boxcox_survey(
  data    = _temp7, 
  subject = adm_rno, 
  var     = otherfoods,  
  covars  = seq2 r24_weekend  agec_71plus ,  
  weight  = bsw0,  
  print   = N, 
  plot    = N, 
  ntitle  = 4 ); 
 

 
*3) Save <best_lambda> macro variable from <boxcox_survey> for current dietary constituent ;
data _x7;
length varname $ 32;
tran_paramindex=7;
varname="otherfoods";
 tran_lambda=&best_lambda ;
 * The value of best_lambda was ... ;
 *
tran_lambda=0.22;
run;
 
 
 /********************************************************************/
 /*   Dietary constituents #8 - pfab: min. amount and lambda   */
 /********************************************************************/
 
* 1) Get minimum non-zero consumption amount ;
proc means data=_tempstratum1 min noprint;
var pfab ;
where pfab > 0;
output out=_min8(keep=minamount) min=minamount ;
run;
 
data _min8;
set _min8;
minamount=minamount/2;
tran_paramindex=8;
length varname $ 32;
varname="pfab";
run;
 
data _temp8;
set inboxcox(keep=adm_rno bsw0 pfab seq2 r24_weekend agec_71plus );
if (pfab > 0) then output;
run;
 
* 2) Call the <boxcox_survey> macro to find best normal transformation ;
 %boxcox_survey(
  data    = _temp8, 
  subject = adm_rno, 
  var     = pfab,  
  covars  = seq2 r24_weekend  agec_71plus ,  
  weight  = bsw0,  
  print   = N, 
  plot    = N, 
  ntitle  = 4 ); 
 

 
*3) Save <best_lambda> macro variable from <boxcox_survey> for current dietary constituent ;
data _x8;
length varname $ 32;
tran_paramindex=8;
varname="pfab";
 tran_lambda=&best_lambda ;
 * The value of best_lambda was ... ;
 *
tran_lambda=0.25;
run;
 
 
 /********************************************************************/
 /*   Dietary constituents #9 - water: min. amount and lambda   */
 /********************************************************************/
 
* 1) Get minimum non-zero consumption amount ;
proc means data=_tempstratum1 min noprint;
var water ;
where water > 0;
output out=_min9(keep=minamount) min=minamount ;
run;
 
data _min9;
set _min9;
minamount=minamount/2;
tran_paramindex=9;
length varname $ 32;
varname="water";
run;
 
data _temp9;
set inboxcox(keep=adm_rno bsw0 water seq2 r24_weekend agec_71plus );
if (water > 0) then output;
run;
 
* 2) Call the <boxcox_survey> macro to find best normal transformation ;
 %boxcox_survey(
  data    = _temp9, 
  subject = adm_rno, 
  var     = water,  
  covars  = seq2 r24_weekend  agec_71plus ,  
  weight  = bsw0,  
  print   = N, 
  plot    = N, 
  ntitle  = 4 ); 
 

 
*3) Save <best_lambda> macro variable from <boxcox_survey> for current dietary constituent ;
data _x9;
length varname $ 32;
tran_paramindex=9;
varname="water";
 tran_lambda=&best_lambda ;
 * The value of best_lambda was ... ;
 *
tran_lambda=0.48;
run;
 
 
 /********************************************************************/
 /*   Dietary constituents #10 - mufa: min. amount and lambda   */
 /********************************************************************/
 
* 1) Get minimum non-zero consumption amount ;
proc means data=_tempstratum1 min noprint;
var mufa ;
where mufa > 0;
output out=_min10(keep=minamount) min=minamount ;
run;
 
data _min10;
set _min10;
minamount=minamount/2;
tran_paramindex=10;
length varname $ 32;
varname="mufa";
run;
 
data _temp10;
set inboxcox(keep=adm_rno bsw0 mufa seq2 r24_weekend agec_71plus );
if (mufa > 0) then output;
run;
 
* 2) Call the <boxcox_survey> macro to find best normal transformation ;
 %boxcox_survey(
  data    = _temp10, 
  subject = adm_rno, 
  var     = mufa,  
  covars  = seq2 r24_weekend  agec_71plus ,  
  weight  = bsw0,  
  print   = N, 
  plot    = N, 
  ntitle  = 4 ); 
 

 
*3) Save <best_lambda> macro variable from <boxcox_survey> for current dietary constituent ;
data _x10;
length varname $ 32;
tran_paramindex=10;
varname="mufa";
 tran_lambda=&best_lambda ;
 * The value of best_lambda was ... ;
 *
tran_lambda=0.16;
run;
 
 
 /********************************************************************/
 /*   Dietary constituents #11 - pufa: min. amount and lambda   */
 /********************************************************************/
 
* 1) Get minimum non-zero consumption amount ;
proc means data=_tempstratum1 min noprint;
var pufa ;
where pufa > 0;
output out=_min11(keep=minamount) min=minamount ;
run;
 
data _min11;
set _min11;
minamount=minamount/2;
tran_paramindex=11;
length varname $ 32;
varname="pufa";
run;
 
data _temp11;
set inboxcox(keep=adm_rno bsw0 pufa seq2 r24_weekend agec_71plus );
if (pufa > 0) then output;
run;
 
* 2) Call the <boxcox_survey> macro to find best normal transformation ;
 %boxcox_survey(
  data    = _temp11, 
  subject = adm_rno, 
  var     = pufa,  
  covars  = seq2 r24_weekend  agec_71plus ,  
  weight  = bsw0,  
  print   = N, 
  plot    = N, 
  ntitle  = 4 ); 
 

 
*3) Save <best_lambda> macro variable from <boxcox_survey> for current dietary constituent ;
data _x11;
length varname $ 32;
tran_paramindex=11;
varname="pufa";
 tran_lambda=&best_lambda ;
 * The value of best_lambda was ... ;
 *
tran_lambda=0.17;
run;
 
 
 /********************************************************************/
 /*   Dietary constituents #12 - sfa: min. amount and lambda         */
 /********************************************************************/
 
* 1) Get minimum non-zero consumption amount ;
proc means data=_tempstratum1 min noprint;
var sfa ;
where sfa > 0;
output out=_min12(keep=minamount) min=minamount ;
run;
 
data _min12;
set _min12;
minamount=minamount/2;
tran_paramindex=12;
length varname $ 32;
varname="sfa";
run;
 
data _temp12;
set inboxcox(keep=adm_rno bsw0 sfa seq2 r24_weekend agec_71plus );
if (sfa > 0) then output;
run;
 
* 2) Call the <boxcox_survey> macro to find best normal transformation ;
 %boxcox_survey(
  data    = _temp12, 
  subject = adm_rno, 
  var     = sfa,  
  covars  = seq2 r24_weekend  agec_71plus ,  
  weight  = bsw0,  
  print   = N, 
  plot    = N, 
  ntitle  = 4 ); 
 

 
*3) Save <best_lambda> macro variable from <boxcox_survey> for current dietary constituent ;
data _x12;
length varname $ 32;
tran_paramindex=12;
varname="sfa";
 tran_lambda=&best_lambda ;
 * The value of best_lambda was ... ;
 *
tran_lambda=0.27;
run;
 
 
 /********************************************************************/
 /*   Dietary constituents #13 - freesugars: min. amount and lambda   */
 /********************************************************************/
 
* 1) Get minimum non-zero consumption amount ;
proc means data=_tempstratum1 min noprint;
var freesugars ;
where freesugars > 0;
output out=_min13(keep=minamount) min=minamount ;
run;
 
data _min13;
set _min13;
minamount=minamount/2;
tran_paramindex=13;
length varname $ 32;
varname="freesugars";
run;
 
data _temp13;
set inboxcox(keep=adm_rno bsw0 freesugars seq2 r24_weekend agec_71plus );
if (freesugars > 0) then output;
run;
 
* 2) Call the <boxcox_survey> macro to find best normal transformation ;
 %boxcox_survey(
  data    = _temp13, 
  subject = adm_rno, 
  var     = freesugars,  
  covars  = seq2 r24_weekend  agec_71plus ,  
  weight  = bsw0,  
  print   = N, 
  plot    = N, 
  ntitle  = 4 ); 
 

 
*3) Save <best_lambda> macro variable from <boxcox_survey> for current dietary constituent ;
data _x13;
length varname $ 32;
tran_paramindex=13;
varname="freesugars";
 tran_lambda=&best_lambda ;
 * The value of best_lambda was ... ;
 *
tran_lambda=0.39;
run;
 
 
 /********************************************************************/
 /*   Dietary constituents #14 - sodium: min. amount and lambda   */
 /********************************************************************/
 
* 1) Get minimum non-zero consumption amount ;
proc means data=_tempstratum1 min noprint;
var sodium ;
where sodium > 0;
output out=_min14(keep=minamount) min=minamount ;
run;
 
data _min14;
set _min14;
minamount=minamount/2;
tran_paramindex=14;
length varname $ 32;
varname="sodium";
run;
 
data _temp14;
set inboxcox(keep=adm_rno bsw0 sodium seq2 r24_weekend agec_71plus );
if (sodium > 0) then output;
run;
 
* 2) Call the <boxcox_survey> macro to find best normal transformation ;
 %boxcox_survey(
  data    = _temp14, 
  subject = adm_rno, 
  var     = sodium,  
  covars  = seq2 r24_weekend  agec_71plus ,  
  weight  = bsw0,  
  print   = N, 
  plot    = N, 
  ntitle  = 4 ); 
 

 
*3) Save <best_lambda> macro variable from <boxcox_survey> for current dietary constituent ;
data _x14;
length varname $ 32;
tran_paramindex=14;
varname="sodium";
 tran_lambda=&best_lambda ;
 * The value of best_lambda was ... ;
 *
tran_lambda=0.41;
run;
 
 
 /********************************************************************/
 /*   Dietary constituents #15 - energy: min. amount and lambda   */
 /********************************************************************/
 
* 1) Get minimum non-zero consumption amount ;
proc means data=_tempstratum1 min noprint;
var energy ;
where energy > 0;
output out=_min15(keep=minamount) min=minamount ;
run;
 
data _min15;
set _min15;
minamount=minamount/2;
tran_paramindex=15;
length varname $ 32;
varname="energy";
run;
 
data _temp15;
set inboxcox(keep=adm_rno bsw0 energy seq2 r24_weekend agec_71plus );
if (energy > 0) then output;
run;
 
* 2) Call the <boxcox_survey> macro to find best normal transformation ;
 %boxcox_survey(
  data    = _temp15, 
  subject = adm_rno, 
  var     = energy,  
  covars  = seq2 r24_weekend  agec_71plus ,  
  weight  = bsw0,  
  print   = N, 
  plot    = N, 
  ntitle  = 4 ); 
 

 
*3) Save <best_lambda> macro variable from <boxcox_survey> for current dietary constituent ;
data _x15;
length varname $ 32;
tran_paramindex=15;
varname="energy";
 tran_lambda=&best_lambda ;
 * The value of best_lambda was ... ;
 *
tran_lambda=0.35;
run;
 
 
 /********************************************************************/
 /*   Dietary constituents #16 - calcium: min. amount and lambda   */
 /********************************************************************/
 
* 1) Get minimum non-zero consumption amount ;
proc means data=_tempstratum1 min noprint;
var calcium ;
where calcium > 0;
output out=_min16(keep=minamount) min=minamount ;
run;
 
data _min16;
set _min16;
minamount=minamount/2;
tran_paramindex=16;
length varname $ 32;
varname="calcium";
run;
 
data _temp16;
set inboxcox(keep=adm_rno bsw0 calcium seq2 r24_weekend agec_71plus );
if (calcium > 0) then output;
run;
 
* 2) Call the <boxcox_survey> macro to find best normal transformation ;
 %boxcox_survey(
  data    = _temp16, 
  subject = adm_rno, 
  var     = calcium,  
  covars  = seq2 r24_weekend  agec_71plus ,  
  weight  = bsw0,  
  print   = N, 
  plot    = N, 
  ntitle  = 4 ); 
 

 
*3) Save <best_lambda> macro variable from <boxcox_survey> for current dietary constituent ;
data _x16;
length varname $ 32;
tran_paramindex=16;
varname="calcium";
 tran_lambda=&best_lambda ;
 * The value of best_lambda was ... ;
 *
tran_lambda=0.19;
run;
 
 
 /********************************************************************/
 /*   Dietary constituents #17 - vit_d: min. amount and lambda       */
 /********************************************************************/
 
* 1) Get minimum non-zero consumption amount ;
proc means data=_tempstratum1 min noprint;
var vit_d ;
where vit_d > 0;
output out=_min17(keep=minamount) min=minamount ;
run;
 
data _min17;
set _min17;
minamount=minamount/2;
tran_paramindex=17;
length varname $ 32;
varname="vit_d";
run;
 
data _temp17;
set inboxcox(keep=adm_rno bsw0 vit_d seq2 r24_weekend agec_71plus );
if (vit_d > 0) then output;
run;
 
* 2) Call the <boxcox_survey> macro to find best normal transformation ;
 %boxcox_survey(
  data    = _temp17, 
  subject = adm_rno, 
  var     = vit_d,  
  covars  = seq2 r24_weekend  agec_71plus ,  
  weight  = bsw0,  
  print   = N, 
  plot    = N, 
  ntitle  = 4 ); 
 

 
*3) Save <best_lambda> macro variable from <boxcox_survey> for current dietary constituent ;
data _x17;
length varname $ 32;
tran_paramindex=17;
varname="vit_d";
 tran_lambda=&best_lambda ;
 * The value of best_lambda was ... ;
 *
tran_lambda=0.21;
run;
 
 /*****************************************************************/
 /*         Append data sets of all dietary constituents          */
 /*****************************************************************/
	 
	data work.xlambdas_s1_0;
	set _x1-_x17;
	run;
	 
	data work.xminamount_s1_0;
	set _min1-_min17;
	label minamount= "Min. non zero amount divided by 2";
	run;
 
 /********************************************************************/
 /*   Standardize data according to lambda and min. amount values    */
 /********************************************************************/
 
  * Call the <std_cov_boxcox24hr_conday_minamt> macro to standardize data and continuous covariates, if any; 
 %std_cov_boxcox24hr_conday_minamt(
  data                       = _tempstratum1, 
  prestand_continuous_covars =  , 
  rec24hr_epis_vars          = wg pfpb otherbevs milk rg, 
  rec24hr_daily_vars         = vf otherfoods pfab water mufa pufa sfa freesugars sodium energy calcium vit_d, 
  boxcox_tran_lambda_data    = xlambdas_s1_0, 
  minamount_data             = xminamount_s1_0, 
  print                      = y, 
  titles                     = 3 ); 
 
 
/*************************************************************************/
/* Fit the multivariate measurement error model (multivar_mcmc)          */
/*************************************************************************/
 
 
  * Call the <multivar_mcmc> macro to fit the measurement error model; 
   title1 "Fit Multivariate Measurement Error Model Using MCMC with 24-Hour Recall as Main Instrument"; 
 %multivar_mcmc(
  data                        = stdcov_stdbc24hr_conday_out, 
  subject                     = adm_rno, 
  weight_var                  = bsw0 , 
  repeat                      = suppid, 
  conday_epis_vars            = conday_wg  conday_pfpb  conday_otherbevs  conday_milk  conday_rg, 
  gst_rec24hr_epis_vars       = stdbc_wg  stdbc_pfpb  stdbc_otherbevs  stdbc_milk  stdbc_rg, 
  gst_rec24hr_daily_vars      = stdbc_vf  stdbc_otherfoods  stdbc_pfab  stdbc_water  stdbc_mufa  stdbc_pufa  stdbc_sfa  
 stdbc_freesugars  stdbc_sodium  stdbc_energy  stdbc_calcium  stdbc_vit_d, 
  covars_epis_prob            = constant1 seq2 r24_weekend agec_71plus   , 
  covars_epis_amt             = constant1 seq2 r24_weekend agec_71plus   , 
  covars_daily_amt            = constant1 seq2 r24_weekend agec_71plus   , 
  set_seed_mcmc               = 42941, 
  set_number_mcmc_iterations  = 8000,
  set_number_burn_iterations  = 3000,
  set_thin                    = 10,
  prior_sigmau_mean_data      = , 
  sigmau_constant             = , 
  gen_inverse                 = y, 
  print                       = y, 
  titles                      = 1, 
  std_print_store             = y, 
  notes_print                 = y, 
  out_lib                     = baselib, 
  out_store_label             = mcmc_s1_rep0, 
  out_save_label_max5char     = s10, 
  set_number_saved_out_data   = , 
  save_mcmc_u_out_data        = y, 
  set_number_post_mcmc_u_out  = , 
  traceplots_method1_gpath    = , 
  traceplots_method2_file_pdf = trace_rep0_s1.pdf, 
  optional_iml_store_data     = backtran_out, 
  optional_iml_store_names    = constant1 seq2 r24_weekend agec_71plus   tran_paramindex tran_lambda tran_center tran_scale 
 minamount 
  ); 
 
 
* Save lambdas (Box-Cox transformation lambda values) and minimum amount data ;
data baselib.backtran_out0_s1;
retain replicate sex;
set work.backtran_out;
* indicate replicate number, current stratum value, variable name;
replicate = 0;
sex = 1;
length varname $ 32 ;
if _N_=1 then varname = "wg" ;
if _N_=2 then varname = "pfpb" ;
if _N_=3 then varname = "otherbevs" ;
if _N_=4 then varname = "milk" ;
if _N_=5 then varname = "rg" ;
if _N_=6 then varname = "vf" ;
if _N_=7 then varname = "otherfoods" ;
if _N_=8 then varname = "pfab" ;
if _N_=9 then varname = "water" ;
if _N_=10 then varname = "mufa" ;
if _N_=11 then varname = "pufa" ;
if _N_=12 then varname = "sfa" ;
if _N_=13 then varname = "freesugars" ;
if _N_=14 then varname = "sodium" ;
if _N_=15 then varname = "energy" ;
if _N_=16 then varname = "calcium" ;
if _N_=17 then varname = "vit_d" ;
run;
 
 
/*************************************************************************/
/* Simulation of pseudo-individual (multivar_distrib)                    */
/*************************************************************************/
 
* Prepare an input data for the <optional_input_data> option in <multivar_distrib>;
	proc sort data=_tempstratum1 nodupkey out=optional_input_data(keep= adm_rno sex bsw0 
	 agec_71plus sex drig );
	by adm_rno ;
	run;

 
* Call the <multivar_distrib> macro to simulate usual intakes for pseudo-individuals; 
 %multivar_distrib(
  multivar_mcmc_out_lib           = baselib ,  
  multivar_mcmc_out_store_label   = mcmc_s1_rep0, 
  t_weightavg_covariates_list1    = constant1 constant0 constant0 agec_71plus   ,  
  t_weightavg_covariates_list2    = constant1 constant0 constant1 agec_71plus  , 
  set_value_for_weight_cov_list1  = 4, 
  set_value_for_weight_cov_list2  = 3, 
  optional_input_data             = optional_input_data , 
  optional_input_data_var_list    = , 
  optional_input_mcmc_u_out_data  = , 
  additional_output_var_list      = sex sex drig  , 
  additional_output_subject_var   = adm_rno , 
  output_mcmc_weight_var          = y  , 
  set_seed_distrib                = 89009890, 
  set_number_monte_carlo_rand_obs = 500,  
  print                           = y 
  ); 
 
 
* Save the Monte Carlo simulation data for current stratum;
data baselib.mc_t_distrib_out0_s1;
set mc_t_distrib_out;
run;
 
* delete temporary data sets ;
proc datasets lib=work nolist nodetails;
delete mc_t_distrib_out optional_input_data stdcov_stdbc24hr_conday_out xlambdas_: 
 xminamount_: backtran_out ;
run;
 
proc datasets lib=baselib nolist nodetails ;
delete multivar_mcmc_samples_u_outs10 ;
run;
 
 
/*************************************************************************/
/* Stratum selection based on STRATA=sex (2/2)                           */
/*************************************************************************/
 
 
 * Output data of participants (adm_rno) in stratum 2 only ;
	
	data _tempstratum2 ;
	set preNCI ;
	if (sex = 2) then output;
	run;
	 
	proc sort ;
	by adm_rno ;
	run;

/*************************************************************************/
/* Box-Cox transformations (loopboxcox2, std_cov_boxcox24hr...)          */
/*************************************************************************/
 
	* Output data for the first 24-h dietary recall completed;
	proc sort data = _tempstratum2 nodupkey out=inboxcox;
	by adm_rno ;
	run;
	 
	* Make sure the macro variable <best_lambda> is available outside <boxcox_survey>;
	 %global best_lambda ;
	 
 /********************************************************************/
 /*   Dietary constituents #1 - wg: min. amount and lambda   */
 /********************************************************************/
 
* 1) Get minimum non-zero consumption amount ;
proc means data=_tempstratum2 min noprint;
var wg ;
where wg > 0;
output out=_min1(keep=minamount) min=minamount ;
run;
 
data _min1;
set _min1;
minamount=minamount/2;
tran_paramindex=1;
length varname $ 32;
varname="wg";
run;
 
data _temp1;
set inboxcox(keep=adm_rno bsw0 wg seq2 r24_weekend agec_71plus );
if (wg > 0) then output;
run;
 
* 2) Call the <boxcox_survey> macro to find best normal transformation ;
 %boxcox_survey(
  data    = _temp1, 
  subject = adm_rno, 
  var     = wg,  
  covars  = seq2 r24_weekend  agec_71plus ,  
  weight  = bsw0,  
  print   = N, 
  plot    = N, 
  ntitle  = 4 ); 
 

 
*3) Save <best_lambda> macro variable from <boxcox_survey> for current dietary constituent ;
data _x1;
length varname $ 32;
tran_paramindex=1;
varname="wg";
 tran_lambda=&best_lambda ;
 * The value of best_lambda was ... ;
 *
tran_lambda=0.26;
run;
 
 
 /********************************************************************/
 /*   Dietary constituents #2 - pfpb: min. amount and lambda   */
 /********************************************************************/
 
* 1) Get minimum non-zero consumption amount ;
proc means data=_tempstratum2 min noprint;
var pfpb ;
where pfpb > 0;
output out=_min2(keep=minamount) min=minamount ;
run;
 
data _min2;
set _min2;
minamount=minamount/2;
tran_paramindex=2;
length varname $ 32;
varname="pfpb";
run;
 
data _temp2;
set inboxcox(keep=adm_rno bsw0 pfpb seq2 r24_weekend agec_71plus );
if (pfpb > 0) then output;
run;
 
* 2) Call the <boxcox_survey> macro to find best normal transformation ;
 %boxcox_survey(
  data    = _temp2, 
  subject = adm_rno, 
  var     = pfpb,  
  covars  = seq2 r24_weekend  agec_71plus ,  
  weight  = bsw0,  
  print   = N, 
  plot    = N, 
  ntitle  = 4 ); 
 

 
*3) Save <best_lambda> macro variable from <boxcox_survey> for current dietary constituent ;
data _x2;
length varname $ 32;
tran_paramindex=2;
varname="pfpb";
 tran_lambda=&best_lambda ;
 * The value of best_lambda was ... ;
 *
tran_lambda=0.23;
run;
 
 
 /********************************************************************/
 /*   Dietary constituents #3 - otherbevs: min. amount and lambda   */
 /********************************************************************/
 
* 1) Get minimum non-zero consumption amount ;
proc means data=_tempstratum2 min noprint;
var otherbevs ;
where otherbevs > 0;
output out=_min3(keep=minamount) min=minamount ;
run;
 
data _min3;
set _min3;
minamount=minamount/2;
tran_paramindex=3;
length varname $ 32;
varname="otherbevs";
run;
 
data _temp3;
set inboxcox(keep=adm_rno bsw0 otherbevs seq2 r24_weekend agec_71plus );
if (otherbevs > 0) then output;
run;
 
* 2) Call the <boxcox_survey> macro to find best normal transformation ;
 %boxcox_survey(
  data    = _temp3, 
  subject = adm_rno, 
  var     = otherbevs,  
  covars  = seq2 r24_weekend  agec_71plus ,  
  weight  = bsw0,  
  print   = N, 
  plot    = N, 
  ntitle  = 4 ); 
 

 
*3) Save <best_lambda> macro variable from <boxcox_survey> for current dietary constituent ;
data _x3;
length varname $ 32;
tran_paramindex=3;
varname="otherbevs";
 tran_lambda=&best_lambda ;
 * The value of best_lambda was ... ;
 *
tran_lambda=0.34;
run;
 
 
 /********************************************************************/
 /*   Dietary constituents #4 - milk: min. amount and lambda   */
 /********************************************************************/
 
* 1) Get minimum non-zero consumption amount ;
proc means data=_tempstratum2 min noprint;
var milk ;
where milk > 0;
output out=_min4(keep=minamount) min=minamount ;
run;
 
data _min4;
set _min4;
minamount=minamount/2;
tran_paramindex=4;
length varname $ 32;
varname="milk";
run;
 
data _temp4;
set inboxcox(keep=adm_rno bsw0 milk seq2 r24_weekend agec_71plus );
if (milk > 0) then output;
run;
 
* 2) Call the <boxcox_survey> macro to find best normal transformation ;
 %boxcox_survey(
  data    = _temp4, 
  subject = adm_rno, 
  var     = milk,  
  covars  = seq2 r24_weekend  agec_71plus ,  
  weight  = bsw0,  
  print   = N, 
  plot    = N, 
  ntitle  = 4 ); 
 

 
*3) Save <best_lambda> macro variable from <boxcox_survey> for current dietary constituent ;
data _x4;
length varname $ 32;
tran_paramindex=4;
varname="milk";
 tran_lambda=&best_lambda ;
 * The value of best_lambda was ... ;
 *
tran_lambda=0.29;
run;
 
 
 /********************************************************************/
 /*   Dietary constituents #5 - rg: min. amount and lambda   */
 /********************************************************************/
 
* 1) Get minimum non-zero consumption amount ;
proc means data=_tempstratum2 min noprint;
var rg ;
where rg > 0;
output out=_min5(keep=minamount) min=minamount ;
run;
 
data _min5;
set _min5;
minamount=minamount/2;
tran_paramindex=5;
length varname $ 32;
varname="rg";
run;
 
data _temp5;
set inboxcox(keep=adm_rno bsw0 rg seq2 r24_weekend agec_71plus );
if (rg > 0) then output;
run;
 
* 2) Call the <boxcox_survey> macro to find best normal transformation ;
 %boxcox_survey(
  data    = _temp5, 
  subject = adm_rno, 
  var     = rg,  
  covars  = seq2 r24_weekend  agec_71plus ,  
  weight  = bsw0,  
  print   = N, 
  plot    = N, 
  ntitle  = 4 ); 
 

 
*3) Save <best_lambda> macro variable from <boxcox_survey> for current dietary constituent ;
data _x5;
length varname $ 32;
tran_paramindex=5;
varname="rg";
 tran_lambda=&best_lambda ;
 * The value of best_lambda was ... ;
 *
tran_lambda=0.25;
run;
 
 
 /********************************************************************/
 /*   Dietary constituents #6 - vf: min. amount and lambda   */
 /********************************************************************/
 
* 1) Get minimum non-zero consumption amount ;
proc means data=_tempstratum2 min noprint;
var vf ;
where vf > 0;
output out=_min6(keep=minamount) min=minamount ;
run;
 
data _min6;
set _min6;
minamount=minamount/2;
tran_paramindex=6;
length varname $ 32;
varname="vf";
run;
 
data _temp6;
set inboxcox(keep=adm_rno bsw0 vf seq2 r24_weekend agec_71plus );
if (vf > 0) then output;
run;
 
* 2) Call the <boxcox_survey> macro to find best normal transformation ;
 %boxcox_survey(
  data    = _temp6, 
  subject = adm_rno, 
  var     = vf,  
  covars  = seq2 r24_weekend  agec_71plus ,  
  weight  = bsw0,  
  print   = N, 
  plot    = N, 
  ntitle  = 4 ); 
 

 
*3) Save <best_lambda> macro variable from <boxcox_survey> for current dietary constituent ;
data _x6;
length varname $ 32;
tran_paramindex=6;
varname="vf";
 tran_lambda=&best_lambda ;
 * The value of best_lambda was ... ;
 *
tran_lambda=0.35;
run;
 
 
 /********************************************************************/
 /*   Dietary constituents #7 - otherfoods: min. amount and lambda   */
 /********************************************************************/
 
* 1) Get minimum non-zero consumption amount ;
proc means data=_tempstratum2 min noprint;
var otherfoods ;
where otherfoods > 0;
output out=_min7(keep=minamount) min=minamount ;
run;
 
data _min7;
set _min7;
minamount=minamount/2;
tran_paramindex=7;
length varname $ 32;
varname="otherfoods";
run;
 
data _temp7;
set inboxcox(keep=adm_rno bsw0 otherfoods seq2 r24_weekend agec_71plus );
if (otherfoods > 0) then output;
run;
 
* 2) Call the <boxcox_survey> macro to find best normal transformation ;
 %boxcox_survey(
  data    = _temp7, 
  subject = adm_rno, 
  var     = otherfoods,  
  covars  = seq2 r24_weekend  agec_71plus ,  
  weight  = bsw0,  
  print   = N, 
  plot    = N, 
  ntitle  = 4 ); 
 

 
*3) Save <best_lambda> macro variable from <boxcox_survey> for current dietary constituent ;
data _x7;
length varname $ 32;
tran_paramindex=7;
varname="otherfoods";
 tran_lambda=&best_lambda ;
 * The value of best_lambda was ... ;
 *
tran_lambda=0.25;
run;
 
 
 /********************************************************************/
 /*   Dietary constituents #8 - pfab: min. amount and lambda   */
 /********************************************************************/
 
* 1) Get minimum non-zero consumption amount ;
proc means data=_tempstratum2 min noprint;
var pfab ;
where pfab > 0;
output out=_min8(keep=minamount) min=minamount ;
run;
 
data _min8;
set _min8;
minamount=minamount/2;
tran_paramindex=8;
length varname $ 32;
varname="pfab";
run;
 
data _temp8;
set inboxcox(keep=adm_rno bsw0 pfab seq2 r24_weekend agec_71plus );
if (pfab > 0) then output;
run;
 
* 2) Call the <boxcox_survey> macro to find best normal transformation ;
 %boxcox_survey(
  data    = _temp8, 
  subject = adm_rno, 
  var     = pfab,  
  covars  = seq2 r24_weekend  agec_71plus ,  
  weight  = bsw0,  
  print   = N, 
  plot    = N, 
  ntitle  = 4 ); 
 

 
*3) Save <best_lambda> macro variable from <boxcox_survey> for current dietary constituent ;
data _x8;
length varname $ 32;
tran_paramindex=8;
varname="pfab";
 tran_lambda=&best_lambda ;
 * The value of best_lambda was ... ;
 *
tran_lambda=0.26;
run;
 
 
 /********************************************************************/
 /*   Dietary constituents #9 - water: min. amount and lambda   */
 /********************************************************************/
 
* 1) Get minimum non-zero consumption amount ;
proc means data=_tempstratum2 min noprint;
var water ;
where water > 0;
output out=_min9(keep=minamount) min=minamount ;
run;
 
data _min9;
set _min9;
minamount=minamount/2;
tran_paramindex=9;
length varname $ 32;
varname="water";
run;
 
data _temp9;
set inboxcox(keep=adm_rno bsw0 water seq2 r24_weekend agec_71plus );
if (water > 0) then output;
run;
 
* 2) Call the <boxcox_survey> macro to find best normal transformation ;
 %boxcox_survey(
  data    = _temp9, 
  subject = adm_rno, 
  var     = water,  
  covars  = seq2 r24_weekend  agec_71plus ,  
  weight  = bsw0,  
  print   = N, 
  plot    = N, 
  ntitle  = 4 ); 
 

 
*3) Save <best_lambda> macro variable from <boxcox_survey> for current dietary constituent ;
data _x9;
length varname $ 32;
tran_paramindex=9;
varname="water";
 tran_lambda=&best_lambda ;
 * The value of best_lambda was ... ;
 *
tran_lambda=0.3;
run;
 
 
 /********************************************************************/
 /*   Dietary constituents #10 - mufa: min. amount and lambda   */
 /********************************************************************/
 
* 1) Get minimum non-zero consumption amount ;
proc means data=_tempstratum2 min noprint;
var mufa ;
where mufa > 0;
output out=_min10(keep=minamount) min=minamount ;
run;
 
data _min10;
set _min10;
minamount=minamount/2;
tran_paramindex=10;
length varname $ 32;
varname="mufa";
run;
 
data _temp10;
set inboxcox(keep=adm_rno bsw0 mufa seq2 r24_weekend agec_71plus );
if (mufa > 0) then output;
run;
 
* 2) Call the <boxcox_survey> macro to find best normal transformation ;
 %boxcox_survey(
  data    = _temp10, 
  subject = adm_rno, 
  var     = mufa,  
  covars  = seq2 r24_weekend  agec_71plus ,  
  weight  = bsw0,  
  print   = N, 
  plot    = N, 
  ntitle  = 4 ); 
 

 
*3) Save <best_lambda> macro variable from <boxcox_survey> for current dietary constituent ;
data _x10;
length varname $ 32;
tran_paramindex=10;
varname="mufa";
 tran_lambda=&best_lambda ;
 * The value of best_lambda was ... ;
 *
tran_lambda=0.29;
run;
 
 
 /********************************************************************/
 /*   Dietary constituents #11 - pufa: min. amount and lambda   */
 /********************************************************************/
 
* 1) Get minimum non-zero consumption amount ;
proc means data=_tempstratum2 min noprint;
var pufa ;
where pufa > 0;
output out=_min11(keep=minamount) min=minamount ;
run;
 
data _min11;
set _min11;
minamount=minamount/2;
tran_paramindex=11;
length varname $ 32;
varname="pufa";
run;
 
data _temp11;
set inboxcox(keep=adm_rno bsw0 pufa seq2 r24_weekend agec_71plus );
if (pufa > 0) then output;
run;
 
* 2) Call the <boxcox_survey> macro to find best normal transformation ;
 %boxcox_survey(
  data    = _temp11, 
  subject = adm_rno, 
  var     = pufa,  
  covars  = seq2 r24_weekend  agec_71plus ,  
  weight  = bsw0,  
  print   = N, 
  plot    = N, 
  ntitle  = 4 ); 
 

 
*3) Save <best_lambda> macro variable from <boxcox_survey> for current dietary constituent ;
data _x11;
length varname $ 32;
tran_paramindex=11;
varname="pufa";
 tran_lambda=&best_lambda ;
 * The value of best_lambda was ... ;
 *
tran_lambda=0.26;
run;
 
 
 /********************************************************************/
 /*   Dietary constituents #12 - sfa: min. amount and lambda   */
 /********************************************************************/
 
* 1) Get minimum non-zero consumption amount ;
proc means data=_tempstratum2 min noprint;
var sfa ;
where sfa > 0;
output out=_min12(keep=minamount) min=minamount ;
run;
 
data _min12;
set _min12;
minamount=minamount/2;
tran_paramindex=12;
length varname $ 32;
varname="sfa";
run;
 
data _temp12;
set inboxcox(keep=adm_rno bsw0 sfa seq2 r24_weekend agec_71plus );
if (sfa > 0) then output;
run;
 
* 2) Call the <boxcox_survey> macro to find best normal transformation ;
 %boxcox_survey(
  data    = _temp12, 
  subject = adm_rno, 
  var     = sfa,  
  covars  = seq2 r24_weekend  agec_71plus ,  
  weight  = bsw0,  
  print   = N, 
  plot    = N, 
  ntitle  = 4 ); 
 

 
*3) Save <best_lambda> macro variable from <boxcox_survey> for current dietary constituent ;
data _x12;
length varname $ 32;
tran_paramindex=12;
varname="sfa";
 tran_lambda=&best_lambda ;
 * The value of best_lambda was ... ;
 *
tran_lambda=0.17;
run;
 
 
 /********************************************************************/
 /*   Dietary constituents #13 - freesugars: min. amount and lambda   */
 /********************************************************************/
 
* 1) Get minimum non-zero consumption amount ;
proc means data=_tempstratum2 min noprint;
var freesugars ;
where freesugars > 0;
output out=_min13(keep=minamount) min=minamount ;
run;
 
data _min13;
set _min13;
minamount=minamount/2;
tran_paramindex=13;
length varname $ 32;
varname="freesugars";
run;
 
data _temp13;
set inboxcox(keep=adm_rno bsw0 freesugars seq2 r24_weekend agec_71plus );
if (freesugars > 0) then output;
run;
 
* 2) Call the <boxcox_survey> macro to find best normal transformation ;
 %boxcox_survey(
  data    = _temp13, 
  subject = adm_rno, 
  var     = freesugars,  
  covars  = seq2 r24_weekend  agec_71plus ,  
  weight  = bsw0,  
  print   = N, 
  plot    = N, 
  ntitle  = 4 ); 
 

 
*3) Save <best_lambda> macro variable from <boxcox_survey> for current dietary constituent ;
data _x13;
length varname $ 32;
tran_paramindex=13;
varname="freesugars";
 tran_lambda=&best_lambda ;
 * The value of best_lambda was ... ;
 *
tran_lambda=0.43;
run;
 
 
 /********************************************************************/
 /*   Dietary constituents #14 - sodium: min. amount and lambda   */
 /********************************************************************/
 
* 1) Get minimum non-zero consumption amount ;
proc means data=_tempstratum2 min noprint;
var sodium ;
where sodium > 0;
output out=_min14(keep=minamount) min=minamount ;
run;
 
data _min14;
set _min14;
minamount=minamount/2;
tran_paramindex=14;
length varname $ 32;
varname="sodium";
run;
 
data _temp14;
set inboxcox(keep=adm_rno bsw0 sodium seq2 r24_weekend agec_71plus );
if (sodium > 0) then output;
run;
 
* 2) Call the <boxcox_survey> macro to find best normal transformation ;
 %boxcox_survey(
  data    = _temp14, 
  subject = adm_rno, 
  var     = sodium,  
  covars  = seq2 r24_weekend  agec_71plus ,  
  weight  = bsw0,  
  print   = N, 
  plot    = N, 
  ntitle  = 4 ); 
 

 
*3) Save <best_lambda> macro variable from <boxcox_survey> for current dietary constituent ;
data _x14;
length varname $ 32;
tran_paramindex=14;
varname="sodium";
 tran_lambda=&best_lambda ;
 * The value of best_lambda was ... ;
 *
tran_lambda=0.24;
run;
 
 
 /********************************************************************/
 /*   Dietary constituents #15 - energy: min. amount and lambda   */
 /********************************************************************/
 
* 1) Get minimum non-zero consumption amount ;
proc means data=_tempstratum2 min noprint;
var energy ;
where energy > 0;
output out=_min15(keep=minamount) min=minamount ;
run;
 
data _min15;
set _min15;
minamount=minamount/2;
tran_paramindex=15;
length varname $ 32;
varname="energy";
run;
 
data _temp15;
set inboxcox(keep=adm_rno bsw0 energy seq2 r24_weekend agec_71plus );
if (energy > 0) then output;
run;
 
* 2) Call the <boxcox_survey> macro to find best normal transformation ;
 %boxcox_survey(
  data    = _temp15, 
  subject = adm_rno, 
  var     = energy,  
  covars  = seq2 r24_weekend  agec_71plus ,  
  weight  = bsw0,  
  print   = N, 
  plot    = N, 
  ntitle  = 4 ); 
 

 
*3) Save <best_lambda> macro variable from <boxcox_survey> for current dietary constituent ;
data _x15;
length varname $ 32;
tran_paramindex=15;
varname="energy";
 tran_lambda=&best_lambda ;
 * The value of best_lambda was ... ;
 *
tran_lambda=0.41;
run;
 
 
 /********************************************************************/
 /*   Dietary constituents #16 - calcium: min. amount and lambda   */
 /********************************************************************/
 
* 1) Get minimum non-zero consumption amount ;
proc means data=_tempstratum2 min noprint;
var calcium ;
where calcium > 0;
output out=_min16(keep=minamount) min=minamount ;
run;
 
data _min16;
set _min16;
minamount=minamount/2;
tran_paramindex=16;
length varname $ 32;
varname="calcium";
run;
 
data _temp16;
set inboxcox(keep=adm_rno bsw0 calcium seq2 r24_weekend agec_71plus );
if (calcium > 0) then output;
run;
 
* 2) Call the <boxcox_survey> macro to find best normal transformation ;
 %boxcox_survey(
  data    = _temp16, 
  subject = adm_rno, 
  var     = calcium,  
  covars  = seq2 r24_weekend  agec_71plus ,  
  weight  = bsw0,  
  print   = N, 
  plot    = N, 
  ntitle  = 4 ); 
 

 
*3) Save <best_lambda> macro variable from <boxcox_survey> for current dietary constituent ;
data _x16;
length varname $ 32;
tran_paramindex=16;
varname="calcium";
 tran_lambda=&best_lambda ;
 * The value of best_lambda was ... ;
 *
tran_lambda=0.23;
run;
 
 
 /********************************************************************/
 /*   Dietary constituents #17 - vit_d: min. amount and lambda   */
 /********************************************************************/
 
* 1) Get minimum non-zero consumption amount ;
proc means data=_tempstratum2 min noprint;
var vit_d ;
where vit_d > 0;
output out=_min17(keep=minamount) min=minamount ;
run;
 
data _min17;
set _min17;
minamount=minamount/2;
tran_paramindex=17;
length varname $ 32;
varname="vit_d";
run;
 
data _temp17;
set inboxcox(keep=adm_rno bsw0 vit_d seq2 r24_weekend agec_71plus );
if (vit_d > 0) then output;
run;
 
* 2) Call the <boxcox_survey> macro to find best normal transformation ;
 %boxcox_survey(
  data    = _temp17, 
  subject = adm_rno, 
  var     = vit_d,  
  covars  = seq2 r24_weekend  agec_71plus ,  
  weight  = bsw0,  
  print   = N, 
  plot    = N, 
  ntitle  = 4 ); 
 

 
*3) Save <best_lambda> macro variable from <boxcox_survey> for current dietary constituent ;
data _x17;
length varname $ 32;
tran_paramindex=17;
varname="vit_d";
 tran_lambda=&best_lambda ;
 * The value of best_lambda was ... ;
 *
tran_lambda=0.27;
run;
 
 /*****************************************************************/
 /*         Append data sets of all dietary constituents          */
 /*****************************************************************/
 
	data work.xlambdas_s2_0;
	set _x1-_x17;
	run;
	 
	data work.xminamount_s2_0;
	set _min1-_min17;
	label minamount= "Min. non zero amount divided by 2";
	run;
	 
 
 /********************************************************************/
 /*   Standardize data according to lambda and min. amount values    */
 /********************************************************************/
 
  * Call the <std_cov_boxcox24hr_conday_minamt> macro to standardize data and continuous covariates, if any; 
 %std_cov_boxcox24hr_conday_minamt(
  data                       = _tempstratum2, 
  prestand_continuous_covars =  , 
  rec24hr_epis_vars          = wg pfpb otherbevs milk rg, 
  rec24hr_daily_vars         = vf otherfoods pfab water mufa pufa sfa freesugars sodium energy calcium vit_d, 
  boxcox_tran_lambda_data    = xlambdas_s2_0, 
  minamount_data             = xminamount_s2_0, 
  print                      = y, 
  titles                     = 3 ); 
 
 
/*************************************************************************/
/* Fit the multivariate measurement error model (multivar_mcmc)          */
/*************************************************************************/
 
 
* Call the <multivar_mcmc> macro to fit the measurement error model; 
   title1 "Fit Multivariate Measurement Error Model Using MCMC with 24-Hour Recall as Main Instrument"; 
 %multivar_mcmc(
  data                        = stdcov_stdbc24hr_conday_out, 
  subject                     = adm_rno, 
  weight_var                  = bsw0 , 
  repeat                      = suppid, 
  conday_epis_vars            = conday_wg  conday_pfpb  conday_otherbevs  conday_milk  conday_rg, 
  gst_rec24hr_epis_vars       = stdbc_wg  stdbc_pfpb  stdbc_otherbevs  stdbc_milk  stdbc_rg, 
  gst_rec24hr_daily_vars      = stdbc_vf  stdbc_otherfoods  stdbc_pfab  stdbc_water  stdbc_mufa  stdbc_pufa  stdbc_sfa  
 stdbc_freesugars  stdbc_sodium  stdbc_energy  stdbc_calcium  stdbc_vit_d, 
  covars_epis_prob            = constant1 seq2 r24_weekend agec_71plus   , 
  covars_epis_amt             = constant1 seq2 r24_weekend agec_71plus   , 
  covars_daily_amt            = constant1 seq2 r24_weekend agec_71plus   , 
  set_seed_mcmc               = 42941, 
  set_number_mcmc_iterations  = 8000,
  set_number_burn_iterations  = 3000,
  set_thin                    = 10,
  prior_sigmau_mean_data      = , 
  sigmau_constant             = , 
  gen_inverse                 = y, 
  print                       = y, 
  titles                      = 1, 
  std_print_store             = y, 
  notes_print                 = y, 
  out_lib                     = baselib, 
  out_store_label             = mcmc_s2_rep0, 
  out_save_label_max5char     = s20, 
  set_number_saved_out_data   = , 
  save_mcmc_u_out_data        = y, 
  set_number_post_mcmc_u_out  = , 
  traceplots_method1_gpath    = , 
  traceplots_method2_file_pdf = trace_rep0_s2.pdf, 
  optional_iml_store_data     = backtran_out, 
  optional_iml_store_names    = constant1 seq2 r24_weekend agec_71plus   tran_paramindex tran_lambda tran_center tran_scale 
 minamount 
  ); 
 
 
* Save lambdas (Box-Cox transformation lambda values) and minimum amount data ;
data baselib.backtran_out0_s2;
retain replicate sex;
set work.backtran_out;
* indicate replicate number, current stratum value, variable name;
replicate = 0;
sex = 2;
length varname $ 32 ;
if _N_=1 then varname = "wg" ;
if _N_=2 then varname = "pfpb" ;
if _N_=3 then varname = "otherbevs" ;
if _N_=4 then varname = "milk" ;
if _N_=5 then varname = "rg" ;
if _N_=6 then varname = "vf" ;
if _N_=7 then varname = "otherfoods" ;
if _N_=8 then varname = "pfab" ;
if _N_=9 then varname = "water" ;
if _N_=10 then varname = "mufa" ;
if _N_=11 then varname = "pufa" ;
if _N_=12 then varname = "sfa" ;
if _N_=13 then varname = "freesugars" ;
if _N_=14 then varname = "sodium" ;
if _N_=15 then varname = "energy" ;
if _N_=16 then varname = "calcium" ;
if _N_=17 then varname = "vit_d" ;
run;
 
 
/*************************************************************************/
/* Simulation of pseudo-individual (multivar_distrib)                    */
/*************************************************************************/
 
* Prepare an input data for the <optional_input_data> option in <multivar_distrib>;
	proc sort data=_tempstratum2 nodupkey out=optional_input_data(keep= adm_rno sex bsw0 
	 agec_71plus sex drig );
	by adm_rno ;
	run;
	 

 
* Call the <multivar_distrib> macro to simulate usual intakes for pseudo-individuals; 
 %multivar_distrib(
  multivar_mcmc_out_lib           = baselib ,  
  multivar_mcmc_out_store_label   = mcmc_s2_rep0, 
  t_weightavg_covariates_list1    = constant1 constant0 constant0 agec_71plus   ,  
  t_weightavg_covariates_list2    = constant1 constant0 constant1 agec_71plus  , 
  set_value_for_weight_cov_list1  = 4, 
  set_value_for_weight_cov_list2  = 3, 
  optional_input_data             = optional_input_data , 
  optional_input_data_var_list    = , 
  optional_input_mcmc_u_out_data  = , 
  additional_output_var_list      = sex sex drig  , 
  additional_output_subject_var   = adm_rno , 
  output_mcmc_weight_var          = y  , 
  set_seed_distrib                = 89009890, 
  set_number_monte_carlo_rand_obs = 500,  
  print                           = y 
  ); 
 
 
* Save the Monte Carlo simulation data for current stratum;
	data baselib.mc_t_distrib_out0_s2;
	set mc_t_distrib_out;
	run;
	 
* delete temporary data sets ;
	proc datasets lib=work nolist nodetails;
	delete mc_t_distrib_out optional_input_data stdcov_stdbc24hr_conday_out xlambdas_: 
	 xminamount_: backtran_out ;
	run;
	 
	proc datasets lib=baselib nolist nodetails ;
	delete multivar_mcmc_samples_u_outs20 ;
	run;
 
/*************************************************************************/
/* Combine stratum-specific data into one (i.e., Monte Carlo simulations)*/
/*************************************************************************/
 
* Append all stratum-specific <mc_t_distrib_out> data ;
	data baselib.mc_t_distrib_out0 ;
	set baselib.mc_t_distrib_out0_s1-baselib.mc_t_distrib_out0_s2;
	run;
	 
* Append all stratum-specific <backtran_out> data ;
	data baselib.backtran_out0;
	set baselib.backtran_out0_s1-baselib.backtran_out0_s2;
	run;
	 
* Delete temporary data;
	proc datasets lib=baselib nolist nodetails;
	delete backtran_out0_s1-backtran_out0_s2 ;
	run;
	 
	proc datasets lib=work nolist nodetails;
	delete _tempstratum1-_tempstratum2 CKNEGATIVE24HR LAMBDA_MINAMT_MULTIREC NOTNEG24HR_DATA 
	 MCMC_SUBJ1RECDATA trace_afterburn_thin_data sigmae_paneldata sigmau_paneldata beta_paneldata PARAMVARLABELS1REC _SGSORT_ ;
	run;
 
 /*********************  End of measurement error model  ******************/

 
/* end of code */
