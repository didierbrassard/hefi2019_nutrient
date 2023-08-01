 /*************************************************************************/
 /*                                                                       */
 /*                     CCHS 2015 - Nutrition (PUMF)                      */
 /*                                                                       */
 /* Output results: prepare data, linear reg, logistic reg, distribution  */
 /*                       Code for outcome: PROTEIN                       */
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
	%let suffix = _pro ;
	
/* NCI macros */
	/* percentiles*/ %include "&path./Macros/percentiles_survey.macro.v1.1.sas";

	/* Available at: https://prevention.cancer.gov/research-groups/biometry/measurement-error-impact/software-measurement-error/several-regularly-consumed-or-0 */

/* HEFI-2019 scoring macro */
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
 /*  Format the output data of <multivar_mcmc>                            */
 /*                                                                       */
 /*************************************************************************/

	data baselib._usintake_mc_t_out0(drop=mc_t1-mc_t16 i );
		set baselib.mc_t_distrib_out0;
		 
	* Divide weights by set_number_monte_carlo_rand_obs value used in the MULTIVAR_DISTRIB macro ;
		weight_nw_sumw_div = weight_nw_sumw / 500 ;
		
	* Assign variable names for usual intakes from the MULTIVAR_DISTRIB macro ;
		array outmc (*) mc_t1-mc_t16 ;
		array clean (*) wg pfpb otherbevs milk rg vf otherfoods pfab water mufa pufa 
			sfa freesugars sodium energy proteins;
		do i=1 to dim(outmc);
			clean(i)=outmc(i);
		end;
		 
	* Additional formatting (derived variables, binary cutpoints);
		proteins_kg = proteins / bodyweight ;
		 
		if not missing(proteins_kg) then do;
		 
			if proteins_kg < 0.66 then binary_06 =1;
			else binary_06 = 0;
		
			if proteins_kg < 0.8 then binary_08 =1;
				else binary_08 = 0;
		
			if proteins_kg < 1.0 then binary_10 =1;
				else binary_10 = 0;
			
			if proteins_kg < 1.2 then binary_12 =1;
				else binary_12 = 0;
			
			if proteins_kg < 1.4 then binary_14 =1;
				else binary_14 = 0;
			 
		end;
		 
	label 
		binary_06 = "Intake < than EAR (0.66g/kg)"
		binary_08 = "Intake < than 0.8g/kg"
		binary_10 = "Intake < than 1.0g/kg"
		binary_12 = "Intake < than 1.2g/kg"
		binary_14 = "Intake < than 1.4g/kg"
		;
	run;

 /*************************************************************************/
 /*                                                                       */
 /*  Calculate `usual` hefi-2019 scores                                   */
 /*                                                                       */
 /*************************************************************************/

/* use scoring algorithm to calculate HEFI19 score */
	%HEFI2019(indata          = baselib._usintake_mc_t_out0,
			  vegfruits       = vf,
			  wholegrfoods	  = wg,
			  nonwholegrfoods = rg,
			  profoodsanimal  = pfab,
			  profoodsplant	  = pfpb,
			  otherfoods	  = otherfoods,
			  waterhealthybev = water,
			  unsweetmilk	  = milk,
			  unsweetplantbevpro = 0,
			  otherbeverages  = otherbevs,
			  mufat			  = mufa,
			  pufat			  = pufa,
			  satfat		  = sfa,
			  freesugars	  = freesugars,
			  sodium		  = sodium,
			  energy		  = energy,
			  outdata         = baselib._usintake_mc_t_out0
		  	);

* indicate comparison percentile of HEFI2019_TOTAL_SCORE for all models ;
	%let t0perc = 50 ;
	%let t1perc = 90 ;

 /*************************************************************************/
 /*                                                                       */
 /*  Linear model: (Nutrient intake | HEFI2019)                           */
 /*                                                                       */
 /*************************************************************************/

  * Call the <percentiles_Survey> macro to output distribution of usual intakes; 
 %percentiles_Survey(
 data      = baselib._usintake_mc_t_out0, 
 byvar     = , 
 var       = HEFI2019_TOTAL_SCORE , 
 weight    = weight_nw_sumw_div, 
 cutpoints = , 
 print     = N, 
 ntitle    = 3 
  ); 
 
 
data Reg.pctHEFI2019_TOTAL_SCORE_w0;
	retain varname;
	set _percentiles(drop=variance);
	varname="HEFI2019_TOTAL_SCORE";
run;

proc datasets lib=work nolist nodetails;
	delete _percentiles _percentiles2;
	run;

data _null_;
	set Reg.pctHEFI2019_TOTAL_SCORE_w0;
	* data saved above ;
	* output value of t1 and t0;
	call symputx("t1", Pctile&t1perc);
	call symputx("t0", Pctile&t0perc);
	* output p1 and p99 to restrict the regression plot to actual values ;
	call symputx("p1", Pctile1);
	call symputx("p99", Pctile99);
	* define knot values based on distribution ;
	call symputx("k1", put(Pctile5, 5.2));
	call symputx("k2", put(Pctile27, 5.2));
	call symputx("k3", put(Pctile50, 5.2));
	call symputx("k4", put(Pctile73, 5.2));
	call symputx("k5", put(Pctile95, 5.2));
run;
 
/*************************************************************************/
/* Linear regression model (proteins = rcs(HEFI2019_TOTAL_SCORE,5)  )    */
/*************************************************************************/

proc surveyreg data=baselib._usintake_mc_t_out0;
	effect spl_x=spline(HEFI2019_TOTAL_SCORE / details naturalcubic 
		basis=tpf(noint) knotmethod=list(31.69 41.48 46.79 51.66 58.82));
	model proteins=spl_x / adjrsq;
	weight weight_nw_sumw_div;
	estimate "proteins|rcs(HEFI2019_TOTAL_SCORE,5), high(p90) vs. low(p50)" 
		spl_x [1 , 56.442513934][-1 , 46.785789336];
	ods output Estimates=_param(keep=Label Estimate) 
		fitstatistics=_fit(where=(Label1 in ('Adjusted R-Square' 'R carré ajusté')));
	store Reg.rc_rcs_hefi_raw;
run;

ods graphics on;

proc plm restore=Reg.rc_rcs_hefi_raw noinfo noclprint;
	effectplot fit(x=HEFI2019_TOTAL_SCORE) / clm;
	ods output FitPlot=Reg.plmrcst_hefi_raw0(keep=_xcont1 _predicted 
		rename=(_xcont1=HEFI2019_TOTAL_SCORE _predicted=predicted) );
run;

ods graphics off;
* Create data for independant variables ;

data _t1t0;
	HEFI2019_TOTAL_SCORE=&t0;
		* HEFI2019_TOTAL_SCORE = 46.785789336;
	tperc=&t0perc;
		* tperc = 50;
	output;
	HEFI2019_TOTAL_SCORE=&t1;
		* HEFI2019_TOTAL_SCORE = 56.442513934;
	tperc=&t1perc;
		* tperc = 90;
	output;
run;

proc plm restore=Reg.rc_rcs_hefi_raw noinfo noclprint;
	score data=_t1t0 out=_t1t0 predicted;
run;

proc transpose data=_t1t0 out=_t1t0w(drop=_: rename=(COL1=t0_score 
		COL2=t1_score)) prefix=COL;
	var predicted;
run;

data _null_;
	set _fit;
	call symputx("r2", nvalue1);
run;

proc transpose data=_param out=paramt(drop=_:) prefix=parm;
	var Estimate;
	copy label;
run;

data Reg.rcs_hefi_raw0;
	retain replicate label z;
	merge paramt(rename=(parm1=rc_beta) ) _t1t0w;
	* save replicate number ;
	replicate=0;
	* Save the comparison values;
	t0=&t0;
		* t0 = 46.785789336;
	t1=&t1;
		* t1 = 56.442513934;
	* save r2 ;
	r2=&r2;
		* r2 = 0.0215047484;
	* save p1, p99 plot limit ;
	p1=&p1;
		* p1=25.284244094;
	p99=&p99;
		* p99=62.880843051;
	* save model covariates, if need be;
	z="NA";
	label z="Model did not include covariates";
	* label variables;
	label t1="90th perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		t0="50th perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		t1_score="proteins value at 90th perc. of usual HEFI2019_TOTAL_SCORE" 
		t0_score="proteins value at 50th perc. of usual HEFI2019_TOTAL_SCORE" 
		p1="1st perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		p99="99th perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		rc_beta="(proteins change | p90-p50 HEFI2019_TOTAL_SCORE diff.)" 
		r2="Adjusted r-square";
run;

proc transpose data=Reg.plmrcst_hefi_raw0 out=plot0 (rename=(_NAME_=name)) 
		prefix=pred;
	var HEFI2019_TOTAL_SCORE: predicted;
run;

data Reg.plm_rcs_hefi_raw0;
	retain replicate;
	set plot0(drop=_:);
	replicate=0;
	label name=" ";
run;

 /*************************************************************************/
 /*                                                                       */
 /*  Log model: (Pr(X<x) | HEFI2019)                                      */
 /*                                                                       */
 /*************************************************************************/

/*************************************************************************/
/* Logistic reg. model (logit(binary_06) = rcs(HEFI2019_TOTAL_SCORE,5)   */
/*************************************************************************/

proc logistic data=baselib._usintake_mc_t_out0 descending;
	effect spl_x=spline(HEFI2019_TOTAL_SCORE / details naturalcubic 
		basis=tpf(noint) knotmethod=list(31.69 41.48 46.79 51.66 58.82));
	model binary_06=spl_x;
	weight weight_nw_sumw_div;
	estimate "binary_06|rcs(HEFI2019_TOTAL_SCORE,5), high(p90) vs. low(p50)" 
		spl_x [1 , 56.442513934][-1 , 46.785789336] / exp;
	ods output Estimates=_param(keep=Label Estimate ExpEstimate) 
		ParameterEstimates=_modelparms(Keep=Variable Estimate);
	store reg.log_rcs_hefi_y06_raw;
run;

* output model intercept ;

data _NULL_;
	set _modelparms;

	if _N_=1 then
		call symputx("log_beta0", Estimate);
run;

* create curve of predicted probabilities at each value of x ;
ods graphics on;

proc plm restore=reg.log_rcs_hefi_y06_raw noinfo noclprint;
	effectplot fit(x=HEFI2019_TOTAL_SCORE) / clm ilink;
	ods output FitPlot=reg.logrrcst_hefi_y06_raw0(keep=_xcont1 _predicted 
		rename=(_xcont1=HEFI2019_TOTAL_SCORE _predicted=pr_y) );
run;

ods graphics off;
* Create data for independent variables ;

data _t1t0;
	HEFI2019_TOTAL_SCORE=&t0;
	*
HEFI2019_TOTAL_SCORE = 46.785789336;
	tperc=&t0perc;
	*
tperc = 50;
	output;
	HEFI2019_TOTAL_SCORE=&t1;
	*
HEFI2019_TOTAL_SCORE = 56.442513934;
	tperc=&t1perc;
	*
tperc = 90;
	output;
run;

* Calculate Pr(Y) for t1 and t0 values ;

proc plm restore=reg.log_rcs_hefi_y06_raw noinfo noclprint;
	score data=_t1t0 out=_t1t0 predicted / ilink;
	*ilink option requests Pr(Y);
run;

proc transpose data=_t1t0 out=_t1t0w(drop=_: rename=(COL1=y0 COL2=y1)) 
		prefix=COL;
	var predicted;
run;

* save model information and result of t1-t0 contrast;

data reg.logr_mod_hefi_y06_raw0;
	retain replicate label z t0 t1 log_beta exp_y0 exp_y1 exp_beta y0 y1 risk_diff 
		risk_ratio;
	merge _param(rename=(estimate=log_beta expestimate=exp_beta) ) _t1t0w;
	* save replicate number ;
	replicate=0;
	* save the comparison values;
	t0=&t0;
	*
t0 = 46.785789336;
	t1=&t1;
	*
t1 = 56.442513934;
	* save p1, p99 plot limit ;
	p1=&p1;
	*
p1=25.284244094;
	p99=&p99;
	*
p99=62.880843051;
	* save model covariates, if need be;
	z="NA";
	label z="Model did not include covariates";
	* calculate odds(y) at the percentiles;
	exp_y1=y1/(1-y1);
	exp_y0=y0/(1-y0);
	* calculate absolute and relative risk;
	risk_diff=y1 - y0;
	risk_ratio=y1 / y0;
	* input model intercept ;
	log_beta0=&log_beta0;
	*
log_beta0 = 0.3437837288 ;
	* label variables;
	label t1="90th perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		t0="50th perc. of usual  HEFI2019_TOTAL_SCORE intake distribution" 
		y1="Pr(binary_06) at 90th perc. of usual HEFI2019_TOTAL_SCORE" 
		y0="Pr(binary_06) at 50th perc. of usual HEFI2019_TOTAL_SCORE" 
		exp_y1="odds(binary_06) at 90th perc. of usual HEFI2019_TOTAL_SCORE" 
		exp_y0="odds(binary_06) at 50th perc. of usual HEFI2019_TOTAL_SCORE" 
		p1="1st perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		p99="99th perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		log_beta="logit(binary_06) | p90-p50 HEFI2019_TOTAL_SCORE  diff." 
		exp_beta="odds(binary_06) | p90-p50 HEFI2019_TOTAL_SCORE diff." 
		risk_diff="risk diff. binary_06 | p90-p50 HEFI2019_TOTAL_SCORE diff." 
		risk_ratio="risk ratio binary_06 | p90-p50 HEFI2019_TOTAL_SCORE diff." 
		log_beta0="Model intercept (ln scale)";
run;

* transpose predicted probability curve from long to wide ;

proc transpose data=reg.logrrcst_hefi_y06_raw0 out=plot0 (rename=(_NAME_=name)) 
		prefix=pred;
	var HEFI2019_TOTAL_SCORE: pr_y;
run;

data reg.logrrcsw_hefi_y06_raw0;
	retain replicate;
	set plot0(drop=_:);
	replicate=0;
	label name=" ";
run;

* note: logrrcsw = LOGistic Regression with Restricted Cubic Spline plot in the Wide format ;
 
/*************************************************************************/
/* Logistic reg. model (logit(binary_08) = rcs(HEFI2019_TOTAL_SCORE,5)   */
/*************************************************************************/

proc logistic data=baselib._usintake_mc_t_out0 descending;
	effect spl_x=spline(HEFI2019_TOTAL_SCORE / details naturalcubic 
		basis=tpf(noint) knotmethod=list(31.69 41.48 46.79 51.66 58.82));
	model binary_08=spl_x;
	weight weight_nw_sumw_div;
	estimate "binary_08|rcs(HEFI2019_TOTAL_SCORE,5), high(p90) vs. low(p50)" 
		spl_x [1 , 56.442513934][-1 , 46.785789336] / exp;
	ods output Estimates=_param(keep=Label Estimate ExpEstimate) 
		ParameterEstimates=_modelparms(Keep=Variable Estimate);
	store reg.log_rcs_hefi_y08_raw;
run;

* output model intercept ;

data _NULL_;
	set _modelparms;

	if _N_=1 then
		call symputx("log_beta0", Estimate);
run;

* create curve of predicted probabilities at each value of x ;
ods graphics on;

proc plm restore=reg.log_rcs_hefi_y08_raw noinfo noclprint;
	effectplot fit(x=HEFI2019_TOTAL_SCORE) / clm ilink;
	ods output FitPlot=reg.logrrcst_hefi_y08_raw0(keep=_xcont1 _predicted 
		rename=(_xcont1=HEFI2019_TOTAL_SCORE _predicted=pr_y) );
run;

ods graphics off;
* Create data for independent variables ;

data _t1t0;
	HEFI2019_TOTAL_SCORE=&t0;
	*
HEFI2019_TOTAL_SCORE = 46.785789336;
	tperc=&t0perc;
	*
tperc = 50;
	output;
	HEFI2019_TOTAL_SCORE=&t1;
	*
HEFI2019_TOTAL_SCORE = 56.442513934;
	tperc=&t1perc;
	*
tperc = 90;
	output;
run;

* Calculate Pr(Y) for t1 and t0 values ;

proc plm restore=reg.log_rcs_hefi_y08_raw noinfo noclprint;
	score data=_t1t0 out=_t1t0 predicted / ilink;
	*ilink option requests Pr(Y);
run;

proc transpose data=_t1t0 out=_t1t0w(drop=_: rename=(COL1=y0 COL2=y1)) 
		prefix=COL;
	var predicted;
run;

* save model information and result of t1-t0 contrast;

data reg.logr_mod_hefi_y08_raw0;
	retain replicate label z t0 t1 log_beta exp_y0 exp_y1 exp_beta y0 y1 risk_diff 
		risk_ratio;
	merge _param(rename=(estimate=log_beta expestimate=exp_beta) ) _t1t0w;
	* save replicate number ;
	replicate=0;
	* save the comparison values;
	t0=&t0;
	*
t0 = 46.785789336;
	t1=&t1;
	*
t1 = 56.442513934;
	* save p1, p99 plot limit ;
	p1=&p1;
	*
p1=25.284244094;
	p99=&p99;
	*
p99=62.880843051;
	* save model covariates, if need be;
	z="NA";
	label z="Model did not include covariates";
	* calculate odds(y) at the percentiles;
	exp_y1=y1/(1-y1);
	exp_y0=y0/(1-y0);
	* calculate absolute and relative risk;
	risk_diff=y1 - y0;
	risk_ratio=y1 / y0;
	* input model intercept ;
	log_beta0=&log_beta0;
	*
log_beta0 = 1.698808694 ;
	* label variables;
	label t1="90th perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		t0="50th perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		y1="Pr(binary_08) at 90th perc. of usual HEFI2019_TOTAL_SCORE" 
		y0="Pr(binary_08) at 50th perc. of usual HEFI2019_TOTAL_SCORE" 
		exp_y1="odds(binary_08) at 90th perc. of usual HEFI2019_TOTAL_SCORE" 
		exp_y0="odds(binary_08) at 50th perc. of usual HEFI2019_TOTAL_SCORE" 
		p1="1st perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		p99="99th perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		log_beta="logit(binary_08) | p90-p50 HEFI2019_TOTAL_SCORE diff." 
		exp_beta="odds(binary_08) | p90-p50 HEFI2019_TOTAL_SCORE diff." 
		risk_diff="risk diff. binary_08 | p90-p50 HEFI2019_TOTAL_SCORE diff." 
		risk_ratio="risk ratio binary_08 | p90-p50 HEFI2019_TOTAL_SCORE diff." 
		log_beta0="Model intercept (ln scale)";
run;

* transpose predicted probability curve from long to wide ;

proc transpose data=reg.logrrcst_hefi_y08_raw0 out=plot0 (rename=(_NAME_=name)) 
		prefix=pred;
	var HEFI2019_TOTAL_SCORE: pr_y;
run;

data reg.logrrcsw_hefi_y08_raw0;
	retain replicate;
	set plot0(drop=_:);
	replicate=0;
	label name=" ";
run;

* note: logrrcsw = LOGistic Regression with Restricted Cubic Spline plot in the Wide format ;

/*************************************************************************/
/* Logistic reg. model (logit(binary_10) = rcs(HEFI2019_TOTAL_SCORE,5)   */
/*************************************************************************/

proc logistic data=baselib._usintake_mc_t_out0 descending;
	effect spl_x=spline(HEFI2019_TOTAL_SCORE / details naturalcubic 
		basis=tpf(noint) knotmethod=list(31.69 41.48 46.79 51.66 58.82));
	model binary_10=spl_x;
	weight weight_nw_sumw_div;
	estimate "binary_10|rcs(HEFI2019_TOTAL_SCORE,5), high(p90) vs. low(p50)" 
		spl_x [1 , 56.442513934][-1 , 46.785789336] / exp;
	ods output Estimates=_param(keep=Label Estimate ExpEstimate) 
		ParameterEstimates=_modelparms(Keep=Variable Estimate);
	store reg.log_rcs_hefi_y10_raw;
run;

* output model intercept ;

data _NULL_;
	set _modelparms;

	if _N_=1 then
		call symputx("log_beta0", Estimate);
run;

* create curve of predicted probabilities at each value of x ;
ods graphics on;

proc plm restore=reg.log_rcs_hefi_y10_raw noinfo noclprint;
	effectplot fit(x=HEFI2019_TOTAL_SCORE) / clm ilink;
	ods output FitPlot=reg.logrrcst_hefi_y10_raw0(keep=_xcont1 _predicted 
		rename=(_xcont1=HEFI2019_TOTAL_SCORE _predicted=pr_y) );
run;

ods graphics off;
* Create data for independent variables ;

data _t1t0;
	HEFI2019_TOTAL_SCORE=&t0;
	*
HEFI2019_TOTAL_SCORE = 46.785789336;
	tperc=&t0perc;
	*
tperc = 50;
	output;
	HEFI2019_TOTAL_SCORE=&t1;
	*
HEFI2019_TOTAL_SCORE = 56.442513934;
	tperc=&t1perc;
	*
tperc = 90;
	output;
run;

* Calculate Pr(Y) for t1 and t0 values ;

proc plm restore=reg.log_rcs_hefi_y10_raw noinfo noclprint;
	score data=_t1t0 out=_t1t0 predicted / ilink;
	*ilink option requests Pr(Y);
run;

proc transpose data=_t1t0 out=_t1t0w(drop=_: rename=(COL1=y0 COL2=y1)) 
		prefix=COL;
	var predicted;
run;

* save model information and result of t1-t0 contrast;

data reg.logr_mod_hefi_y10_raw0;
	retain replicate label z t0 t1 log_beta exp_y0 exp_y1 exp_beta y0 y1 risk_diff 
		risk_ratio;
	merge _param(rename=(estimate=log_beta expestimate=exp_beta) ) _t1t0w;
	* save replicate number ;
	replicate=0;
	* save the comparison values;
	t0=&t0;
	* t0 = 46.785789336;
	t1=&t1;
	* t1 = 56.442513934;
	* save p1, p99 plot limit ;
	p1=&p1;
	* p1=25.284244094;
	p99=&p99;
	* p99=62.880843051;
	* save model covariates, if need be;
	z="NA";
	label z="Model did not include covariates";
	* calculate odds(y) at the percentiles;
	exp_y1=y1/(1-y1);
	exp_y0=y0/(1-y0);
	* calculate absolute and relative risk;
	risk_diff=y1 - y0;
	risk_ratio=y1 / y0;
	* input model intercept ;
	log_beta0=&log_beta0;
	*
log_beta0 = 3.3635134759 ;
	* label variables;
	label t1="90th perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		t0="50th perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		y1="Pr(binary_10) at 90th perc. of usual HEFI2019_TOTAL_SCORE" 
		y0="Pr(binary_10) at 50th perc. of usual HEFI2019_TOTAL_SCORE" 
		exp_y1="odds(binary_10) at 90th perc. of usual HEFI2019_TOTAL_SCORE" 
		exp_y0="odds(binary_10) at 50th perc. of usual HEFI2019_TOTAL_SCORE" 
		p1="1st perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		p99="99th perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		log_beta="logit(binary_10) | p90-p50 HEFI2019_TOTAL_SCORE diff." 
		exp_beta="odds(binary_10) | p90-p50 HEFI2019_TOTAL_SCORE diff." 
		risk_diff="risk diff. binary_10 | p90-p50 HEFI2019_TOTAL_SCORE diff." 
		risk_ratio="risk ratio binary_10 | p90-p50 HEFI2019_TOTAL_SCORE diff." 
		log_beta0="Model intercept (ln scale)";
run;

* transpose predicted probability curve from long to wide ;

proc transpose data=reg.logrrcst_hefi_y10_raw0 out=plot0 (rename=(_NAME_=name)) 
		prefix=pred;
	var HEFI2019_TOTAL_SCORE: pr_y;
run;

data reg.logrrcsw_hefi_y10_raw0;
	retain replicate;
	set plot0(drop=_:);
	replicate=0;
	label name=" ";
run;

* note: logrrcsw = LOGistic Regression with Restricted Cubic Spline plot in the Wide format ;
 
 /*************************************************************************/
 /*                                                                       */
 /*  Distribution: marginal Pr(X<x) with percentile svy                   */
 /*                                                                       */
 /*************************************************************************/
	
/* 1) ABSOLUTE assessment: grams of protein (raw) */

	/* List of cutpoints of interest for outcome variable */
	%let cutpoints = 60 80 100 120 140 ;
	
	%percentiles_Survey(data      = baselib._usintake_mc_t_out0 ,
	                    byvar     = ,
	                    var       = proteins ,
	                    weight    = weight_nw_sumw_div,
	                    cutpoints = &cutpoints.,
	                    print     = N,
	                    ntitle    = 0
	                    );
	                    
	data mrg_all_w ;
		set _percentiles;
	run;
	
	data mrg_all_t;
		set _percentiles2;
	run;

	%percentiles_Survey(data      = baselib._usintake_mc_t_out0 ,
	                    byvar     = drig,
	                    var       = proteins ,
	                    weight    = weight_nw_sumw_div,
	                    cutpoints = &cutpoints.,
	                    print     = N,
	                    ntitle    = 0
	                    );

	data reslib.distrib_y_abs_w0 ;
	retain replicate drig varname ;
	length varname $ 32;
		set mrg_all_w _percentiles;
	* label current replicate, outcome and subgroup ;
	replicate=0;
	varname="proteins" ;
	if missing(drig) then drig=0;
	run;
	
	data reslib.distrib_y_abs_t0;
	retain replicate drig varname ;
	length varname $ 32;
		set mrg_all_t _percentiles2;
	* label current replicate, outcome and subgroup ;
	replicate=0;
	varname="&nut_y" ;
	if missing(drig) then drig=0;
	* add numerical percentile count ;
	if index(Statistic,"Pctile")>0 then p = input(compress(Statistic,,'a'),10.) ;
	run;
	
	/* clean temporary data */
	proc datasets lib=work nodetails nolist ;
	delete mrg_all_w mrg_all_t _prob _percentiles _percentiles2 ;
	run;

	
/* 2) RELATIVE assessment: grams of protein per kg */

	/* List of cutpoints of interest for outcome variable */
	%let cutpoints = 0.66 0.8 1.0 1.2 1.4 ;

	%percentiles_Survey(data      = baselib._usintake_mc_t_out0 ,
	                    byvar     = ,
	                    var       = proteins_kg ,
	                    weight    = weight_nw_sumw_div,
	                    cutpoints = &cutpoints. ,
	                    print     = N,
	                    ntitle    = 0
	                    );
	                    
	data mrg_all_w ;
		set _percentiles;
	run;
	
	data mrg_all_t;
		set _percentiles2;
	run;

	%percentiles_Survey(data      = baselib._usintake_mc_t_out0 ,
	                    byvar     = drig,
	                    var       = proteins_kg ,
	                    weight    = weight_nw_sumw_div,
	                    cutpoints = &cutpoints.,
	                    print     = N,
	                    ntitle    = 0
	                    );

	data reslib.distrib_y_rel_w0 ;
	retain replicate drig varname ;
	length varname $ 32;
		set mrg_all_w _percentiles;
	* label current replicate, outcome and subgroup ;
	replicate=0;
	varname="proteins_kg" ;
	if missing(drig) then drig=0;
	run;
	
	data reslib.distrib_y_rel_t0;
	retain replicate drig varname ;
	length varname $ 32;
		set mrg_all_t _percentiles2;
	* label current replicate, outcome and subgroup ;
	replicate=0;
	varname="proteins_kg" ;
	if missing(drig) then drig=0;
	* add numerical percentile count ;
	if index(Statistic,"Pctile")>0 then p = input(compress(Statistic,,'a'),10.) ;
	run;
	
	/* clean temporary files */
	proc datasets lib=work nodetails nolist ;
	delete mrg_all_w mrg_all_t _prob _percentiles _percentiles2 ;
	run;

/* end of code */