 /*************************************************************************/
 /*                                                                       */
 /*                     CCHS 2015 - Nutrition (PUMF)                      */
 /*                                                                       */
 /* Output results: prepare data, linear reg, logistic reg, distribution  */
 /*                 Code for outcome: DFE, Mg, Fibers, Potassium          */
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
	%let suffix =  _miscB; 
	
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

	data baselib._usintake_mc_t_out0(
			drop=mc_t1-mc_t19 i mg_EAR dfe_EAR pota_AI fibers_AI
			rename=(folate=dfe magnesium=mg potassium=pota)
			);
		set baselib.mc_t_distrib_out0;
		 
	* Divide weights by set_number_monte_carlo_rand_obs value used in the MULTIVAR_DISTRIB macro ;
		weight_nw_sumw_div = weight_nw_sumw / 500 ;
		
	* Assign variable names for usual intakes from the MULTIVAR_DISTRIB macro ;
		array outmc (*) mc_t1-mc_t19 ;
		array clean (*) wg pfpb otherbevs milk rg vf otherfoods pfab water mufa pufa 
			sfa freesugars sodium energy folate magnesium fibers potassium  ;
		do i=1 to dim(outmc);
			clean(i)=outmc(i);
		end;
		 
	* Additional formatting (derived variables, binary cutpoints);

	* note: magnesium, fibers, potassium AI or EAR varies according to sex only ;
	
	if (sex=1) then do;
	
		mg_EAR = 350 ;
		fibers_AI = 30 ;
		pota_AI = 3400 ;

	end;
	else if (sex=2) then do ;
	
		mg_EAR = 265 ;
		fibers_AI = 21 ;
		pota_AI = 2600 ;
		
	end;
	
	dfe_EAR = 320 ;
	
	
	* Make outcome variables ;

	if not missing(folate) then do;
		if folate < dfe_EAR then binary_dfe = 1;
			else binary_dfe =0;
	end;
	
	if not missing(magnesium) then do;
		if magnesium < mg_EAR then binary_mg =1;
			else binary_mg = 0;
	end;

	if not missing(fibers) then do;
		if fibers < fibers_AI then binary_fib =1;
			else binary_fib = 0;
	end;
	
	if not missing(potassium) then do;
		if potassium < pota_AI then binary_pota = 1;
			else binary_pota =0;
	end;
	
	label 
		binary_dfe = "Intake < than Diet.Folate.Eq. EAR (320mg)"
		binary_mg = "Intake < than Magnesium EAR (265-350mg)"
		binary_fib = "Intake < than Fibers AI (21-30mg)"
		binary_pota = "Intake < than Potassium AI (2.6-3.4g)"
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
/* Linear regression model (dfe=rcs(HEFI2019_TOTAL_SCORE, 5) )           */
/*************************************************************************/

proc surveyreg data=baselib._usintake_mc_t_out0;
	effect spl_x=spline(HEFI2019_TOTAL_SCORE / details naturalcubic 
		basis=tpf(noint) knotmethod=list(30.69 40.73 46.31 51.51 59.31));
	model dfe=spl_x / adjrsq;
	weight weight_nw_sumw_div;
	estimate "dfe|rcs(HEFI2019_TOTAL_SCORE,5), high(p90) vs. low(p50)" spl_x 
	[1 , 56.699478711][-1 , 46.310059657];
	ods output Estimates=_param(keep=Label Estimate) 
		fitstatistics=_fit(where=(Label1 in ('Adjusted R-Square' 'R carré ajusté')));
	store Reg.rc_rcs_dfe_raw;
run;

ods graphics on;
	
	proc plm restore=Reg.rc_rcs_dfe_raw noinfo noclprint;
		effectplot fit(x=HEFI2019_TOTAL_SCORE) / clm;
		ods output FitPlot=Reg.plmrcst_dfe_raw0(keep=_xcont1 _predicted 
			rename=(_xcont1=HEFI2019_TOTAL_SCORE _predicted=predicted) );
	run;

ods graphics off;

* Create data for independant variables ;
data _t1t0;
	HEFI2019_TOTAL_SCORE=&t0;
		*HEFI2019_TOTAL_SCORE = 46.310059657;
	tperc=&t0perc;
		*tperc = 50;
	output;
	HEFI2019_TOTAL_SCORE=&t1;
		*HEFI2019_TOTAL_SCORE = 56.699478711;
	tperc=&t1perc;
		*tperc = 90;
	output;
run;

* Calculate predicted values (Y) for t1 and t0 (X));
	proc plm restore=Reg.rc_rcs_dfe_raw noinfo noclprint;
		score data=_t1t0 out=_t1t0 predicted;
	run;
	
	proc transpose data=_t1t0 out=_t1t0w(drop=_: rename=(COL1=t0_score 
			COL2=t1_score)) prefix=COL;
		var predicted;
	run;

* Save regression output ;
data _null_;
	set _fit;
	call symputx("r2", nvalue1);
run;

proc transpose data=_param out=paramt(drop=_:) prefix=parm;
	var Estimate;
	copy label;
run;

data Reg.rcs_dfe_raw0;
	retain replicate label z;
	merge paramt(rename=(parm1=rc_beta) ) _t1t0w;
	* save replicate number ;
	replicate=0;
	* Save the comparison values;
	t0=&t0;
		* t0 = 46.310059657;
	t1=&t1;
		*t1 = 56.699478711;
	* save r2 ;
	r2=&r2;
		*r2 = 0.0140768229;
	* save p1, p99 plot limit ;
	p1=&p1;
		*p1=24.335132654;
	p99=&p99;
		*p99=63.622905466;
	* save model covariates, if need be;
	z="NA";
	label z="Model did not include covariates";
	* label variables;
	label t1="90th perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		t0="50th perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		t1_score="dfe value at 90th perc. of usual HEFI2019_TOTAL_SCORE" 
		t0_score="dfe value  at 50th perc. of usual HEFI2019_TOTAL_SCORE" 
		p1="1st perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		p99="99th perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		rc_beta="(dfe change | p90-p50 HEFI2019_TOTAL_SCORE diff.)" 
		r2="Adjusted r-square";
run;

proc transpose data=Reg.plmrcst_dfe_raw0 out=plot0 (rename=(_NAME_=name)) 
		prefix=pred;
	var HEFI2019_TOTAL_SCORE: predicted;
run;

data Reg.plm_rcs_dfe_raw0;
	retain replicate;
	set plot0(drop=_:);
	replicate=0;
	label name=" ";
run;

/*************************************************************************/
/* Linear regression model (mg = rcs(HEFI2019_TOTAL_SCORE,5)  )          */
/*************************************************************************/

proc surveyreg data=baselib._usintake_mc_t_out0;
	effect spl_x=spline(HEFI2019_TOTAL_SCORE / details naturalcubic 
		basis=tpf(noint) knotmethod=list(30.69 40.73 46.31 51.51 59.31));
	model mg=spl_x / adjrsq;
	weight weight_nw_sumw_div;
	estimate "mg|rcs(HEFI2019_TOTAL_SCORE,5), high(p90) vs. low(p50)" spl_x
	[1 , 56.699478711][-1 , 46.310059657];
	ods output Estimates=_param(keep=Label Estimate) 
		fitstatistics=_fit(where=(Label1 in ('Adjusted R-Square' 'R carré ajusté')));
	store Reg.rc_rcs_mg_raw;
run;

ods graphics on;

proc plm restore=Reg.rc_rcs_mg_raw noinfo noclprint;
	effectplot fit(x=HEFI2019_TOTAL_SCORE) / clm;
	ods output FitPlot=Reg.plmrcst_mg_raw0(keep=_xcont1 _predicted 
		rename=(_xcont1=HEFI2019_TOTAL_SCORE _predicted=predicted) );
run;

ods graphics off;

* Create data for independant variables ;
data _t1t0;
	HEFI2019_TOTAL_SCORE=&t0;
	*HEFI2019_TOTAL_SCORE = 46.310059657;
	tperc=&t0perc;
	*tperc = 50;
	output;
	HEFI2019_TOTAL_SCORE=&t1;
	*HEFI2019_TOTAL_SCORE = 56.699478711;
	tperc=&t1perc;
	*tperc = 90;
	output;
run;

* Calculate predicted values (Y) for t1 and t0 (X));
proc plm restore=Reg.rc_rcs_mg_raw noinfo noclprint;
	score data=_t1t0 out=_t1t0 predicted;
run;

proc transpose data=_t1t0 out=_t1t0w(drop=_: rename=(COL1=t0_score 
		COL2=t1_score)) prefix=COL;
	var predicted;
run;

* Save regression output ;
data _null_;
	set _fit;
	call symputx("r2", nvalue1);
run;

proc transpose data=_param out=paramt(drop=_:) prefix=parm;
	var Estimate;
	copy label;
run;

data Reg.rcs_mg_raw0;
	retain replicate label z;
	merge paramt(rename=(parm1=rc_beta) ) _t1t0w;
	* save replicate number ;
	replicate=0;
	* Save the comparison values;
	t0=&t0;
	*t0 = 46.310059657;
	t1=&t1;
	*t1 = 56.699478711;
	* save r2 ;
	r2=&r2;
	*r2 = 0.0922226353;
	* save p1, p99 plot limit ;
	p1=&p1;
	*p1=24.335132654;
	p99=&p99;
	*p99=63.622905466;
	* save model covariates, if need be;
	z="NA";
	label z="Model did not include covariates";
	* label variables;
	label t1="90th perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		t0="50th perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		t1_score="mg value at 90th perc. of usual HEFI2019_TOTAL_SCORE" 
		t0_score="mg value at 50th perc. of usual HEFI2019_TOTAL_SCORE" 
		p1="1st perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		p99="99th perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		rc_beta="(mg change | p90-p50 HEFI2019_TOTAL_SCORE diff.)" 
		r2="Adjusted r-square";
run;

proc transpose data=Reg.plmrcst_mg_raw0 out=plot0 (rename=(_NAME_=name)) 
		prefix=pred;
	var HEFI2019_TOTAL_SCORE: predicted;
run;

data Reg.plm_rcs_mg_raw0;
	retain replicate;
	set plot0(drop=_:);
	replicate=0;
	label name=" ";
run;


/*************************************************************************/
/* Linear regression model (fibers = rcs(HEFI2019_TOTAL_SCORE,5)  )      */
/*************************************************************************/

proc surveyreg data=baselib._usintake_mc_t_out0;
	effect spl_x=spline(HEFI2019_TOTAL_SCORE / details naturalcubic 
		basis=tpf(noint) knotmethod=list(30.69 40.73 46.31 51.51 59.31));
	model fibers=spl_x / adjrsq;
	weight weight_nw_sumw_div;
	estimate "fibers|rcs(HEFI2019_TOTAL_SCORE,5), high(p90) vs. low(p50)" 
		spl_x [1 , 56.699478711][-1 , 46.310059657];
	ods output Estimates=_param(keep=Label Estimate) 
		fitstatistics=_fit(where=(Label1 in ('Adjusted R-Square' 'R carré ajusté')));
	store Reg.rc_rcs_fibers_raw;
run;

ods graphics on;

	proc plm restore=Reg.rc_rcs_fibers_raw noinfo noclprint;
		effectplot fit(x=HEFI2019_TOTAL_SCORE) / clm;
		ods output FitPlot=Reg.plmrcst_fibers_raw0(keep=_xcont1 _predicted 
			rename=(_xcont1=HEFI2019_TOTAL_SCORE _predicted=predicted) );
	run;

ods graphics off;

* Create data for independant variables ;
data _t1t0;
	HEFI2019_TOTAL_SCORE=&t0;
	* HEFI2019_TOTAL_SCORE = 46.310059657;
	tperc=&t0perc;
	* tperc = 50;
	output;
	HEFI2019_TOTAL_SCORE=&t1;
	* HEFI2019_TOTAL_SCORE = 56.699478711;
	tperc=&t1perc;
	* tperc = 90;
	output;
run;

* Calculate predicted values (Y) for t1 and t0 (X));
proc plm restore=Reg.rc_rcs_fibers_raw noinfo noclprint;
	score data=_t1t0 out=_t1t0 predicted;
run;

proc transpose data=_t1t0 out=_t1t0w(drop=_: rename=(COL1=t0_score 
		COL2=t1_score)) prefix=COL;
	var predicted;
run;

* Save regression output ;
data _null_;
	set _fit;
	call symputx("r2", nvalue1);
run;

proc transpose data=_param out=paramt(drop=_:) prefix=parm;
	var Estimate;
	copy label;
run;

data Reg.rcs_fibers_raw0;
	retain replicate label z;
	merge paramt(rename=(parm1=rc_beta) ) _t1t0w;
	* save replicate number ;
	replicate=0;
	* Save the comparison values;
	t0=&t0;
	* t0 = 46.310059657;
	t1=&t1;
	*t1 = 56.699478711;
	* save r2 ;
	r2=&r2;
	*r2 = 0.1987851648;
	* save p1, p99 plot limit ;
	p1=&p1;
	*p1=24.335132654;
	p99=&p99;
	*p99=63.622905466;
	* save model covariates, if need be;
	z="NA";
	label z="Model did not include covariates";
	* label variables;
	label t1="90th perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		t0="50th perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		t1_score="fibers value at 90th perc. of usual HEFI2019_TOTAL_SCORE" 
		t0_score="fibers value at 50th perc. of usual HEFI2019_TOTAL_SCORE" 
		p1="1st perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		p99="99th perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		rc_beta="(fibers change | p90-p50 HEFI2019_TOTAL_SCORE diff.)" 
		r2="Adjusted r-square";
run;

proc transpose data=Reg.plmrcst_fibers_raw0 out=plot0 (rename=(_NAME_=name)) 
		prefix=pred;
	var HEFI2019_TOTAL_SCORE: predicted;
run;

data Reg.plm_rcs_fibers_raw0;
	retain replicate;
	set plot0(drop=_:);
	replicate=0;
	label name=" ";
run;

/*************************************************************************/
/* Linear regression model (pota = rcs(HEFI2019_TOTAL_SCORE,5)  )        */
/*************************************************************************/

proc surveyreg data=baselib._usintake_mc_t_out0;
	effect spl_x=spline(HEFI2019_TOTAL_SCORE / details naturalcubic 
		basis=tpf(noint) knotmethod=list(30.69 40.73 46.31 51.51 59.31));
	model pota=spl_x / adjrsq;
	weight weight_nw_sumw_div;
	estimate "pota|rcs(HEFI2019_TOTAL_SCORE,5), high(p90) vs. low(p50)" spl_x [1 , 
		56.699478711][-1 , 46.310059657];
	ods output Estimates=_param(keep=Label Estimate) 
		fitstatistics=_fit(where=(Label1 in ('Adjusted R-Square' 'R carré ajusté')));
	store Reg.rc_rcs_pota_raw;
run;


ods graphics on;

	proc plm restore=Reg.rc_rcs_pota_raw noinfo noclprint;
		effectplot fit(x=HEFI2019_TOTAL_SCORE) / clm;
		ods output FitPlot=Reg.plmrcst_pota_raw0(keep=_xcont1 _predicted 
			rename=(_xcont1=HEFI2019_TOTAL_SCORE _predicted=predicted) );
	run;

ods graphics off;

* Create data for independant variables ;
data _t1t0;
	HEFI2019_TOTAL_SCORE=&t0;
	*HEFI2019_TOTAL_SCORE = 46.310059657;
	tperc=&t0perc;
	*tperc = 50;
	output;
	HEFI2019_TOTAL_SCORE=&t1;
	*HEFI2019_TOTAL_SCORE = 56.699478711;
	tperc=&t1perc;
	*tperc = 90;
	output;
run;

* Calculate predicted values (Y) for t1 and t0 (X));
proc plm restore=Reg.rc_rcs_pota_raw noinfo noclprint;
	score data=_t1t0 out=_t1t0 predicted;
run;

proc transpose data=_t1t0 out=_t1t0w(drop=_: rename=(COL1=t0_score 
		COL2=t1_score)) prefix=COL;
	var predicted;
run;

* Save regression output ;
data _null_;
	set _fit;
	call symputx("r2", nvalue1);
run;

proc transpose data=_param out=paramt(drop=_:) prefix=parm;
	var Estimate;
	copy label;
run;

data Reg.rcs_pota_raw0;
	retain replicate label z;
	merge paramt(rename=(parm1=rc_beta) ) _t1t0w;
	* save replicate number ;
	replicate=0;
	* Save the comparison values;
	t0=&t0;
	*t0 = 46.310059657;
	t1=&t1;
	*t1 = 56.699478711;
	* save r2 ;
	r2=&r2;
	*r2 = 0.0266528219;
	* save p1, p99 plot limit ;
	p1=&p1;
	*p1=24.335132654;
	p99=&p99;
	*p99=63.622905466;
	* save model covariates, if need be;
	z="NA";
	label z="Model did not include covariates";
	* label variables;
	label t1="90th perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		t0="50th perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		t1_score="pota value at 90th perc. of usual HEFI2019_TOTAL_SCORE" 
		t0_score="pota value at 50th perc. of usual HEFI2019_TOTAL_SCORE" 
		p1="1st perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		p99="99th perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		rc_beta="(pota change | p90-p50 HEFI2019_TOTAL_SCORE diff.)" 
		r2="Adjusted r-square";
run;

proc transpose data=Reg.plmrcst_pota_raw0 out=plot0 (rename=(_NAME_=name)) 
		prefix=pred;
	var HEFI2019_TOTAL_SCORE: predicted;
run;

data Reg.plm_rcs_pota_raw0;
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
/* Logistic reg. model (logit(binary_dfe) = rcs(HEFI2019_TOTAL_SCORE,5) )*/
/*************************************************************************/

proc logistic data=baselib._usintake_mc_t_out0 descending;
	effect spl_x=spline(HEFI2019_TOTAL_SCORE / details naturalcubic 
		basis=tpf(noint) knotmethod=list(30.69 40.73 46.31 51.51 59.31));
	model binary_dfe=spl_x;
	weight weight_nw_sumw_div;
	estimate "binary_dfe|rcs(HEFI2019_TOTAL_SCORE,5), high(p90) vs. low(p50)" 
		spl_x [1 , 56.699478711][-1 , 46.310059657] / exp;
	ods output Estimates=_param(keep=Label Estimate ExpEstimate) 
		ParameterEstimates=_modelparms(Keep=Variable Estimate);
	store reg.log_rcs_y_dfe_raw;
run;

* output model intercept ;
data _NULL_;
	set _modelparms;

	if _N_=1 then
		call symputx("log_beta0", Estimate);
run;

* create curve of predicted probabilities at each value of x ;
ods graphics on;
	
	proc plm restore=reg.log_rcs_y_dfe_raw noinfo noclprint;
		effectplot fit(x=HEFI2019_TOTAL_SCORE) / clm ilink;
		ods output FitPlot=reg.logrrcst_y_dfe_raw0(keep=_xcont1 _predicted 
			rename=(_xcont1=HEFI2019_TOTAL_SCORE _predicted=pr_y) );
	run;

ods graphics off;

* Create data for independent variables ;
data _t1t0;
	HEFI2019_TOTAL_SCORE=&t0;
	*HEFI2019_TOTAL_SCORE = 46.310059657;
	tperc=&t0perc;
	*tperc = 50;
	output;
	HEFI2019_TOTAL_SCORE=&t1;
	*HEFI2019_TOTAL_SCORE = 56.699478711;
	tperc=&t1perc;
	*tperc = 90;
	output;
run;

* Calculate Pr(Y) for t1 and t0 values ;
proc plm restore=reg.log_rcs_y_dfe_raw noinfo noclprint;
	score data=_t1t0 out=_t1t0 predicted / ilink;
	*ilink option requests Pr(Y);
run;

proc transpose data=_t1t0 out=_t1t0w(drop=_: rename=(COL1=y0 COL2=y1)) 
		prefix=COL;
	var predicted;
run;

* Save regression output and result of t1-t0 contrast;
data reg.logr_mod_y_dfe_raw0;
	retain replicate label z t0 t1 log_beta exp_y0 exp_y1 exp_beta y0 y1 risk_diff 
		risk_ratio;
	merge _param(rename=(estimate=log_beta expestimate=exp_beta) ) _t1t0w;
	* save replicate number ;
	replicate=0;
	* save the comparison values;
	t0=&t0;
	*t0 = 46.310059657;
	t1=&t1;
	*t1 = 56.699478711;
	* save p1, p99 plot limit ;
	p1=&p1;
	*p1=24.335132654;
	p99=&p99;
	*p99=63.622905466;
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
	* log_beta0 = -1.73773143614928 ;
	* label variables;
	label t1="90th perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		t0="50th perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		y1="Pr(binary_dfe) at 90th perc. of usual HEFI2019_TOTAL_SCORE" 
		y0="Pr(binary_dfe) at 50th perc. of usual HEFI2019_TOTAL_SCORE" 
		exp_y1="odds(binary_dfe) at 90th perc. of usual HEFI2019_TOTAL_SCORE" 
		exp_y0="odds(binary_dfe) at 50th perc. of usual HEFI2019_TOTAL_SCORE" 
		p1="1st perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		p99="99th perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		log_beta="logit(binary_dfe) | p90-p50 HEFI2019_TOTAL_SCORE diff." 
		exp_beta="odds(binary_dfe) | p90-p50 HEFI2019_TOTAL_SCORE diff." 
		risk_diff="risk diff. binary_dfe | p90-p50 HEFI2019_TOTAL_SCORE diff." 
		risk_ratio="risk ratio binary_dfe | p90-p50 HEFI2019_TOTAL_SCORE diff." 
		log_beta0="Model intercept (ln scale)";
run;

* transpose predicted probability curve from long to wide ;
proc transpose data=reg.logrrcst_y_dfe_raw0 out=plot0 (rename=(_NAME_=name)) 
		prefix=pred;
	var HEFI2019_TOTAL_SCORE: pr_y;
run;

data reg.logrrcsw_y_dfe_raw0;
	retain replicate;
	set plot0(drop=_:);
	replicate=0;
	label name=" ";
run;

* note: logrrcsw = LOGistic Regression with Restricted Cubic Spline plot in the Wide format ;

/*************************************************************************/
/* Logistic reg. model (logit(binary_mg) = rcs(HEFI2019_TOTAL_SCORE,5) ) */
/*************************************************************************/

proc logistic data=baselib._usintake_mc_t_out0 descending;
	effect spl_x=spline(HEFI2019_TOTAL_SCORE / details naturalcubic 
		basis=tpf(noint) knotmethod=list(30.69 40.73 46.31 51.51 59.31));
	model binary_mg=spl_x;
	weight weight_nw_sumw_div;
	estimate "binary_mg|rcs(HEFI2019_TOTAL_SCORE,5), high(p90) vs. low(p50)" 
		spl_x [1 , 56.699478711][-1 , 46.310059657] / exp;
	ods output Estimates=_param(keep=Label Estimate ExpEstimate) 
		ParameterEstimates=_modelparms(Keep=Variable Estimate);
	store reg.log_rcs_y_mg_raw;
run;

* output model intercept ;
data _NULL_;
	set _modelparms;

	if _N_=1 then
		call symputx("log_beta0", Estimate);
run;

* create curve of predicted probabilities at each value of x ;
ods graphics on;

	proc plm restore=reg.log_rcs_y_mg_raw noinfo noclprint;
		effectplot fit(x=HEFI2019_TOTAL_SCORE) / clm ilink;
		ods output FitPlot=reg.logrrcst_y_mg_raw0(keep=_xcont1 _predicted 
			rename=(_xcont1=HEFI2019_TOTAL_SCORE _predicted=pr_y) );
	run;

ods graphics off;

* Create data for independent variables ;
data _t1t0;
	HEFI2019_TOTAL_SCORE=&t0;
	*HEFI2019_TOTAL_SCORE = 46.310059657;
	tperc=&t0perc;
	*tperc = 50;
	output;
	HEFI2019_TOTAL_SCORE=&t1;
	*HEFI2019_TOTAL_SCORE = 56.699478711;
	tperc=&t1perc;
	*tperc = 90;
	output;
run;

* Calculate Pr(Y) for t1 and t0 values ;

proc plm restore=reg.log_rcs_y_mg_raw noinfo noclprint;
	score data=_t1t0 out=_t1t0 predicted / ilink;
	*ilink option requests Pr(Y);
run;

proc transpose data=_t1t0 out=_t1t0w(drop=_: rename=(COL1=y0 COL2=y1)) 
		prefix=COL;
	var predicted;
run;

* Save regression output and result of t1-t0 contrast;

data reg.logr_mod_y_mg_raw0;
	retain replicate label z t0 t1 log_beta exp_y0 exp_y1 exp_beta y0 y1 risk_diff 
		risk_ratio;
	merge _param(rename=(estimate=log_beta expestimate=exp_beta) ) _t1t0w;
	* save replicate number ;
	replicate=0;
	* save the comparison values;
	t0=&t0;
	*t0 = 46.310059657;
	t1=&t1;
	*t1 = 56.699478711;
	* save p1, p99 plot limit ;
	p1=&p1;
	*
p1=24.335132654;
	p99=&p99;
	*p99=63.622905466;
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
	*log_beta0 = 4.4026930772 ;
	* label variables;
	label t1="90th perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		t0="50th perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		y1="Pr(binary_mg) at 90th perc. of usual HEFI2019_TOTAL_SCORE" 
		y0="Pr(binary_mg) at 50th perc. of usual HEFI2019_TOTAL_SCORE" 
		exp_y1="odds(binary_mg) at 90th perc. of usual HEFI2019_TOTAL_SCORE" 
		exp_y0="odds(binary_mg) at 50th perc. of usual HEFI2019_TOTAL_SCORE" 
		p1="1st perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		p99="99th perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		log_beta="logit(binary_mg) | p90-p50 HEFI2019_TOTAL_SCORE diff."
		exp_beta="odds(binary_mg) | p90-p50 HEFI2019_TOTAL_SCORE diff." 
		risk_diff="risk diff. binary_mg | p90-p50 HEFI2019_TOTAL_SCORE diff." 
		risk_ratio="risk ratio binary_mg | p90-p50 HEFI2019_TOTAL_SCORE diff." 
		log_beta0="Model intercept (ln scale)";
run;

* transpose predicted probability curve from long to wide ;
proc transpose data=reg.logrrcst_y_mg_raw0 out=plot0 (rename=(_NAME_=name)) 
		prefix=pred;
	var HEFI2019_TOTAL_SCORE: pr_y;
run;

data reg.logrrcsw_y_mg_raw0;
	retain replicate;
	set plot0(drop=_:);
	replicate=0;
	label name=" ";
run;

* note: logrrcsw = LOGistic Regression with Restricted Cubic Spline plot in the Wide format ;

/*************************************************************************/
/* Logistic reg. model (logit(binary_fib) = rcs(HEFI2019_TOTAL_SCORE,5) )*/
/*************************************************************************/

proc logistic data=baselib._usintake_mc_t_out0 descending;
	effect spl_x=spline(HEFI2019_TOTAL_SCORE / details naturalcubic 
		basis=tpf(noint) knotmethod=list(30.69 40.73 46.31 51.51 59.31));
	model binary_fib=spl_x;
	weight weight_nw_sumw_div;
	estimate "binary_fib|rcs(HEFI2019_TOTAL_SCORE,5), high(p90) vs. low(p50)" 
		spl_x [1 , 56.699478711][-1 , 46.310059657] / exp;
	ods output Estimates=_param(keep=Label Estimate ExpEstimate) 
		ParameterEstimates=_modelparms(Keep=Variable Estimate);
	store reg.log_rcs_y_fib_raw;
run;

* output model intercept ;
data _NULL_;
	set _modelparms;

	if _N_=1 then
		call symputx("log_beta0", Estimate);
run;

* create curve of predicted probabilities at each value of x ;
ods graphics on;

	proc plm restore=reg.log_rcs_y_fib_raw noinfo noclprint;
		effectplot fit(x=HEFI2019_TOTAL_SCORE) / clm ilink;
		ods output FitPlot=reg.logrrcst_y_fib_raw0(keep=_xcont1 _predicted 
			rename=(_xcont1=HEFI2019_TOTAL_SCORE _predicted=pr_y) );
	run;

ods graphics off;

* Create data for independent variables ;
data _t1t0;
	HEFI2019_TOTAL_SCORE=&t0;
	* HEFI2019_TOTAL_SCORE = 46.310059657;
	tperc=&t0perc;
	* tperc = 50;
	output;
	HEFI2019_TOTAL_SCORE=&t1;
	* HEFI2019_TOTAL_SCORE = 56.699478711;
	tperc=&t1perc;
	*tperc = 90;
	output;
run;

* Calculate Pr(Y) for t1 and t0 values ;
proc plm restore=reg.log_rcs_y_fib_raw noinfo noclprint;
	score data=_t1t0 out=_t1t0 predicted / ilink;
	*ilink option requests Pr(Y);
run;

proc transpose data=_t1t0 out=_t1t0w(drop=_: rename=(COL1=y0 COL2=y1)) 
		prefix=COL;
	var predicted;
run;

* Save regression output and result of t1-t0 contrast;

data reg.logr_mod_y_fib_raw0;
	retain replicate label z t0 t1 log_beta exp_y0 exp_y1 exp_beta y0 y1 risk_diff 
		risk_ratio;
	merge _param(rename=(estimate=log_beta expestimate=exp_beta) ) _t1t0w;
	* save replicate number ;
	replicate=0;
	* save the comparison values;
	t0=&t0;
	* t0 = 46.310059657;
	t1=&t1;
	*t1 = 56.699478711;
	* save p1, p99 plot limit ;
	p1=&p1;
	*p1=24.335132654;
	p99=&p99;
	*p99=63.622905466;
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
log_beta0 = 9.7529845345 ;
	* label variables;
	label t1="90th perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		t0="50th perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		y1="Pr(binary_fib) at 90th perc. of usual HEFI2019_TOTAL_SCORE" 
		y0="Pr(binary_fib) at 50th perc. of usual HEFI2019_TOTAL_SCORE" 
		exp_y1="odds(binary_fib) at 90th perc. of usual HEFI2019_TOTAL_SCORE" 
		exp_y0="odds(binary_fib) at 50th perc. of usual HEFI2019_TOTAL_SCORE" 
		p1="1st perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		p99="99th perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		log_beta="logit(binary_fib) | p90-p50 HEFI2019_TOTAL_SCORE diff."
		exp_beta="odds(binary_fib) | p90-p50 HEFI2019_TOTAL_SCORE diff." 
		risk_diff="risk diff. binary_fib | p90-p50 HEFI2019_TOTAL_SCORE diff." 
		risk_ratio="risk ratio binary_fib | p90-p50 HEFI2019_TOTAL_SCORE diff." 
		log_beta0="Model intercept (ln scale)";
run;

* transpose predicted probability curve from long to wide ;
proc transpose data=reg.logrrcst_y_fib_raw0 out=plot0 (rename=(_NAME_=name)) 
		prefix=pred;
	var HEFI2019_TOTAL_SCORE: pr_y;
run;

data reg.logrrcsw_y_fib_raw0;
	retain replicate;
	set plot0(drop=_:);
	replicate=0;
	label name=" ";
run;

* note: logrrcsw = LOGistic Regression with Restricted Cubic Spline plot in the Wide format ;

/*************************************************************************/
/* Logistic reg. model (logit(binary_pota) = rcs(HEFI2019_TOTAL_SCORE,5))*/
/*************************************************************************/

proc logistic data=baselib._usintake_mc_t_out0 descending;
	effect spl_x=spline(HEFI2019_TOTAL_SCORE / details naturalcubic 
		basis=tpf(noint) knotmethod=list(30.69 40.73 46.31 51.51 59.31));
	model binary_pota=spl_x;
	weight weight_nw_sumw_div;
	estimate "binary_pota|rcs(HEFI2019_TOTAL_SCORE,5), high(p90) vs. low(p50)" 
		spl_x [1 , 56.699478711][-1 , 46.310059657] / exp;
	ods output Estimates=_param(keep=Label Estimate ExpEstimate) 
		ParameterEstimates=_modelparms(Keep=Variable Estimate);
	store reg.log_rcs_y_pota_raw;
run;

* output model intercept ;

data _NULL_;
	set _modelparms;

	if _N_=1 then
		call symputx("log_beta0", Estimate);
run;

* create curve of predicted probabilities at each value of x ;

ods graphics on;
	
	proc plm restore=reg.log_rcs_y_pota_raw noinfo noclprint;
		effectplot fit(x=HEFI2019_TOTAL_SCORE) / clm ilink;
		ods output FitPlot=reg.logrrcst_y_pota_raw0(keep=_xcont1 _predicted 
			rename=(_xcont1=HEFI2019_TOTAL_SCORE _predicted=pr_y) );
	run;

ods graphics off;
* Create data for independent variables ;

data _t1t0;
	HEFI2019_TOTAL_SCORE=&t0;
	*HEFI2019_TOTAL_SCORE = 46.310059657;
	tperc=&t0perc;
	*tperc = 50;
	output;
	HEFI2019_TOTAL_SCORE=&t1;
	*HEFI2019_TOTAL_SCORE = 56.699478711;
	tperc=&t1perc;
	*tperc = 90;
	output;
run;

* Calculate Pr(Y) for t1 and t0 values ;
proc plm restore=reg.log_rcs_y_pota_raw noinfo noclprint;
	score data=_t1t0 out=_t1t0 predicted / ilink;
	*ilink option requests Pr(Y);
run;

proc transpose data=_t1t0 out=_t1t0w(drop=_: rename=(COL1=y0 COL2=y1)) 
		prefix=COL;
	var predicted;
run;

* Save regression output and result of t1-t0 contrast;
data reg.logr_mod_y_pota_raw0;
	retain replicate label z t0 t1 log_beta exp_y0 exp_y1 exp_beta y0 y1 risk_diff 
		risk_ratio;
	merge _param(rename=(estimate=log_beta expestimate=exp_beta) ) _t1t0w;
	* save replicate number ;
	replicate=0;
	* save the comparison values;
	t0=&t0;
	*t0 = 46.310059657;
	t1=&t1;
	*t1 = 56.699478711;
	* save p1, p99 plot limit ;
	p1=&p1;
	*p1=24.335132654;
	p99=&p99;
	*p99=63.622905466;
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
log_beta0 = 3.1541152726 ;
	* label variables;
	label t1="90th perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		t0="50th perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		y1="Pr(binary_pota) at 90th perc. of usual HEFI2019_TOTAL_SCORE" 
		y0="Pr(binary_pota) at 50th perc. of usual HEFI2019_TOTAL_SCORE" 
		exp_y1="odds(binary_pota) at 90th perc. of usual HEFI2019_TOTAL_SCORE" 
		exp_y0="odds(binary_pota) at 50th perc. of usual HEFI2019_TOTAL_SCORE" 
		p1="1st perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		p99="99th perc. of usual HEFI2019_TOTAL_SCORE intake distribution" 
		log_beta="logit(binary_pota) | p90-p50 HEFI2019_TOTAL_SCORE diff." 
		exp_beta="odds(binary_pota) | p90-p50 HEFI2019_TOTAL_SCORE diff." 
		risk_diff="risk diff. binary_pota | p90-p50 HEFI2019_TOTAL_SCORE diff." 
		risk_ratio="risk ratio binary_pota | p90-p50 HEFI2019_TOTAL_SCORE diff." 
		log_beta0="Model intercept (ln scale)";
run;

* transpose predicted probability curve from long to wide ;
proc transpose data=reg.logrrcst_y_pota_raw0 out=plot0 (rename=(_NAME_=name)) 
		prefix=pred;
	var HEFI2019_TOTAL_SCORE: pr_y;
run;

data reg.logrrcsw_y_pota_raw0;
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
 
/* Pr(X < EAR) */
	proc means data = baselib._usintake_mc_t_out0 noprint;
	class drig; 
	var  binary_dfe binary_mg binary_fib binary_pota ;
	weight weight_nw_sumw_div ;
	output out=cutprob_manual mean= ; 
	run;
	
	/* some formatting */
	data cutprob_manual;
		set cutprob_manual;
	samplesize = _FREQ_/500;
	run;
	
/* 1) Folate (DFE) */
	%percentiles_Survey(data      = baselib._usintake_mc_t_out0 ,
	                    byvar     = ,
	                    var       = dfe ,
	                    weight    = weight_nw_sumw_div,
	                    cutpoints = ,
	                    print     = N,
	                    ntitle    = 0
	                    );
	
	data mrg_all_w ;
		merge _percentiles cutprob_manual(keep= _TYPE_ samplesize binary_dfe where=(_TYPE_=0));
		/* no <by> since both data have 1 row */
	run;
	
	data mrg_all_t;
		set _percentiles2;
	run;
	
	%percentiles_Survey(data      = baselib._usintake_mc_t_out0 ,
	                    byvar     = drig,
	                    var       = dfe ,
	                    weight    = weight_nw_sumw_div,
	                    cutpoints = ,
	                    print     = N,
	                    ntitle    = 0
	                    );

/* merge cutprob by drig */
	data _percentiles;
		merge _percentiles cutprob_manual(keep=_TYPE_ samplesize binary_dfe drig where=(_TYPE_=1)) ;
		by drig;
	run;
		
/* Combine overall estimates with drig-specific estimates */
	data reslib.distrib_y_dfe_w0 ;
	retain replicate drig varname ;
	length varname $ 32;
		set mrg_all_w _percentiles;
	* label current replicate, outcome and subgroup ;
	replicate=0;
	varname="dfe" ;
	if missing(drig) then drig=0;
	drop _TYPE_ ;
	run;
	
	data reslib.distrib_y_dfe_t0;
	retain replicate drig varname ;
	length varname $ 32;
		set mrg_all_t _percentiles2;
	* label current replicate, outcome and subgroup ;
	replicate=0;
	varname="dfe" ;
	if missing(drig) then drig=0;
	* add numerical percentile count ;
	if index(Statistic,"Pctile")>0 then p = input(compress(Statistic,,'a'),10.) ;
	run;
	
	proc datasets lib=work nodetails nolist ;
	delete mrg_all_w mrg_all_t _prob _percentiles _percentiles2 ;
	run;

/* 2) Magnesium */
	%percentiles_Survey(data      = baselib._usintake_mc_t_out0 ,
	                    byvar     = ,
	                    var       = mg ,
	                    weight    = weight_nw_sumw_div,
	                    cutpoints = ,
	                    print     = N,
	                    ntitle    = 0
	                    );
	
	data mrg_all_w ;
		merge _percentiles cutprob_manual(keep= _TYPE_ samplesize binary_mg where=(_TYPE_=0));
		/* no <by> since both data have 1 row */
	run;
	
	data mrg_all_t;
		set _percentiles2;
	run;

	%percentiles_Survey(data      = baselib._usintake_mc_t_out0 ,
	                    byvar     = drig,
	                    var       = mg ,
	                    weight    = weight_nw_sumw_div,
	                    cutpoints = ,
	                    print     = N,
	                    ntitle    = 0
	                    );

/* merge cutprob by drig */
	data _percentiles;
		merge _percentiles cutprob_manual(keep=_TYPE_ samplesize binary_mg drig where=(_TYPE_=1)) ;
		by drig;
	run;
		
/* Combine overall estimates with drig-specific estimates */
	data  reslib.distrib_y_mg_w0 ;
	retain replicate drig varname ;
	length varname $ 32;
		set mrg_all_w _percentiles;
	* label current replicate, outcome and subgroup ;
	replicate=0;
	varname="mg" ;
	if missing(drig) then drig=0;
	drop _TYPE_ ;
	run;
	
	data reslib.distrib_y_mg_t0;
	retain replicate drig varname ;
	length varname $ 32;
		set mrg_all_t _percentiles2;
	* label current replicate, outcome and subgroup ;
	replicate=0;
	varname="mg" ;
	if missing(drig) then drig=0;
	* add numerical percentile count ;
	if index(Statistic,"Pctile")>0 then p = input(compress(Statistic,,'a'),10.) ;
	run;
	
	proc datasets lib=work nodetails nolist ;
	delete mrg_all_w mrg_all_t _prob _percentiles _percentiles2 ;
	run;

/* 3) Fibers */
	%percentiles_Survey(data      = baselib._usintake_mc_t_out0 ,
	                    byvar     = ,
	                    var       = fibers ,
	                    weight    = weight_nw_sumw_div,
	                    cutpoints = ,
	                    print     = N,
	                    ntitle    = 0
	                    );
	
	data mrg_all_w ;
		merge _percentiles cutprob_manual(keep= _TYPE_ samplesize binary_fib where=(_TYPE_=0));
		/* no <by> since both data have 1 row */
	run;
	
	data mrg_all_t;
		set _percentiles2;
	run;
	
	%percentiles_Survey(data      = baselib._usintake_mc_t_out0 ,
	                    byvar     = drig,
	                    var       = fibers ,
	                    weight    = weight_nw_sumw_div,
	                    cutpoints = ,
	                    print     = N,
	                    ntitle    = 0
	                    );

/* merge cutprob by drig */
	data _percentiles;
		merge _percentiles cutprob_manual(keep=_TYPE_ samplesize binary_fib drig where=(_TYPE_=1)) ;
		by drig;
	run;
		
/* Combine overall estimates with drig-specific estimates */
	data  reslib.distrib_y_fib_w0 ;
	retain replicate drig varname ;
	length varname $ 32;
		set mrg_all_w _percentiles;
	* label current replicate, outcome and subgroup ;
	replicate=0;
	varname="fibers" ;
	if missing(drig) then drig=0;
	drop _TYPE_ ;
	run;
	
	data reslib.distrib_y_fib_t0;
	retain replicate drig varname ;
	length varname $ 32;
		set mrg_all_t _percentiles2;
	* label current replicate, outcome and subgroup ;
	replicate=0;
	varname="fibers" ;
	if missing(drig) then drig=0;
	* add numerical percentile count ;
	if index(Statistic,"Pctile")>0 then p = input(compress(Statistic,,'a'),10.) ;
	run;

	proc datasets lib=work nodetails nolist ;
	delete mrg_all_w mrg_all_t _prob _percentiles _percentiles2 ;
	run;

/* 4) Potassium */
	%percentiles_Survey(data      = baselib._usintake_mc_t_out0 ,
	                    byvar     = ,
	                    var       = pota ,
	                    weight    = weight_nw_sumw_div,
	                    cutpoints = ,
	                    print     = N,
	                    ntitle    = 0
	                    );
	
	data mrg_all_w ;
		merge _percentiles cutprob_manual(keep= _TYPE_ samplesize binary_pota where=(_TYPE_=0));
		/* no <by> since both data have 1 row */
	run;
	
	
	data mrg_all_t;
		set _percentiles2;
	run;

	%percentiles_Survey(data      = baselib._usintake_mc_t_out0 ,
	                    byvar     = drig,
	                    var       = pota ,
	                    weight    = weight_nw_sumw_div,
	                    cutpoints = ,
	                    print     = N,
	                    ntitle    = 0
	                    );

/* merge cutprob by drig */
	data _percentiles;
		merge _percentiles cutprob_manual(keep=_TYPE_ samplesize binary_pota drig where=(_TYPE_=1)) ;
		by drig;
	run;
		
/* Combine overall estimates with drig-specific estimates */
	data  reslib.distrib_y_pota_w0 ;
	retain replicate drig varname ;
	length varname $ 32;
		set mrg_all_w _percentiles;
	* label current replicate, outcome and subgroup ;
	replicate=0;
	varname="pota" ;
	if missing(drig) then drig=0;
	drop _TYPE_ ;
	run;
	
	data reslib.distrib_y_pota_t0;
	retain replicate drig varname ;
	length varname $ 32;
		set mrg_all_t _percentiles2;
	* label current replicate, outcome and subgroup ;
	replicate=0;
	varname="pota" ;
	if missing(drig) then drig=0;
	* add numerical percentile count ;
	if index(Statistic,"Pctile")>0 then p = input(compress(Statistic,,'a'),10.) ;
	run;

	proc datasets lib=work nodetails nolist ;
	delete mrg_all_w mrg_all_t _prob _percentiles _percentiles2 cutprob_manual;
	run;
 
/* end of code */
