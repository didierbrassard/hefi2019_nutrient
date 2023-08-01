 /*************************************************************************/
 /*                                                                       */
 /*                     CCHS 2015 - Nutrition (PUMF)                      */
 /*                                                                       */
 /*                     Bootstrap variance estimation                     */
 /*                                                                       */
 /*                        Author: Didier Brassard                        */
 /*                                                                       */
 /*                               Version 2                               */
 /*                                MAY2023                                */
 /*                                                                       */
 /*************************************************************************/
 /*                                                                       */
 /* 1 - (CODE 2) : Variance estimation for PROTEIN                        */
 /* 2 - (CODE 3) : Variance estimation for CA/VITD                        */
 /* 3 - (CODE 4) : Variance estimation for MISC. A                        */
 /* 4 - (CODE 5) : Variance estimation for MISC. B                        */
 /* 5 - (CODE 6) : Variance estimation for VIT. A                         */
 /*                                                                       */
 /*************************************************************************/


 /*************************************************************************/
 /*                    SET MACRO VARIABLES AND LIBRARIES                  */
 /*************************************************************************/

/* indicate file location - warning: case sensitive */
%global path suffix;
	%let path = /home/DIBRA22/hefi2019_nutrients/ ;
	
/* include homemade macros */
	/* boot aux.  */ %include "&path./Macros/boot_auxiliary.sas";

/* auto assign proper libraries */
	%macro mcmclib(suffix);
	
	libname MCMC "&path./NCI/MCMC&suffix./";		  		/* create MCMC folder, if need be */
	libname chain "&path./NCI/MCMC&suffix./Markovchains/";	/* folder for trace plots */
	libname reslib "&path./NCI/MCMC&suffix./Results/";		/* results */
	libname baselib "&path./NCI/MCMC&suffix./Model/";		/* base estimates */
	libname bootlib "&path./NCI/MCMC&suffix./Bootlib/";		/* bootstrap replicates */
	libname reg "&path./NCI/MCMC&suffix./Reg/";			  	 		/* base estimate for linear reg */
	libname regb "&path./NCI/MCMC&suffix./Reg/Bootstraps/";  		/* bootstrap data for linear reg */

	%mend mcmclib;
	
	%macro MakeMCMCFolder(NCIPath=,prefix=,suffix=);
	
	 /************************************************/
	 /*                                              */
	 /*             MakeMCMCFolder macro             */
	 /*                                              */
	 /*        Create subfolders to save data        */
	 /*                                              */
	 /************************************************/
	 /*                                              */
	 /* NCIPath = Location of NCI files              */
	 /* prefix  = Name of folder prefix (eg, reg)    */
	 /* suffix  = 1 or more folder suffix (eg, age)  */
	 /*                                              */
	 /************************************************/
	
		%put ;
		%put - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -;
		%put &SYSMACRONAME: Create subfolder(s) in: &NCIPath. ;
		%put - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -;
		%put ;
	
	 %if (%sysevalf(%superq(NCIPath)=,boolean)=1) %then %do;
	 %put ERROR: &sysmacroname. <NCIPath> is not defined. Please indicate value.;
	 %return;
	 %end;
	 %if (%sysevalf(%superq(prefix)=,boolean)=1) %then %do;
	 %local prefix;
	 %let prefix = _;
	 %end;
	 %if (%sysevalf(%superq(suffix)=,boolean)=1) %then %do;
	 %put ERROR: &sysmacroname. <suffix> is not defined. Please indicate value.;
	 %return;
	 %end;
	
	/*%local nSuffix ;
	%let nSuffix = %sysfunc(countw(&suffix,' '));*/
	
	/* initialize numbering */
		%local jth; %let jth=1;
	
	/* loop through suffix, ie create as many outpub library as there are suffix */
		options mprint;
		%do %while (%scan(&Suffix, &jth) ne);
			%let CurrentSuffix = %scan(&Suffix, &jth);
	
		/* create subgroup-specific folder */
		options dlcreatedir;
		libname &prefix.&jth. "&NCIPath./&prefix.&CurrentSuffix./" ;
		libname &prefix.B&jth. "&NCIPath./&prefix.&CurrentSuffix./Bootstraps/" ;
	
			%let jth = %eval(&jth + 1);
		%end; /* end of suffix loop */
		options nomprint;
		
		/* reset value, just in case */
		%let jth=;
	
	%mend MakeMCMCFolder;

/* macro variable common to all analyses */
	%let subgrp = sex drig ;

/* Percentile values to look at when distribution are examined */

%let pctile_list = Pctile1 Pctile5 Pctile10 Pctile25 Pctile50 Pctile75 Pctile90 Pctile95 Pctile99 ;
	
 /*************************************************************************/
 /*                                                                       */
 /*         REF. CODE 2: Variance estimation for PROTEIN outcome          */
 /*                                                                       */
 /*************************************************************************/

 /************************************************/
 /*            Files and folder set-up           */
 /************************************************/

/* 0) define suffix of current outcome */
	
	%let suffix = _pro ;

/* 1) assign proper library using the same suffix */
	%mcmclib(&suffix.);
	%MakeMCMCFolder(NCIPath = &path./NCI/MCMC&suffix.,
					prefix  = Reg_,
					suffix  = &subgrp.
					);

 /*************************************************************************/
 /* ALL: Linear regression results, model estimates                       */
 /*************************************************************************/

/* 1.1) Output linear regression estimates data */
	
	data lin_hefi_raw0 lin_hefi_raw_boot;
		set regb._rcs_hefi_rawBS;
	if replicate = 0 then output lin_hefi_raw0;
	else output lin_hefi_raw_boot;
	run;
	
	/* note: rcs_hefi_raw = (Linear) regression with Restricted Cubic Spline transformation
		of the HEFI, raw association (i.e., unadjusted) */

/* 1.2) Use <boot_converge.sas> to look a normal-theory bootstrap variance estimation assumption */
%boot_converge(indata    = lin_hefi_raw_boot, /* bootstrap estimates data */
			   checkconv = rc_beta, /* variable for which we want to look at convergence */
			   where_sub = , /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );
	
/* 1.3) Use <boot_variance.sas> macro for efficient variance estimation */
%boot_variance(inboot    = lin_hefi_raw_boot, /* bootstrap estimates */
			   inbase    = lin_hefi_raw0,  /* original estimatee */
			   byvar     = , /* If <by> variable used */
			   varlist   = rc_beta r2,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.1, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-2), 
			   outdata   = reslib.lin_beta_allf /* name of output dataset */
			   ); 

	proc print data=reslib.lin_beta_allf label;
	id name nboot;
	var  estimate se lcl95 ucl95 pvalue svalue cv;
	run;

 /*************************************************************************/
 /* ALL: Linear regression results, predicted values for plot             */
 /*************************************************************************/

/* 2.1) Output linear regression plot data */
	data lin_hefi_plot0 lin_hefi_plot_boot;
		set regb._plm_rcs_hefi_rawBS;
	if replicate>0 then output lin_hefi_plot_boot;
		else output lin_hefi_plot0;
	run;
	
	/* note: plm_rcs_hefi_raw = Plot of Linear Model with RCS for HEFI,
		raw association (i.e., unadjusted) */
	
/* 2.2) Use <boot_converge.sas> to look a normal-theory bootstrap variance estimation assumption */
%boot_converge(indata    = lin_hefi_plot_boot, /* bootstrap estimates data */
			   checkconv = pred25 pred100 pred175, /* variable for which we want to look at convergence */
			   where_sub = , /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );

/* 2.3) Use <boot_variance.sas> macro for efficient variance estimation */
%boot_variance(inboot    = lin_hefi_plot_boot, /* bootstrap estimates */
			   inbase    = lin_hefi_plot0 ,  /* original estimatee */
			   where     = %str(name = "predicted"), /* variance for predicted values only */
			   byvar     = , /* If <by> variable used */
			   varlist   = pred1-pred200,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.1, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-1),
			   outdata   = lin_plot, /* name of output dataset */
			   echo =1
			   ); 

	/* 2.3.1) Add X values to the data */
		data reslib.lin_plot_allf;
		retain name HEFI2019_TOTAL_SCORE estimate ;
			merge reg.plmrcst_hefi_raw0(keep=HEFI2019_TOTAL_SCORE) lin_plot  ;
			/* no need for <by var;>, both data have same nb of obs */
		run;
	
	/* 2.3.2) Draft plot */
		proc sgplot data=reslib.lin_plotf ;
		title1 "Usual protein intake (grams)  | HEFI-2019 total score in adults 65y+ from CCHS 2015 - Nutrition";
		series x=HEFI2019_TOTAL_SCORE y=estimate / lineattrs=(color=black thickness=3);
		series x=HEFI2019_TOTAL_SCORE y=lcl95 / lineattrs=(color=black pattern=longdash thickness=3);
		series x=HEFI2019_TOTAL_SCORE y=ucl95 / lineattrs=(color=black pattern=longdash thickness=3);
		run;
		title1;

 /*************************************************************************/
 /* ALL: Logistic regression results, model estimates                     */
 /*************************************************************************/

/* 3.1) Output logistic regression estimates data */
		data log_hefi_mod0 log_hefi_mod_boot; ;
			set regb._logr_mod_hefi_y_rawBS;
		* split bootstrap and original;
		if replicate= 0 then output log_hefi_mod0; 
		else output log_hefi_mod_boot ;
		run;

	/* note: depending on the proportion of outcome, using log_beta or risk_ratio may be more stable.
		Hence, variance of both are estimated. */
	
/* 3.2) Use <boot_converge.sas> to look a normal-theory bootstrap variance estimation assumption */
%boot_converge(indata    = log_hefi_mod_boot, /* bootstrap estimates data */
			   checkconv = log_beta risk_ratio risk_diff, /* variable for which we want to look at convergence */
			   where_sub = %str(where (outcome=('binary_06'))), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );

/* 3.3) Use <boot_variance.sas> macro for efficient variance estimation */
%boot_variance(inboot    = log_hefi_mod_boot, /* bootstrap estimates */
			   inbase    = log_hefi_mod0 ,  /* original estimatee */
			   byvar     = outcome, /* If <by> variable used */
			   varlist   = risk_ratio,
			   null      = 1, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.2, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-2),
			   outdata   = log_beta_all1 /* name of output dataset */
			   ); 
			   
%boot_variance(inboot    = log_hefi_mod_boot, /* bootstrap estimates */
			   inbase    = log_hefi_mod0 ,  /* original estimatee */
			   byvar     = outcome, /* If <by> variable used */
			   varlist   = log_beta ,
			   exp       = 1 , /* exponentiate lcl + ucl once estimated */
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.2, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-1),
			   outdata   = log_beta_all2 /* name of output dataset */
			   );
			   
%boot_variance(inboot    = log_hefi_mod_boot, /* bootstrap estimates */
			   inbase    = log_hefi_mod0 ,  /* original estimatee */
			   byvar     = outcome, /* If <by> variable used */
			   varlist   = risk_diff ,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.2, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-2),
			   outdata   = log_beta_all3 /* name of output dataset */
			   ); 

	data reslib.log_beta_allf(rename=(_LABEL_=label));
		set log_beta_all1 log_beta_all2 log_beta_all3;
	* rename to facilitate interpretation ;
	if name = "log_beta" then name ="odds_ratio";
	if missing(exp) then exp=0;
	* change erroneous label ;
	_LABEL_ = compress(_LABEL_,'_6');
	label _LABEL_ = " ";
	run;
	
	proc sort data=reslib.log_beta_allf;
		by outcome;
	run;

/* Clean temp */
	proc datasets lib=work nolist nodetails ;
	delete _: log_beta: log_hefi: ;
	run;

 /*************************************************************************/
 /* BY DRI GROUP: Linear regression results, model estimates              */
 /*************************************************************************/

/* 1.1) Output linear regression estimates data */
	data lin_hefi_raw0 lin_hefi_raw_boot;
		set reg_b2._rcs_hefi_rawBS;
		if replicate>0 then output lin_hefi_raw_boot ;
		else output lin_hefi_raw0;
	run;

/* 1.2) Use <boot_converge.sas> to look a normal-theory bootstrap variance estimation assumption */
%boot_converge(indata    = lin_hefi_raw_boot, /* bootstrap estimates data */
			   checkconv = rc_beta, /* variable for which we want to look at convergence */
			   where_sub = %str(where (drig=12)), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );

	
/* 1.3) Use <boot_variance.sas> macro for efficient variance estimation */
%boot_variance(inboot    = lin_hefi_raw_boot, /* bootstrap estimates */
			   inbase    = lin_hefi_raw0,  /* original estimatee */
			   byvar     = drig, /* If <by> variable used */
			   varlist   = rc_beta ,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.1, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-2), 
			   outdata   = reslib.lin_beta_drigf /* name of output dataset */
			   ); 

	proc print data=reslib.lin_beta_drigf label;
	id name nboot drig;
	var  estimate se lcl95 ucl95 pvalue svalue cv;
	run;

 /*************************************************************************/
 /* BY DRI GROUP: Linear regression results, predicted values for plot    */
 /*************************************************************************/

/* 2.1) Output linear regression plot data */
	data lin_hefi_plot0 lin_hefi_plot_boot;	
		set reg_b2._plm_rcs_hefi_rawBS;
	if replicate>0 then output lin_hefi_plot_boot;
	else output lin_hefi_plot0 ;
	run;
	
/* 2.2) Use <boot_converge.sas> to look a normal-theory bootstrap variance estimation assumption */
%boot_converge(indata    = lin_hefi_plot_boot, /* bootstrap estimates data */
			   checkconv = pred25 pred100 pred175, /* variable for which we want to look at convergence */
			   where_sub = %str(where drig='12'), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );

/* 2.3) Use <boot_variance.sas> macro for efficient variance estimation */
%boot_variance(inboot    = lin_hefi_plot_boot, /* bootstrap estimates */
			   inbase    = lin_hefi_plot0 ,  /* original estimatee */
			   where     = %str(name = "predicted"), /* variance for predicted values only */
			   byvar     = drig, /* If <by> variable used */
			   varlist   = pred1-pred200,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.1, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-1),
			   outdata   = lin_plot /* name of output dataset */
			   ); 

	/* 2.3.1) Add X values to the data */
		data reslib.lin_plot_drigf(rename=(_drig=drig));
		retain _drig name pred_id HEFI2019_TOTAL_SCORE estimate ;
			merge reg_2.plmrcst_hefi_raw0(keep=drig HEFI2019_TOTAL_SCORE) lin_plot  ;
			by drig;
		* change char -> num ;
		_drig =input(drig,10.);
		drop drig;
		* make numerical index for predicted values;
		pred_id = input(compress(name,,'a'),10.);
		run;
	
	/* 2.3.2) Draft plot */
		proc sgplot data=reslib.lin_plot_drigf ;
		title1 "Usual protein intake (grams)  | HEFI-2019 total score in adults 65y+ from CCHS 2015 - Nutrition, by DRIG";
		series x=HEFI2019_TOTAL_SCORE y=estimate / lineattrs=( thickness=3) group=drig;
		*series x=HEFI2019_TOTAL_SCORE y=lcl95 / lineattrs=(color=black pattern=longdash thickness=3) ;
		*series x=HEFI2019_TOTAL_SCORE y=ucl95 / lineattrs=(color=black pattern=longdash thickness=3) ;
		run;
		title1;

 /*************************************************************************/
 /* BY DRI GROUP: Logistic regression results, model estimates            */
 /*************************************************************************/

/* 3.1) Output logistic regression estimates data */
	proc sort data=reg_b2._logr_mod_hefi_y_rawBS;
		by replicate drig outcome ;
	run;

	/* split sample for the <boot_variance> macro */
		data log_hefi_mod0 log_hefi_mod_boot;
			set reg_b2._logr_mod_hefi_y_rawBS;
		if replicate= 0 then output log_hefi_mod0; 
		else output log_hefi_mod_boot ;
		run;


/* 3.2) Use <boot_converge.sas> to look a normal-theory bootstrap variance estimation assumption */
%boot_converge(indata    = log_hefi_mod_boot, /* bootstrap estimates data */
			   checkconv = exp_beta risk_diff risk_ratio, /* variable for which we want to look at convergence */
			   where_sub = %str(where (outcome='binary_06') & drig=15), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );

/* 3.3) Use <boot_variance.sas> macro for efficient variance estimation */
%boot_variance(inboot    = log_hefi_mod_boot, /* bootstrap estimates */
			   inbase    = log_hefi_mod0 ,  /* original estimatee */
			   byvar     = drig outcome, /* If <by> variable used */
			   varlist   = risk_ratio,
			   null      = 1, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.2, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-2),
			   outdata   = log_beta_all1 /* name of output dataset */
			   ); 

%boot_variance(inboot    = log_hefi_mod_boot, /* bootstrap estimates */
			   inbase    = log_hefi_mod0 ,  /* original estimatee */
			   byvar     = drig outcome, /* If <by> variable used */
			   varlist   = log_beta ,
			   exp       = 1,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.2, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-1),
			   outdata   = log_beta_all2 /* name of output dataset */
			   ); 

%boot_variance(inboot    = log_hefi_mod_boot, /* bootstrap estimates */
			   inbase    = log_hefi_mod0 ,  /* original estimatee */
			   byvar     = drig outcome, /* If <by> variable used */
			   varlist   = risk_diff ,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.2, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-2),
			   outdata   = log_beta_all3 /* name of output dataset */
			   ); 

	data reslib.log_beta_drigf(rename=(_LABEL_=label));
		set log_beta_all1 log_beta_all2 log_beta_all3;
	* rename to facilitate interpretation ;
	if name = "log_beta" then name ="odds_ratio";
	if missing(exp) then exp=0;
	* change erroneous label ;
	_LABEL_ = compress(_LABEL_,'_6');
	label _LABEL_ = " ";
	run;
	
	proc sort data=reslib.log_beta_drigf;
		by outcome drig;
	run;

/* Clean temp */
	proc datasets lib=work nolist nodetails ;
	delete _: log_beta: log_hefi: ;
	run;

 /*************************************************************************/
 /* ALL: Distribution of usual intake, marginal Pr(X<x)                   */
 /*************************************************************************/

/* 1.1) Output marginal Y distribution data */
	
	/* note: distrib_y_rel_w_boot = Distribution of Y, Relative intakes, Wide format */

	proc sort data=bootlib._distrib_y_rel_wBS;
	by replicate DRIg;
	run;
	
	/* split sample for the <boot_variance> macro */
	data distrib_y_rel_w0 distrib_y_rel_w_boot;
		set bootlib._distrib_y_rel_wBS;
	* For PREVALENCE: rescale prob to percentage point ;
	array prob(*) Prob1-Prob5;
	do i=1 to dim(prob);
	prob(i) = prob(i) * 100 ;
	end;
	drop i;
	if replicate>0 then output distrib_y_rel_w_boot;
	else output distrib_y_rel_w0;
	run;

/* 1.2) Use <boot_converge.sas> to look a normal-theory bootstrap variance estimation assumption */
%boot_converge(indata    = distrib_y_rel_w_boot, /* bootstrap estimates data */
			   checkconv = Mean Pctile10 Pctile90, /* variable for which we want to look at convergence */
			   where_sub = %str(where (drig=0)), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );

	
/* 1.3) Use <boot_variance.sas> macro for efficient variance estimation */

/* 1.3.2) actual variance estimation */
%boot_variance(inboot    = distrib_y_rel_w_boot, /* bootstrap estimates */
			   inbase    = distrib_y_rel_w0,  /* original estimatee */
			   byvar     = drig, /* If <by> variable used */
			   varlist   = Mean &pctile_list ,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 4.2, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-1), 
			   outdata   = distrib_y_rel_wf1 /* name of output dataset */
			   );

%boot_variance(inboot    = distrib_y_rel_w_boot, /* bootstrap estimates */
			   inbase    = distrib_y_rel_w0,  /* original estimatee */
			   byvar     = drig, /* If <by> variable used */
			   varlist   = Prob1-Prob5,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 4.1, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-1), 
			   outdata   = distrib_y_rel_wf2 /* name of output dataset */
			   ); 
			  
	data reslib.distrib_y_rel_wf (rename=(_LABEL_ = label));
	retain drig name _LABEL_ ;
		set distrib_y_rel_wf1(in=a) distrib_y_rel_wf2(in=b) ;
	label _LABEL_ = " ";
	run;
	
	proc sort;
	by drig;
	run;

	proc print data=reslib.distrib_y_rel_wf label;
	id name label nboot drig;
	var  estimate se lcl95 ucl95 pvalue svalue cv;
	run;

/* Clean temp */
	proc datasets lib=work nolist nodetails ;
	delete _: distrib_: ;
	run;

 /*************************************************************************/
 /*                                                                       */
 /*         REF. CODE 3: Variance estimation for CA/VITD outcome          */
 /*                                                                       */
 /*************************************************************************/

 /************************************************/
 /*            Files and folder set-up           */
 /************************************************/

/* 0) define suffix of current outcome */
	
	%let suffix = _ca_vit_d ;

/* 1) assign proper library using the same suffix */
	%mcmclib(&suffix.);
	%MakeMCMCFolder(NCIPath = &path./NCI/MCMC&suffix.,
					prefix  = Reg_,
					suffix  = &subgrp.
					);

 /*************************************************************************/
 /* ALL: Linear regression results, model estimates                       */
 /*************************************************************************/

/* 1.1) Output linear regression estimates data */
	
	proc sort data=regb._rcs_hefi_rawBS;
	by replicate;
	run;

	/* note: rcs_hefi_raw = (Linear) regression with Restricted Cubic Spline transformation
		of the HEFI, raw association (i.e., unadjusted) */

	/* split sample for the <boot_variance> macro */
	data lin_hefi_raw0 lin_hefi_raw_boot ;
		set regb._rcs_hefi_rawBS ;
	if replicate=0 then output lin_hefi_raw0; 
	else output lin_hefi_raw_boot;
	run;

/* 1.2) Use <boot_converge.sas> to look a normal-theory bootstrap variance estimation assumption */
%boot_converge(indata    = regb._rcs_hefi_rawBS, /* bootstrap estimates data */
			   checkconv = rc_beta, /* variable for which we want to look at convergence */
			   where_sub = %str(where varname = 'calcium'), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );
%boot_converge(indata    = regb._rcs_hefi_rawBS, /* bootstrap estimates data */
			   checkconv = rc_beta, /* variable for which we want to look at convergence */
			   where_sub = %str(where varname = 'vit_d'), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );

/* 1.3) Use <boot_variance.sas> macro for efficient variance estimation */
%boot_variance(inboot    = lin_hefi_raw_boot, /* bootstrap estimates */
			   inbase    = lin_hefi_raw0,  /* original estimatee */
			   byvar     = varname ,
			   varlist   = rc_beta r2,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.1, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-2), 
			   outdata   = reslib.lin_beta_allf /* name of output dataset */
			   ); 

	proc print data=reslib.lin_beta_allf label;
	id varname name nboot;
	var  estimate se lcl95 ucl95 pvalue svalue cv;
	run;

 /*************************************************************************/
 /* ALL: Linear regression results, predicted values for plot             */
 /*************************************************************************/

/* 2.1) Output linear regression plot data */
	proc sort data=regb._lin_hefi_plotBS;
	by replicate varname name;
	run;

	/* note: plm_rcs_hefi_raw = Plot of Linear Model with RCS for HEFI,
		raw association (i.e., unadjusted) */
	
	/* split sample for the <boot_variance> macro */
	data  lin_hefi_plot0 lin_hefi_plot_boot;
		set regb._lin_hefi_plotBS;
	* remove x values;
	if (name = "HEFI2019_TOTAL_SCORE") then delete;
	if replicate =0 then output lin_hefi_plot0;
	else output lin_hefi_plot_boot;
	run;
	
	
/* 2.2) Use <boot_converge.sas> to look a normal-theory bootstrap variance estimation assumption */
%boot_converge(indata    = lin_hefi_plot_boot, /* bootstrap estimates data */
			   checkconv = pred25 pred100 pred175, /* variable for which we want to look at convergence */
			   where_sub = %str(where varname = 'calcium'), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );
%boot_converge(indata    = lin_hefi_plot_boot, /* bootstrap estimates data */
			   checkconv = pred25 pred100 pred175, /* variable for which we want to look at convergence */
			   where_sub = %str(where varname = 'vit_d'), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );

/* 2.3) Use <boot_variance.sas> macro for efficient variance estimation */
%boot_variance(inboot    = lin_hefi_plot_boot, /* bootstrap estimates */
			   inbase    = lin_hefi_plot0 ,  /* original estimatee */
			   where     = %str(name = "predicted"), /* variance for predicted values only */	   
			   byvar     = varname ,
			   varlist   = pred1-pred200,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.1, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-1),
			   outdata   = lin_plot /* name of output dataset */
			   ); 

	/* 2.3.1) Add X values to the data */
	
	proc transpose data=regb._lin_hefi_plotBS out=_xvalues prefix=HEFI2019_TOTAL_SCORE;
	by varname ;
	var pred1-pred200;
	where name = "HEFI2019_TOTAL_SCORE";
	run;
	
	proc sort data=_xvalues;
	by varname _NAME_ ;
	run;
	
	proc sort data=lin_plot;
	by varname name ;
	run;
	
	data reslib.lin_plot_allf;
	retain name HEFI2019_TOTAL_SCORE estimate ;
		merge _xvalues(rename=(HEFI2019_TOTAL_SCORE1=HEFI2019_TOTAL_SCORE _NAME_=name))
			  lin_plot  ;
		by varname name ;
	* sort from 1 to 200 ;
	pred_id = input(compress(name,,'a'),10.);
	run;
	
	proc sort data=reslib.lin_plot_allf;
	by varname pred_id ;
	run;
	
	/* 2.3.2) Draft plot */
		proc sgplot data=reslib.lin_plot_allf ;
		title1 "Usual Calcium intake (mg)  | HEFI-2019 total score in adults 65y+ from CCHS 2015 - Nutrition";
		series x=HEFI2019_TOTAL_SCORE y=estimate / lineattrs=(color=black thickness=3);
		series x=HEFI2019_TOTAL_SCORE y=lcl95 / lineattrs=(color=black pattern=longdash thickness=3);
		series x=HEFI2019_TOTAL_SCORE y=ucl95 / lineattrs=(color=black pattern=longdash thickness=3);
		where varname = "calcium" ;
		run;
		title1;
		proc sgplot data=reslib.lin_plot_allf ;
		title1 "Usual Vit. D intake (mcg)  | HEFI-2019 total score in adults 65y+ from CCHS 2015 - Nutrition";
		series x=HEFI2019_TOTAL_SCORE y=estimate / lineattrs=(color=black thickness=3);
		series x=HEFI2019_TOTAL_SCORE y=lcl95 / lineattrs=(color=black pattern=longdash thickness=3);
		series x=HEFI2019_TOTAL_SCORE y=ucl95 / lineattrs=(color=black pattern=longdash thickness=3);
		where varname = "vit_d" ;
		run;
		title1;

 /*************************************************************************/
 /* ALL: Logistic regression results, model estimates                     */
 /*************************************************************************/

/* 3.1) Output logistic regression estimates data */
	
	proc sort data=regb._logr_mod_hefi_y_rawBS;
	by replicate varname ;
	run;

	/* note: logr_mod_hefi_y = Logistic Regression Model for HEFI with outcome Y, 
		raw association (i.e., unadjusted) */

	data log_hefi_mod0 log_hefi_mod_boot;
		set regb._logr_mod_hefi_y_rawBS ;
	* split bootstrap and original;
		if replicate= 0 then output log_hefi_mod0; 
		else output log_hefi_mod_boot ;
	run;

	/* note: depending on the proportion of outcome, using log_beta or risk_ratio may be more stable.
		Hence, variance of both are estimated. */
	
/* 3.2) Use <boot_converge.sas> to look a normal-theory bootstrap variance estimation assumption */
%boot_converge(indata    = log_hefi_mod_boot, /* bootstrap estimates data */
			   checkconv = log_beta risk_ratio, /* variable for which we want to look at convergence */
			   where_sub = %str(where (varname=('binary_ca'))), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );
%boot_converge(indata    = log_hefi_mod_boot, /* bootstrap estimates data */
			   checkconv = log_beta risk_ratio, /* variable for which we want to look at convergence */
			   where_sub = %str(where (varname=('binary_vitd'))), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );


/* 3.3) Use <boot_variance.sas> macro for efficient variance estimation */
%boot_variance(inboot    = log_hefi_mod_boot, /* bootstrap estimates */
			   inbase    = log_hefi_mod0 ,  /* original estimatee */
			   byvar     = varname, /* If <by> variable used */
			   varlist   = risk_ratio,
			   null      = 1, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.2, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-2),
			   outdata   = log_beta_all1 /* name of output dataset */
			   );

%boot_variance(inboot    = log_hefi_mod_boot, /* bootstrap estimates */
			   inbase    = log_hefi_mod0 ,  /* original estimatee */
			   byvar     = varname, /* If <by> variable used */
			   varlist   = log_beta ,
			   exp       = 1 , /* exponentiate lcl + ucl once estimated */
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.2, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-1),
			   outdata   = log_beta_all2 /* name of output dataset */
			   ); 

%boot_variance(inboot    = log_hefi_mod_boot, /* bootstrap estimates */
			   inbase    = log_hefi_mod0 ,  /* original estimatee */
			   byvar     = varname, /* If <by> variable used */
			   varlist   = risk_diff ,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.2, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-2),
			   outdata   = log_beta_all3 /* name of output dataset */
			   ); 

	data reslib.log_beta_allf(rename=(_LABEL_=label));
		set log_beta_all1 log_beta_all2 log_beta_all3;
	* rename to facilitate interpretation ;
		if name = "log_beta" then name ="odds_ratio";
		if missing(exp) then exp=0;
	* change erroneous label ;
		_LABEL_ = substr(_LABEL_,index(_LABEL_,'|')+2);
		label _LABEL_ = " ";
	run;
	
	proc sort data=reslib.log_beta_allf;
		by varname;
	run;

/* Clean temp */
	proc datasets lib=work nolist nodetails ;
	delete _: lin: log_: ;
	run;

 /*************************************************************************/
 /* BY DRI GROUP: Linear regression results, model estimates              */
 /*************************************************************************/

/* 1.1) Output linear regression estimates data */	
	proc sort data=reg_b2._rcs_hefi_rawBS;
	by replicate varname drig;
	run;

	/* note: rcs_hefi_raw = (Linear) regression with Restricted Cubic Spline transformation
		of the HEFI, raw association (i.e., unadjusted) */

	/* split sample for the <boot_variance> macro */
	data lin_hefi_raw0 lin_hefi_raw_boot ;
		set reg_b2._rcs_hefi_rawBS ;
	if replicate=0 then output lin_hefi_raw0; 
	else output lin_hefi_raw_boot;
	run;

/* 1.2) Use <boot_converge.sas> to look a normal-theory bootstrap variance estimation assumption */
%boot_converge(indata    = lin_hefi_raw_boot, /* bootstrap estimates data */
			   checkconv = rc_beta, /* variable for which we want to look at convergence */
			   where_sub = %str(where (drig=12 and varname='calcium')), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );
%boot_converge(indata    = lin_hefi_raw_boot, /* bootstrap estimates data */
			   checkconv = rc_beta, /* variable for which we want to look at convergence */
			   where_sub = %str(where (drig=15 and varname='calcium')), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );

/* 1.3) Use <boot_variance.sas> macro for efficient variance estimation */
%boot_variance(inboot    = lin_hefi_raw_boot, /* bootstrap estimates */
			   inbase    = lin_hefi_raw0,  /* original estimatee */
			   byvar     = varname drig, /* If <by> variable used */
			   varlist   = rc_beta ,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.1, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-2), 
			   outdata   = reslib.lin_beta_drigf /* name of output dataset */
			   ); 

	proc print data=reslib.lin_beta_drigf label;
	id varname drig name nboot ;
	var  estimate se lcl95 ucl95 pvalue svalue cv;
	run;

 /*************************************************************************/
 /* BY DRI GROUP: Linear regression results, predicted values for plot    */
 /*************************************************************************/

/* 2.1) Output linear regression plot data */	
	proc sort data=reg_b2._lin_hefi_plotBS;
	by replicate varname drig name;
	run;

	/* note: plm_rcs_hefi_raw = Plot of Linear Model with RCS for HEFI,
		raw association (i.e., unadjusted) */
	
	/* split sample for the <boot_variance> macro */
	data  lin_hefi_plot0 lin_hefi_plot_boot;
		set reg_b2._lin_hefi_plotBS;
	* remove x values;
	if (name = "HEFI2019_TOTAL_SCORE") then delete;
	if replicate =0 then output lin_hefi_plot0;
	else output lin_hefi_plot_boot;
	run;
	
/* 2.2) Use <boot_converge.sas> to look a normal-theory bootstrap variance estimation assumption */
%boot_converge(indata    = lin_hefi_plot_boot, /* bootstrap estimates data */
			   checkconv = pred25 pred100 pred175, /* variable for which we want to look at convergence */
			   where_sub = %str(where drig='12' and varname='calcium'), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );
%boot_converge(indata    = lin_hefi_plot_boot, /* bootstrap estimates data */
			   checkconv = pred25 pred100 pred175, /* variable for which we want to look at convergence */
			   where_sub = %str(where drig='15' and varname='calcium'), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );
			   
/* 2.3) Use <boot_variance.sas> macro for efficient variance estimation */
%boot_variance(inboot    = lin_hefi_plot_boot, /* bootstrap estimates */
			   inbase    = lin_hefi_plot0 ,  /* original estimatee */
			   where     = %str(name = "predicted"), /* variance for predicted values only */
			   byvar     = varname drig,
			   varlist   = pred1-pred200,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.1, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-1),
			   outdata   = lin_plot /* name of output dataset */
			   ); 

	/* 2.3.1) Add X values to the data */
		
	data _xvalues;
		set reg_2.plmrcst_calcium_raw0(in=a)
			reg_2.plmrcst_vit_d_raw0 (in=b);
	length varname $ 32;
	if a then varname = "calcium";
	else if b then varname = "vit_d";
	run;
		
		/* some manipulation to output <name> variable*/
		proc sort;
			by varname drig;
		run;
		
		data _xvalues;
		set _xvalues;
		by varname drig;
		retain pred_id ;
		if first.drig then pred_id=1;
		else pred_id+1;
		length name $  8;
		name = cats("pred",pred_id);
		run;
		
	proc sort data=_xvalues;
	by varname drig name ;
	run;
	
	proc sort data=lin_plot;
	by varname drig name ;
	run;
	
	data reslib.lin_plot_drigf(rename=(_drig=drig));
	retain varname _drig pred_id name HEFI2019_TOTAL_SCORE estimate ;
		merge _xvalues
			  lin_plot  ;
		by varname drig name ;
	* change drig to num ;
		_drig = input(drig, 10.);
	* add merge error flag for extra safety ;
		if predicted ne estimate then _MergeError_=1;
		else _MergeError_=0;
		label _MergeError_ = "Flag for any merge error";
	*clean; 
		drop drig predicted;
	run;
	
	proc sort data=reslib.lin_plot_drigf;
	by varname drig pred_id ;
	run;	

	
	/* 2.3.2) Draft plot */
		proc sgplot data=reslib.lin_plot_drigf ;
		title1 "Usual calcium (mg)  | HEFI-2019 total score in adults 65y+ from CCHS 2015 - Nutrition, by DRIG";
		series x=HEFI2019_TOTAL_SCORE y=estimate / lineattrs=( thickness=3) group=drig;
		*series x=HEFI2019_TOTAL_SCORE y=lcl95 / lineattrs=(color=black pattern=longdash thickness=3) ;
		*series x=HEFI2019_TOTAL_SCORE y=ucl95 / lineattrs=(color=black pattern=longdash thickness=3) ;
		where varname = "calcium";
		run;
		title1;
		proc sgplot data=reslib.lin_plot_drigf ;
		title1 "Usual vit_d (mcg)  | HEFI-2019 total score in adults 65y+ from CCHS 2015 - Nutrition, by DRIG";
		series x=HEFI2019_TOTAL_SCORE y=estimate / lineattrs=( thickness=3) group=drig;
		*series x=HEFI2019_TOTAL_SCORE y=lcl95 / lineattrs=(color=black pattern=longdash thickness=3) ;
		*series x=HEFI2019_TOTAL_SCORE y=ucl95 / lineattrs=(color=black pattern=longdash thickness=3) ;
		where varname = "vit_d";
		run;
		title1;


/* Clean temp */
	proc datasets lib=work nolist nodetails ;
	delete _: lin: ;
	run;

 /*************************************************************************/
 /* BY DRI GROUP: Logistic regression results, model estimates            */
 /*************************************************************************/

/* 3.1) Output logistic regression estimates data */	
	proc sort data=reg_b2._logr_mod_hefi_y_rawBS;
	by replicate varname drig;
	run;

	/* note: logr_mod_hefi_y = Logistic Regression Model for HEFI with outcome Y, 
		raw association (i.e., unadjusted) */

	data log_hefi_mod0 log_hefi_mod_boot;
		set reg_b2._logr_mod_hefi_y_rawBS ;
	* split bootstrap and original;
		if replicate= 0 then output log_hefi_mod0; 
		else output log_hefi_mod_boot ;
	run;

/* 3.2) Use <boot_converge.sas> to look a normal-theory bootstrap variance estimation assumption */
*%boot_converge(indata    = log_hefi_mod_boot, /* bootstrap estimates data */
			   checkconv = exp_beta risk_diff risk_ratio, /* variable for which we want to look at convergence */
			   where_sub = %str(where (outcome=('binary_06'))), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );

/* 3.3) Use <boot_variance.sas> macro for efficient variance estimation */
%boot_variance(inboot    = log_hefi_mod_boot, /* bootstrap estimates */
			   inbase    = log_hefi_mod0 ,  /* original estimatee */
			   byvar     = varname drig, /* If <by> variable used */
			   varlist   = risk_ratio,
			   null      = 1, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.2, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-2),
			   outdata   = log_beta_all1 /* name of output dataset */
			   ); 

%boot_variance(inboot    = log_hefi_mod_boot, /* bootstrap estimates */
			   inbase    = log_hefi_mod0 ,  /* original estimatee */
			   byvar     = varname drig, /* If <by> variable used */
			   varlist   = log_beta ,
			   exp       = 1,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.2, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-1),
			   outdata   = log_beta_all2 /* name of output dataset */
			   );

%boot_variance(inboot    = log_hefi_mod_boot, /* bootstrap estimates */
			   inbase    = log_hefi_mod0 ,  /* original estimatee */
			   byvar     = varname drig, /* If <by> variable used */
			   varlist   = risk_diff ,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.2, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-2),
			   outdata   = log_beta_all3 /* name of output dataset */
			   ); 

	data reslib.log_beta_drigf(rename=(_LABEL_=label));
		set log_beta_all1 log_beta_all2 log_beta_all3;
	* rename to facilitate interpretation ;
		if name = "log_beta" then name ="odds_ratio";
		if missing(exp) then exp=0;
	* change erroneous label ;
		_LABEL_ = substr(_LABEL_,index(_LABEL_,'|')+2);
		label _LABEL_ = " ";
	run;
	
	proc sort data=reslib.log_beta_drigf;
		by varname drig;
	run;

/* Clean temp */
	proc datasets lib=work nolist nodetails ;
	delete _: log_beta: log_hefi: ;
	run;

 /*************************************************************************/
 /* ALL: Distribution of usual intake, marginal Pr(X<x)                   */
 /*************************************************************************/

/* 1.1) Output marginal Y distribution data */
	data bootlib._distrib_y_wBS;
		set reslib.distrib_y_ca_w0
			reslib.distrib_y_vitd_w0
			bootlib.distrib_y_: ;
	* For PREVALENCE: rescale prob to percentage point ;
	array prob(*) binary_ca binary_vitd;
	do i=1 to dim(prob);
	prob(i) = prob(i) * 100 ;
	end;
	drop i;
	* Combine both variables for easy variance estimation;
	if varname = "calcium" then binary_EAR = binary_ca;
	else if varname = "vit_d" then binary_EAR = binary_vitd;
	label binary_EAR = "Pr(X<EAR)";
	run;

	proc sort data=bootlib._distrib_y_wBS;
	by replicate varname DRIg;
	run;
	
	/* note: distrib_y_rel_w_boot = Distribution of Y, Relative intakes, Wide format */

	/* 1.2) split sample for the <boot_variance> macro */
	data distrib_y_w0 distrib_y_w_boot ;
		set bootlib._distrib_y_wBS;
	if replicate= 0 then output distrib_y_w0;
	else output distrib_y_w_boot;
	run;
	

/* 1.2) Use <boot_converge.sas> to look a normal-theory bootstrap variance estimation assumption */
%boot_converge(indata    = distrib_y_w_boot, /* bootstrap estimates data */
			   checkconv = Mean Pctile10 Pctile90, /* variable for which we want to look at convergence */
			   where_sub = %str(where (drig=0 and varname='calcium')), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );

	
/* 1.3) Use <boot_variance.sas> macro for efficient variance estimation */

/* 1.3.2) actual variance estimation */
%boot_variance(inboot    = distrib_y_w_boot, /* bootstrap estimates */
			   inbase    = distrib_y_w0,  /* original estimatee */
			   byvar     = varname drig, /* If <by> variable used */
			   varlist   = Mean &pctile_list ,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 4.2, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-1), 
			   outdata   = distrib_y_wf1 /* name of output dataset */
			   );

%boot_variance(inboot    = distrib_y_w_boot, /* bootstrap estimates */
			   inbase    = distrib_y_w0,  /* original estimatee */
			   byvar     = varname drig, /* If <by> variable used */
			   varlist   = binary_ear ,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 4.1, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-1), 
			   outdata   = distrib_y_wf2 /* name of output dataset */
			   ); 
			  
	data reslib.distrib_y_wf (rename=(_LABEL_ = label));
	retain drig name _LABEL_ ;
	length name $ 32;
		set distrib_y_wf1(in=a) distrib_y_wf2(in=b) ;
	label _LABEL_ = " ";
	run;
	
	proc sort;
	by varname drig;
	run;

	proc print data=reslib.distrib_y_wf label;
	id varname name label nboot drig;
	var  estimate se lcl95 ucl95 pvalue svalue cv;
	run;

 /*************************************************************************/
 /* BY CATEGORY OF SCORE: Distribution of usual intake                    */
 /*************************************************************************/

/* 2.1) Append marginal Y distribution data */
	data bootlib._distrib_catscore_wBS;
		set reslib.distrib_w_catscore0
			bootlib.distrib_w_catscore: 
			;
	* PREVALENCE: rescale prob to percentage point;
	if catscore < 9000 then binary_ear = binary_ear*100;
	run;
	
	/* note: distrib_y_rel_w_boot = Distribution of Y, Relative intakes, Wide format */

	proc sort data=bootlib._distrib_catscore_wBS;
	by replicate varname catscore;
	run;
	
	/* 2.1.2) split sample for the <boot_variance> macro */
	data distrib_catscore0 distrib_catscore_boot ;
		set bootlib._distrib_catscore_wBS;
	if replicate= 0 then output distrib_catscore0;
	else output distrib_catscore_boot;
	run;


/* 2.2) Use <boot_converge.sas> to look a normal-theory bootstrap variance estimation assumption */
%boot_converge(indata    = distrib_catscore_boot, /* bootstrap estimates data */
			   checkconv = Mean Pctile10 Pctile90, /* variable for which we want to look at convergence */
			   where_sub = %str(where (catscore=4 and varname='calcium')), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );

	
/* 2.3) Use <boot_variance.sas> macro for efficient variance estimation */

	/* 2.3.2) Variance estimation */
	%boot_variance(inboot    = distrib_catscore_boot, /* bootstrap estimates */
				   inbase    = distrib_catscore0,  /* original estimatee */
				   byvar     = varname catscore, /* If <by> variable used */
				   varlist   = Mean &pctile_list ,
				   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
				   set_confidence_level=95, /* 95 = 95%CI */
				   numberfmt = 5.2, /* formatting of estimates */
				   set_degrees_freedom = %str(nboot-1), 
				   outdata   = distrib_catscoref1 /* name of output dataset */
				   ); 
	
	%boot_variance(inboot    = distrib_catscore_boot, /* bootstrap estimates */
				   inbase    = distrib_catscore0,  /* original estimatee */
				   byvar     = varname catscore, /* If <by> variable used */
				   varlist   = binary_ear ,
				   null      = 1, /* Null H for pvalues (i.e., 0 for a difference) */
				   set_confidence_level=95, /* 95 = 95%CI */
				   numberfmt = 4.2, /* formatting of estimates */
				   set_degrees_freedom = %str(nboot-2), 
				   outdata   = distrib_catscoref2 /* name of output dataset */
				   ); 
				   
 	/* 2.3.3) Check and save output */
	data reslib.distrib_w_catscoref(rename=(_LABEL_=label));
	retain name _LABEL_ catscore;
	length name $ 32;
		set distrib_catscoref1 distrib_catscoref2 ;
	label _LABEL_ = " ";
	run;
	
	proc sort;
	by varname catscore;
	run;
	
	title1 "Nutrient intake by quarter of total HEFI-2019 score";
	proc print data=reslib.distrib_w_catscoref label;
	id varname name nboot catscore;
	var  estimate se lcl95 ucl95 pvalue svalue cv;
	where index(name,"Prob")=0;
	run;

	title2 "Pr(X<x) by quarter of total HEFI-2019 score";
	proc print data=reslib.distrib_w_catscoref label;
	id varname name label nboot catscore;
	var  estimate se lcl95 ucl95 pvalue svalue cv;
	where index(name,"Prob")>0;
	run;
	title1;
	title2;

/* Clean temp */
	proc datasets lib=work nolist nodetails ;
	delete _: distrib_: ;
	run;

 /*************************************************************************/
 /*                                                                       */
 /*     REF. CODE 4: Variance estimation for MISC.A NUTRIENT outcome      */
 /*                                                                       */
 /*************************************************************************/

data _null_;
	file print;
	put #3 @10 "Bootstrap variance estimation for miscellaneous A (&SYSDATE)";
run;

/* nutrient intakes: iron, zinc, vitamin B6 and B12 */
/* binary outcomes: binary_iron binary_zinc binary_b6 binary_b12 */

 /************************************************/
 /*            Files and folder set-up           */
 /************************************************/

/* 0) define suffix of current outcome */
	
	%let suffix = _miscA ;

/* 1) assign proper library using the same suffix */
	%mcmclib(&suffix.);
	%MakeMCMCFolder(NCIPath = &path./NCI/MCMC&suffix.,
					prefix  = Reg_,
					suffix  = &subgrp.
					);

 /*************************************************************************/
 /* ALL: Linear regression results, model estimates                       */
 /*************************************************************************/

/* 1.1) Output linear regression estimates data */
	data regb._rcs_hefi_rawBS;
	retain replicate varname ;
	length label $ 65;
		set reg.rcs_iron_raw0  reg.rcs_zinc_raw0 
			regb.rcs_iron_raw: regb.rcs_zinc_raw: 
			reg.rcs_vit_b6_raw0  reg.rcs_vit_b12_raw0 
			regb.rcs_vit_b6_raw: regb.rcs_vit_b12_raw: 
			;
	* use the model label to make shorter nutrient identifier (ie, varname); 
	length varname $ 32;
	varname = substr(label,1,index(label,'|')-1);
	run;
	
	/* note: rcs_hefi_raw = (Linear) regression with Restricted Cubic Spline transformation
		of the HEFI, raw association (i.e., unadjusted) */
	
	proc sort data=regb._rcs_hefi_rawBS;
	by replicate;
	run;

	/* split sample for the <boot_variance> macro */
	data lin_hefi_raw0 lin_hefi_raw_boot ;
		set regb._rcs_hefi_rawBS ;
	if replicate=0 then output lin_hefi_raw0; 
	else output lin_hefi_raw_boot;
	run;

/* 1.2) Use <boot_converge.sas> to look a normal-theory bootstrap variance estimation assumption */
%boot_converge(indata    = regb._rcs_hefi_rawBS, /* bootstrap estimates data */
			   checkconv = rc_beta, /* variable for which we want to look at convergence */
			   where_sub = %str(where varname = 'iron'), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );
%boot_converge(indata    = regb._rcs_hefi_rawBS, /* bootstrap estimates data */
			   checkconv = rc_beta, /* variable for which we want to look at convergence */
			   where_sub = %str(where varname = 'vit_b12'), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );

/* 1.3) Use <boot_variance.sas> macro for efficient variance estimation */
%boot_variance(inboot    = lin_hefi_raw_boot, /* bootstrap estimates */
			   inbase    = lin_hefi_raw0,  /* original estimatee */
			   byvar     = varname ,
			   varlist   = rc_beta r2,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.1, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-2), 
			   outdata   = reslib.lin_beta_allf /* name of output dataset */
			   ); 

	proc print data=reslib.lin_beta_allf label;
	id varname name nboot;
	var  estimate se lcl95 ucl95 pvalue svalue cv;
	run;

 /*************************************************************************/
 /* ALL: Linear regression results, predicted values for plot             */
 /*************************************************************************/

/* 2.1) Output linear regression plot data */
	data regb._lin_hefi_plotBS;
	retain replicate varname ;
		set reg.plm_rcs_iron_raw0(in=a1)
			regb.plm_rcs_iron_raw:(in=a2)
			reg.plm_rcs_zinc_raw0(in=b1)
			regb.plm_rcs_zinc_raw:(in=b2)
			reg.plm_rcs_vit_b6_raw0(in=c1)
			regb.plm_rcs_vit_b6_raw:(in=c2)
			reg.plm_rcs_vit_b12_raw0(in=d1)
			regb.plm_rcs_vit_b12_raw:(in=d2)
			;
	*make shorter nutrient identifier (ie, varname);
	length varname $ 32;
		if (a1 or a2) then varname = "iron";
		else if (b1 or b2) then varname ="zinc";
		else if (c1 or c2) then varname = "vit_b6";
		else if (d1 or d2) then varname = "vit_b12";
		
	* remove x values (for which we dont calculate variance estimates); 
	if (replicate>0) and (name = "HEFI2019_TOTAL_SCORE") then delete; 
	run;

	/* note: plm_rcs_hefi_raw = Plot of Linear Model with RCS for HEFI,
		raw association (i.e., unadjusted) */
	
	proc sort data=regb._lin_hefi_plotBS;
	by replicate varname name;
	run;
	
	/* split sample for the <boot_variance> macro */
	data  lin_hefi_plot0 lin_hefi_plot_boot;
		set regb._lin_hefi_plotBS;
	* remove x values;
	if (name = "HEFI2019_TOTAL_SCORE") then delete;
	if replicate =0 then output lin_hefi_plot0;
	else output lin_hefi_plot_boot;
	run;
	
	
/* 2.2) Use <boot_converge.sas> to look a normal-theory bootstrap variance estimation assumption */
%boot_converge(indata    = lin_hefi_plot_boot, /* bootstrap estimates data */
			   checkconv = pred25 pred100 pred175, /* variable for which we want to look at convergence */
			   where_sub = %str(where varname = 'iron'), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );
%boot_converge(indata    = lin_hefi_plot_boot, /* bootstrap estimates data */
			   checkconv = pred25 pred100 pred175, /* variable for which we want to look at convergence */
			   where_sub = %str(where varname = 'zinc'), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );

/* 2.3) Use <boot_variance.sas> macro for efficient variance estimation */
%boot_variance(inboot    = lin_hefi_plot_boot, /* bootstrap estimates */
			   inbase    = lin_hefi_plot0 ,  /* original estimatee */
			   where     = %str(name = "predicted"), /* variance for predicted values only */	   
			   byvar     = varname ,
			   varlist   = pred1-pred200,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.1, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-1),
			   outdata   = lin_plot /* name of output dataset */
			   ); 

	/* 2.3.1) Add X values to the data */
	
	proc transpose data=regb._lin_hefi_plotBS out=_xvalues prefix=HEFI2019_TOTAL_SCORE;
	by varname ;
	var pred1-pred200;
	where name = "HEFI2019_TOTAL_SCORE";
	run;
	
	proc sort data=_xvalues;
	by varname _NAME_ ;
	run;
	
	proc sort data=lin_plot;
	by varname name ;
	run;
	
	data reslib.lin_plot_allf;
	retain name HEFI2019_TOTAL_SCORE estimate ;
		merge _xvalues(rename=(HEFI2019_TOTAL_SCORE1=HEFI2019_TOTAL_SCORE _NAME_=name))
			  lin_plot  ;
		by varname name ;
	* sort from 1 to 200 ;
	pred_id = input(compress(name,,'a'),10.);
	run;
	
	proc sort data=reslib.lin_plot_allf;
	by varname pred_id ;
	run;
	
	/* 2.3.2) Draft plot */
		proc sgplot data=reslib.lin_plot_allf ;
		title1 "Usual iron intake (mg)  | HEFI-2019 total score in adults 65y+ from CCHS 2015 - Nutrition";
		series x=HEFI2019_TOTAL_SCORE y=estimate / lineattrs=(color=black thickness=3);
		series x=HEFI2019_TOTAL_SCORE y=lcl95 / lineattrs=(color=black pattern=longdash thickness=3);
		series x=HEFI2019_TOTAL_SCORE y=ucl95 / lineattrs=(color=black pattern=longdash thickness=3);
		where varname = "iron" ;
		run;
		title1;
		proc sgplot data=reslib.lin_plot_allf ;
		title1 "Usual zinc intake (mcg)  | HEFI-2019 total score in adults 65y+ from CCHS 2015 - Nutrition";
		series x=HEFI2019_TOTAL_SCORE y=estimate / lineattrs=(color=black thickness=3);
		series x=HEFI2019_TOTAL_SCORE y=lcl95 / lineattrs=(color=black pattern=longdash thickness=3);
		series x=HEFI2019_TOTAL_SCORE y=ucl95 / lineattrs=(color=black pattern=longdash thickness=3);
		where varname = "zinc" ;
		run;
		title1;
		proc sgplot data=reslib.lin_plot_allf ;
		title1 "Usual vit_b6 intake (mcg)  | HEFI-2019 total score in adults 65y+ from CCHS 2015 - Nutrition";
		series x=HEFI2019_TOTAL_SCORE y=estimate / lineattrs=(color=black thickness=3);
		series x=HEFI2019_TOTAL_SCORE y=lcl95 / lineattrs=(color=black pattern=longdash thickness=3);
		series x=HEFI2019_TOTAL_SCORE y=ucl95 / lineattrs=(color=black pattern=longdash thickness=3);
		where varname = "vit_b6" ;
		run;
		title1;
		proc sgplot data=reslib.lin_plot_allf ;
		title1 "Usual vit_b12 intake (mcg)  | HEFI-2019 total score in adults 65y+ from CCHS 2015 - Nutrition";
		series x=HEFI2019_TOTAL_SCORE y=estimate / lineattrs=(color=black thickness=3);
		series x=HEFI2019_TOTAL_SCORE y=lcl95 / lineattrs=(color=black pattern=longdash thickness=3);
		series x=HEFI2019_TOTAL_SCORE y=ucl95 / lineattrs=(color=black pattern=longdash thickness=3);
		where varname = "vit_b12" ;
		run;
		title1;

 /*************************************************************************/
 /* ALL: Logistic regression results, model estimates                     */
 /*************************************************************************/

/* 3.1) Output logistic regression estimates data */
	proc sort data=regb._logr_mod_hefi_y_rawBS;
	by replicate varname ;
	run;
	
	/* note: logr_mod_hefi_y = Logistic Regression Model for HEFI with outcome Y, 
		raw association (i.e., unadjusted) */

	data log_hefi_mod0 log_hefi_mod_boot;
		set regb._logr_mod_hefi_y_rawBS ;
	* split bootstrap and original;
		if replicate= 0 then output log_hefi_mod0; 
		else output log_hefi_mod_boot ;
	run;

	/* note: depending on the proportion of outcome, using log_beta or risk_ratio may be more stable.
		Hence, variance of both are estimated. */
	
/* 3.2) Use <boot_converge.sas> to look a normal-theory bootstrap variance estimation assumption */
%boot_converge(indata    = log_hefi_mod_boot, /* bootstrap estimates data */
			   checkconv = log_beta risk_diff risk_ratio, /* variable for which we want to look at convergence */
			   where_sub = %str(where (varname=('binary_iron'))), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );
%boot_converge(indata    = log_hefi_mod_boot, /* bootstrap estimates data */
			   checkconv = log_beta risk_ratio, /* variable for which we want to look at convergence */
			   where_sub = %str(where (varname=('binary_b12'))), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );

/* 3.3) Use <boot_variance.sas> macro for efficient variance estimation */
%boot_variance(inboot    = log_hefi_mod_boot, /* bootstrap estimates */
			   inbase    = log_hefi_mod0 ,  /* original estimatee */
			   byvar     = varname, /* If <by> variable used */
			   varlist   = risk_ratio,
			   null      = 1, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.2, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-2),
			   outdata   = log_beta_all1 /* name of output dataset */
			   ); 
			   
%boot_variance(inboot    = log_hefi_mod_boot, /* bootstrap estimates */
			   inbase    = log_hefi_mod0 ,  /* original estimatee */
			   byvar     = varname, /* If <by> variable used */
			   varlist   = log_beta ,
			   exp       = 1 , /* exponentiate lcl + ucl once estimated */
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.2, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-1),
			   outdata   = log_beta_all2 /* name of output dataset */
			   ); 
			   
%boot_variance(inboot    = log_hefi_mod_boot, /* bootstrap estimates */
			   inbase    = log_hefi_mod0 ,  /* original estimatee */
			   byvar     = varname , /* If <by> variable used */
			   varlist   = risk_diff ,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.2, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-2),
			   outdata   = log_beta_all3 /* name of output dataset */
			   ); 


	data reslib.log_beta_allf(rename=(_LABEL_=label));
		set log_beta_all1 log_beta_all2 log_beta_all3;
	* rename to facilitate interpretation ;
		if name = "log_beta" then name ="odds_ratio";
		if missing(exp) then exp=0;
	* change erroneous label ;
		_LABEL_ = substr(_LABEL_,index(_LABEL_,'|')+2);
		label _LABEL_ = " ";
	run;
	
	proc sort data=reslib.log_beta_allf;
		by varname;
	run;

/* Clean temp */
	proc datasets lib=work nolist nodetails ;
	delete _: log_beta: lin: ;
	run;

 /*************************************************************************/
 /* BY DRI GROUP: Linear regression results, model estimates              */
 /*************************************************************************/

/* 1.1) Output linear regression estimates data */
	proc sort data=reg_b2._rcs_hefi_rawBS;
	by replicate varname drig;
	run;
	
	/* note: rcs_hefi_raw = (Linear) regression with Restricted Cubic Spline transformation
		of the HEFI, raw association (i.e., unadjusted) */

	/* split sample for the <boot_variance> macro */
	data lin_hefi_raw0 lin_hefi_raw_boot ;
		set reg_b2._rcs_hefi_rawBS ;
	if replicate=0 then output lin_hefi_raw0; 
	else output lin_hefi_raw_boot;
	run;

/* 1.2) Use <boot_converge.sas> to look a normal-theory bootstrap variance estimation assumption */
%boot_converge(indata    = lin_hefi_raw_boot, /* bootstrap estimates data */
			   checkconv = rc_beta, /* variable for which we want to look at convergence */
			   where_sub = %str(where (drig=12 and varname='iron')), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );
%boot_converge(indata    = lin_hefi_raw_boot, /* bootstrap estimates data */
			   checkconv = rc_beta, /* variable for which we want to look at convergence */
			   where_sub = %str(where (drig=15 and varname='iron')), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );

/* 1.3) Use <boot_variance.sas> macro for efficient variance estimation */
%boot_variance(inboot    = lin_hefi_raw_boot, /* bootstrap estimates */
			   inbase    = lin_hefi_raw0,  /* original estimatee */
			   byvar     = varname drig, /* If <by> variable used */
			   varlist   = rc_beta ,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.1, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-2), 
			   outdata   = reslib.lin_beta_drigf /* name of output dataset */
			   ); 
	
	title2 "Nutrient intake difference between respondents with high (90th) vs. average (50th) HEFI-2019 score";
	proc print data=reslib.lin_beta_drigf label;
	id varname drig name nboot ;
	var  estimate se lcl95 ucl95 pvalue svalue cv;
	run;
	title2;

 /*************************************************************************/
 /* BY DRI GROUP: Linear regression results, predicted values for plot    */
 /*************************************************************************/

/* 2.1) Output linear regression plot data */
	proc sort data=reg_b2._lin_hefi_plotBS;
	by replicate varname drig name;
	run;
	
	/* note: plm_rcs_hefi_raw = Plot of Linear Model with RCS for HEFI,
		raw association (i.e., unadjusted) */
		
	data  lin_hefi_plot0 lin_hefi_plot_boot;
		set reg_b2._lin_hefi_plotBS;
	* remove x values;
	if (name = "HEFI2019_TOTAL_SCORE") then delete;
	if replicate =0 then output lin_hefi_plot0;
	else output lin_hefi_plot_boot;
	run;
	
/* 2.2) Use <boot_converge.sas> to look a normal-theory bootstrap variance estimation assumption */
%boot_converge(indata    = lin_hefi_plot_boot, /* bootstrap estimates data */
			   checkconv = pred25 pred100 pred175, /* variable for which we want to look at convergence */
			   where_sub = %str(where drig='12' and varname='iron'), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );
%boot_converge(indata    = lin_hefi_plot_boot, /* bootstrap estimates data */
			   checkconv = pred25 pred100 pred175, /* variable for which we want to look at convergence */
			   where_sub = %str(where drig='15' and varname='iron'), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );
			   
/* 2.3) Use <boot_variance.sas> macro for efficient variance estimation */
%boot_variance(inboot    = lin_hefi_plot_boot, /* bootstrap estimates */
			   inbase    = lin_hefi_plot0 ,  /* original estimatee */
			   where     = %str(name = "predicted"), /* variance for predicted values only */
			   byvar     = varname drig,
			   varlist   = pred1-pred200,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.1, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-1),
			   outdata   = lin_plot /* name of output dataset */
			   ); 

	/* 2.3.1) Add X values to the data */
		
	data _xvalues;
		set reg_2.plmrcst_iron_raw0(in=a)
			reg_2.plmrcst_zinc_raw0 (in=b)
			reg_2.plmrcst_vit_b6_raw0(in=c)
			reg_2.plmrcst_vit_b12_raw0 (in=d)
			;
	length varname $ 32;
	if a then varname = "iron";
	else if b then varname = "zinc";
	else if c then varname = "vit_b6";
	else if d then varname = "vit_b12";
	run;
		
		/* some manipulation to output <name> variable*/
		proc sort;
			by varname drig;
		run;
		
		data _xvalues;
		set _xvalues;
		by varname drig;
		retain pred_id ;
		if first.drig then pred_id=1;
		else pred_id+1;
		length name $  8;
		name = cats("pred",pred_id);
		run;
		
	proc sort data=_xvalues;
	by varname drig name ;
	run;
	
	proc sort data=lin_plot;
	by varname drig name ;
	run;
	
	data reslib.lin_plot_drigf(rename=(_drig=drig));
	retain varname _drig pred_id name HEFI2019_TOTAL_SCORE estimate ;
		merge _xvalues
			  lin_plot  ;
		by varname drig name ;
	* change drig to num ;
		_drig = input(drig, 10.);
	* add merge error flag for extra safety ;
		if predicted ne estimate then _MergeError_=1;
		else _MergeError_=0;
		label _MergeError_ = "Flag for any merge error";
	*clean; 
		drop drig predicted;
	run;
	
	proc sort data=reslib.lin_plot_drigf;
	by varname drig pred_id ;
	run;	

	
	/* 2.3.2) Draft plot */
		proc sgplot data=reslib.lin_plot_drigf ;
		title1 "Usual iron (mg)  | HEFI-2019 total score in adults 65y+ from CCHS 2015 - Nutrition, by DRIG";
		series x=HEFI2019_TOTAL_SCORE y=estimate / lineattrs=( thickness=3) group=drig;
		*series x=HEFI2019_TOTAL_SCORE y=lcl95 / lineattrs=(color=black pattern=longdash thickness=3) ;
		*series x=HEFI2019_TOTAL_SCORE y=ucl95 / lineattrs=(color=black pattern=longdash thickness=3) ;
		where varname = "iron";
		run;
		title1;
		proc sgplot data=reslib.lin_plot_drigf ;
		title1 "Usual zinc (mcg)  | HEFI-2019 total score in adults 65y+ from CCHS 2015 - Nutrition, by DRIG";
		series x=HEFI2019_TOTAL_SCORE y=estimate / lineattrs=( thickness=3) group=drig;
		*series x=HEFI2019_TOTAL_SCORE y=lcl95 / lineattrs=(color=black pattern=longdash thickness=3) ;
		*series x=HEFI2019_TOTAL_SCORE y=ucl95 / lineattrs=(color=black pattern=longdash thickness=3) ;
		where varname = "zinc";
		run;
		title1;
		proc sgplot data=reslib.lin_plot_drigf ;
		title1 "Usual vit_b6  | HEFI-2019 total score in adults 65y+ from CCHS 2015 - Nutrition, by DRIG";
		series x=HEFI2019_TOTAL_SCORE y=estimate / lineattrs=( thickness=3) group=drig;
		*series x=HEFI2019_TOTAL_SCORE y=lcl95 / lineattrs=(color=black pattern=longdash thickness=3) ;
		*series x=HEFI2019_TOTAL_SCORE y=ucl95 / lineattrs=(color=black pattern=longdash thickness=3) ;
		where varname = "vit_b6";
		run;
		title1;
		proc sgplot data=reslib.lin_plot_drigf ;
		title1 "Usual vit_b12  | HEFI-2019 total score in adults 65y+ from CCHS 2015 - Nutrition, by DRIG";
		series x=HEFI2019_TOTAL_SCORE y=estimate / lineattrs=( thickness=3) group=drig;
		*series x=HEFI2019_TOTAL_SCORE y=lcl95 / lineattrs=(color=black pattern=longdash thickness=3) ;
		*series x=HEFI2019_TOTAL_SCORE y=ucl95 / lineattrs=(color=black pattern=longdash thickness=3) ;
		where varname = "vit_b12";
		run;
		title1;

 /*************************************************************************/
 /* BY DRI GROUP: Logistic regression results, model estimates            */
 /*************************************************************************/

/* 3.1) Output logistic regression estimates data */
	proc sort data=reg_b2._logr_mod_hefi_y_rawBS;
	by replicate varname drig;
	run;
	
	/* note: logr_mod_hefi_y = Logistic Regression Model for HEFI with outcome Y, 
		raw association (i.e., unadjusted) */

	data log_hefi_mod0 log_hefi_mod_boot;
		set reg_b2._logr_mod_hefi_y_rawBS ;
	* split bootstrap and original;
		if replicate= 0 then output log_hefi_mod0; 
		else output log_hefi_mod_boot ;
	run;

/* 3.2) Use <boot_converge.sas> to look a normal-theory bootstrap variance estimation assumption */
%boot_converge(indata    = log_hefi_mod_boot, /* bootstrap estimates data */
			   checkconv = exp_beta risk_diff risk_ratio, /* variable for which we want to look at convergence */
			   where_sub = %str(where (drig=12 & varname='binary_iron')), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );

/* 3.3) Use <boot_variance.sas> macro for efficient variance estimation */
%boot_variance(inboot    = log_hefi_mod_boot, /* bootstrap estimates */
			   inbase    = log_hefi_mod0 ,  /* original estimatee */
			   byvar     = varname drig, /* If <by> variable used */
			   varlist   = risk_ratio,
			   null      = 1, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.2, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-2),
			   outdata   = log_beta_all1 /* name of output dataset */
			   ); 

%boot_variance(inboot    = log_hefi_mod_boot, /* bootstrap estimates */
			   inbase    = log_hefi_mod0 ,  /* original estimatee */
			   byvar     = varname drig, /* If <by> variable used */
			   varlist   = log_beta ,
			   exp       = 1,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.2, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-1),
			   outdata   = log_beta_all2 /* name of output dataset */
			   ); 

%boot_variance(inboot    = log_hefi_mod_boot, /* bootstrap estimates */
			   inbase    = log_hefi_mod0 ,  /* original estimatee */
			   byvar     = varname drig, /* If <by> variable used */
			   varlist   = risk_diff ,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.2, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-2),
			   outdata   = log_beta_all3 /* name of output dataset */
			   ); 
			   		   
	data reslib.log_beta_drigf(rename=(_LABEL_=label));
		set log_beta_all1 log_beta_all2 log_beta_all3;
	* rename to facilitate interpretation ;
		if name = "log_beta" then name ="odds_ratio";
		if missing(exp) then exp=0;
	* change erroneous label ;
		_LABEL_ = substr(_LABEL_,index(_LABEL_,'|')+2);
		label _LABEL_ = " ";
	run;
	
	proc sort data=reslib.log_beta_drigf;
		by varname drig;
	run;

/* Clean temp */
	proc datasets lib=work nolist nodetails ;
	delete _: log_beta: log_hefi: lin_: ;
	run;

 /*************************************************************************/
 /* ALL: Distribution of usual intake, marginal Pr(X<x)                   */
 /*************************************************************************/

/* 1.1) Output marginal Y distribution data */
	data bootlib._distrib_y_wBS;
		set reslib.distrib_y_iron_w0
			reslib.distrib_y_zinc_w0
			reslib.distrib_y_b6_w0
			reslib.distrib_y_b12_w0
			bootlib.distrib_y_: ;
	* For PREVALENCE: rescale prob to percentage point ;
	array prob(*) binary_iron binary_zinc binary_b6 binary_b12;
	do i=1 to dim(prob);
	prob(i) = prob(i) * 100 ;
	end;
	drop i;
	* Combine both variables for easy variance estimation;
	if varname = "iron" then binary_EAR = binary_iron;
	else if varname = "zinc" then binary_EAR = binary_zinc;
	else if varname = "vit_b6" then binary_EAR = binary_b6;
	else if varname = "vit_b12" then binary_EAR = binary_b12;
	label binary_EAR = "Pr(X<EAR)";
	run;
	
	/* note: distrib_y_rel_w_boot = Distribution of Y, Relative intakes, Wide format */

	proc sort data=bootlib._distrib_y_wBS;
	by replicate varname DRIg;
	run;

	/* 1.2) split sample for the <boot_variance> macro */
	data distrib_y_w0 distrib_y_w_boot ;
		set bootlib._distrib_y_wBS;
	if replicate= 0 then output distrib_y_w0;
	else output distrib_y_w_boot;
	run;
	

/* 1.2) Use <boot_converge.sas> to look a normal-theory bootstrap variance estimation assumption */
%boot_converge(indata    = distrib_y_w_boot, /* bootstrap estimates data */
			   checkconv = Mean Pctile10 Pctile90, /* variable for which we want to look at convergence */
			   where_sub = %str(where (drig=0 and varname='iron')), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );
%boot_converge(indata    = distrib_y_w_boot, /* bootstrap estimates data */
			   checkconv = Mean Pctile10 Pctile90, /* variable for which we want to look at convergence */
			   where_sub = %str(where (drig=0 and varname='vit_b12')), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );

	
/* 1.3) Use <boot_variance.sas> macro for efficient variance estimation */

/* 1.3.2) actual variance estimation */
%boot_variance(inboot    = distrib_y_w_boot, /* bootstrap estimates */
			   inbase    = distrib_y_w0,  /* original estimatee */
			   byvar     = varname drig, /* If <by> variable used */
			   varlist   = Mean &pctile_list ,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 4.2, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-1), 
			   outdata   = distrib_y_wf1 /* name of output dataset */
			   );

%boot_variance(inboot    = distrib_y_w_boot, /* bootstrap estimates */
			   inbase    = distrib_y_w0,  /* original estimatee */
			   byvar     = varname drig, /* If <by> variable used */
			   varlist   = binary_ear ,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 4.1, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-1), 
			   outdata   = distrib_y_wf2 /* name of output dataset */
			   ); 
			  
	data reslib.distrib_y_wf (rename=(_LABEL_ = label));
	retain drig name _LABEL_ ;
	length name $ 32;
		set distrib_y_wf1(in=a) distrib_y_wf2(in=b) ;
	label _LABEL_ = " ";
	run;
	
	proc sort;
	by varname drig;
	run;

	proc print data=reslib.distrib_y_wf label;
	id varname name label nboot drig;
	var  estimate se lcl95 ucl95 pvalue svalue cv;
	run;

 /*************************************************************************/
 /* BY CATEGORY OF SCORE: Distribution of usual intake                    */
 /*************************************************************************/

/* 2.1) Append marginal Y distribution data */
	data bootlib._distrib_catscore_wBS;
		set reslib.distrib_w_catscore0
			bootlib.distrib_w_catscore: 
			;
	* PREVALENCE: rescale prob to percentage point;
	if catscore < 9000 then binary_ear = binary_ear*100;
	run;
	
	/* note: distrib_y_rel_w_boot = Distribution of Y, Relative intakes, Wide format */

	proc sort data=bootlib._distrib_catscore_wBS;
	by replicate varname catscore;
	run;
	
	/* 2.1.2) split sample for the <boot_variance> macro */
	data distrib_catscore0 distrib_catscore_boot ;
		set bootlib._distrib_catscore_wBS;
	if replicate= 0 then output distrib_catscore0;
	else output distrib_catscore_boot;
	run;


/* 2.2) Use <boot_converge.sas> to look a normal-theory bootstrap variance estimation assumption */
%boot_converge(indata    = distrib_catscore_boot, /* bootstrap estimates data */
			   checkconv = Mean Pctile10 Pctile90, /* variable for which we want to look at convergence */
			   where_sub = %str(where (catscore=4 and varname='iron')), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );

	
/* 2.3) Use <boot_variance.sas> macro for efficient variance estimation */

	/* 2.3.2) Variance estimation */
	%boot_variance(inboot    = distrib_catscore_boot, /* bootstrap estimates */
				   inbase    = distrib_catscore0,  /* original estimatee */
				   byvar     = varname catscore, /* If <by> variable used */
				   varlist   = Mean &pctile_list ,
				   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
				   set_confidence_level=95, /* 95 = 95%CI */
				   numberfmt = 5.2, /* formatting of estimates */
				   set_degrees_freedom = %str(nboot-1), 
				   outdata   = distrib_catscoref1 /* name of output dataset */
				   ); 
	
	%boot_variance(inboot    = distrib_catscore_boot, /* bootstrap estimates */
				   inbase    = distrib_catscore0,  /* original estimatee */
				   byvar     = varname catscore, /* If <by> variable used */
				   varlist   = binary_ear ,
				   null      = 1, /* Null H for pvalues (i.e., 0 for a difference) */
				   set_confidence_level=95, /* 95 = 95%CI */
				   numberfmt = 4.2, /* formatting of estimates */
				   set_degrees_freedom = %str(nboot-2), 
				   outdata   = distrib_catscoref2 /* name of output dataset */
				   ); 
				   
 	/* 2.3.3) Check and save output */
	data reslib.distrib_w_catscoref(rename=(_LABEL_=label));
	retain name _LABEL_ catscore;
	length name $ 32;
		set distrib_catscoref1 distrib_catscoref2 ;
	label _LABEL_ = " ";
	run;
	
	proc sort;
	by varname catscore;
	run;
	
	title1 "Nutrient intake by quarter of total HEFI-2019 score";
	proc print data=reslib.distrib_w_catscoref label;
	id varname name nboot catscore;
	var  estimate se lcl95 ucl95 pvalue svalue cv;
	where index(name,"Prob")=0;
	run;

	title2 "Pr(X<x) by quarter of total HEFI-2019 score";
	proc print data=reslib.distrib_w_catscoref label;
	id varname name label nboot catscore;
	var  estimate se lcl95 ucl95 pvalue svalue cv;
	where index(name,"Prob")>0;
	run;
	title1;
	title2;

/* Clean temp */
	proc datasets lib=work nolist nodetails ;
	delete _: distrib_: ;
	run;


 /*************************************************************************/
 /*                                                                       */
 /*     REF. CODE 5: Variance estimation for MISC.B NUTRIENT outcome      */
 /*                                                                       */
 /*************************************************************************/

data _null_;
	file print;
	put #3 @10 "Bootstrap variance estimation for miscellaneous B (&SYSDATE)";
run;

/* nutrient intakes: dfe mg fibers pota */
/* binary outcomes: binary_dfe binary_mg binary_fib binary_pota */

/* note: convergence plot revealed extreme/implausible observations for replicate 484,
	suggesting failed run. Replicate removed. */

 /************************************************/
 /*            Files and folder set-up           */
 /************************************************/

/* 0) define suffix of current outcome */
	
	%let suffix = _miscB ;

/* 1) assign proper library using the same suffix */
	%mcmclib(&suffix.);
	%MakeMCMCFolder(NCIPath = &path./NCI/MCMC&suffix.,
					prefix  = Reg_,
					suffix  = &subgrp.
					);

 /*************************************************************************/
 /* ALL: Linear regression results, model estimates                       */
 /*************************************************************************/

/* 1.1) Output linear regression estimates data */	
	proc sort data=regb._rcs_hefi_rawBS;
	by replicate;
	run;

	/* note: rcs_hefi_raw = (Linear) regression with Restricted Cubic Spline transformation
		of the HEFI, raw association (i.e., unadjusted) */

	/* check extreme replicate */
	ods select moments quantiles extremeobs ;
	title1 color=red "Extreme replicates for Misc. B";
	proc univariate data=regb._rcs_hefi_rawBS;
	class varname;
	var rc_beta ;
	id replicate;
	where replicate ne 0 ;
	run;
	title1;

	/* split sample for the <boot_variance> macro */
	data lin_hefi_raw0 lin_hefi_raw_boot ;
		set regb._rcs_hefi_rawBS ;
	* remove replicate with extreme values ;
	if replicate = 484 then delete;
	if replicate=0 then output lin_hefi_raw0; 
	else output lin_hefi_raw_boot;
	run;

/* 1.2) Use <boot_converge.sas> to look a normal-theory bootstrap variance estimation assumption */
%boot_converge(indata    = regb._rcs_hefi_rawBS, /* bootstrap estimates data */
			   checkconv = rc_beta, /* variable for which we want to look at convergence */
			   where_sub = %str(where varname = 'dfe' & replicate ne 484), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );
%boot_converge(indata    = regb._rcs_hefi_rawBS, /* bootstrap estimates data */
			   checkconv = rc_beta, /* variable for which we want to look at convergence */
			   where_sub = %str(where varname = 'pota' & replicate ne 484 ), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );
	%macro notrun; /* extra checks */
	%boot_converge(indata    = regb._rcs_hefi_rawBS, /* bootstrap estimates data */
				   checkconv = rc_beta, /* variable for which we want to look at convergence */
				   where_sub = %str(where varname = 'mg' & replicate ne 484 ), /* if a subgroup is used - NA for this example */
				   sefmt     = 5.2, /* formatting for plot */
				   normality = 1 /* if 1, checks normality */
				   );
	%boot_converge(indata    = regb._rcs_hefi_rawBS, /* bootstrap estimates data */
				   checkconv = rc_beta, /* variable for which we want to look at convergence */
				   where_sub = %str(where varname = 'fibers' & replicate ne 484 ), /* if a subgroup is used - NA for this example */
				   sefmt     = 5.2, /* formatting for plot */
				   normality = 1 /* if 1, checks normality */
				   );
	%mend notrun;

/* 1.3) Use <boot_variance.sas> macro for efficient variance estimation */
%boot_variance(inboot    = lin_hefi_raw_boot, /* bootstrap estimates */
			   inbase    = lin_hefi_raw0,  /* original estimatee */
			   byvar     = varname ,
			   varlist   = rc_beta r2,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.1, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-2), 
			   outdata   = reslib.lin_beta_allf /* name of output dataset */
			   ); 

	proc print data=reslib.lin_beta_allf label;
	id varname name nboot;
	var  estimate se lcl95 ucl95 pvalue svalue cv;
	run;

 /*************************************************************************/
 /* ALL: Linear regression results, predicted values for plot             */
 /*************************************************************************/

/* 2.1) Output linear regression plot data */
	proc sort data=regb._lin_hefi_plotBS;
	by replicate varname name;
	run;

	/* note: plm_rcs_hefi_raw = Plot of Linear Model with RCS for HEFI,
		raw association (i.e., unadjusted) */

	/* check extreme replicate */
		ods select moments quantiles extremeobs ;
		title1 color=red "Extreme replicates for Misc. B";
		proc univariate data=regb._lin_hefi_plotBS;
		class varname;
		var pred100 ;
		id replicate;
		where replicate ne 0 ;
		run;
		title1;

	/* split sample for the <boot_variance> macro */
	data  lin_hefi_plot0 lin_hefi_plot_boot;
		set regb._lin_hefi_plotBS;
	* remove replicate with extreme values ;
		if replicate = 484 then delete;
	* remove x values;
		if (name = "HEFI2019_TOTAL_SCORE") then delete;
	if replicate =0 then output lin_hefi_plot0;
	else output lin_hefi_plot_boot;
	run;
	
/* 2.2) Use <boot_converge.sas> to look a normal-theory bootstrap variance estimation assumption */
%boot_converge(indata    = lin_hefi_plot_boot, /* bootstrap estimates data */
			   checkconv = pred25 pred100 pred175, /* variable for which we want to look at convergence */
			   where_sub = %str(where varname = 'dfe' & replicate ne 484), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );
%boot_converge(indata    = lin_hefi_plot_boot, /* bootstrap estimates data */
			   checkconv = pred25 pred100 pred175, /* variable for which we want to look at convergence */
			   where_sub = %str(where varname = 'mg' & replicate ne 484), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );

/* 2.3) Use <boot_variance.sas> macro for efficient variance estimation */
%boot_variance(inboot    = lin_hefi_plot_boot, /* bootstrap estimates */
			   inbase    = lin_hefi_plot0 ,  /* original estimatee */
			   where     = %str(name = "predicted"), /* variance for predicted values only */	   
			   byvar     = varname ,
			   varlist   = pred1-pred200,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.1, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-1),
			   outdata   = lin_plot /* name of output dataset */
			   ); 

	/* 2.3.1) Add X values to the data */
	
	proc transpose data=regb._lin_hefi_plotBS out=_xvalues prefix=HEFI2019_TOTAL_SCORE;
	by varname ;
	var pred1-pred200;
	where name = "HEFI2019_TOTAL_SCORE";
	run;
	
	proc sort data=_xvalues;
	by varname _NAME_ ;
	run;
	
	proc sort data=lin_plot;
	by varname name ;
	run;
	
	data reslib.lin_plot_allf;
	retain name HEFI2019_TOTAL_SCORE estimate ;
		merge _xvalues(rename=(HEFI2019_TOTAL_SCORE1=HEFI2019_TOTAL_SCORE _NAME_=name))
			  lin_plot  ;
		by varname name ;
	* sort from 1 to 200 ;
	pred_id = input(compress(name,,'a'),10.);
	run;
	
	proc sort data=reslib.lin_plot_allf;
	by varname pred_id ;
	run;
	
	/* 2.3.2) Draft plot */
		proc sgplot data=reslib.lin_plot_allf ;
		title1 "Usual dfe intake (mg)  | HEFI-2019 total score in adults 65y+ from CCHS 2015 - Nutrition";
		series x=HEFI2019_TOTAL_SCORE y=estimate / lineattrs=(color=black thickness=3);
		series x=HEFI2019_TOTAL_SCORE y=lcl95 / lineattrs=(color=black pattern=longdash thickness=3);
		series x=HEFI2019_TOTAL_SCORE y=ucl95 / lineattrs=(color=black pattern=longdash thickness=3);
		where varname = "dfe" ;
		run;
		title1;
		proc sgplot data=reslib.lin_plot_allf ;
		title1 "Usual mg intake (mcg)  | HEFI-2019 total score in adults 65y+ from CCHS 2015 - Nutrition";
		series x=HEFI2019_TOTAL_SCORE y=estimate / lineattrs=(color=black thickness=3);
		series x=HEFI2019_TOTAL_SCORE y=lcl95 / lineattrs=(color=black pattern=longdash thickness=3);
		series x=HEFI2019_TOTAL_SCORE y=ucl95 / lineattrs=(color=black pattern=longdash thickness=3);
		where varname = "mg" ;
		run;
		title1;
		proc sgplot data=reslib.lin_plot_allf ;
		title1 "Usual fib intake (mcg)  | HEFI-2019 total score in adults 65y+ from CCHS 2015 - Nutrition";
		series x=HEFI2019_TOTAL_SCORE y=estimate / lineattrs=(color=black thickness=3);
		series x=HEFI2019_TOTAL_SCORE y=lcl95 / lineattrs=(color=black pattern=longdash thickness=3);
		series x=HEFI2019_TOTAL_SCORE y=ucl95 / lineattrs=(color=black pattern=longdash thickness=3);
		where varname = "fibers" ;
		run;
		title1;
		proc sgplot data=reslib.lin_plot_allf ;
		title1 "Usual potassium intake (mcg)  | HEFI-2019 total score in adults 65y+ from CCHS 2015 - Nutrition";
		series x=HEFI2019_TOTAL_SCORE y=estimate / lineattrs=(color=black thickness=3);
		series x=HEFI2019_TOTAL_SCORE y=lcl95 / lineattrs=(color=black pattern=longdash thickness=3);
		series x=HEFI2019_TOTAL_SCORE y=ucl95 / lineattrs=(color=black pattern=longdash thickness=3);
		where varname = "pota" ;
		run;
		title1;

 /*************************************************************************/
 /* ALL: Logistic regression results, model estimates                     */
 /*************************************************************************/

/* 3.1) Output logistic regression estimates data */
	proc sort data=regb._logr_mod_hefi_y_rawBS;
	by replicate varname ;
	run;

	/* note: logr_mod_hefi_y = Logistic Regression Model for HEFI with outcome Y, 
		raw association (i.e., unadjusted) */

	/* check extreme replicate */
		ods select moments quantiles extremeobs ;
		title1 color=red "Extreme replicates for Misc. B";
		proc univariate data=regb._logr_mod_hefi_y_rawBS;
		class varname;
		var risk_diff ;
		id replicate;
		where replicate ne 0 ;
		run;
		title1;

	data log_hefi_mod0 log_hefi_mod_boot;
		set regb._logr_mod_hefi_y_rawBS ;
	* split bootstrap and original;
		if replicate= 0 then output log_hefi_mod0; 
		else output log_hefi_mod_boot ;
	run;

	/* note: depending on the proportion of outcome, using log_beta or risk_ratio may be more stable.
		Hence, variance of both are estimated. */
	
/* 3.2) Use <boot_converge.sas> to look a normal-theory bootstrap variance estimation assumption */
%boot_converge(indata    = log_hefi_mod_boot, /* bootstrap estimates data */
			   checkconv = log_beta risk_ratio, /* variable for which we want to look at convergence */
			   where_sub = %str(where (varname=('binary_dfe' ))), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );
%boot_converge(indata    = log_hefi_mod_boot, /* bootstrap estimates data */
			   checkconv = log_beta risk_ratio, /* variable for which we want to look at convergence */
			   where_sub = %str(where (varname=('binary_pota' ))), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );
%boot_converge(indata    = log_hefi_mod_boot, /* bootstrap estimates data */
			   checkconv = log_beta risk_ratio, /* variable for which we want to look at convergence */
			   where_sub = %str(where (varname=('binary_fib'))), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );
%boot_converge(indata    = log_hefi_mod_boot, /* bootstrap estimates data */
			   checkconv = log_beta risk_ratio, /* variable for which we want to look at convergence */
			   where_sub = %str(where (varname=('binary_mg' ))), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );

/* 3.3) Use <boot_variance.sas> macro for efficient variance estimation */
%boot_variance(inboot    = log_hefi_mod_boot, /* bootstrap estimates */
			   inbase    = log_hefi_mod0 ,  /* original estimatee */
			   byvar     = varname, /* If <by> variable used */
			   varlist   = risk_ratio,
			   null      = 1, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.2, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-2),
			   outdata   = log_beta_all1 /* name of output dataset */
			   ); 
			   
%boot_variance(inboot    = log_hefi_mod_boot, /* bootstrap estimates */
			   inbase    = log_hefi_mod0 ,  /* original estimatee */
			   byvar     = varname, /* If <by> variable used */
			   varlist   = log_beta ,
			   exp       = 1 , /* exponentiate lcl + ucl once estimated */
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.2, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-1),
			   outdata   = log_beta_all2 /* name of output dataset */
			   ); 
			   
%boot_variance(inboot    = log_hefi_mod_boot, /* bootstrap estimates */
			   inbase    = log_hefi_mod0 ,  /* original estimatee */
			   byvar     = varname, /* If <by> variable used */
			   varlist   = risk_diff ,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.2, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-2),
			   outdata   = log_beta_all3 /* name of output dataset */
			   ); 	   

	data reslib.log_beta_allf(rename=(_LABEL_=label));
		set log_beta_all1 log_beta_all2 log_beta_all3;
	* rename to facilitate interpretation ;
		if name = "log_beta" then name ="odds_ratio";
		if missing(exp) then exp=0;
	* change erroneous label ;
		_LABEL_ = substr(_LABEL_,index(_LABEL_,'|')+2);
		label _LABEL_ = " ";
	run;
	
	proc sort data=reslib.log_beta_allf;
		by varname;
	run;

/* Clean temp */
	proc datasets lib=work nolist nodetails ;
	delete _: lin: log_beta: log_hefi: ;
	run;

 /*************************************************************************/
 /* BY DRI GROUP: Linear regression results, model estimates              */
 /*************************************************************************/

/* 1.1) Output linear regression estimates data */
	proc sort data=reg_b2._rcs_hefi_rawBS;
	by replicate varname drig;
	run;

	/* note: rcs_hefi_raw = (Linear) regression with Restricted Cubic Spline transformation
		of the HEFI, raw association (i.e., unadjusted) */


	/* check extreme replicate */
		ods select moments quantiles extremeobs ;
		title1 color=red "Extreme replicates for Misc. B";
		proc univariate data=reg_b2._rcs_hefi_rawBS;
		class drig varname;
		var rc_beta ;
		id replicate;
		where replicate ne 0 ;
		run;
		title1;

	/* split sample for the <boot_variance> macro */
	data lin_hefi_raw0 lin_hefi_raw_boot ;
		set reg_b2._rcs_hefi_rawBS ;
	* remove replicate with extreme values ;
		if replicate = 484 then delete;
	if replicate=0 then output lin_hefi_raw0; 
	else output lin_hefi_raw_boot;
	run;

/* 1.2) Use <boot_converge.sas> to look a normal-theory bootstrap variance estimation assumption */
%boot_converge(indata    = lin_hefi_raw_boot, /* bootstrap estimates data */
			   checkconv = rc_beta, /* variable for which we want to look at convergence */
			   where_sub = %str(where (drig=12 and varname='dfe' & replicate ne 484)), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );
%boot_converge(indata    = lin_hefi_raw_boot, /* bootstrap estimates data */
			   checkconv = rc_beta, /* variable for which we want to look at convergence */
			   where_sub = %str(where (drig=15 and varname='dfe' & replicate ne 484)), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );

/* 1.3) Use <boot_variance.sas> macro for efficient variance estimation */
%boot_variance(inboot    = lin_hefi_raw_boot, /* bootstrap estimates */
			   inbase    = lin_hefi_raw0,  /* original estimatee */
			   byvar     = varname drig, /* If <by> variable used */
			   varlist   = rc_beta ,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.1, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-2), 
			   outdata   = reslib.lin_beta_drigf /* name of output dataset */
			   ); 
	
	title2 "Nutrient intake difference between respondents with high (90th) vs. average (50th) HEFI-2019 score";
	proc print data=reslib.lin_beta_drigf label;
	id varname drig name nboot ;
	var  estimate se lcl95 ucl95 pvalue svalue cv;
	run;
	title2;

 /*************************************************************************/
 /* BY DRI GROUP: Linear regression results, predicted values for plot    */
 /*************************************************************************/

/* 2.1) Output linear regression plot data */
	proc sort data=reg_b2._lin_hefi_plotBS;
	by replicate varname drig name;
	run;
	
	/* note: plm_rcs_hefi_raw = Plot of Linear Model with RCS for HEFI,
		raw association (i.e., unadjusted) */

	/* check extreme replicate */
		ods select moments quantiles extremeobs ;
		title1 color=red "Extreme replicates for Misc. B";
		proc univariate data=reg_b2._lin_hefi_plotBS;
		class drig varname;
		var pred100 ;
		id replicate;
		where replicate ne 0 ;
		run;
		title1;	

	/* split sample for the <boot_variance> macro */
	data  lin_hefi_plot0 lin_hefi_plot_boot;
		set reg_b2._lin_hefi_plotBS;
	* remove replicate with extreme values ;
		if replicate = 484 then delete;
	* remove x values;
		if (name = "HEFI2019_TOTAL_SCORE") then delete;
	if replicate =0 then output lin_hefi_plot0;
	else output lin_hefi_plot_boot;
	run;
	
/* 2.2) Use <boot_converge.sas> to look a normal-theory bootstrap variance estimation assumption */
%boot_converge(indata    = lin_hefi_plot_boot, /* bootstrap estimates data */
			   checkconv = pred25 pred100 pred175, /* variable for which we want to look at convergence */
			   where_sub = %str(where drig='12' and varname='dfe' & replicate ne 484), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );
%boot_converge(indata    = lin_hefi_plot_boot, /* bootstrap estimates data */
			   checkconv = pred25 pred100 pred175, /* variable for which we want to look at convergence */
			   where_sub = %str(where drig='15' and varname='dfe' & replicate ne 484), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );
			   
/* 2.3) Use <boot_variance.sas> macro for efficient variance estimation */
%boot_variance(inboot    = lin_hefi_plot_boot, /* bootstrap estimates */
			   inbase    = lin_hefi_plot0 ,  /* original estimatee */
			   where     = %str(name = "predicted"), /* variance for predicted values only */
			   byvar     = varname drig,
			   varlist   = pred1-pred200,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.1, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-1),
			   outdata   = lin_plot /* name of output dataset */
			   ); 

	/* 2.3.1) Add X values to the data */
		
	data _xvalues;
		set reg_2.plmrcst_dfe_raw0(in=a)
			reg_2.plmrcst_mg_raw0 (in=b)
			reg_2.plmrcst_fibers_raw0(in=c)
			reg_2.plmrcst_pota_raw0 (in=d)
			;
	length varname $ 32;
	if a then varname = "dfe";
	else if b then varname = "mg";
	else if c then varname = "fibers";
	else if d then varname = "pota";
	run;
		
		/* some manipulation to output <name> variable*/
		proc sort;
			by varname drig;
		run;
		
		data _xvalues;
		set _xvalues;
		by varname drig;
		retain pred_id ;
		if first.drig then pred_id=1;
		else pred_id+1;
		length name $  8;
		name = cats("pred",pred_id);
		run;
		
	proc sort data=_xvalues;
	by varname drig name ;
	run;
	
	proc sort data=lin_plot;
	by varname drig name ;
	run;
	
	data reslib.lin_plot_drigf(rename=(_drig=drig));
	retain varname _drig pred_id name HEFI2019_TOTAL_SCORE estimate ;
		merge _xvalues
			  lin_plot  ;
		by varname drig name ;
	* change drig to num ;
		_drig = input(drig, 10.);
	* add merge error flag for extra safety ;
		if predicted ne estimate then _MergeError_=1;
		else _MergeError_=0;
		label _MergeError_ = "Flag for any merge error";
	*clean; 
		drop drig predicted;
	run;
	
	proc sort data=reslib.lin_plot_drigf;
	by varname drig pred_id ;
	run;	

	
	/* 2.3.2) Draft plot */
		proc sgplot data=reslib.lin_plot_drigf ;
		title1 "Usual dfe | HEFI-2019 total score in adults 65y+ from CCHS 2015 - Nutrition, by DRIG";
		series x=HEFI2019_TOTAL_SCORE y=estimate / lineattrs=( thickness=3) group=drig;
		*series x=HEFI2019_TOTAL_SCORE y=lcl95 / lineattrs=(color=black pattern=longdash thickness=3) ;
		*series x=HEFI2019_TOTAL_SCORE y=ucl95 / lineattrs=(color=black pattern=longdash thickness=3) ;
		where varname = "dfe";
		run;
		title1;
		proc sgplot data=reslib.lin_plot_drigf ;
		title1 "Usual mg  | HEFI-2019 total score in adults 65y+ from CCHS 2015 - Nutrition, by DRIG";
		series x=HEFI2019_TOTAL_SCORE y=estimate / lineattrs=( thickness=3) group=drig;
		*series x=HEFI2019_TOTAL_SCORE y=lcl95 / lineattrs=(color=black pattern=longdash thickness=3) ;
		*series x=HEFI2019_TOTAL_SCORE y=ucl95 / lineattrs=(color=black pattern=longdash thickness=3) ;
		where varname = "mg";
		run;
		title1;
		proc sgplot data=reslib.lin_plot_drigf ;
		title1 "Usual fibers | HEFI-2019 total score in adults 65y+ from CCHS 2015 - Nutrition, by DRIG";
		series x=HEFI2019_TOTAL_SCORE y=estimate / lineattrs=( thickness=3) group=drig;
		*series x=HEFI2019_TOTAL_SCORE y=lcl95 / lineattrs=(color=black pattern=longdash thickness=3) ;
		*series x=HEFI2019_TOTAL_SCORE y=ucl95 / lineattrs=(color=black pattern=longdash thickness=3) ;
		where varname = "fibers";
		run;
		title1;
		proc sgplot data=reslib.lin_plot_drigf ;
		title1 "Usual pota | HEFI-2019 total score in adults 65y+ from CCHS 2015 - Nutrition, by DRIG";
		series x=HEFI2019_TOTAL_SCORE y=estimate / lineattrs=( thickness=3) group=drig;
		*series x=HEFI2019_TOTAL_SCORE y=lcl95 / lineattrs=(color=black pattern=longdash thickness=3) ;
		*series x=HEFI2019_TOTAL_SCORE y=ucl95 / lineattrs=(color=black pattern=longdash thickness=3) ;
		where varname = "pota";
		run;
		title1;

 /*************************************************************************/
 /* BY DRI GROUP: Logistic regression results, model estimates            */
 /*************************************************************************/

/* 3.1) Output logistic regression estimates data */
	proc sort data=reg_b2._logr_mod_hefi_y_rawBS;
	by replicate varname drig;
	run;

	/* note: logr_mod_hefi_y = Logistic Regression Model for HEFI with outcome Y, 
		raw association (i.e., unadjusted) */

	/* check extreme replicate */
		ods select moments quantiles extremeobs ;
		title1 color=red "Extreme replicates for Misc. B";
		proc univariate data=reg_b2._logr_mod_hefi_y_rawBS;
		class drig varname;
		var risk_diff ;
		id replicate;
		where replicate ne 0 ;
		run;
		title1;	

	/* split sample for the <boot_variance> macro */
	data log_hefi_mod0 log_hefi_mod_boot;
		set reg_b2._logr_mod_hefi_y_rawBS ;
	* split bootstrap and original;
		if replicate= 0 then output log_hefi_mod0; 
		else output log_hefi_mod_boot ;
	run;

/* 3.2) Use <boot_converge.sas> to look a normal-theory bootstrap variance estimation assumption */
%boot_converge(indata    = log_hefi_mod_boot, /* bootstrap estimates data */
			   checkconv = exp_beta risk_diff risk_ratio, /* variable for which we want to look at convergence */
			   where_sub = %str(where (varname='binary_dfe' & drig=15 & replicate ne 484)), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );

/* 3.3) Use <boot_variance.sas> macro for efficient variance estimation */
%boot_variance(inboot    = log_hefi_mod_boot, /* bootstrap estimates */
			   inbase    = log_hefi_mod0 ,  /* original estimatee */
			   byvar     = varname drig, /* If <by> variable used */
			   varlist   = risk_ratio,
			   null      = 1, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.2, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-2),
			   outdata   = log_beta_all1 /* name of output dataset */
			   ); 
			   
%boot_variance(inboot    = log_hefi_mod_boot, /* bootstrap estimates */
			   inbase    = log_hefi_mod0 ,  /* original estimatee */
			   byvar     = varname drig, /* If <by> variable used */
			   varlist   = log_beta ,
			   exp       = 1,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.2, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-1),
			   outdata   = log_beta_all2 /* name of output dataset */
			   ); 

%boot_variance(inboot    = log_hefi_mod_boot, /* bootstrap estimates */
			   inbase    = log_hefi_mod0 ,  /* original estimatee */
			   byvar     = varname drig, /* If <by> variable used */
			   varlist   = risk_diff ,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.2, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-2),
			   outdata   = log_beta_all3 /* name of output dataset */
			   ); 
			   
	data reslib.log_beta_drigf(rename=(_LABEL_=label));
		set log_beta_all1 log_beta_all2 log_beta_all3;
	* rename to facilitate interpretation ;
		if name = "log_beta" then name ="odds_ratio";
		if missing(exp) then exp=0;
	* change erroneous label ;
		_LABEL_ = substr(_LABEL_,index(_LABEL_,'|')+2);
		label _LABEL_ = " ";
	run;
	
	proc sort data=reslib.log_beta_drigf;
		by varname drig;
	run;

/* Clean temp */
	proc datasets lib=work nolist nodetails ;
	delete _: log_beta: log_hefi: lin_: ;
	run;

 /*************************************************************************/
 /* ALL: Distribution of usual intake, marginal Pr(X<x)                   */
 /*************************************************************************/

/* 1.1) Output marginal Y distribution data */
	proc sort data=bootlib._distrib_y_wBS;
	by replicate varname DRIg;
	run;

	/* note: distrib_y_rel_w_boot = Distribution of Y, Relative intakes, Wide format */

	/* check extreme replicate */
		ods select moments quantiles extremeobs ;
		title1 color=red "Extreme replicates for Misc. B";
		proc univariate data=bootlib._distrib_y_wBS;
		class  varname;
		var Mean ;
		id replicate;
		where replicate ne 0 and drig =0;
		run;
		title1;	

	/* 1.2) split sample for the <boot_variance> macro */
	data distrib_y_w0 distrib_y_w_boot ;
		set bootlib._distrib_y_wBS;
	* remove replicate with extreme values ;
		if replicate = 484 then delete;
	if replicate= 0 then output distrib_y_w0;
	else output distrib_y_w_boot;
	run;
	

/* 1.2) Use <boot_converge.sas> to look a normal-theory bootstrap variance estimation assumption */
%boot_converge(indata    = distrib_y_w_boot, /* bootstrap estimates data */
			   checkconv = Mean Pctile10 Pctile90, /* variable for which we want to look at convergence */
			   where_sub = %str(where (drig=0 and varname='dfe' & replicate ne 484)), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );
%boot_converge(indata    = distrib_y_w_boot, /* bootstrap estimates data */
			   checkconv = Mean Pctile10 Pctile90, /* variable for which we want to look at convergence */
			   where_sub = %str(where (drig=0 and varname='pota' & replicate ne 484)), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );
%boot_converge(indata    = distrib_y_w_boot, /* bootstrap estimates data */
			   checkconv = Mean Pctile10 Pctile90, /* variable for which we want to look at convergence */
			   where_sub = %str(where (drig=0 and varname='fibers' & replicate ne 484)), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );
%boot_converge(indata    = distrib_y_w_boot, /* bootstrap estimates data */
			   checkconv = Mean Pctile10 Pctile90, /* variable for which we want to look at convergence */
			   where_sub = %str(where (drig=0 and varname='mg' & replicate ne 484)), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );
	
/* 1.3) Use <boot_variance.sas> macro for efficient variance estimation */

/* 1.3.2) actual variance estimation */
%boot_variance(inboot    = distrib_y_w_boot, /* bootstrap estimates */
			   inbase    = distrib_y_w0,  /* original estimatee */
			   byvar     = varname drig, /* If <by> variable used */
			   varlist   = Mean &pctile_list ,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 4.2, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-1), 
			   outdata   = distrib_y_wf1 /* name of output dataset */
			   );

%boot_variance(inboot    = distrib_y_w_boot, /* bootstrap estimates */
			   inbase    = distrib_y_w0,  /* original estimatee */
			   byvar     = varname drig, /* If <by> variable used */
			   varlist   = binary_ear ,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 4.1, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-1), 
			   outdata   = distrib_y_wf2 /* name of output dataset */
			   ); 
			  
	data reslib.distrib_y_wf (rename=(_LABEL_ = label));
	retain drig name _LABEL_ ;
	length name $ 32;
		set distrib_y_wf1(in=a) distrib_y_wf2(in=b) ;
	label _LABEL_ = " ";
	run;
	
	proc sort;
	by varname drig;
	run;

	proc print data=reslib.distrib_y_wf label;
	id varname name label nboot drig;
	var  estimate se lcl95 ucl95 pvalue svalue cv;
	run;

 /*************************************************************************/
 /* BY CATEGORY OF SCORE: Distribution of usual intake                    */
 /*************************************************************************/

/* 2.1) Append marginal Y distribution data */
	proc sort data=bootlib._distrib_catscore_wBS;
	by replicate varname catscore;
	run;

	/* note: distrib_y_rel_w_boot = Distribution of Y, Relative intakes, Wide format */

	/* check extreme replicate */
		ods select moments quantiles extremeobs ;
		title1 color=red "Extreme replicates for Misc. B";
		proc univariate data=bootlib._distrib_catscore_wBS;
		class catscore;
		var Mean ;
		id replicate;
		where replicate ne 0 and catscore <9000;
		run;
		title1;	

	/* split sample for the <boot_variance> macro */
	data distrib_catscore0 distrib_catscore_boot ;
		set bootlib._distrib_catscore_wBS;
	* remove replicate with extreme values ;
		if replicate = 484 then delete;
	if replicate= 0 then output distrib_catscore0;
	else output distrib_catscore_boot;
	run;


/* 2.2) Use <boot_converge.sas> to look a normal-theory bootstrap variance estimation assumption */
%boot_converge(indata    = distrib_catscore_boot, /* bootstrap estimates data */
			   checkconv = Mean Pctile10 Pctile90, /* variable for which we want to look at convergence */
			   where_sub = %str(where (catscore=4 and varname='dfe' & replicate ne 484)), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );

	
/* 2.3) Use <boot_variance.sas> macro for efficient variance estimation */

	/* 2.3.2) Variance estimation */
	%boot_variance(inboot    = distrib_catscore_boot, /* bootstrap estimates */
				   inbase    = distrib_catscore0,  /* original estimatee */
				   byvar     = varname catscore, /* If <by> variable used */
				   varlist   = Mean &pctile_list ,
				   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
				   set_confidence_level=95, /* 95 = 95%CI */
				   numberfmt = 5.2, /* formatting of estimates */
				   set_degrees_freedom = %str(nboot-1), 
				   outdata   = distrib_catscoref1 /* name of output dataset */
				   ); 
	
	%boot_variance(inboot    = distrib_catscore_boot, /* bootstrap estimates */
				   inbase    = distrib_catscore0,  /* original estimatee */
				   byvar     = varname catscore, /* If <by> variable used */
				   varlist   = binary_ear ,
				   null      = 1, /* Null H for pvalues (i.e., 0 for a difference) */
				   set_confidence_level=95, /* 95 = 95%CI */
				   numberfmt = 4.2, /* formatting of estimates */
				   set_degrees_freedom = %str(nboot-2), 
				   outdata   = distrib_catscoref2 /* name of output dataset */
				   ); 
				   
 	/* 2.3.3) Check and save output */
	data reslib.distrib_w_catscoref(rename=(_LABEL_=label));
	retain name _LABEL_ catscore;
	length name $ 32;
		set distrib_catscoref1 distrib_catscoref2 ;
	label _LABEL_ = " ";
	run;
	
	proc sort;
	by varname catscore;
	run;
	
	title1 "Nutrient intake by quarter of total HEFI-2019 score";
	proc print data=reslib.distrib_w_catscoref label;
	id varname name nboot catscore;
	var  estimate se lcl95 ucl95 pvalue svalue cv;
	where index(name,"Prob")=0;
	run;

	title2 "Pr(X<x) by quarter of total HEFI-2019 score";
	proc print data=reslib.distrib_w_catscoref label;
	id varname name label nboot catscore;
	var  estimate se lcl95 ucl95 pvalue svalue cv;
	where index(name,"Prob")>0;
	run;
	title1;
	title2;

/* Clean temp */
	proc datasets lib=work nolist nodetails ;
	delete _: distrib_: ;
	run;

 /*************************************************************************/
 /*                                                                       */
 /*     REF. CODE 6: Variance estimation for VIT. A NUTRIENT outcome      */
 /*                                                                       */
 /*************************************************************************/

 /************************************************/
 /*            Files and folder set-up           */
 /************************************************/

/* 0) define suffix of current outcome */
	
	%let suffix = _vit_a ;

/* 1) assign proper library using the same suffix */
	%mcmclib(&suffix.);
	%MakeMCMCFolder(NCIPath = &path./NCI/MCMC&suffix.,
					prefix  = Reg_,
					suffix  = &subgrp.
					);

 /*************************************************************************/
 /* ALL: Linear regression results, model estimates                       */
 /*************************************************************************/

/* 1.1) Output linear regression estimates data */
	proc sort data=regb._rcs_hefi_rawBS;
	by replicate;
	run;

	/* note: rcs_hefi_raw = (Linear) regression with Restricted Cubic Spline transformation
		of the HEFI, raw association (i.e., unadjusted) */

	/* split sample for the <boot_variance> macro */
	data lin_hefi_raw0 lin_hefi_raw_boot ;
		set regb._rcs_hefi_rawBS ;
	if replicate=0 then output lin_hefi_raw0; 
	else output lin_hefi_raw_boot;
	run;

/* 1.2) Use <boot_converge.sas> to look a normal-theory bootstrap variance estimation assumption */
%boot_converge(indata    = regb._rcs_hefi_rawBS, /* bootstrap estimates data */
			   checkconv = rc_beta, /* variable for which we want to look at convergence */
			   where_sub = %str(where varname = 'rae'), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );

/* 1.3) Use <boot_variance.sas> macro for efficient variance estimation */
%boot_variance(inboot    = lin_hefi_raw_boot, /* bootstrap estimates */
			   inbase    = lin_hefi_raw0,  /* original estimatee */
			   byvar     = varname ,
			   varlist   = rc_beta r2,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.1, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-2), 
			   outdata   = reslib.lin_beta_allf /* name of output dataset */
			   ); 

	proc print data=reslib.lin_beta_allf label;
	id name nboot;
	var  estimate se lcl95 ucl95 pvalue svalue cv;
	run;

 /*************************************************************************/
 /* ALL: Linear regression results, predicted values for plot             */
 /*************************************************************************/

/* 2.1) Output linear regression plot data */
	proc sort data=regb._lin_hefi_plotBS;
	by replicate varname name;
	run;

	/* note: plm_rcs_hefi_raw = Plot of Linear Model with RCS for HEFI,
		raw association (i.e., unadjusted) */
	
	/* split sample for the <boot_variance> macro */
	data  lin_hefi_plot0 lin_hefi_plot_boot;
		set regb._lin_hefi_plotBS;
	* remove x values;
	if (name = "HEFI2019_TOTAL_SCORE") then delete;
	if replicate =0 then output lin_hefi_plot0;
	else output lin_hefi_plot_boot;
	run;
	
	
/* 2.2) Use <boot_converge.sas> to look a normal-theory bootstrap variance estimation assumption */
%boot_converge(indata    = lin_hefi_plot_boot, /* bootstrap estimates data */
			   checkconv = pred25 pred100 pred175, /* variable for which we want to look at convergence */
			   where_sub = %str(where varname = 'rae'), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );

/* 2.3) Use <boot_variance.sas> macro for efficient variance estimation */
%boot_variance(inboot    = lin_hefi_plot_boot, /* bootstrap estimates */
			   inbase    = lin_hefi_plot0 ,  /* original estimatee */
			   where     = %str(name = "predicted"), /* variance for predicted values only */	   
			   byvar     = varname ,
			   varlist   = pred1-pred200,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.1, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-1),
			   outdata   = lin_plot /* name of output dataset */
			   ); 

	/* 2.3.1) Add X values to the data */
	
	proc transpose data=regb._lin_hefi_plotBS out=_xvalues prefix=HEFI2019_TOTAL_SCORE;
	by varname ;
	var pred1-pred200;
	where name = "HEFI2019_TOTAL_SCORE";
	run;
	
	proc sort data=_xvalues;
	by varname _NAME_ ;
	run;
	
	proc sort data=lin_plot;
	by varname name ;
	run;
	
	data reslib.lin_plot_allf;
	retain name HEFI2019_TOTAL_SCORE estimate ;
		merge _xvalues(rename=(HEFI2019_TOTAL_SCORE1=HEFI2019_TOTAL_SCORE _NAME_=name))
			  lin_plot  ;
		by varname name ;
	* sort from 1 to 200 ;
	pred_id = input(compress(name,,'a'),10.);
	run;
	
	proc sort data=reslib.lin_plot_allf;
	by varname pred_id ;
	run;
	
	/* 2.3.2) Draft plot */
		proc sgplot data=reslib.lin_plot_allf ;
		title1 "Usual rae intake (mg)  | HEFI-2019 total score in adults 65y+ from CCHS 2015 - Nutrition";
		series x=HEFI2019_TOTAL_SCORE y=estimate / lineattrs=(color=black thickness=3);
		series x=HEFI2019_TOTAL_SCORE y=lcl95 / lineattrs=(color=black pattern=longdash thickness=3);
		series x=HEFI2019_TOTAL_SCORE y=ucl95 / lineattrs=(color=black pattern=longdash thickness=3);
		where varname = "rae" ;
		run;
		title1;

 /*************************************************************************/
 /* ALL: Logistic regression results, model estimates                     */
 /*************************************************************************/

/* 3.1) Output logistic regression estimates data */
	proc sort data=regb._logr_mod_hefi_y_rawBS;
	by replicate varname ;
	run;

	/* note: logr_mod_hefi_y = Logistic Regression Model for HEFI with outcome Y, 
		raw association (i.e., unadjusted) */

	data log_hefi_mod0 log_hefi_mod_boot;
		set regb._logr_mod_hefi_y_rawBS ;
	* split bootstrap and original;
		if replicate= 0 then output log_hefi_mod0; 
		else output log_hefi_mod_boot ;
	run;

	/* note: depending on the proportion of outcome, using log_beta or risk_ratio may be more stable.
		Hence, variance of both are estimated. */
	
/* 3.2) Use <boot_converge.sas> to look a normal-theory bootstrap variance estimation assumption */
%boot_converge(indata    = log_hefi_mod_boot, /* bootstrap estimates data */
			   checkconv = log_beta risk_ratio, /* variable for which we want to look at convergence */
			   where_sub = %str(where (varname=('binary_rae'))), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );


/* 3.3) Use <boot_variance.sas> macro for efficient variance estimation */
%boot_variance(inboot    = log_hefi_mod_boot, /* bootstrap estimates */
			   inbase    = log_hefi_mod0 ,  /* original estimatee */
			   byvar     = varname, /* If <by> variable used */
			   varlist   = risk_ratio,
			   null      = 1, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.2, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-2),
			   outdata   = log_beta_all1 /* name of output dataset */
			   );

%boot_variance(inboot    = log_hefi_mod_boot, /* bootstrap estimates */
			   inbase    = log_hefi_mod0 ,  /* original estimatee */
			   byvar     = varname, /* If <by> variable used */
			   varlist   = log_beta ,
			   exp       = 1 , /* exponentiate lcl + ucl once estimated */
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.2, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-1),
			   outdata   = log_beta_all2 /* name of output dataset */
			   ); 

%boot_variance(inboot    = log_hefi_mod_boot, /* bootstrap estimates */
			   inbase    = log_hefi_mod0 ,  /* original estimatee */
			   byvar     = varname, /* If <by> variable used */
			   varlist   = risk_diff ,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.2, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-2),
			   outdata   = log_beta_all3 /* name of output dataset */
			   ); 

	data reslib.log_beta_allf(rename=(_LABEL_=label));
		set log_beta_all1 log_beta_all2 log_beta_all3;
	* rename to facilitate interpretation ;
		if name = "log_beta" then name ="odds_ratio";
		if missing(exp) then exp=0;
	* change erroneous label ;
		_LABEL_ = substr(_LABEL_,index(_LABEL_,'|')+2);
		label _LABEL_ = " ";
	run;
	
	proc sort data=reslib.log_beta_allf;
		by varname;
	run;

/* Clean temp */
	proc datasets lib=work nolist nodetails ;
	delete _: lin: ;
	run;

 /*************************************************************************/
 /* BY DRI GROUP: Linear regression results, model estimates              */
 /*************************************************************************/

/* 1.1) Output linear regression estimates data */
	proc sort data=reg_b2._rcs_hefi_rawBS;
	by replicate varname drig;
	run;
	
	/* note: rcs_hefi_raw = (Linear) regression with Restricted Cubic Spline transformation
		of the HEFI, raw association (i.e., unadjusted) */

	/* split sample for the <boot_variance> macro */
	data lin_hefi_raw0 lin_hefi_raw_boot ;
		set reg_b2._rcs_hefi_rawBS ;
	if replicate=0 then output lin_hefi_raw0; 
	else output lin_hefi_raw_boot;
	run;

/* 1.2) Use <boot_converge.sas> to look a normal-theory bootstrap variance estimation assumption */
%boot_converge(indata    = lin_hefi_raw_boot, /* bootstrap estimates data */
			   checkconv = rc_beta, /* variable for which we want to look at convergence */
			   where_sub = %str(where (drig=12 and varname='rae')), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );
%boot_converge(indata    = lin_hefi_raw_boot, /* bootstrap estimates data */
			   checkconv = rc_beta, /* variable for which we want to look at convergence */
			   where_sub = %str(where (drig=15 and varname='rae')), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );

/* 1.3) Use <boot_variance.sas> macro for efficient variance estimation */
%boot_variance(inboot    = lin_hefi_raw_boot, /* bootstrap estimates */
			   inbase    = lin_hefi_raw0,  /* original estimatee */
			   byvar     = varname drig, /* If <by> variable used */
			   varlist   = rc_beta ,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.1, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-2), 
			   outdata   = reslib.lin_beta_drigf /* name of output dataset */
			   ); 

	proc print data=reslib.lin_beta_drigf label;
	id varname drig name nboot ;
	var  estimate se lcl95 ucl95 pvalue svalue cv;
	run;

 /*************************************************************************/
 /* BY DRI GROUP: Linear regression results, predicted values for plot    */
 /*************************************************************************/

/* 2.1) Output linear regression plot data */
	proc sort data=reg_b2._lin_hefi_plotBS;
	by replicate varname drig name;
	run;

	/* note: plm_rcs_hefi_raw = Plot of Linear Model with RCS for HEFI,
		raw association (i.e., unadjusted) */
	
	/* split sample for the <boot_variance> macro */
	data  lin_hefi_plot0 lin_hefi_plot_boot;
		set reg_b2._lin_hefi_plotBS;
	* remove x values;
	if (name = "HEFI2019_TOTAL_SCORE") then delete;
	if replicate =0 then output lin_hefi_plot0;
	else output lin_hefi_plot_boot;
	run;
	
/* 2.2) Use <boot_converge.sas> to look a normal-theory bootstrap variance estimation assumption */
%boot_converge(indata    = lin_hefi_plot_boot, /* bootstrap estimates data */
			   checkconv = pred25 pred100 pred175, /* variable for which we want to look at convergence */
			   where_sub = %str(where drig='12' and varname='rae'), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );
%boot_converge(indata    = lin_hefi_plot_boot, /* bootstrap estimates data */
			   checkconv = pred25 pred100 pred175, /* variable for which we want to look at convergence */
			   where_sub = %str(where drig='15' and varname='rae'), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );
			   
/* 2.3) Use <boot_variance.sas> macro for efficient variance estimation */
%boot_variance(inboot    = lin_hefi_plot_boot, /* bootstrap estimates */
			   inbase    = lin_hefi_plot0 ,  /* original estimatee */
			   where     = %str(name = "predicted"), /* variance for predicted values only */
			   byvar     = varname drig,
			   varlist   = pred1-pred200,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.1, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-1),
			   outdata   = lin_plot /* name of output dataset */
			   ); 

	/* 2.3.1) Add X values to the data */
		
	data _xvalues;
		set reg_2.plmrcst_rae_raw0(in=a) ;
	length varname $ 32;
	if a then varname = "rae";
	run;
		
		/* some manipulation to output <name> variable*/
		proc sort;
			by varname drig;
		run;
		
		data _xvalues;
		set _xvalues;
		by varname drig;
		retain pred_id ;
		if first.drig then pred_id=1;
		else pred_id+1;
		length name $  8;
		name = cats("pred",pred_id);
		run;
		
	proc sort data=_xvalues;
	by varname drig name ;
	run;
	
	proc sort data=lin_plot;
	by varname drig name ;
	run;
	
	data reslib.lin_plot_drigf(rename=(_drig=drig));
	retain varname _drig pred_id name HEFI2019_TOTAL_SCORE estimate ;
		merge _xvalues
			  lin_plot  ;
		by varname drig name ;
	* change drig to num ;
		_drig = input(drig, 10.);
	* add merge error flag for extra safety ;
		if predicted ne estimate then _MergeError_=1;
		else _MergeError_=0;
		label _MergeError_ = "Flag for any merge error";
	*clean; 
		drop drig predicted;
	run;
	
	proc sort data=reslib.lin_plot_drigf;
	by varname drig pred_id ;
	run;	

	
	/* 2.3.2) Draft plot */
		proc sgplot data=reslib.lin_plot_drigf ;
		title1 "Usual rae (mg)  | HEFI-2019 total score in adults 65y+ from CCHS 2015 - Nutrition, by DRIG";
		series x=HEFI2019_TOTAL_SCORE y=estimate / lineattrs=( thickness=3) group=drig;
		*series x=HEFI2019_TOTAL_SCORE y=lcl95 / lineattrs=(color=black pattern=longdash thickness=3) ;
		*series x=HEFI2019_TOTAL_SCORE y=ucl95 / lineattrs=(color=black pattern=longdash thickness=3) ;
		where varname = "rae";
		run;
		title1;


/* Clean temp */
	proc datasets lib=work nolist nodetails ;
	delete _: lin: ;
	run;

 /*************************************************************************/
 /* BY DRI GROUP: Logistic regression results, model estimates            */
 /*************************************************************************/

/* 3.1) Output logistic regression estimates data */
	proc sort data=reg_b2._logr_mod_hefi_y_rawBS;
	by replicate varname drig;
	run;

	/* note: logr_mod_hefi_y = Logistic Regression Model for HEFI with outcome Y, 
		raw association (i.e., unadjusted) */

	data log_hefi_mod0 log_hefi_mod_boot;
		set reg_b2._logr_mod_hefi_y_rawBS ;
	* split bootstrap and original;
		if replicate= 0 then output log_hefi_mod0; 
		else output log_hefi_mod_boot ;
	run;

/* 3.2) Use <boot_converge.sas> to look a normal-theory bootstrap variance estimation assumption */
*%boot_converge(indata    = log_hefi_mod_boot, /* bootstrap estimates data */
			   checkconv = exp_beta risk_diff risk_ratio, /* variable for which we want to look at convergence */
			   where_sub = %str(where (outcome=('binary_06'))), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );

/* 3.3) Use <boot_variance.sas> macro for efficient variance estimation */
%boot_variance(inboot    = log_hefi_mod_boot, /* bootstrap estimates */
			   inbase    = log_hefi_mod0 ,  /* original estimatee */
			   byvar     = varname drig, /* If <by> variable used */
			   varlist   = risk_ratio,
			   null      = 1, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.2, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-2),
			   outdata   = log_beta_all1 /* name of output dataset */
			   ); 

%boot_variance(inboot    = log_hefi_mod_boot, /* bootstrap estimates */
			   inbase    = log_hefi_mod0 ,  /* original estimatee */
			   byvar     = varname drig, /* If <by> variable used */
			   varlist   = log_beta ,
			   exp       = 1,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.2, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-1),
			   outdata   = log_beta_all2 /* name of output dataset */
			   );

%boot_variance(inboot    = log_hefi_mod_boot, /* bootstrap estimates */
			   inbase    = log_hefi_mod0 ,  /* original estimatee */
			   byvar     = varname drig, /* If <by> variable used */
			   varlist   = risk_diff ,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 5.2, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-2),
			   outdata   = log_beta_all3 /* name of output dataset */
			   ); 

	data reslib.log_beta_drigf(rename=(_LABEL_=label));
		set log_beta_all1 log_beta_all2 log_beta_all3;
	* rename to facilitate interpretation ;
		if name = "log_beta" then name ="odds_ratio";
		if missing(exp) then exp=0;
	* change erroneous label ;
		_LABEL_ = substr(_LABEL_,index(_LABEL_,'|')+2);
		label _LABEL_ = " ";
	run;
	
	proc sort data=reslib.log_beta_drigf;
		by varname drig;
	run;

/* Clean temp */
	proc datasets lib=work nolist nodetails ;
	delete _: log_beta: log_hefi: ;
	run;

 /*************************************************************************/
 /* ALL: Distribution of usual intake, marginal Pr(X<x)                   */
 /*************************************************************************/

/* 1.1) Output marginal Y distribution data */
	proc sort data=bootlib._distrib_y_wBS;
	by replicate varname DRIg;
	run;

	/* note: distrib_y_rel_w_boot = Distribution of Y, Relative intakes, Wide format */

	/* split sample for the <boot_variance> macro */
	data distrib_y_w0 distrib_y_w_boot ;
		set bootlib._distrib_y_wBS;
	if replicate= 0 then output distrib_y_w0;
	else output distrib_y_w_boot;
	run;
	

/* 1.2) Use <boot_converge.sas> to look a normal-theory bootstrap variance estimation assumption */
%boot_converge(indata    = distrib_y_w_boot, /* bootstrap estimates data */
			   checkconv = Mean Pctile10 Pctile90, /* variable for which we want to look at convergence */
			   where_sub = %str(where (drig=0 and varname='rae')), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );

	
/* 1.3) Use <boot_variance.sas> macro for efficient variance estimation */

/* 1.3.2) actual variance estimation */
%boot_variance(inboot    = distrib_y_w_boot, /* bootstrap estimates */
			   inbase    = distrib_y_w0,  /* original estimatee */
			   byvar     = varname drig, /* If <by> variable used */
			   varlist   = Mean &pctile_list ,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 4.2, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-1), 
			   outdata   = distrib_y_wf1 /* name of output dataset */
			   );

%boot_variance(inboot    = distrib_y_w_boot, /* bootstrap estimates */
			   inbase    = distrib_y_w0,  /* original estimatee */
			   byvar     = varname drig, /* If <by> variable used */
			   varlist   = binary_ear ,
			   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
			   set_confidence_level=95, /* 95 = 95%CI */
			   numberfmt = 4.1, /* formatting of estimates */
			   set_degrees_freedom = %str(nboot-1), 
			   outdata   = distrib_y_wf2 /* name of output dataset */
			   ); 
			  
	data reslib.distrib_y_wf (rename=(_LABEL_ = label));
	retain drig name _LABEL_ ;
	length name $ 32;
		set distrib_y_wf1(in=a) distrib_y_wf2(in=b) ;
	label _LABEL_ = " ";
	run;
	
	proc sort;
	by varname drig;
	run;

	proc print data=reslib.distrib_y_wf label;
	id varname name label nboot drig;
	var  estimate se lcl95 ucl95 pvalue svalue cv;
	run;

 /*************************************************************************/
 /* BY CATEGORY OF SCORE: Distribution of usual intake                    */
 /*************************************************************************/

/* 2.1) Append marginal Y distribution data */
	data bootlib._distrib_catscore_wBS;
		set reslib.distrib_w_catscore0
			bootlib.distrib_w_catscore: 
			;
	* PREVALENCE: rescale prob to percentage point;
	if catscore < 9000 then binary_ear = binary_ear*100;
	run;
	
	/* note: distrib_y_rel_w_boot = Distribution of Y, Relative intakes, Wide format */

	proc sort data=bootlib._distrib_catscore_wBS;
	by replicate varname catscore;
	run;
	
	/* 2.1.2) split sample for the <boot_variance> macro */
	data distrib_catscore0 distrib_catscore_boot ;
		set bootlib._distrib_catscore_wBS;
	if replicate= 0 then output distrib_catscore0;
	else output distrib_catscore_boot;
	run;


/* 2.2) Use <boot_converge.sas> to look a normal-theory bootstrap variance estimation assumption */
%boot_converge(indata    = distrib_catscore_boot, /* bootstrap estimates data */
			   checkconv = Mean Pctile10 Pctile90, /* variable for which we want to look at convergence */
			   where_sub = %str(where (catscore=4 and varname='rae')), /* if a subgroup is used - NA for this example */
			   sefmt     = 5.2, /* formatting for plot */
			   normality = 1 /* if 1, checks normality */
			   );

	
/* 2.3) Use <boot_variance.sas> macro for efficient variance estimation */

	/* 2.3.2) Variance estimation */
	%boot_variance(inboot    = distrib_catscore_boot, /* bootstrap estimates */
				   inbase    = distrib_catscore0,  /* original estimatee */
				   byvar     = varname catscore, /* If <by> variable used */
				   varlist   = Mean &pctile_list ,
				   null      = 0, /* Null H for pvalues (i.e., 0 for a difference) */
				   set_confidence_level=95, /* 95 = 95%CI */
				   numberfmt = 5.2, /* formatting of estimates */
				   set_degrees_freedom = %str(nboot-1), 
				   outdata   = distrib_catscoref1 /* name of output dataset */
				   ); 
	
	%boot_variance(inboot    = distrib_catscore_boot, /* bootstrap estimates */
				   inbase    = distrib_catscore0,  /* original estimatee */
				   byvar     = varname catscore, /* If <by> variable used */
				   varlist   = binary_ear ,
				   null      = 1, /* Null H for pvalues (i.e., 0 for a difference) */
				   set_confidence_level=95, /* 95 = 95%CI */
				   numberfmt = 4.2, /* formatting of estimates */
				   set_degrees_freedom = %str(nboot-2), 
				   outdata   = distrib_catscoref2 /* name of output dataset */
				   ); 
				   
 	/* 2.3.3) Check and save output */
	data reslib.distrib_w_catscoref(rename=(_LABEL_=label));
	retain name _LABEL_ catscore;
	length name $ 32;
		set distrib_catscoref1 distrib_catscoref2 ;
	label _LABEL_ = " ";
	run;
	
	proc sort;
	by varname catscore;
	run;
	
	title1 "Nutrient intake by quarter of total HEFI-2019 score";
	proc print data=reslib.distrib_w_catscoref label;
	id varname name nboot catscore;
	var  estimate se lcl95 ucl95 pvalue svalue cv;
	where index(name,"Prob")=0;
	run;

	title2 "Pr(X<x) by quarter of total HEFI-2019 score";
	proc print data=reslib.distrib_w_catscoref label;
	id varname name label nboot catscore;
	var  estimate se lcl95 ucl95 pvalue svalue cv;
	where index(name,"Prob")>0;
	run;
	title1;
	title2;

/* Clean temp */
	proc datasets lib=work nolist nodetails ;
	delete _: distrib_: ;
	run;
