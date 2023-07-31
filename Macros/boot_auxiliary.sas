/********************************************************************************/
/*                        Auxilirary bootstrap macros                           */
/*               Suit of macros to perform bootstrap related tasks              */
/*                          Didier Brassard  2020-10                            */
/*                                                                              */
/* %boot_convergence; --> checks bootstrap convergence                          */
/* %ci; --> caclulate bootstrap confidence intervals                            */
/* %boot_variance; --> merge files and calculate 95%ci based on %ci             */
/*                                                                              */
/********************************************************************************/


	/********************************************************************************/
	/***                       The BOOT_CONVERGE macro                            ***/
	/*                           Adapted from Stat Can                              */
	/*                        Didier Brassard  2020-09 v2.2                         */
	/*                                                                              */
	/*  This macro looks graphically at the convergence of one/many estimate(s)     */
	/*                                                                              */
	/********************************************************************************/
	/*                                                                              */
	/* indata    = name of a data with all bootstrap estimates in tall format       */
	/* checkconv = name of the variable(s) to check                                 */
	/* where_sub = where statement (eg, where visit=1, where sex="male", ...)       */
	/* b         = number of bootstrap resamples                                    */
	/* showb     = add stderr across b 1 to showb (eg, 100 add stderr of b1-b100)   */
	/* normality = if 1, will assess normality of the bootstrap estimates           */
	/* clean     = if 1, delete temporary datasets (default option)                 */
	/*                                                                              */
	/********************************************************************************/
	
	%macro boot_converge(indata=, checkconv=, where_sub=, b=, showb=,sefmt=4.1,normality=0,clean=1);
	
	/* check showb (showb cannot be b)*/
		%if (&showb ne %str()) & (&showb = &b) %then %do;
		%local showb;
			%let showb = %eval(&b-1);
		%end;
	
	/* make sure user call doesnt erase defaults */
		%if (&sefmt = %str()) %then %do;
			%local sefmt;
			%let sefmt=4.1;
		%end;
		
		%if (&normality = %str()) %then %do;
			%local normality;
			%let normality=0;
		%end;

		%if (&clean = %str()) %then %do;
			%local clean;
			%let clean=1;
		%end;

	/* create temporary data for analysis */	
		data temp;
			set &indata;
			&where_sub;
		run;
				
			/* count number of bootstraps */
			%local nboot;
			data _NULL_;
				if 0 then set work.temp nobs=n;
					call symputx('nboot',n);
				stop;
			run;
			
			/* show number */
			%put # &sysmacroname. &=nboot for &=indata ; 
			
			/* use observed if inconsistent */
			%if &nboot ne &b %then %do;
			%if (&b ne %str()) %then %put WARNING: &sysmacroname. Number of observed bootstrap (&nboot) inconsistent with user call (&b);
			%local b;
			%let b=&nboot;
			%end;
	
		%do kth=1 %to %sysfunc(countw(&checkconv));
			%let next&kth = %scan(&checkconv, &kth, %str( ));
	
			/* START OF LOOP */
			proc transpose data=TEMP out=variance_t(keep=_NAME_ COL1-COL&b ) prefix=COL;
				variable &&next&kth;
			RUN;
	
			%local i;
			
			/* calculate variance across bootstraps */
			data variance_t1;
				set variance_t;
				array col(&b) col1-col&b;
				array variance(%eval(&b-1)) variance2-variance&b;
	
				%do i=1 %to %eval(&b-1);
					variance(&i)=var(of col1-col%eval(&i+1));
				%end;
			run;
			
			/* transpose */
			proc transpose data=variance_t1 out=variance_t2;
				variable variance:;
			run;
			
			/* assign obs n */
			%local showSE ;
			data variance_t2;
				set variance_t2;
				retain OBS 0;
				OBS+1;
				se=sqrt(&&next&kth);
				%if (&showb ne %str()) %then %do;
						if (OBS = &showb) then call symputx ("showSE", put(se, &sefmt. ));
				%end;
			run;
			
			%if (&showb ne %str()) %then %do;
				%put # b1-b&showb resamples stderr = &showSE ;
			%end;
				
				/* output max stderr */
				%local max;
					proc means data=variance_t2 min max noprint;
					var Se;
					output out=_mean min=min max=max;
					run;
					
					data _null_;
						set _mean;
					call symputx ("max", put(max, &sefmt. ));
					run;
				
				%put # Max bootstrap variance = &max ;
				
			%if (%sysevalf(%superq(where_sub)=,boolean)=1) %then %do;
				title2 "&SYSMACRONAME assessment, variable = &&next&kth.";
			%end;
			%else %do;
				title2 "&SYSMACRONAME assessment, variable = &&next&kth., &where_sub.";
			%end;
	
			proc sgplot data=variance_t2;
				scatter y=se x=OBS / markerattrs=(color=black);
				xaxis label="Number of bootstraps" labelattrs=(weight=bold);
				yaxis label="Standard errors across bootstraps" labelattrs=(weight=bold);
				
				%if (&showb ne %str()) %then %do;
					refline &showb / axis=x label=("b1-b&showb = &showSE") lineattrs=(pattern=longdash) ;
				%end;
				
				refline &max / axis=y label=("Max=&max") lineattrs=(pattern=longdash) ;
			run;
	
			title2;
			
				/* Assess normality of the bootstrap distribution*/
				%if (&normality=1) %then %do;
					%if (%sysevalf(%superq(where_sub)=,boolean)=1) %then %do;
						title2 "Bootstrap distribution of parameters estimate, variable = &&next&kth.";
					%end;
					%else %do;
						title2 "Bootstrap distribution of parameters estimate, variable = &&next&kth., &where_sub.";
					%end;
					proc sgplot data=temp;
					histogram &&next&kth;
					density &&next&kth / type=kernel;
					density &&next&kth / type=normal;
					run;
					
					title3 italic "Normality check - should follow the linear curve";
					PROC UNIVARIATE DATA=temp noprint;
						probplot &&next&kth / normal(mu=est sigma=est);
						VAR &&next&kth;
					RUN;
					title3;
					title2;
				%end;
	
			/* END OF LOOP */
	
	%end;
	
	/* clean temporary datasets */
	%if (&clean=1) %then %do;
		proc datasets lib=work nolist nodetails ;
			delete temp _mean variance_t variance_t1 variance_t2 ;
		run;
	%end;
	%mend;


	/********************************************************************************/
	/***                       The BOOT_VARIANCE macro                            ***/
	/*                        Didier Brassard  2021-05 v2.4                         */
	/*                                                                              */
	/*      This macro estimates bootstrap variance from bootstrap samples          */
	/*                                                                              */
	/********************************************************************************/
	/*                                                                              */
	/* inboot    = name of a data with bootstrap estimates (1 row per boot.)        */
	/* inbase    = name of the data with base estimates (1 row OR 1 row per byvar)  */
	/* where_estimate = where statement (eg, where visit=1, where sex="male", ...)  */
	/* byvar     = if specified, bootstrap std are calculated indendently by byvar  */
	/* outdata   = name of the ouput dataset with estimate and variance             */
	/* varlist   = stderr are estimated on those variables (can be many)            */
	/* set_degrees_freedom = DF for standard errors (default=number of bootstrap-1) */
	/* null      = value for null hypothesis testing (default=0)                    */
	/* numberfmt = format for estimate (e.g., 4.1)                                  */
	/* exp       = if 1, estimate, lcl and ucl are exponentiated (default=0)        */
	/*             To use when estimates are expressed on the ln scale              */
	/*                                                                              */
	/********************************************************************************/

	/********************************************************************************/
	/***   Sub-macro for 95ci and auto-formatting of base and bootstrap estimates ***/
	/********************************************************************************/
	
		%macro ci95(estimate=,set_confidence_level=95,set_degrees_freedom=,null=0,ptail=2);
		
		/* input parameters */
		    alpha = 1 - (&set_confidence_level/100);
		    one_minus_half_alpha = 1 - alpha/2;
		    
		/* to avoid passing invalid arguments to <probt> */
		if &set_degrees_freedom >0 then do;
			    t_quant = quantile('t', one_minus_half_alpha, &set_degrees_freedom);
			    
			/* confidence limits at set confidence level */
			    lcl&set_confidence_level = &estimate - t_quant * se;
			    ucl&set_confidence_level = &estimate + t_quant * se;
			  
			/* tvalue + pvalue - different than null */
			if (se ne 0) then do;
				tvalue =abs( ( &estimate - &null ) / se );
			    format pvalue PVALUE6.4;
				pvalue = &ptail * (1 - probt(tvalue, &set_degrees_freedom ) );
			end;
			else do;
				tvalue = .;
				pvalue = .;
			end; /* end of conditional on se ne 0 */
		end; /* end of conditional on set_degrees_freedom >0 */
		else do;
			t_quant=.;
			lcl&set_confidence_level=.;
			ucl&set_confidence_level=.;
			tvalue = .;
			pvalue = .;
			
		end; /* end of conditional on set_degrees_freedom <=0  */
		
		label pvalue = "&ptail.-sided pvalue (null H value=&null.)";
		
		%mend ci95;
	
	/* sub-macro for autoimport of base and bootstrap estimates	*/
	
		/* ### both input data should have within-bootstrap information 
			on a single row and between-bootstrap using tall format */
			
	%macro boot_variance(inboot=,inbase=,where_estimate=,byvar=,outdata=bootstrap,varlist=,set_confidence_level=95,set_degrees_freedom=%str(nboot-1),null=0,exp=0,echo=0,numberfmt=4.1);

	/*********************************************/
	/*        confirm correct macro call         */
	/*********************************************/
	
	%if (%sysevalf(%superq(inboot)=,boolean)) OR (%sysevalf(%superq(inbase)=,boolean)) %then %do;
		%put ERROR: Missing <inboot> and/or <inbase> data ;
		%return;
	%end;
	
	/* avoid writing over null value */
		%if (&echo = %str()) %then %do;
			%local echo;
			%let echo=0;
		%end;
		
	/* ECHO ON: hide preliminary steps */
	%if (&echo=1) %then %do;
		options nomprint;
		options nosource nonotes;
	%end;
	
		%if (&set_confidence_level = %str()) %then %do;
			%local set_confidence_level;
			%let set_confidence_level = 95;
			%put # &sysmacroname. Confidence level = &set_confidence_level;
		%end;
		%if (&null = %str()) %then %do;
			%local null;
			%let null=0;
		%end;
		%if (&numberfmt = %str()) %then %do;
			%local numberfmt;
			%let numberfmt=4.1;
		%end;
		%if (&exp = %str()) %then %do;
			%local exp ;
			%let exp=0;
		%end;
		%else %if (&exp ne 0) %then %do;
			%put # &sysmacroname. Estimates, LCL and UCL are exponentiated (StdErr expressed on ln scale);

			%if (&null = 1) %then %do;
				%put WARNING: &sysmacroname. The null hypothesis is usually zero (0) for ln(b), but could be 1 for exp(b). Please revise. ;
			%end;
		%end;
		%if %sysfunc(countw(&set_confidence_level,' '))>1 %then %do;
			%put # &sysmacroname. Confidence levels = &set_confidence_level;
			%if (&exp ne 0) %then %put WARNING: Currently not possible to use the exp option with multiple confidence levels;
		%end;
		%if %sysevalf(%superq(set_degrees_freedom)=,boolean) %then %do;
			%local set_degrees_freedom;
			%let set_degrees_freedom = %str(nboot-1);
		%end;
		%else %if (&set_degrees_freedom ne %str(nboot-1)) %then %do;
			%put WARNING: &sysmacroname. Using a number of DF different than the number of (non-missing) bootstrap estimates is not recommended;
		%end;
	
		%put # &sysmacroname. Degrees of freedom = &set_degrees_freedom ;
	
	/* ECHO ON: show data manipulation */
	%if (&echo=1) %then %do;
		options mprint;
		options nosource nonotes;
		%put ;
		%put - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -;
		%put # &SYSmacroname.: PARAMETRIC BOOTSTRAP VARIANCE ESTIMATION ;
		%put - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -;
		%put ;
	%end;
	
	/* sort both datasets by byvar, if need be */
		%if (&byvar ne %str()) %then %do; 
			proc sort data=&inboot;
				by &byvar;
			run;
			
			%if (&echo=1) %then %put ;
	
			proc sort data=&inbase;
				by &byvar;
			run;
			
			%if (&echo=1) %then %put ;
		
		%end;
		
		/*********************************************/
		/* 1) output number of non-missing boostraps */
		/*********************************************/
			
			/* 1.1) Get number of missing values */
			proc means data=&inboot %if (&where_estimate ne %str()) %then %do; (where=(&where_estimate))%end; nmiss noprint;
			%if (&byvar ne %str()) %then %do; by &byvar ; %end;
				var &varlist ;
				output out=_nboot(drop=_TYPE_) nmiss=;
			run;

			%if (&echo=1) %then %put ;
		
			/* 1.2) Get number of non-missing bootstrap estimates */
			data _nboot;
				set _nboot;
				array varlist (*) &varlist ;
				* Number of non-missing bootstrap = total obs. (_FREQ_) - nmiss ;
				do i=1 to dim(varlist);
				varlist(i) = _FREQ_ - (varlist(i));
				end;
			run;
			
			%if (&echo=1) %then %put ;
			
			/* 1.3) Transpose nboot values */
			proc transpose data=_nboot out=_nboott(rename=(COL1=nboot)) prefix=COL ;
			%if (&byvar ne %str()) %then %do; by &byvar ; %end;
			var &varlist ;
			run;
		
			%if (&echo=1) %then %put ;
			
		/*********************************************/
		/* 2) calculate bootstrap std                */
		/*********************************************/
		
			proc means data=&inboot %if (&where_estimate ne %str()) %then %do; (where=(&where_estimate))%end; std noprint;
			%if (&byvar ne %str()) %then %do; by &byvar ; %end;
			var &varlist ;
			output out=_std(drop=_TYPE_) std=;
			run;
			
			%if (&echo=1) %then %put ;
						
			/* 2.1) transpose stderr */
			proc transpose data=_std out=_stdt(rename=(COL1=se)) prefix=COL ;
			%if (&byvar ne %str()) %then %do; by &byvar ; %end;
			var &varlist ;
			run;
			
			%if (&echo=1) %then %put ;
			
		/*********************************************/
		/* 3) transpose base estimates               */
		/*********************************************/
		
			proc transpose data=&inbase%if (&where_estimate ne %str()) %then %do; (where=(&where_estimate))%end;  out=_estimate (rename=(estimate1=estimate)) prefix=estimate;
			%if (&byvar ne %str()) %then %do; by &byvar ; %end;			
			var &varlist ;
			run;

			%if (&echo=1) %then %put ;
			
		/*********************************************/
		/* 4) merge base with stderr,
			calculate 95ci, test statistics ; */
		/*********************************************/
		
			data &outdata. (rename=(_NAME_=name));
				merge _estimate _stdt _nboott;
				/* no need for by - both data ordered, unless byvar specified */
				%if (&byvar ne %str()) %then %do; by &byvar ; %end;
				
			%if (&echo=1) %then %put ;
			
			/* Classic scenario, calculating only 1 confidence interval (eg 95ci) */
			%if %sysfunc(countw(&set_confidence_level,' '))=1 %then %do;
			
			* calculate 95ci, pvalue, coefficient of variation;	
			%ci95(estimate               = estimate, /* ### hardcoded name from proc means above */
				  set_confidence_level   = &set_confidence_level.,
				  set_degrees_freedom    = &set_degrees_freedom.,
				  null                   = &null.,
				  ptail                  = 2
				  );
			
			%if (&echo=1) %then %put ;
			
			* svalue - ignores missing and truncates at 50 ;
				if not missing(pvalue) then do;
					if pvalue > 8.88178419700125e-16 then svalue=round(-log2(pvalue),1);
						else svalue=round(-log2(8.88178419700125e-16),1);
				end;
				else svalue=.;
				
				label svalue = "Suprisal value (bits of information against H0)";
				
					/* The S-value is a logarithmic transformation of the p-value that 
						rescales it on an additive scale and tells us how much information is embedded 
						in the test statistic and can be used as evidence against the test hypothesis. 
						
						Rafi, Z. P-Values Are Tough And S-Values Can Help. Less Likely. 
						/statistics/s-values/. Published November 11, 2018.  */
			
			%if (&echo=1) %then %put ;
			* coefficient of variation ;
				format cv PERCENT6.1;
				if not missing(estimate) & (estimate ne 0) then cv = abs(se / estimate);
					else cv = .;
			
				%if (&exp>0) %then %do;
				* exponentiate data ;
				%if (&echo=1) %then %put ;
					array _estimate_ (*) estimate lcl&set_confidence_level ucl&set_confidence_level ;
					do kth=1 to dim(_estimate_);
						_estimate_(kth) = exp(_estimate_(kth));
					end;
					drop kth cv;
					exp=1;
					label exp = "Estimate, LCL and UCL were exponentiated"
						  se = "Bootstrap standard errors (ln scale)";
				%end;

			%if (&echo=1) %then %put ;
			* make formatted estimates ;
				length _ci $ 15 _se $ 8;
				_ci = cats('(',put(lcl&set_confidence_level,&numberfmt.),',',put(ucl&set_confidence_level,&numberfmt.),')');
				_se = cats('(',put(se, &numberfmt.),')');
			
				length estimate_ci $ 24 estimate_se $ 18;
				estimate_ci = cat(put(estimate,&numberfmt.),' ',_ci);
				estimate_se = cat(put(estimate,&numberfmt.),' ',_se);
			
				%if (&exp>0) %then %do;
				* se are expressed on ln scale, drop to avoid confusion ;
				%if (&echo=1) %then %put ;
					drop estimate_se _se ;
				%end;

			%end; /* end of conditional on classical scenario */
			%else %if %sysfunc(countw(&set_confidence_level,' '))>1 %then %do; /* special scenario, multiple ci levels */
				
				%do ith=1 %to %sysfunc(countw(&set_confidence_level));
					%let set_confidence_level&ith = %scan(&set_confidence_level,&ith,%str( ));
					
					%if (&echo=1) %then %put ;
					%if (&echo=0) %then %put # &sysmacroname. #&ith Estimating data for confidence level &&set_confidence_level&ith;
					%else %put %str(* #&ith Estimating data for confidence level &&set_confidence_level&ith ;);
					
					%ci95(estimate               = estimate,
						  set_confidence_level   = &&set_confidence_level&ith,
						  set_degrees_freedom    = &set_degrees_freedom.,
						  null                   = &null.,
						  ptail                  = 2
						  );
						  
						/* ### of note, macro is named ci95, but can estimate
							other confidence level than the classical alpha=0.05 */
							
					%if (&echo=1) %then %put ;
					pvalue&&set_confidence_level&ith = pvalue;
					tvalue&&set_confidence_level&ith = tvalue;
					t_quant&&set_confidence_level&ith = t_quant;

					label pvalue&&set_confidence_level&ith = "2-sided pvalue (null H value=&null.)";
					
				%end; /* end of confidence level loop */
				
				%if (&echo=1) %then %put ;
				
				drop alpha one_minus_half_alpha pvalue: tvalue: t_quant;
					
					/* of note, pvalue and tvalue are not corrected for different alpha and hence, cant be used*/
				
			%end; /* end of conditional on > 1 confidence level */

			%if (&echo=1) %then %put ;
			* clear ugly labels;
				label _NAME_=" " nboot =" ";
			run;
			
	/* ECHO OFF */
	%if (&echo=1) %then %do;
		options nomprint;
		options source notes;
	%end;
	
	/*********************************************/
	/*          clean temporary data             */
	/*********************************************/
	
		proc datasets lib=work nolist nodetails;
		delete _nboot _nboott ;
		run;
	
	%mend boot_variance;
	
	/* example */
	/*
	%boot_variance(inboot         = bootstrap,
				   inbase         = base,
				   where_estimate = ,
				   outdata        = boot_energy,
				   varlist        =intercept rc_beta_unit adjr2 ,
				   set_degrees_freedom = 300 
				   );
				   */

