/**	Tests no attenuation bias
	The full model is CKSVAR
**/

#include <oxstd.oxh>
#include <oxdraw.oxh>
#include <oxprob.oxh>
#import <maximize>
#include <cksvar.ox>

main()
{
	decl R = 1000;							// # of MC simulations for importance sampling to integrate out latent lags
	decl fLatLags = TRUE;						// if TRUE (default), model has lags of latent variable Y2* (shadow rate) on the RHS
	
	for(decl iJP = 0; iJP < 2; ++iJP){			// do US and JP

	decl fUseKSVARinit = iJP ? FALSE : TRUE;	// if TRUE, estimate KSVAR first and use as initial values for estimation of CKSVAR. KSVAR initialisation performs better than without it for US data but worse for Japan

	decl pmax = iJP ? 2 : 3;
	decl pmin = iJP ? 2 : 3;;
	decl b = iJP ? 0.07 : 0.2;					// lower bound = IOR+0.07 for JP, 0.2 for US
	decl mresults = constant(.NaN, pmax-pmin+1,8);	// store p, lik, cp, pv-p, aic, dlr, df, pval

	decl lrs_us = {"LONGRATE", "LR7Y", "LR5Y", "LR3Y", "LR2Y", "LR1Y"};
	decl slrs_us = sizeof(lrs_us);
	decl lrs_jp = {"LR9Y", "LR7Y", "LR5Y", "LR3Y", "LR2Y", "LR1Y"};
	decl slrs_jp = sizeof(lrs_jp);
	
	for(decl ilr = 0; ilr < (iJP? slrs_jp : slrs_us); ++ilr){

	for(decl ip = 0; ip < pmax-pmin+1; ++ip){									// VAR order
		decl p = pmax-ip;
		if(iJP == 0){
		 	if(ilr == 0 || ilr == 4)
				p = 4;
		 }
	
	///////////////////////////////////////////////////////////////////////
	//	Setup the model	
		decl model = new CKSVAR();					// create object, holds all functions needed to estimated the model
		model.SetBound(b);
		model.SetNumMCreps(R);
		model.SetLatentLags(fLatLags);
		decl sdata = sprint(iJP ? "JPN" : "US", "data.csv");
		model.LoadCsv(sdata);	// load data
		model.Deterministic(FALSE);					// create constant term in database (change to "TRUE" to create seasonal dummies, too)
		model.Select(CKSVAR::Y1_VAR, {"INFL", 0, p, iJP ? "YGAP_BOJ" : "YGAP", 0, p, iJP? lrs_jp[ilr] : lrs_us[ilr], 0, p});	  // Select Y1 (unconstrained) variables in the SVAR
		decl sSHORTRATE = iJP ? "SR" : "FEDFUNDS";
		model.Select(CKSVAR::Y2_VAR, {sSHORTRATE, 0, p});	   // Selects the variable that is subject to a lower bound constraint, Y2 (short term interest rate).
		model.Select(CKSVAR::X_VAR, {"Constant", 0, 0});	// Select any other deterministic or exogenous regressors
		if(iJP){
			model.Select(CKSVAR::X_VAR, {"PGR", 0, 0});
			model.Select(CKSVAR::Y2_BOUNDS, {"IOR", 0, 0}); // Variable ELB as in Hayashi and Koeda
		}
			
		decl iYear0 = iJP ? 1985 : 1960;
		decl iQuarter0 = iJP ? 3 : 1;
		decl iYear1 = 2019;
		decl iQuarter1 = iJP ? 1 : 4;
		if (iJP == 0){
		switch(ilr){
			case 0:
				iYear0 = 1960;
				iQuarter0 = 1;
				break;
			case 1:
				iYear0 = 1970;
				iQuarter0 = 4;
				break;
			case 4:
				iYear0 = 1977;
				iQuarter0 = 4;
				break;
			default:
				iYear0 = 1963;
				iQuarter0 = 2;
				break;
		}
		}
		model.SetSelSample(/*Start Year*/iYear0,/*Start Period*/ iQuarter0, /*End Year*/ 2019, /*End Period*/ 1);		// -1 defaults to earliest (for start) or latest (for end) available dates in the sample 

		ranseed(-1);
		
	//  Estimate unsrestricted model and record results
		decl vp;
		if(fUseKSVARinit){
			decl ksvar = clone(model);
			ksvar.SetLatentLags(FALSE);
			ksvar.SetPrint(FALSE);
			ksvar.Estimate();
			decl vpK = ksvar.GetPar();
			decl asNamesK = ksvar.GetParNames();
			delete ksvar;
			model.InitPar();
			vp = model.GetPar();
			vp[find(model.GetParNames(),asNamesK)] = vpK;			
			model.SetStartPar(vp);				// initialize at KSVAR estimates
		}
	
		decl time = timer();
			
		model.Estimate();								// estimates the model by ML and stores results inside the model object using BFGS
		decl aic = (-2 * model.GetLogLik() + 2 * model.GetcDfLoss())/model.GetcT();
		mresults[ip][:2] = p~model.GetLogLik()~model.GetParCount();
		if(ip)
			mresults[ip][3] = tailchi(2*(mresults[0][1]-mresults[ip][1]), mresults[0][2]-mresults[ip][2]);		// p-value of LR test of SVAR(p) against SVAR(pmax)
		mresults[ip][4] = aic;
		println("\nlapsed time: ", timespan(time),"\n");
	    model.SetUnrestrictedModel(clone(model));		// store settings of unrestricted CKSVAR model, to test against

	//  Impose exclusion restrictions and estimate
		decl cY1 = columns(model.GetGroupLag(CKSVAR::Y1_VAR, 0,0));
//		decl m_asbeta = new array[2];				// names of beta parameters in Y1 equation(s) 
//		for(decl j = 0; j < sizeof(m_asbeta); ++j)
		decl m_asbeta = sprint("$\\tilde{\\beta}_3$");		// index refers to equation
		decl vidx = strfind(model.GetParNames(), m_asbeta);			// locate the betas within the parameter vector
		for(decl i = 0; i < sizerc(vidx);++i)
			model.FixPar(vidx[i], 0);
//		if(p){
//			decl vIdxExcludeLags = <>;
//			decl asY2lags;
//			model.GetGroupLagNames(CKSVAR::Y2_VAR, 1, p, &asY2lags);
//			for(decl i = 0; i < cY1-1; ++i){
//				for(decl j = 0; j < sizeof(asY2lags); ++j)
//					vIdxExcludeLags |= find(model.GetParNames(), sprint("Eq.", i+1, " ", asY2lags[j]));
//				for(decl k = 0; k < sizeof(asY2lags); ++k)
//					vIdxExcludeLags |= find(model.GetParNames(), sprint("Eq.", i+1, " ", "l", asY2lags[k]));
//			}
//			model.FixPar(vIdxExcludeLags, zeros(sizerc(vIdxExcludeLags), 1));
//		}
		println("\n*****   Started estimation of the restricted model");	
		model.Estimate();
		decl atest = model.TestAgainstUnrestr();
		if(atest)
			mresults[ip][5:] = atest[0]~atest[1]~atest[2];
		delete model;
		println("\nlapsed time: ", timespan(time),"\n");
		}
		decl asnames = {"p", "lik", "par","pv-p", "aic", "LR", "df", "pval"};
		println("Results", "%c", asnames, mresults);
	}
	}
}

