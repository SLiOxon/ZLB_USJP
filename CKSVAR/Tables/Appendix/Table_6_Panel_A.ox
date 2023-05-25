/**	Tests ZLB irrelevance hypothesis in US and Japan by testing exclusion of short rate (and no kink) in Y1 equations, with credit spread
	The full model is KSVAR
**/

#include <oxstd.oxh>
#include <oxdraw.oxh>
#include <oxprob.oxh>
#import <maximize>
#include <cksvar.ox>

main()
{
	decl R = 1000;							// # of MC simulations for importance sampling to integrate out latent lags
	decl fLatLags = FALSE;						// if TRUE (default), model has lags of latent variable Y2* (shadow rate) on the RHS
	
	for(decl iJP = 0; iJP < 2; ++iJP){			// do US and JP

	decl pmax = iJP ? 5 : 5;
	decl pmin = 1;
	decl b = iJP ? 0.07 : 0.2;					// lower bound = IOR+0.07 for JP, 0.2 for US
	decl mresults = constant(.NaN, pmax-pmin+1,8);	// store p, lik, cp, pv-p, aic, dlr, df, pval

	decl year0idx; // record the latest starting period for the longest lags

	for(decl ip = 0; ip < pmax-pmin+1; ++ip){									// VAR order
		decl p = pmax-ip;
	
	///////////////////////////////////////////////////////////////////////
	//	Setup the model	
		decl model = new CKSVAR();					// create object, holds all functions needed to estimated the model
		model.SetBound(b);
		model.SetNumMCreps(R);
		model.SetLatentLags(fLatLags);
		decl sdata = sprint(iJP ? "JPN" : "US", "data.csv");
		model.LoadCsv(sdata);	// load data
		model.Deterministic(FALSE);					// create constant term in database
		model.Select(CKSVAR::Y1_VAR, {"INFL", 0, p, iJP ? "YGAP_BOJ" : "YGAP", 0, p, iJP? "LR9Y" : "LONGRATE", 0, p, iJP ? "EBP" : "BAA10Y", 0, p});	  // Select Y1 (unconstrained) variables in the SVAR
		decl sSHORTRATE = iJP ? "SR" : "FEDFUNDS";
		model.Select(CKSVAR::Y2_VAR, {sSHORTRATE, 0, p});	   // Selects the variable that is subject to a lower bound constraint, Y2 (short term interest rate).
		model.Select(CKSVAR::X_VAR, {"Constant", 0, 0});	// Select any other deterministic or exogenous regressors
		if(iJP){
			model.Select(CKSVAR::X_VAR, {"PGR", 0, 0});	
			model.Select(CKSVAR::Y2_BOUNDS, {"IOR", 0, 0}); // Variable ELB as in Hayashi and Koeda
		}
		
		decl iYear0 = iJP ? 1985 : 1960;
		decl iQuarter0 = iJP ? 3 : 1;
		model.SetSelSample(/*Start Year*/iYear0,/*Start Period*/ iQuarter0, /*End Year*/ 2019, /*End Period*/ 1);		// -1 defaults to earliest (for start) or latest (for end) available dates in the sample 
		if(ip == 0)	
			year0idx = model.GetSelStart();
		model.SetSelSampleByIndex(year0idx, model.GetSelEnd());  // if the sample for given money supply is availabla later than 1960Q1, keep selected sample fixed
		
		ranseed(-1);
		
		model.InitData();
		model.InitPar();
	
		decl time = timer();
			
	//  Estimate unsrestricted model and record results
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
		model.Fixbeta(zeros(cY1,1));
		if(p){
			decl vIdxExcludeLags = <>;
			decl asY2lags;
			model.GetGroupLagNames(CKSVAR::Y2_VAR, 1, p, &asY2lags);
			for(decl i = 0; i < cY1; ++i){
				for(decl j = 0; j < sizeof(asY2lags); ++j)
					vIdxExcludeLags |= find(model.GetParNames(), sprint("Eq.", i+1, " ", asY2lags[j]));
			}
			model.FixPar(vIdxExcludeLags, zeros(sizerc(vIdxExcludeLags), 1));
		}
		println("\n*****   Started estimation of the restricted model");	
		model.Estimate();
		decl atest = model.TestAgainstUnrestr();
		if(atest)
			mresults[ip][5:] = atest[0]~atest[1]~atest[2];
		//mresults_grand[ilr][][] = mresults;
		delete model;
		println("\nlapsed time: ", timespan(time),"\n");
		}
		decl asnames = {"p", "lik", "par","pv-p", "aic", "LR", "df", "pval"};
		println("Results", "%c", asnames, mresults);
	}
}

