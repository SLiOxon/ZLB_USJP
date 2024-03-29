/**	Sensitivity: Tests with alternative measures of long term rates and thus shorter sample periods
	Tests ZLB irrelevance hypothesis in Japan by testing CSVAR against CKSVAR
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

	decl fUseKSVARinit = TRUE;					// if TRUE, estimate KSVAR first and use as initial values for estimation of CKSVAR

	decl p = 2;
	decl b = 0.07;					// lower bound = IOR+0.07 for JP
	decl mresults = constant(.NaN,1,8);	// store p, lik, cp, pv-p, aic, dlr, df, pval
	
	///////////////////////////////////////////////////////////////////////
	//	Setup the model	
		decl model = new CKSVAR();					// create object, holds all functions needed to estimated the model
		model.SetBound(b);
		model.SetNumMCreps(R);
		model.SetLatentLags(fLatLags);
		decl sdata = sprint("JPNdata.csv");
		model.LoadCsv(sdata);	// load data
		model.Deterministic(FALSE);					// create constant term in database (change to "TRUE" to create seasonal dummies, too)
		model.Select(CKSVAR::Y1_VAR, {"INFL", 0, p, "YGAP_BOJ", 0, p, "LR10Y", 0, p});	  // Select Y1 (unconstrained) variables in the SVAR
		decl sSHORTRATE = "SR";
		model.Select(CKSVAR::Y2_VAR, {sSHORTRATE, 0, p});	   // Selects the variable that is subject to a lower bound constraint, Y2 (here, the call rate).
		model.Select(CKSVAR::X_VAR, {"Constant", 0, 0, "PGR", 0, 0});	// Select any other deterministic or exogenous regressors
		model.Select(CKSVAR::Y2_BOUNDS, {"IOR", 0, 0}); // Variable ELB as in Hayashi and Koeda
			
		decl iYear0 = 1987;
		model.SetSelSample(/*Start Year*/iYear0,/*Start Period*/ 4, /*End Year*/ -1, /*End Period*/ -1);		// -1 defaults to earliest (for start) or latest (for end) available dates in the sample 

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
		mresults[0][:2] = p~model.GetLogLik()~model.GetParCount();
		mresults[0][4] = aic;
		println("\nlapsed time: ", timespan(time),"\n");
	    model.SetUnrestrictedModel(clone(model));		// store settings of unrestricted CKSVAR model, to test against

	//  Impose CSVAR restrictions and estimate	
		decl asNamesCK = model.GetParNames();
		vp = model.GetPar();							// store final estimates
		model.SetCSVAR(TRUE);							// impose CSVAR restrictions
		model.InitPar();					
		model.SetStartPar(vp[find(asNamesCK,model.GetParNames())]);		// use estimates of CKSVAR as initial conditions
		
		println("\n*****   Started estimation of the restricted model");	
		model.Estimate();
		decl atest = model.TestAgainstUnrestr();
		if(atest)
			mresults[0][5:] = atest[0]~atest[1]~atest[2];
		delete model;
		println("\nlapsed time: ", timespan(time),"\n");
		decl asnames = {"p", "lik", "par","pv-p", "aic", "LR", "df", "pval"};
		println("Results", "%c", asnames, mresults);
}