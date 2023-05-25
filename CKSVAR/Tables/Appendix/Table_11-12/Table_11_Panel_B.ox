/**	Tests ZLB irrelevance hypothesis in US and Japan by testing exclusion of short rate (and no kink) in Y1 equations
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
	decl fUseKSVARinit = TRUE;
	
	decl p = 1;
	decl b = 0.002481393;					// lower bound = IOR+0.07 for JP, 0.2 for US
	decl v_xi = <70,75,80,85,90,95,99>;		// possible value of xi
	decl mcreject = constant(.NaN,7,4);
	decl time0 = timer();

	for(decl indxi = 0; indxi < 7; ++indxi){		// for each value of xi

	decl mresults = constant(.NaN,100,6);	// store inds, lik, aic, dlr, df, pval
	decl sxi = v_xi[indxi];
	decl cReject1 = 0;			// rejection count with 0.01 significance
	decl cReject2 = 0;			// rejection count with 0.05 significance
	decl cReject3 = 0;			// rejection count with 0.1 significance
	println("Starting estimation for xi=", sxi/100);

	for(decl inds = 1; inds < 101; ++inds){		// for each MC simulation of dataset (total of 100 datasets for each xi)
	
	///////////////////////////////////////////////////////////////////////
	//	Setup the model	
		decl model = new CKSVAR();					// create object, holds all functions needed to estimated the model
		model.SetBound(b);
		model.SetNumMCreps(R);
		model.SetLatentLags(fLatLags);
		decl sdata = sprint("Data_l", sxi, ".csv");
		model.LoadCsv(sdata);	// load data
		model.Deterministic(FALSE);					// create constant term in database (change to "TRUE" to create seasonal dummies, too)
		model.Select(CKSVAR::Y1_VAR, {sprint("INFL", inds), 0, p, sprint("YGAP", inds), 0, p, sprint("LR",inds), 0, p});	  // Select Y1 (unconstrained) variables in the SVAR
		model.Select(CKSVAR::Y2_VAR, {sprint("RATE", inds), 0, p});	   // Selects the variable that is subject to a lower bound constraint, Y2 (short term interest rate).
		model.Select(CKSVAR::X_VAR, {"Constant", 0, 0});	// Select any other deterministic or exogenous regressors
		
		model.SetSelSample(/*Start Year*/-1,/*Start Period*/-1, /*End Year*/ -1, /*End Period*/ -1);		// -1 defaults to earliest (for start) or latest (for end) available dates in the sample 

//		model.SetFAPF(TRUE);
//		model.SetSimAnn(TRUE, {/*iNS*/20,/*iNT*/5,/*dRt*/.5,/*vM*/1,/*vC*/2,/*dT*/5});
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
		mresults[inds-1][:1] = inds~model.GetLogLik();
		mresults[inds-1][2] = aic;
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
				for(decl k = 0; k < sizeof(asY2lags); ++k)
					vIdxExcludeLags |= find(model.GetParNames(), sprint("Eq.", i+1, " ", "l", asY2lags[k]));
			}
			model.FixPar(vIdxExcludeLags, zeros(sizerc(vIdxExcludeLags), 1));
		}
		println("\n*****   Started estimation of the restricted model");	
		model.Estimate();
		decl atest = model.TestAgainstUnrestr();
		if(atest)
			mresults[inds-1][3:] = atest[0]~atest[1]~atest[2];

		cReject1 = cReject1 + (atest[2] <= 0.01);	// count of rejections with 0.01 significance level
		cReject2 = cReject2 + (atest[2] <= 0.05);	// count of rejections with 0.05 significance level
		cReject3 = cReject3 + (atest[2] <= 0.1);	// count of rejections with 0.1 significance level		
		delete model;
		println("\nlapsed time: ", timespan(time),"\n");
		}

//		decl asnames = {"simulation", "lik", "aic", "LR", "df", "pval"};
//		println("Results", "%c", asnames, mresults);
		println("No. of Rejections for xi = ", sxi/100, " with 0.01 significance level: ", cReject1);
		println("No. of Rejections for xi = ", sxi/100, " with 0.05 significance level: ", cReject2);
		println("No. of Rejections for xi = ", sxi/100, " with 0.1 significance level: ", cReject3);
		mcreject[indxi][] = sxi/100 ~ cReject1 ~ cReject2 ~ cReject3;
		println("\nlapsed time: ", timespan(time0),"\n");

	}

	decl asnames2 = {"xi", "p <= 0.01", "p <= 0.05", "p <= 0.1"}; 
	println("Results", "%c", asnames2, mcreject);
	println("\nlapsed time: ", timespan(time0),"\n");
}

