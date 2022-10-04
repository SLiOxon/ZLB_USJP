//	Estimates CKSVAR models for Japan and USA and computes Choleski-idnetified IRFs


#include <oxstd.oxh>
#include <oxdraw.oxh>
#include <oxprob.oxh>
#import <maximize>
#include <cksvar.ox>

class CKSVARshadow : CKSVAR
{
	CKSVARshadow();
	GetShadow(const vY2bar, const xi);
	GetIRFbounds(const dcoverage, const ximin, const ximax);
	funcIRF0(const avIRF, const vP);
	fcritIM(const avF, const vX);
	fJacoIM(const avF, const vX);
	decl m_dterm1, m_dterm2;		// constants in equation of Imbens and Manski critical value
};
CKSVARshadow::CKSVARshadow()
{
	CKSVAR();
}
CKSVARshadow::fcritIM(const avF, const vX)
/** equation that determines the critical value in Imbens and Manski (2004)
**/
{
	avF[0] = probn(vX+m_dterm1) - probn(-vX) - m_dterm2;
	return TRUE;
}
CKSVARshadow::fJacoIM(const avF, const vX)
/** Jacobian of equation for Imbens and Manski critical value
**/
{
	avF[0] = densn(vX+m_dterm1) + densn(vX);
	return TRUE;
}
CKSVARshadow::GetShadow(const vY2bar, const xi)
{
	decl vpall;
	if(!Reparameterize(&vpall, this.GetFreePar()))						// obtain full parameter vector
		return 0;
	decl tau, vC2obs, vC2lat, vbetatilde, mC1obs, mC1lat, vdelta, mOm1g2Ch;			// obtain the parameters first
	[tau, vC2obs, vC2lat, vbetatilde, mC1obs, mC1lat, vdelta, mOm1g2Ch] = this.GetModelPar(vpall);
	decl mOmega;

	decl mOm1g2 = mOm1g2Ch*mOm1g2Ch';					// variance of u_1|u_2
	decl mOm11 = mOm1g2 + tau^2*vdelta*vdelta';	// unconditional variance matrix for Y1
	mOmega = mOm11 ~ vdelta*tau^2
			| vdelta'tau^2 ~ tau^2;			// variance matrix of reduced form errors for all Y
		
	decl betabar = CKSVAR::GetBetabar(vbetatilde, xi, mOmega);
	betabar = betabar[][0];
	decl gamma = CKSVAR::GetGamma(betabar, mOmega);
	decl arg = (1-gamma*betabar)/(1-xi*gamma*betabar);
	return vY2bar .> m_vb .? vY2bar .: vY2bar*arg + (1-arg)*m_vb;
}
CKSVARshadow::funcIRF0(const avIRF, const vP)
/** like funcIRF, but dropping the ones corresponding to second betabar solution when there are two
**/
{
	decl mIRFs, ires;
	ires = this.funcIRF(&mIRFs, vP);
	mIRFs = shape(mIRFs,m_iSol*(m_cY1+1), m_cHorizon+1);
	if(m_iSol>1)
		mIRFs = mIRFs[:m_cY1][];			// drop second solution
	avIRF[0] = vec(mIRFs);
	return ires;
}
CKSVARshadow::GetIRFbounds(const dcoverage, const ximin, const ximax)
/** Computes error bands on the bounds of Impulse Response function to structural shock to Y2*
	@param dcoverage, double, in (0,1).
	@param ximin, double, minimum value of xi
	@param ximax, double, maximum value of xi
	returns array of [m_cY][m_cHorizon+1] matrices of lower and upper confidence bands for the IRF
**/
{
	decl lambda = m_dlambda;
	decl aresult = new array[2];
	decl msterr, msterrhi = <>, msterrlo = <>;				// standard errors corresponding to boundaries of IRFs
	decl amIRFs = new array[2];
	for(decl i = 0; i < 2; ++i){
		// set xi at the specified end points
		if(i == 0)			
			this.Setlambda(ximin);
		else
			this.Setlambda(ximax);
		decl mIRFs;
		funcIRF0(&mIRFs, this.GetFreePar());
		mIRFs = shape(mIRFs,(m_cY1+1), m_cHorizon+1);
		amIRFs[i] = mIRFs;
		// compute error bands at each endpoint
		if(!m_fBootstrap){			// if you don't have bootstrapped samples, use delta method (prior estimation suggests Delta is accurate for Japan, i.e., ~ bootstrap, but not for US -- too big s.e.s)
			decl vJac;
			decl ires = NumJacobian(funcIRF0, this.GetFreePar(), &vJac);
			decl mIRFvar = vJac*this.GetCovar()*vJac';			// variance matrix of IRFs by the delta method
			msterr = ires ? shape(sqrt(diagonal(mIRFvar)), m_cY1+1, m_cHorizon+1) : constant(.NaN, m_cY1+1, m_cHorizon+1);	// to hold matrix of asymptotic standard errors
		}
		if(m_fBootstrap){
			decl mBootPar = m_mBootData[][m_fTestUnrestr+range(1,this.GetFreeParCount())];		// retrieve bootstrapped parameter estimates
			decl amIRF = new array[(m_cY1+1)];													// array to hold results per variable
			for(decl j = 0; j < sizeof(amIRF); ++j)
				amIRF[j] = constant(.NaN, sizeof(mBootPar),m_cHorizon+1);	// bootstrap reps in rows for IRF of var i at horizon in columns. 
			for(decl j = 0; j < sizeof(mBootPar); ++j){
				decl vParBoot = mBootPar[j][]';								// extract parameter estimates from jth bootstrap replication
				decl aroots = this.GetCompanionRoots(vParBoot);
				if(any(cabs(aroots[0]).>1) || any(cabs(aroots[1]).>1))		// if explosive roots, skip
					continue;
				decl arg;
				if(!funcIRF0(&arg, vParBoot))								// if IRF computation fails, skip
					continue;									
				arg = shape(arg,(m_cY1+1), m_cHorizon+1);
				for(decl h = 0; h < sizeof(amIRF); ++h)
					amIRF[h][j][] = arg[h][];
			}		
			msterr = new matrix[m_cY1+1][m_cHorizon+1];
			for(decl h = 0; h < sizeof(amIRF); ++h)
				msterr[h][] = sqrt(varc(deleter(amIRF[h])));
		}
		if(i == 0)
			msterrhi = msterrlo = msterr; 
		else{
			msterrhi = mIRFs .> amIRFs[0] .? msterr .: msterrhi;
			msterrlo = mIRFs .< amIRFs[0] .? msterr .: msterrlo;
		}
	}
	decl mIRFlo = msterr, mIRFup = msterr;					// initialize dimensions
	for(decl i = 0; i < m_cY1+1; ++i){
		for(decl j = 0; j < m_cHorizon+1; ++j){
			decl Delta = fabs(amIRFs[1][i][j]-amIRFs[0][i][j]);			// Imbens Manski p. 1850
			m_dterm1 = Delta/max(msterrhi[i][j],msterrlo[i][j]);
			m_dterm2 = dcoverage;
			decl dcritIM=quann(dcoverage); 							// initialize at z_alpha
			SolveNLE(this.fcritIM, &dcritIM, -1, fJacoIM);			// compute Imbens-Manski cv
			mIRFlo[i][j] = min(amIRFs[0][i][j],amIRFs[1][i][j])-dcritIM*msterrlo[i][j];
			mIRFup[i][j] = max(amIRFs[0][i][j],amIRFs[1][i][j])+dcritIM*msterrhi[i][j];
		}
	}
	aresult = {mIRFlo,mIRFup};
	this.Setlambda(lambda);
	return aresult;
}

main()
{
	decl time = timer();
	decl R = 1000;								// # of MC simulations for importance sampling to integrate out latent lags
	decl fLatLags = TRUE;						// if TRUE (default), model has lags of latent variable Y2* (shadow rate) on the RHS
	decl fUseKSVARinit = TRUE;					// if TRUE, estimate KSVAR first and use as initial values for estimation of CKSVAR
	decl fboot = FALSE;
	
	for(decl iJP = 0; iJP < 2; ++iJP){			// do US and JP
  		decl p = iJP ? 2 : 3;

		decl sbootfile = sprint("cksvar_", iJP ? "JP_" : "US_", "p_",p,"_R_",R, "_chol_bootstrap");		// file name where bootstrap replications are to be stored (or have been stored already)
		decl sestimfile = sprint("cksvar_", iJP ? "JP_" : "US_", "p_",p,"_R_",R, "_chol");

	///////////////////////////////////////////////////////////////////////
///////	Setup the model	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		decl model = new CKSVARshadow();					// create object, holds all functions needed to estimated the model
		decl b = iJP ? 0.07 : 0.2;					// lower bound = IOR + 0.07 for JP, 0.2 for US
		model.SetBound(b);
		model.SetNumMCreps(R);
		model.SetLatentLags(fLatLags);
		decl sdata = sprint(iJP ? "JPN" : "US", "data.csv");
		model.LoadCsv(sdata);	// load data
		model.Deterministic(FALSE);					// create constant term in database (change to "TRUE" to create seasonal dummies, too)
		model.Select(CKSVAR::Y1_VAR, {"INFL", 0, p, iJP ? "YGAP_BOJ" : "YGAP", 0, p});	  // Select Y1 (unconstrained) variables in the SVAR
		decl sSHORTRATE = iJP ? "SR" : "FEDFUNDS";
		model.Select(CKSVAR::Y2_VAR, {sSHORTRATE, 0, p});	   // Selects the variable that is subject to a lower bound constraint, Y2 (here, the short term interest rate).
		model.Select(CKSVAR::X_VAR, {"Constant", 0, 0});	// Select any other deterministic or exogenous regressors
		if(iJP){
			model.Select(CKSVAR::X_VAR, {"PGR", 0, 0});	// Select any other deterministic or exogenous regressors
			model.Select(CKSVAR::Y2_BOUNDS, {"IOR", 0, 0}); // Variable ELB as in Hayashi and Koeda
		}
			
		decl iYear0 = iJP ? 1985 : 1960;
		decl iQuarter0 = iJP ? 3 : 1;
		model.SetSelSample(/*Start Year*/iYear0,/*Start Period*/ iQuarter0, /*End Year*/ 2019, /*End Period*/ 1);		// -1 defaults to earliest (for start) or latest (for end) available dates in the sample 

		decl asY1;									// array of strings, to hold variable names
		model.GetGroupLagNames(CKSVAR::Y1_VAR,0,0,&asY1);			// get Y1 variable names
		
		ranseed(-1);
		
		decl vp;
		decl fusesaved = FALSE;
		if(fLatLags){
			model.InitData();
			decl festim = fopen(sprint(sbootfile),"v");
			if(!isfile(festim))
				festim = fopen(sprint(sestimfile),"v");
		
			decl asEstimNames = 0;		// to store parameter names from previous estimation
			if(isfile(festim)){
					fscan(festim, "%v", &model);
			        asEstimNames = model.GetParNames();		
					model.Output();
					fclose(festim);
            }			   
			if(!any(find(model.GetParNames(), asEstimNames).==-1)){
				    fusesaved = TRUE;
			}
			
			else{
				if(fUseKSVARinit){
					decl ksvar = clone(model);
					ksvar.SetLatentLags(FALSE);
					ksvar.SetPrint(FALSE);
					ksvar.Fixbeta(zeros(sizeof(asY1),1));				// impose recursive ID (policy rate ordered last)
					ksvar.Estimate();
					decl vpK = ksvar.GetPar();
					decl asNamesK = ksvar.GetParNames();
					delete ksvar;
					model.InitPar();
					vp = model.GetPar();
					vp[find(model.GetParNames(),asNamesK)] = vpK;			
					model.SetStartPar(vp);				// initialize at KSVAR estimates
				}
			}
		}
		decl time = timer();

		if(!fusesaved){
			model.Fixbeta(zeros(sizeof(asY1),1));				// impose recursive ID (policy rate ordered last)
			model.Estimate();								// estimates the model by ML and stores results inside the model object using BFGS
			decl fsaved = fopen(sestimfile, "w");
	        fprintln(fsaved, "%v", model);
	        fclose(fsaved);
		}

		if(fboot){
//         Perform bootstrap and store bootstrap replications, for computing CIs later
	       model.SetBootstrap(/*Boot reps*/999, /*InitTrue*/TRUE,/*#parallel*/15,/*StoreIntermediate*/TRUE);	   // I use around 15 parallel computations because it seems to be the fastest (more slows it down)
	   	   model.SetBootstrapSeed(13);
		   model.Bootstrap();
	       decl fo = fopen(sbootfile, "w");
	       fprintln(fo, "%v", model);
	       fclose(fo);
		   println("\nTime to complete Bootstraps ", timespan(time),"\n");	
	     } 
   


////////////////////////////////////////////////////////////////////////////////////////////////////
///////	Settings for IRF /////////////////////////////////////////////////////////////////////////////

		decl asY = {"Inflation", "Output gap", iJP ? "Call rate" : "Fed Funds rate"};
		
		decl t0, iYear, iquarter, iT_irf, cIRFreps;
		decl cHorizon = 20;									// # of periods for plotting IRFs
		cIRFreps = 1000;									// # of replications for MC integration wrt to other shocks to compute IRF. Use 0 to set other shocks to zero.
		t0 = model.GetSelStart();							// index over entire db of first obs in estimation sample
		iYear = 1985;										// starting year
		iquarter = 3;								// starting quarter
		iT_irf = (iYear - iYear0)*4+iquarter-iQuarter0;		// origin of IRF expressed as observation index 
		decl SRper = range(0,3);					// range of periods over which you impose sign restrictions
		
///////// get IRFs /////////////////////////////////////////////////////////////////////////////////////////////

		println("Now computing IRFs at all dates...");
		model.Setlambda(0.72);
		decl vxi = range(0,0,0.01);
		decl amIRF, amIRFlo, amIRFhi, aamIRF = {}, aamIRFlo = {}, aamIRFhi = {};
		// first collect all IRFs and sign-restrictions across all dates
		for(decl iori = 0; iori < model.GetSelEnd()-iT_irf-t0+1; ++iori){		    
			model.SetIRF(/*horizon*/cHorizon, /*dImpulse*/-0.25, /*#MC reps*/cIRFreps, /*Origin t*/iT_irf+iori);
			amIRF = model.GetIRFset(vxi, &vxi);			// computes IRFs over range of values for xi 
 			aamIRF ~= {amIRF};									// store IRFs to avoid repeating later				
			if(iJP){
				  aamIRF[iori][0][] = 4*aamIRF[iori][0][];		// rescale the IRFs for inflation as annualised	for Japanese data
			}
		}

///////	plot IRFs for particular origins ////////////////////////////////////////////////////////////////////////////////////////////////////////////		
		decl iYear_zlb = iJP ? 2010 : 2009;
		decl iquarter_zlb = 1;
		decl iYear_nonzlb = iJP ? 1990 : 1999;
		decl iquarter_nonzlb = 1;
		decl iT_zlb = (iYear_zlb - iYear)*4+iquarter_zlb-iquarter;
		decl iT_nonzlb = (iYear_nonzlb - iYear)*4+iquarter_nonzlb-iquarter;
		decl sirf_origin_zlb = sprint(model.ObsYear(iT_zlb+iT_irf+t0), "q", model.ObsPeriod(iT_zlb+iT_irf+t0));		// date at which to compute IRF
		decl sirf_origin_nonzlb = sprint(model.ObsYear(iT_nonzlb+iT_irf+t0), "q", model.ObsPeriod(iT_nonzlb+iT_irf+t0));		// date at which to compute IRF

		
		SetDrawWindow(sprint(iJP ? "JP" : "US", " Choleski IRF CKSVAR(",p,")"));
		decl mlimits = new matrix[3][2];
		decl mIRF, mIRFlo, mIRFup, dcoverage;
		dcoverage = 0.67;
		decl mlimits_min = new matrix[3][2];	// record the scale of each origin to adjust in the end
		decl mlimits_max = new matrix[3][2];	// record the scale of each origin to adjust in the end
		for(decl i = 0; i < 2; ++i){
			decl iTreg = i ? iT_zlb :iT_nonzlb;
			decl sorig = i ? sirf_origin_zlb : sirf_origin_nonzlb;
			model.SetIRF(/*horizon*/cHorizon, /*dImpulse*/-0.25, /*#MC reps*/cIRFreps, /*Origin t*/iTreg+iT_irf);
			[aamIRFlo, aamIRFhi] = model.GetIRFbounds(dcoverage, 0, 0);
			if(iJP){
				  aamIRFlo[0][] = 4*aamIRFlo[0][];
				  aamIRFhi[0][] = 4*aamIRFhi[0][];
			}
			
			for(decl j = 0; j < 3; ++j){
				mlimits[j][] = min(aamIRFlo[j][])~max(aamIRFhi[j][]);							// store limits for next plot
				mlimits[j][0] = mlimits[j][0].>0 .? mlimits[j][0]/1.1 .: mlimits[j][0]*1.1;		// set range for plots
				mlimits[j][1] = mlimits[j][1].>0 .? mlimits[j][1]*1.1 .: mlimits[j][1]/1.1;		// set range for plots
				if(i == 1 && j==2)
					mlimits[2][0] = min(-.25, .1);
				if(j==0)
					DrawTitle(2*j+i, sprint("IRF to monetary policy shock in ", sorig));	// column label
				if(i==0)
					DrawText(2*j, asY[j], 0, 0, /*iFontNo*/-1, /*iFontSize*/360, TEXT_YLABEL, 90);		// Y-label
				DrawMatrix(2*j+i, aamIRF[iTreg][j][], "", 0, 1, 0, 2*ones(sizeof(aamIRF[iTreg][j][]),1));					 // IRFs at all other values of xi
				DrawMatrix(2*j+i, aamIRFlo[j][], sprint(100*dcoverage,"\% error bands"), 0, 1, 0, 5);							
				DrawMatrix(2*j+i, aamIRFhi[j][], "", 0, 1, 0, 5);							 
				DrawAxis(2*j+i, TRUE, .NaN, 0, cHorizon, 2, 2, 0, 0);
				DrawLegend(2*j+i, 400,0, j);
				DrawAdjust(ADJ_AREA_Y, 2*j+i, mlimits[j][0], mlimits[j][1]);
			}
			mlimits_min[][i] = mlimits[][0];
			mlimits_max[][i] = mlimits[][1];			
		}
		decl mlimits_draw = new matrix[3][2];
		for(decl j = 0; j < 3; ++j){
				mlimits_draw[j][] = min(mlimits_min[j][])~max(mlimits_max[j][]);							// store limits for next plot
				DrawAdjust(ADJ_AREA_Y, 2*j, mlimits_draw[j][0], mlimits_draw[j][1]);
				DrawAdjust(ADJ_AREA_Y, 2*j+1, mlimits_draw[j][0], mlimits_draw[j][1]);
		}		
		ShowDrawWindow();		
		SaveDrawWindow(sprint("Figures\\", iJP ? sprint("JP_p_",p) : "US","_Chol_IRF_", 100*dcoverage, ".eps"));
		
		
		delete model;

	}
}


