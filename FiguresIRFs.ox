//	Estimates CKSVAR models for Japan and USA and computes set-identified IRFs with sign restrictions


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
	
	for(decl iJP = 0; iJP < 2; ++iJP){			// do US and JP
  		decl p = iJP ? 2 : 3;

		decl sbootfile = sprint("cksvar_", iJP ? "JP_" : "US_", "p_",p,"_R_",R, "_bootstrap");		// file name where bootstrap replications are to be stored (or have been stored already)
		decl sestimfile = sprint("cksvar_", iJP ? "JP_" : "US_", "p_",p,"_R_",R);

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
		model.SetSelSample(/*Start Year*/iYear0,/*Start Period*/ iQuarter0, /*End Year*/ -1, /*End Period*/ -1);		// -1 defaults to earliest (for start) or latest (for end) available dates in the sample 

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
			model.Estimate();								// estimates the model by ML and stores results inside the model object using BFGS
			decl fsaved = fopen(sestimfile, "w");
	        fprintln(fsaved, "%v", model);
	        fclose(fsaved);
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
		
///////// get the narrowest identified IRF set and the corresponding IRFs /////////////////////////////////////////////////////////////////////////////////////////////
		decl dxi = 0.01, vxi = range(0,1,dxi);
		decl vxirange, vxirestrange;	    				// hold identified range of xi
		decl amIRF;
		vxirestrange = vxi;
		for(decl iori = 0; iori < model.GetSelEnd()-iT_irf-t0+1; ++iori){		    
			model.SetIRF(/*horizon*/cHorizon, /*dImpulse*/-0.25, /*#MC reps*/cIRFreps, /*Origin t*/iT_irf+iori);
			model.SetPrint(FALSE);
			amIRF = model.GetIRFset(vxirestrange, &vxirange);			// computes IRFs over range of values for xi
			if(iori == 0)
				println("Identified xi without sign restrictions: [", min(vxirange), ", ", max(vxirange), "]");
			decl msignrestri = ((!isdotfeq(amIRF[2][][SRper],0) .&& amIRF[2][][SRper] .> 0) .|| (!isdotfeq(amIRF[0][][SRper],0) .&& amIRF[0][][SRper] .< 0) .|| (!isdotfeq(amIRF[1][][SRper],0) .&& amIRF[1][][SRper] .< 0));
			vxirestrange = deleteifr(vxirange,msignrestri);
			vxirestrange = range(min(vxirestrange),max(vxirestrange),dxi);		// do this to remove duplicates
						
			println(model.ObsYear(iT_irf+t0+iori), "q", model.ObsPeriod(iT_irf+t0+iori),
			": Identified set for lambda with sign restriction: [", min(vxirestrange), ", ", max(vxirestrange),"]");
		}
		vxi = range(min(vxirestrange),max(vxirestrange),iJP ? dxi : dxi/10);			// make the grid finer for US to fill the space
		println("The narrowest set of xi : [", min(vxirestrange), ", ", max(vxirestrange), "]");
		println("\nlapsed time: ", timespan(time),"\n");

		println("Now computing IRFs at all dates over the above range of xi...");

		decl msignrestr = 0, aamIRF = {};
		// first collect all IRFs and sign-restrictions across all dates
		for(decl iori = 0; iori < model.GetSelEnd()-iT_irf-t0+1; ++iori){		    
			model.SetIRF(/*horizon*/cHorizon, /*dImpulse*/-0.25, /*#MC reps*/cIRFreps, /*Origin t*/iT_irf+iori);
			amIRF = model.GetIRFset(vxi, &vxirange);			// computes IRFs over range of values for xi 
 			aamIRF ~= {amIRF};									// store IRFs to avoid repeating later
			// evaluate sign restrictions
			decl msignrestri = ((!isdotfeq(amIRF[2][][SRper],0) .&& amIRF[2][][SRper] .> 0) .|| (!isdotfeq(amIRF[0][][SRper],0) .&& amIRF[0][][SRper] .< 0) .|| (!isdotfeq(amIRF[1][][SRper],0) .&& amIRF[1][][SRper] .< 0));
			msignrestr = msignrestri .|| msignrestr;  			// combine restrictions from previous cases
		}
		// then remove IRFs that fail the joint sign restrictions
		for(decl iori = 0; iori < model.GetSelEnd()-iT_irf-t0+1; ++iori)	    
		    for(decl j = 0; j < 3; ++j){					
		    	aamIRF[iori][j] = deleteifr(aamIRF[iori][j], msignrestr);
				if(j==0 && iJP)
					  aamIRF[iori][j] = 4*aamIRF[iori][j];		// rescale the IRFs for inflation as annualised	for Japanese data
			}

//////// Compute Shadow Rate /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		decl mshadowlimits = <>;		// for Matlab plots
		SetDrawWindow(sprint(iJP ? "JP " : "US ", "shadow rate"));
		SetDraw(SET_LINE, 2, 0, 30);
		decl vy2bar = model.GetY2star(<>);
		decl mY2eff;
		for(decl ix = 0; ix < sizerc(vxirestrange); ++ix){
			decl xi = vxirestrange[ix];
			mY2eff = model.GetShadow(vy2bar, xi);
			DrawTMatrix(0, iJP ? mY2eff' : mY2eff[iT_irf:][]', "", iYear, iquarter, 4, 0, 2);
			mshadowlimits |= iJP ? mY2eff' : mY2eff[iT_irf:][]';
		}
		mshadowlimits = minc(mshadowlimits)'~maxc(mshadowlimits)';
		savemat(sprint("shadowlimits", iJP ? "_JP" : "_US", ".csv"), mshadowlimits);
		ShowDrawWindow();
	   	SaveDrawWindow(sprint("Figures\\", iJP ? sprint("JP_p_",p) : "US","_shadow.eps"));

///////	plot IRFs for particular origins ////////////////////////////////////////////////////////////////////////////////////////////////////////////		
		decl mIRFlimits = <>;		// for Matlab plots
		decl iYear_zlb = iJP ? 2010 : 2009;
		decl iquarter_zlb = 1;
		decl iYear_nonzlb = iJP ? 1990 : 1999;
		decl iquarter_nonzlb = 1;
		decl iT_zlb = (iYear_zlb - iYear)*4+iquarter_zlb-iquarter;
		decl iT_nonzlb = (iYear_nonzlb - iYear)*4+iquarter_nonzlb-iquarter;
		decl sirf_origin_zlb = sprint(model.ObsYear(iT_zlb+iT_irf+t0), "q", model.ObsPeriod(iT_zlb+iT_irf+t0));		// date at which to compute IRF
		decl sirf_origin_nonzlb = sprint(model.ObsYear(iT_nonzlb+iT_irf+t0), "q", model.ObsPeriod(iT_nonzlb+iT_irf+t0));		// date at which to compute IRF

		SetDrawWindow(sprint(iJP ? "JP" : "US", " SetID IRF CKSVAR(",p,")"));
		decl mlimits = new matrix[3][2];
		decl mIRF, mIRFlo, mIRFup, dcoverage;
		dcoverage = 0.67;
		decl mlimits_min = new matrix[3][2];	// record the scale of each origin to adjust in the end
		decl mlimits_max = new matrix[3][2];	// record the scale of each origin to adjust in the end
		for(decl i = 0; i < 2; ++i){
			decl iTreg = i ? iT_zlb :iT_nonzlb;
			decl sorig = i ? sirf_origin_zlb : sirf_origin_nonzlb;
			model.SetIRF(/*horizon*/cHorizon, /*dImpulse*/-0.25, /*#MC reps*/cIRFreps, /*Origin t*/iTreg+iT_irf);
			[mIRFlo, mIRFup] = model.GetIRFbounds(dcoverage, min(vxi), max(vxi));
			mIRFlo[:1][SRper] = mIRFlo[:1][SRper] .> 0 .? mIRFlo[:1][SRper] .: 0;			// impose sign restrictions on C.I. as in Granziera et al step 2 on p.1096
			mIRFup[2][SRper] = mIRFup[2][SRper] .< 0 .? mIRFup[2][SRper] .: 0;				// impose sign restrictions on C.I. as in Granziera et al step 2 on p.1096
			if(iJP){
				mIRFlo[0][] = 4*mIRFlo[0][];
				mIRFup[0][] = 4*mIRFup[0][];
			}
			for(decl j = 0; j < 3; ++j){
				mlimits[j][] = min(mIRFlo[j][])~max(mIRFup[j][]);							// store limits for next plot
				mlimits[j][0] = mlimits[j][0].>0 .? mlimits[j][0]/1.1 .: mlimits[j][0]*1.1;		// set range for plots
				mlimits[j][1] = mlimits[j][1].>0 .? mlimits[j][1]*1.1 .: mlimits[j][1]/1.1;		// set range for plots
				if(j==2)
					mlimits[2][0] = min(-.25, .1);
				if(j==0)
					DrawTitle(2*j+i, sprint("IRF to monetary policy shock in ", sorig));	// column label
				if(i==0)
					DrawText(2*j, asY[j], 0, 0, /*iFontNo*/-1, /*iFontSize*/360, TEXT_YLABEL, 90);		// Y-label
				DrawMatrix(2*j+i, aamIRF[iTreg][j], "", 0, 1, 0, 2*ones(sizeof(aamIRF[iTreg][j]),1));					 // IRFs at all other values of xi
				DrawMatrix(2*j+i, mIRFlo[j][], sprint(100*dcoverage,"\% error bands"), 0, 1, 0, 5);							
				DrawMatrix(2*j+i, mIRFup[j][], "", 0, 1, 0, 5);							 
				DrawAxis(2*j+i, TRUE, .NaN, 0, cHorizon, 2, 2, 0, 0);
				DrawLegend(2*j+i, 400,0, j);
				DrawAdjust(ADJ_AREA_Y, 2*j+i, mlimits[j][0], mlimits[j][1]);
				mIRFlimits ~= minc(aamIRF[iTreg][j])'~maxc(aamIRF[iTreg][j])'~mIRFlo[j][]'~mIRFup[j][]';	// horizon by 24 matrix, after all loops here, each 4-column-unit are upper and lower limits, upper and lower C.I., there are 3 such 4-units for each starting point, 2 starting points
			}
			mlimits_min[][i] = mlimits[][0];
			mlimits_max[][i] = mlimits[][1];			
		}
		savemat(sprint("IRFlimits", iJP ? "_JP" : "_US", ".csv"), range(0,cHorizon)'~mIRFlimits);
		decl mlimits_draw = new matrix[3][2];
		for(decl j = 0; j < 3; ++j){
				mlimits_draw[j][] = min(mlimits_min[j][])~max(mlimits_max[j][]);							// store limits for next plot
				DrawAdjust(ADJ_AREA_Y, 2*j, mlimits_draw[j][0], mlimits_draw[j][1]);
				DrawAdjust(ADJ_AREA_Y, 2*j+1, mlimits_draw[j][0], mlimits_draw[j][1]);
		}		
		ShowDrawWindow();		
		SaveDrawWindow(sprint("Figures\\", iJP ? sprint("JP_p_",p) : "US","_Set_IRF_", 100*dcoverage, ".eps"));
		
/////// Plot IRFs across origins for different horizons ///////////////////////////////////////////////////////////////////////////////////////////////////
		decl mIRFlimits_horz = <>;		// for Matlab plots
		decl horz1 = 4;
		decl horz2 = 8;			// additional horizons to look at

		decl amIRFimp = new array[3];					   // initialise matrices for holding results
		decl amIRFhorz1 = new array[3];
		decl amIRFhorz2 = new array[3];
		for(decl i = 0; i < 3; ++i){
		   amIRFimp[i] = <>;
		   amIRFhorz1[i] = <>;
		   amIRFhorz2[i] = <>;
		}
				
		for(decl iori = 0; iori < model.GetSelEnd()-iT_irf-t0+1; ++iori){	    
			for(decl i = 0; i < 3; ++i){
		    	amIRFimp[i] ~= aamIRF[iori][i][][0];
				amIRFhorz1[i] ~= aamIRF[iori][i][][horz1];
				amIRFhorz2[i] ~= aamIRF[iori][i][][horz2];
			}
		}

		decl cumIRFhorz1 = new array[3];
		decl cumIRFhorz2 = new array[3];
		for(decl i = 0; i < 3; ++i){
		   cumIRFhorz1[i] = new matrix[sizer(amIRFimp[0])][sizec(amIRFimp[0])-horz1-1];
		   cumIRFhorz2[i] = new matrix[sizer(amIRFimp[0])][sizec(amIRFimp[0])-horz2-1];
		}
		for(decl k = 0; k < sizec(amIRFimp[0])-horz1-1; ++k)
			for(decl ii = 0; ii < horz1+1; ++ii)
				for(decl jj =0; jj < 3; ++jj)
					cumIRFhorz1[jj][][k] += aamIRF[k+ii][jj][][horz1-ii];

		for(decl k = 0; k < sizec(amIRFimp[0])-horz2-1; ++k)
			for(decl ii = 0; ii < horz2+1; ++ii)
				for(decl jj =0; jj < 3; ++jj)
					cumIRFhorz2[jj][][k] += aamIRF[k+ii][jj][][horz2-ii];

		for(decl cum = 1; cum < 2; ++cum){
			SetDrawWindow(sprint(iJP ? "JP" : "US", " SetID IRF CKSVAR(",p,") over time"));
			SetDraw(SET_LINE, 2, 0, 30);
			for(decl j = 0; j < 3; ++j){
				for(decl h = 0; h < 3; ++h){
					decl mirf = h == 0 ? amIRFimp[j] : h == 1 ? (cum ? cumIRFhorz1[j] : amIRFhorz1[j]) :
							(cum ? cumIRFhorz2[j] : amIRFhorz2[j]);
					mlimits[j][] = min(mirf)~max(mirf);					// store limits for next plot
					if(j==0)
						DrawTitle(3*j+h, h==0 ? "On impact" : sprint(h == 1 ? horz1 : horz2, " quarters ahead"));
					DrawTMatrix(3*j+h, mirf, "", iYear, iquarter, 4, 0, 2*ones(sizeof(mirf),1));		// IRFs at all other values of lambda
					if(h==0 && isstring(asY[j]))
						DrawText(3*j, asY[j], 0, 0, /*iFontNo*/-1, /*iFontSize*/360, TEXT_YLABEL, 90);		// Y-label
					DrawLegend(3*j+h, 400,0, TRUE);
					DrawAdjust(ADJ_AREA_Y, 3*j+h, mlimits[j][0]>0 ? mlimits[j][0]/1.1 : mlimits[j][0]*1.1,
												mlimits[j][1]<0 ? mlimits[j][1]/1.1 : mlimits[j][1]*1.1);
					mIRFlimits_horz ~= minc(mirf)'~maxc(mirf)';	// horizon by 18 matrix, after all loops here, each 2-column-unit are upper and lower limits, there are 3 such 2-units for each horizon, 3 variables
				}
			}
			ShowDrawWindow();
			SaveDrawWindow(sprint("Figures\\", iJP ? sprint("JP_p_",p) : "US", "_Set_IRF_", cum ? "cum" : "one_off",".eps"));
		}
		savemat(sprint("IRFlimits_horz", iJP ? "_JP" : "_US", ".csv"), mIRFlimits_horz);
		delete model;

	}
}
