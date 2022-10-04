/**to add bootstrap to CKSVAR estimation
**/

#include <oxstd.oxh>
#include <oxfloat.oxh>
#include <oxprob.h>

#import <modelbase>
#import <maxsqp>
#import <solvenle>
#import <maximize>
#import <packages/maxsa/maxsa>

#include "cksvar.oxh"

CKSVAR::GetPackageName()
{
    return "CKSVAR";
}
CKSVAR::CKSVAR()
{
//	ShowBanner(fShowBanner);
	Modelbase();
	Init();
	m_fdiffY2 = FALSE;							// default is that the constrained series is in levels
	m_fLatLags = TRUE;
	m_fUseOLS = TRUE;							// default is to compute OLS estimates of VAR coefficients to use as starting values, unless SetStartPar is used
	m_fCholeski = FALSE;						// default is not to compute Choleski factorization
	m_db = 0.0;
	m_cR = 100;
	m_dkappa = 1;								// default is 1 -- this parameter is not identified seperately from the coefficients on the lags.
	m_dlambda = 0;								// default is 0, unconv policy completely ineffective
	m_dzeta = 1;								// default is 1, no change in reaction function across regimes
	m_mxi = .NaN;								// initialize at NaN. To be set by SetNumMCreps or InitData()
	m_fFAPF = FALSE;
	m_fSimAnn = FALSE;							// Default optimization method is BFGS
	m_fCSVAR = FALSE;							// if TRUE, impose CSVAR restrictions: betatilde = 0 and Clat = (Cobs corresponding to X2)
	m_cOrigin = -1;								// initialize origin of IRF outside estimation sample
	m_cBootRep = 0;								// initialize # of bootstrap replications for inference
	m_iSol = 1;									// number of solutions for betabar as a function of betatilde, xi=lambda*zeta, and Omega
	m_fInitTrue = FALSE;						// if TRUE, use "true" values to initialize bootstrap MLE, otherwise, use OLS
	m_fStoreBoot = FALSE;						// if TRUE, store bootstrap replications
	m_fStoreBootCovar = FALSE;					// if TRUE, store Covariance matrix for each bootstrap replication
	m_sStoreBootFile = "bootstrap.in7";			// file name to store bootstrap repls if desired, can be altered
	m_fBootstrap = FALSE;						// flag, if TRUE, bootstrap has been performed, results stored in m_mBootData;
	m_mBootData = .NaN;							// to hold bootstrap data (first row holds estimates on real data)
	m_dBootSeed = 0;							// holds seed for the bootstrap random numbers (optional), default = 0;
	m_oUnrestricted = .NaN;						// object to hold unrestricted model for LR testing, initialize empty
	m_fTestUnrestr = FALSE;						// flag, if TRUE, compute bootstrap simulations for m_oUnrestricted model under H0 (restricted)
	m_cParallel = 1;							// number of iterations in each parallel for loop
	m_fSavePerSet = FALSE; 						// if TRUE, saves bootstrap replications after the end of each of the parallel runs
	m_dImpulse = 1;								// default value of impulse -- can be changed using SetIRF
	m_fOlsen = FALSE;							// default is to use standard parametrization of likelihood. If true, use Olsen's (1978) param of the Tobit part
	m_fbounded = FALSE;							// flag to indicate whether to do bounded estimation
	m_fUseBetabar = FALSE;						// if TRUE, use the parameterization betabar instead of betatilde in GetModelPar and estimation -- can be used to impose sign restrictions for given value of xi = lambda*zetA
	m_mrannIRForigin = m_mSigns = <>;
	m_vBreakDates = <>;							// vector of dates when parameters are allowed to change -- default is empty
	m_aBreakParIdx = {};						// array of vectors of indices of parameters are allowed to change at m_vBreakDates 
	m_cBreakPar = 0;							// # of parameters that are allowed to change
	m_cBreaks = 0;								// # of breaks
}
CKSVAR::SetOlsen(const fOlsen)
/**	@param fOlsen, flag
**/
{
	m_fOlsen = fOlsen;							// default is to use standard parametrization of likelihood. If true, use Olsen's (1978) param of the Tobit part
}
CKSVAR::GetOlsen()
{
	return m_fOlsen;							
}
CKSVAR::SetUseBetabar(const fUseBetabar) 		// sets flag m_fUseBetabar
{
	m_fUseBetabar = fUseBetabar;
}
CKSVAR::GetUseBetabar() 						// returns flag m_fUseBetabar
{
	return m_fUseBetabar;
}
CKSVAR::SetParBounds(const viPar, const vLo, const vHi)
/**	@param viPar, index of vector of indices for bounded parameters
	@param vLo,vHi vectors of lower/upper bounds, vLo < vHi	**/
{
	m_fbounded = TRUE;
	m_vLo = constant(-.Inf,m_cPar,1);
	m_vHi = -m_vLo;
	m_vLo[viPar] = vLo;
	m_vHi[viPar] = vHi;
}
CKSVAR::GetParBounds()
{
	if(m_fbounded)
		return m_vLo ~ m_vHi;
	return ones(m_cPar,1)*<-.Inf,.Inf>;
}
CKSVAR::ClearParBounds()
{
	m_fbounded = FALSE;
}
CKSVAR::SetBound(const b)
/**	(Re)sets the bound					  
	@param b, double 
**/
{
	m_db = b;
}
CKSVAR::SetFAPF(const fFAPF)
/**	(Re)sets whether to use Fully adapted particle filter or not
	@param fFAPF, flag
**/
{
	m_fFAPF = fFAPF;
}
CKSVAR::SetBreaks(const vBreakDates, const aBreakParIdx)
/**	sets the dates when breaks are allowed and the coefficients that change at those dates
	@param vBreakDates, vector [cb] of indices of break dates, must be between 1 and m_cT-2
	@param aBreakParIdx, array [cb] of indices of parameters that change at those dates -- if -1, all change
**/
{
	if(sizeof(vec(vBreakDates)) != sizeof(aBreakParIdx)){
		println("Break dates set incorrectly: arguments must have the same size");
		m_vBreakDates = <>; m_aBreakParIdx = {};
		return FALSE;
	}
	decl vsort=sortcindex(vec(vBreakDates));
	m_vBreakDates = matrix(vec(vBreakDates)[vsort]);
	m_aBreakParIdx = array(aBreakParIdx[vsort]);
	m_cBreaks = sizeof(m_vBreakDates);
}
CKSVAR::SetStoreBoot(const fStoreBoot, ...)
/**	sets flag to determine whether to store boostrap estimates
	@param fStoreBoot, flag, default = FALSE
	(optional) flag, to store covariance matrix for each bootstrap (default is false as it's expensive)
**/
{
	m_fStoreBoot = fStoreBoot;
	decl arg = va_arglist();
	if(sizeof(arg)>0)							// if nonzero, set flag to compute covariance matrix for each bootstrap
		m_fStoreBootCovar = m_fStoreBoot == FALSE ? FALSE : arg[0];		// never store covariances if you haven't stored parameter estimates!
		
}
CKSVAR::SetStoreBootFile(const sStoreBootFileName)
/**	sets flag to determine whether to store boostrap estimates
	@param sStoreBootFileName, string, name of file to store bootstrap reps, must have valid extension, .in7, csv, xls, dat
**/
{
	m_sStoreBootFile = sStoreBootFileName;
	m_fStoreBoot = TRUE;				// assume that if you give a file name to store, you'll want to store
}
CKSVAR::SetBootstrap(const cBootRep, const fInitTrue, const cParallel, const fSavePerSet)
/**	Sets # of bootstrap replications, initialization of estimation per boostrap, and parallelization settings
	@param cBootRep, positive integer
	@param fInitTrue, flag, if TRUE, use "true" values to initialize bootstrap MLE, otherwise, use OLS
	@param cParallel, positive integer, number of iterations in each run of parallel for loop
	@param fSavePerSet, flag, if TRUE, save file after each set is completed
**/
{
	m_cBootRep = max(0,int(cBootRep));			   // make sure it's a non-negative integer 
	m_fInitTrue = fInitTrue;
	m_cParallel = int(max(1,min(cParallel,m_cBootRep)));		// make sure it's an integer between 1 and m_cBootRep
	m_fSavePerSet = fSavePerSet;
}
CKSVAR::GetBootstrapData(const fSave)
/**	returns the bootstrap data
	@param fSave, flag, if TRUE, save to m_sBootFileName
	returns array [2] {m_mBootData, m_asBootNames}
**/
{
	if(fSave)
		savemat(m_sStoreBootFile, m_mBootData, m_asBootNames);

	return {m_mBootData, m_asBootNames};
}
CKSVAR::SetBootstrapData(const arg)
/**	sets bootstrap data computed in a previous run
	@param arg, either data file name or array {mBootData, asBootNames};
**/
{
	if(isfile(arg))
		m_mBootData = loadmat(arg, &m_asBootNames);
	else{
		m_mBootData = arg[0]; m_asBootNames = arg[1];
	}
}
CKSVAR::SetBootstrapSeed(const dBootSeed)
/**	Sets seed for bootstrap random numbers
	@param dBootSeed, double
**/
{
	m_dBootSeed = dBootSeed;
}
CKSVAR::ClearBootstrap()
/**	Sets m_fBootstrap to FALSE and clears bootstrap data
**/
{
	m_fBootstrap = FALSE;
	m_mBootData = .NaN;
}
CKSVAR::GetZLBcount()
/**	return # of observations at the ZLB
**/
{
	if(!InitData())				// make sure data has been initialized, so that m_cZLB is set
		return .NaN;
	return double(m_cZLB);
}
CKSVAR::SetUnrestrictedModel(const oUnrestricted)
/**	Stores unrestricted model object for LR testing in memory
	@param oUnrestricted, object, typically cloned version of the current object with some restrictions relaxed
**/
{
	m_oUnrestricted = oUnrestricted;
}
CKSVAR::DeleteUnrestrictedModel()
/**	Deletes unrestricted model object from memory
**/
{
	delete m_oUnrestricted;
}
CKSVAR::SimData(const mrann)
/**	Simulates data from the model at current parameter values
	Up to this point (1/2/19), I have not checked whether it works when m_fdiffY2=TRUE
	@param mrann [T][k] matrix of standard normals (to be used in case you don't want to redraw) or <> (draw a new set of normals)
	@return if mrann = <>, [m_cT][k] matrix of data on observables, Y, if mrann == matrix[T][k], then {mX, mXst} of [T][k*p], [T][k] matrices
**/
{
	decl tau, vC2obs, vC2lat, vbetatil, mC1obs, mC1lat, vdelta, mOm1g2Ch;			// obtain the parameters first
	[tau, vC2obs, vC2lat, vbetatil, mC1obs, mC1lat, vdelta, mOm1g2Ch] = this.GetModelPar(this.GetPar());
	decl T, mrannum;
	if(mrann == <>){
		T = m_cT;
		mrannum = rann(T, m_cY1+1);
	}
	else if(!ismatrix(mrann) || columns(mrann) != m_cY1+1){
		eprint("CKSVAR::SimData expects a matrix with ". m_cY1+1, " columns.");
		return FALSE;
	}
	else{
		T = rows(mrann);
		mrannum = mrann;
	}
	if(rows(mrann) > m_cT && !isdotfeq(this.GetGroup(X_VAR),1))
		eprint("WARNING: CKSVAR::SimData keeps exogenous variables fixed at last observed value.");

	decl mu = new matrix[T][m_cY1+1];
	mu[][m_cY1] = mrannum[][m_cY1]*tau;		// u2
	if(m_cY1)
		mu[][:m_cY1-1] = mu[][m_cY1]*vec(vdelta)'+mrannum[][:m_cY1-1]*mOm1g2Ch';

	decl k = sizec(mu);
	decl p = max(this.GetMaxGroupLag(Y1_VAR),
			this.GetMaxGroupLag(Y2_VAR)+m_fdiffY2);		// VAR order, assumes all Y variables have same number of lags
	decl mY = new matrix[T+p][k];						// simulated sample must include p initial values to maintain size T. Simulations conditional on observed initial values.
	decl vy2star = mY[][k-1];
	decl Y0 = (this.GetGroupLag(Y1_VAR, 1, this.GetMaxGroupLag(Y1_VAR))
			~ this.GetGroupLag(Y2_VAR, 1, this.GetMaxGroupLag(Y2_VAR)))[0][];
	decl mX0 = this.GetGroup(X_VAR)[:min(T,m_cT)-1][];					// exogenous (typically deterministic) regressors, to be treated as fixed in the simulation, not valid unless strictly exogenous
	if(T>m_cT)
		mX0 |= ones(T-m_cT,1)*mX0[m_cT-1][];							// if you need bigger sample than observed, fix X0 to its last value -- not good for trends/seasonals
	decl cX0 = columns(this.GetGroup(X_VAR));							// # of exogenous and deterministic regressors
	decl mX = new matrix[mrann == <> ? 1 : T][columns(Y0)];				// initialize observed lags
	decl mXst = m_fLatLags ? new matrix[mrann == <> ? 1 : T][p] : 0;	// initialize latent lags (lags of Y2*) Xt*, only necessary if m_cLatLags == TRUE
	decl vb = constant(m_vb[0], p, 1) | m_vb;			// extrapolate backwards using the same bound as the first observation i the sample
	for(decl t = -p; t < T; ++t){
		if(t<0){									// first, store initial values in first p rows of mY
			mY[t+p][] = Y0[][range(0,k-1)*p-t-1];	// recall that Y0 has Y1_-1,...,Y1_-p, Y2_-1,...,Y2_-p, ..., Yk_-1,...,Yk_-p
			vy2star[t+p] = m_vY2stinit[][-t-1]; 	// and similarly for Y2star
			continue;
		}
		decl it = mrann == <> ? 0 : t;
		for(decl ik = 0; ik < k; ++ik){
			for(decl j = 0; j < p; ++j)
				mX[it][ik*p+j] = mY[t+p-j-1][ik];
		}
		if(m_fLatLags){
			for(decl j = 0; j < p; ++j)					// not executed if m_cY2lags = 0
				mXst[it][j] = vy2star[t+p-j-1];		
		}
		vy2star[t+p] = (cX0 ? mX0[t][]*vC2obs[:cX0-1] : 0) + (sizeof(vC2obs)>cX0 ? mX[it][]*vC2obs[cX0:] : 0) + ( (vC2lat != <>) ? mXst[it][]*vC2lat : 0) + mu[t][k-1];

		mY[t+p][k-1] = max(vy2star[t+p],vb[t+p]);
		if(k>1)
			mY[t+p][:k-2] = (cX0 ? mX0[t][]*mC1obs[:cX0-1][] : 0) + (sizeof(mC1obs)>cX0 ?  mX[it][]*mC1obs[cX0:][] : 0) + ( (mC1lat != <>) ? mXst[it][]*mC1lat : 0)
					- min(vy2star[t+p] - vb[t+p],0) * vbetatil' +mu[t][:k-2];
	}
	return mrann == <> ? mY : {mX, mXst};
}				
CKSVAR::InsideLoop(const ib)
/** gives the contents of for loop in ParBoot
**/
{
	decl oboot = clone(this);			// clone the current object, to copy all model settings (NOTE: this will also keep the uniform random draws the same, m_mxi)
//	if(m_dBootSeed)
//		ranseed(m_dBootSeed+ib*3+iset*);		// choose an arbitrary seed sequence (cannot guarantee it won't interfere with MC replications)
	decl mY = this.SimData(<>);			// simulate new data at estimated parameter values
	oboot.Renew(mY, m_asY1~m_asY2, this.GetSelStart()-m_cY2lags + m_fdiffY2);
	oboot.ClearModel();					// you need this in order to renew the data in the estimation (otherwise it uses the original data!)
	oboot.SetPrint(FALSE);				// no printing of bootstrap estimations!
	if(m_fInitTrue)	oboot.SetStartPar(this.GetPar()');			// initialize at true value of bootstrap DGP = MLE
	if(oboot.Estimate()){				// if estimation is successful, store data, otherwise skip
		decl vidx = m_fTestUnrestr+range(1,this.GetFreeParCount());
		m_mIntermBoot[ib][vidx] = oboot.GetFreePar()';
		m_mIntermBoot[ib][0] = oboot.GetLogLik();
		if(m_fStoreBootCovar)
			m_mIntermBoot[ib][max(vidx)+1:] = vech(oboot.GetCovar())';
		delete oboot;					// clear memory
		if(m_fTestUnrestr){				// if Testing against unrestricted model is required...
			decl obootUnr = clone(m_oUnrestricted);
			obootUnr.Renew(mY, m_asY1~m_asY2, this.GetSelStart()-m_cY2lags + m_fdiffY2);
			obootUnr.ClearModel();					// you need this in order to renew the data in the estimation (otherwise it uses the original data!)
			obootUnr.SetPrint(FALSE);				// no printing of bootstrap estimations!
			if(obootUnr.Estimate())
				m_mIntermBoot[ib][1] = 2*(obootUnr.GetLogLik() - m_mIntermBoot[ib][0]);		// compute LR statistic
			delete obootUnr;
		}
	}
}
CKSVAR::ParBoot(const csize, const mEstimData)
/**	performs parallel for loop for bootstraps -- auxiliary function used twice in Bootstrap()
**/
{
	decl cFreePar = this.GetFreeParCount();
	m_mIntermBoot = constant(.NaN, csize,1+m_fTestUnrestr+cFreePar+m_fStoreBootCovar*cFreePar*(cFreePar+1)/2);		// store bootstrap replications of loglik, LR (if needed), free parameters, and vech(Covar) (if needed)
	if(csize == 1)
		InsideLoop(0);							// no parallelization
	else
		parallel for(decl ib = 0; ib < csize; ++ib)
			InsideLoop(ib);

	decl arg = deleter(m_mIntermBoot);
	m_mBootData |= arg;
	if(m_fSavePerSet && m_fStoreBoot)			 								// store every csize bootstraps
		savemat(m_sStoreBootFile, mEstimData|m_mBootData, m_asBootNames);
	
	return sizeof(arg);						// return number of successful bootstrap runs
}
CKSVAR::Bootstrap()
/**	performs parametric bootstrap, conditioning on initial values of Y_VAR and all values of X_VAR
	NOTE: valid only if model exogenous regressors (X_VAR) are strictly exogenous or determinstic
	returns 1 if successful
**/
{
	if(m_iModelStatus < MS_ESTIMATED)
		return;										// Only to be done after model has been estimated
	println("\nBootstrapping...\n");
	m_mBootData = <>;	
	m_asBootNames = {"LogLik"} ~ ( m_fTestUnrestr ? {"LR"} : {}) ~ this.GetFreeParNames();
	decl mEstimData = this.GetLogLik() ~
				( m_fTestUnrestr ? (m_oUnrestricted && m_oUnrestricted.GetModelStatus()==MS_ESTIMATED ? 2*(m_oUnrestricted.GetLogLik()-this.GetLogLik()) : .NaN) : <>)		
				~ this.GetFreePar()' ~ (m_fStoreBootCovar ? vech(this.GetCovar())' : <>);	// to store values from actual data estimation together with bootstrap samples below
	decl cRemainingBootRep = m_cBootRep-sizeof(m_mBootData);		// number of bootstraps remaining
	if(m_dBootSeed)
		ranseed(m_dBootSeed);							   // set seed for bootstrap ran num generator
	decl cBlocks = ceil(m_cBootRep/m_cParallel);			// # of blocks to force parallel computation over
	decl cextra = 10, iter=0			;				// cap the number of attempts to reach desired number
	// Keep working until you hit the required number (allow some extra iterations in case some bootstraps fail)
	while(cRemainingBootRep>0 && iter<cBlocks+cextra){
		decl time = timer();
		oxprintlevel(-1);			// switch printing off to avoid printing error messages 
		decl csize = min(m_cParallel,cRemainingBootRep);
		decl csuccess = ParBoot(csize, mEstimData);			
		oxprintlevel(1);			// switch printing back on
		println("Performed ", csize, " bootstraps, ", csuccess, " successful, ", csize-csuccess, " failed. Lapsed time: ", timespan(time));
		cRemainingBootRep -= csuccess;				// reduce the remaining bootstrap replications required
		++iter;
	}
	// append the actual estimates into the bootstrap sample (at the top)
	m_fBootstrap = TRUE; 							// set the flag to avoid repeating later, if unnecessary
	m_mBootData = mEstimData|m_mBootData;	 		// add estimates at the top of the BootData matrix
	// finally store everything one last time
	if(m_fStoreBoot)
		savemat(m_sStoreBootFile, m_mBootData, m_asBootNames);
	return m_fBootstrap;
}
CKSVAR::SetCSVAR(const fCSVAR)
/**	Sets whether to impose CSVAR restrictions or not
	@param fCSVAR, flag if TRUE impose restrictions
**/
{
	if(m_fCSVAR != fCSVAR && m_iModelStatus >= MS_DATA){
		this.ClearModel();				// if you're switching from a previously unrestricted model, you need to re-initialize the data
		m_fUseOLS = TRUE;				// reset OLS as default initial values (otherwise code may fail)
	}
	
	m_fCSVAR = fCSVAR;
	if(fCSVAR && !m_fLatLags)			// if Latent lags have not been selected, add them
		m_fLatLags = TRUE;				// this will not re-initialize the parameters, you'll need to do that manually
	if(m_iModelStatus >= MS_DATA)
		this.InitPar();
	
}
CKSVAR::SetNumMCreps(const cR)
/** Sets number of MC reps in Importance sampling to integrate out latent lags
	If not called, the default is 100
**/
{
	if(m_iModelStatus < MS_DATA)		// if function is called before InitData(), just set m_cR (you don't know m_cT yet to do the rest)
		m_cR = cR;
		
	else{								// draw random numbers for Y2* and keep them fixed (for computation of likelihood and moments of Y2*, IRFs)
		m_mxi = (cR == m_cR || cR < 0) ? m_mxi : cR < m_cR ? m_mxi[][:cR-1] : (m_mxi ~ ranu(m_cT, cR-m_cR));	// select subset or append
		m_cR = cR <= 0 ? m_cR : cR;		// number of MC reps in Importance sampling to integrate out latent lags, if <=0, use default value
	}
}
CKSVAR::SetLatentLags(const flatlags)
/** Sets whether the model contains latent lags or not
	If not, we avoid importance sampling, so code is much faster
**/
{
	if(m_fLatLags != flatlags && m_iModelStatus >= MS_DATA){
		m_fLatLags = flatlags;				// flag: if TRUE (default), the model contains latent lags
		this.InitData();					// reinitialize data to drop coefficients on latent lags from parameter count (sets m_cPar)
		if(m_fPrint)
			println("CKSVAR::SetLatentLags: Parameters have been set to initial values.");
		this.InitPar();						// reinitialize parameter to get the dimension of m_vPar right
	}
	m_fLatLags = flatlags;					// flag: if TRUE (default), the model contains latent lags
}
CKSVAR::SetDiffY2(const fdiffY2)
/** Sets flag m_fdiffY2, to indicate whether the constrained variable appears in first differences or not
**/
{
	m_fdiffY2 = fdiffY2;
}
CKSVAR::Setkappa(const kappa)
/** Sets kink in Y2* eqn -- not identified seperately from the coefficients on the lags.
**/
{
	m_dkappa = kappa;
}
CKSVAR::Setlambda(const dlambda)
/** Sets coef on shadow rate, efficacy of unconventional policy
**/
{
	m_dlambda = dlambda;
}
CKSVAR::Setzeta(const dzeta)
/** Sets ratio of vol of eps_2 in ZLB vs nonZLB regim
**/
{
	m_dzeta = dzeta;
}
CKSVAR::InitData()
/**	Modelbase function that extracts the data from the underlying database and if successful sets model status to MS_DATA
**/
{
    ClearModel();
    if (m_iT1sel < 0)
        SetSample(-1, 1, -1, 1);

    m_iT1est = m_iT1sel;  m_iT2est = m_iT2sel;

	m_cT = m_iT2est - m_iT1est + 1;

    m_mY1 = GetGroupLag(Y1_VAR, 0, 0);			// Y1 (unconstrained) variables
    m_mY2 = GetGroupLag(Y2_VAR, 0, 0);			// Y2 (constrained) variable, must be exactly 1
    m_mY = m_mY1~m_mY2; 
    decl k = columns(m_mY);
	m_cY1 = columns(m_mY1);
	if (k == 0)
    {   println("No Y variable selected");
        return FALSE;
    }
	if (columns(m_mY2) != 1)
    {   println("Only 1 constrained variable allowed.");
        return FALSE;
    }
	if (m_fdiffY2 == TRUE)
	{	m_vY2laglev = GetGroupLag(Y2_LAGLEV,1,1);		// load the level of Y2 lagged once for use in latent regression
		if (columns(m_vY2laglev) != 1){
			println("Lagged level of Y2 is needed when Y2 is in first difference.");
			return FALSE;
		}
	}
	decl mY1lags = GetGroupLag(Y1_VAR, 1, this.GetMaxGroupLag(Y1_VAR));	  	// lags of Y1 on RHS
	decl mY2lags = GetGroupLag(Y2_VAR, 1, this.GetMaxGroupLag(Y2_VAR));		// lags of Y2 on RHS
	m_cY2lags = columns(mY2lags);											// # of lags (they don't have to be the same as Y1 lags)
	m_mX = GetGroup(X_VAR)~mY1lags~mY2lags;									// all exogenous regressors
    m_cX = columns(m_mX);													// # of observed exogenous regressors

	if(any(GetGroup(Y2_BOUNDS)))
		m_vb = m_db+GetGroup(Y2_BOUNDS);
	else //if(isdouble(m_db) || (ismatrix(m_db) && sizerc(m_db)==1))
		m_vb = constant(m_db, m_cT, 1);
	m_vb = m_fdiffY2 ? m_vb - m_vY2laglev : m_vb;
	decl vZLB = !(m_mY2 .> m_vb);						// indicator: 1 if at the ZLB, 0 otherwise
	m_cZLB = int(sumc(vZLB));

	decl cY2lat = this.GetMaxGroupLag(Y2_VAR) + m_fdiffY2;					// # of latent lags in unrestricted specification. Equal to lags of Y2 if the latter is in levels, plus 1 if it is in first differences
	if (m_cT-m_cZLB < m_cX+cY2lat+k){										// assuming #Y1lags = #Y2lags, you need at least as many obs above b as coefficients (m_cX+cY2lat) + k for the variance
    	if(m_fPrint)
			println("Only ", m_cT-m_cZLB, " unconstrained observations. Not enough for estimation, you need ", m_cX+cY2lat+k, ".");
        return FALSE;
    }//else

	if(m_cZLB == 0){									// of there are no observations at the ZLB, perform OLS
		if(m_fPrint)
			println("No observations at the ZLB.");
		if(m_fCSVAR)
			if(m_fPrint)
				println("CSVAR reduced-form estimates are equal to OLS.");
		else{
			if(m_fPrint)
				println("Model is under-identified.");
			return FALSE;
		}
	}//else
	if (!m_fCSVAR && m_cZLB < 1+m_fLatLags*cY2lat){			// For KSVAR, you need one observation for betatilde, and for CKSVAR you also need one for each latent lag
    	if(m_fPrint)
			println("Too few observations at the ZLB. Model is underidentified.");
		return FALSE;
    }

	// Set initial values of Y2*
	if(!m_fdiffY2 && m_cY2lags)
		m_vY2stinit = mY2lags[0][];										// set equal to the obseved Y2lags
	else if(m_fdiffY2){
		if(cY2lat == 1)
			m_vY2stinit = m_vY2laglev[0];									// if VAR order p is 1, then returns Y2*_0
		else
			m_vY2stinit = mY2lags[0][]~m_vY2laglev[0][];					// returns DY2*_0 ~ ... DY2*_p-2 and Y2*_0
	}

	m_cPar = m_fCSVAR ?	m_cX * k + k*(k+1)/2								// total # of parameters in CSVAR, same as in standard VAR
			: (m_cX + m_fLatLags*m_cY2lags)*m_cY1 + (m_cX + m_fLatLags*cY2lat)+m_cY1+k*(k+1)/2;			// total # of parameters in CKSVAR: (coefs on obs&lat regr in Y1 eqs) + (coefs on obs&lat regr in Y2 eq) + betatildes + variance
	m_cBreakPar = 0;
	m_aRegimeDates = new array[m_cBreaks+1];
	decl idx = 0;
	for(decl i = 0; i < m_cBreaks; ++i){
		if((m_vBreakDates[i] < 1) || (m_vBreakDates[i] > m_cT-2) || (i>0 && (m_vBreakDates[i]-m_vBreakDates[i-1]<1))){
			println("Break dates must not be too close to each other or endpoints of sample. Reset breaks to none.");
			m_vBreakDates = <>; m_aBreakParIdx = {}; m_cBreaks = 0;				// return all break settings to default
			break;
		}
		if(any(m_aBreakParIdx[i].<-1) || any(m_aBreakParIdx[i].==<>)){
			println("Break parameter indices set incorrectly. Reset breaks to none.");
			m_vBreakDates = <>; m_aBreakParIdx = {}; m_cBreaks = 0;
			break;
		}
		m_aBreakParIdx[i] = (m_aBreakParIdx[i] == -1) ? range(0,m_cPar-1)' : vec(m_aBreakParIdx[i]);	// -1 means all parameters allowed to change, force column vector of indices
		if(any(m_aBreakParIdx[i] .> m_cPar-1)){
			println("Break date ", m_vBreakDates[i], " parameter indices set outside parameter range. Reset breaks to none.");
			m_vBreakDates = <>; m_aBreakParIdx = {}; m_cBreaks = 0;
			break;
		}
		m_cBreakPar += sizerc(m_aBreakParIdx[i]);					// increase the number of parameters by breaks (-1 corresponds to all)
		m_aRegimeDates[i] = range(idx,m_vBreakDates[i]-1)';			// date range for ith regime
		idx = m_vBreakDates[i];
	}
	m_aRegimeDates[m_cBreaks] = range(idx, m_cT-1);					// final regime (only regime if no breaks);
	m_cPar += m_cBreakPar;
	m_mxi = sizeof(m_mxi) != m_cT ? ranu(m_cT, m_cR) : m_mxi;		// draw uniform random numbers to integrate out Y2* from the likelihood using importance sampling (or do nothing if already drawn earlier)
	// get names of regressors
	GetGroupLagNames(Y1_VAR, 0, 0, &m_asY1);			// get Y1 variable names
	GetGroupLagNames(Y2_VAR, 0, 0, &m_asY2);			// get Y2 variable name

	GetGroupNames(X_VAR, &m_asX);
	decl aux;
	GetGroupLagNames(Y1_VAR, 1, 1000, &aux);
	m_asX ~= aux;											
	GetGroupLagNames(Y2_VAR, 1, 1000, &aux);
	m_asX ~= aux;											// names of all observed exogenous regressors
	
	m_iModelStatus = MS_DATA;

	return TRUE;
}
CKSVAR::InitPar()
/**	Modification of the Modelbase::InitPar() to allow for OLS estimates to be used as initial parameters for estimation.
**/
{
    if (m_iModelStatus < MS_DATA)
        InitData();
    if (m_iModelStatus < MS_DATA)                              /* need data */
        return FALSE;

	this.SetParCount(m_cPar);
	if(m_fUseOLS){
		// Start with simple OLS estimation ignoring the ZLB to get starting values
		decl mC, Omega;							// OLS estimates of RF VAR coefficients and error variance
		[mC, Omega] = this.GetOLS();

		// organise the parameters as expected, eq Y2 first, then Y1
		decl tau = sqrt(Omega[m_cY1][m_cY1]);						// sd of error in Y2, u2
		decl vdelta = m_cY1 ? Omega[:m_cY1-1][m_cY1]/tau^2 : <>;	// regression coefficient of u1 on u2
		decl mOm1g2 = m_cY1 ? Omega[:m_cY1-1][:m_cY1-1]-tau^2*vdelta*vdelta' : <>;	// conditional variance of u1|u2
	
		// collect all parameters
		decl vC2 = mC[][m_cY1];										// coeffs on observed regressors in Y2 eq
		decl mC1 = m_cY1 ? mC[][:m_cY1-1]' : <>;					// coefficients on observed regressors in Y1 eqs
		decl vbetatilde = !m_fCSVAR ? zeros(m_cY1,1) : <>;			// contemp impact coefficients, initialize at zero
		decl vC2lat = !m_fCSVAR*m_fLatLags ? zeros(1,m_cY2lags+m_fdiffY2) : <>;			// coefficients on latent lags in Y2 eq, initialize at zero
		decl mC1lat = !m_fCSVAR*m_fLatLags ? zeros(m_cY1,m_cY2lags) : <>;				// coefficients on latent lags in Y1 eqs, initialize at zero
		decl vpsi2 = tau|vecr(vC2)|vecr(vC2lat);										// gather all parameters in Y2 eqn
		decl vpsi1 = vbetatilde|vecr(mC1)|vecr(mC1lat)|vdelta|vech(choleski(mOm1g2));	// and Y1 eqns
		m_vPar = m_fOlsen ? (1|vpsi2[1:])/tau|vpsi1 : vpsi2|vpsi1;
		if(m_cBreakPar)
			m_vPar |= zeros(m_cBreakPar,1);			// initialize changes in coefficients that are allowed to vary at zero
		
	}

	ClearEstimation();
    ResetFixedPar();
    m_iModelStatus = MS_PARAMS;

    return TRUE;
}
CKSVAR::SetStartPar(const ParFree)
/**	Sets starting values for the parameters
	@param: ParFree: either vector of free parameters or string "OLS", if neither, return FALSE
	@return TRUE if successfully set starting values
**/
{
	if(isstring(ParFree) && ParFree == "OLS"){
		m_fUseOLS = TRUE;
		return this.InitPar();
	}
	//else
//	m_fUseOLS = FALSE;							// change the default value so as not to compute OLS estimates
	// the rest is as in the original Modelbase::SetStartPar
	if (!this.InitPar())
        return FALSE;
    SetFreePar(ParFree);
    
    return TRUE;
}
CKSVAR::GetAuxPars(const mOm1g2Ch, const vbetatil, const vdelta, const tau)
/**	computes parameters Omega_{1|2}^{-1} (mOm1g2inv), tau2, Xi1inv, arg2, to be used in SimLik
	@param: mOm1g2Ch [m_cY1][m_cY1] Choleski factor of red form conditional variance matrix Omega_{1|2}
	@param: vbetatil [m_cY1][1]	betatilde vector in CKSVAR
	@param: vdelta [m_cY1][1] regression coef of u_1 on u_2 (controlling for exogenous regressors)
	returns array [4] {mOm1g2inv, Xi1inv, tau2, arg2}
**/
{
	decl mOm1g2, mOm1g2inv = <>, vdeltil, mXi1inv = <>, tau2, arg2 = <>;
	if(m_cY1){
		decl flold = oxwarning(0); 			// disable warnings
		mOm1g2 = mOm1g2Ch*mOm1g2Ch';
		mOm1g2inv = invertsym(mOm1g2);
	
		vdeltil = vdelta - vbetatil;						// delta tilde in the notes
		mXi1inv = invertsym(mOm1g2+vdeltil*vdeltil'*tau^2);
		tau2 = tau*sqrt(1-tau^2*vdeltil'mXi1inv*vdeltil);
		arg2 = tau^2*vdeltil'mXi1inv;
	    oxwarning(flold);        			// reset the previous warning flags.	
	}
	else
		tau2 = tau;
	return {mOm1g2inv, mXi1inv, tau2, arg2};
}	
/** Prepare and estimate the model -- exactly as in Modelbase with eprint replaced by print (so it can be disabled during bootstrap)
@returns TRUE if successful
@comments Calls `Modelbase::InitData` and `Modelbase::InitPar` if not yet called.
Then uses `Modelbase::InitData` to get the free parameters, calls `Modelbase::DoEstimation`
`Modelbase::SetFreePar` and `Modelbase::Output` (unless switched off).
*/
CKSVAR::Estimate()
{
    decl  vpstart, vpfree, estout;

    ClearEstimation();

    if (m_iModelStatus < MS_DATA && !InitData())
    {   print("Failed to get data.\n");
        return FALSE;
    }
    if (m_iModelStatus < MS_PARAMS && !InitPar())
    {   print("Failed to get starting values.\n");
        return FALSE;
    }
    m_iModelStatus = MS_EST_FAILED;
    vpstart = GetFreePar();         // map parameters to format used in iterative procedure
    estout = DoEstimation(vpstart); // do the estimation
    vpfree = isarray(estout) ? estout[0] : estout;

    SetFreePar(vpfree);             // map estimated parameters to normal format

    if (m_iResult >= MAX_CONV && m_iResult < MAX_MAXIT)
        m_iModelStatus = MS_ESTIMATED;
    else
        m_iModelStatus = MS_EST_FAILED;

    if (m_fPrint)
    {   Output();
        if (isarray(estout))
            OutputMax(estout[1], m_iResult, vpstart, estout[2]);
    }
    return m_iModelStatus == MS_ESTIMATED;
}
CKSVAR::LogLik(const vp, const adFunc, const avScore, const amHessian)
/**	log-likelihood function of the model
	vp: unrestricted parameters
**/
{
	decl vpall, vS, mY2star, mw;
	if(!Reparameterize(&vpall, vp))
		return 0;
	[vS, mY2star, mw] = this.SimLik(vpall, /*fSimulate*/m_fLatLags && m_cZLB);		// if m_cZLB==0, you just need to compute Normal likelihood without any simulation
	adFunc[0] = double(sumc(log(vS)));
	if(m_fSimAnn && ismissing(adFunc[0]))
		adFunc[0] = -1e10;							  			// modification for MaxSA because it dislikes missing values
	return !isnan(adFunc[0]);								
}
CKSVAR::GetOmega(const tau, const vdelta, const mOm1g2Ch)
/** Returns reduced form variance matrix
	@param tau: double, st dev of u2
	@param vdelta: vector[m_cY1][1]: coef of regression of u1 on u2
	@param mOm1g2Ch: choleski factor of var(u1|u2)
	@returns matrix [m_cY1+1][m_cY1+1]
**/
{
	decl mOm1g2 = mOm1g2Ch*mOm1g2Ch';					// variance of u_1|u_2
	decl mOm11 = mOm1g2 + tau^2*vdelta*vdelta';	// unconditional variance matrix for Y1
	decl mOmega = mOm11 ~ vdelta*tau^2
			| vdelta'tau^2 ~ tau^2;			// variance matrix of reduced form errors for all Y
	return mOmega;
}
CKSVAR::GetBreakPar(const vpall)
/** Returns set of full parameter vectors one for each regime
	vpall: full (unrestricted) parameter vector
	returns array[m_cBreaks] of vectors of parameters
**/
{
	decl apars = new array[m_cBreaks];
	decl vpbreakall = vpall[m_cPar-m_cBreakPar:];			// extract the break coefficients only
	decl vpfull = vpall[:m_cPar-m_cBreakPar-1];							// initialize coeffcients for all parameters @ default regime
	decl idx0 = 0;
	for(decl i = 0; i < m_cBreaks; ++i){
		decl vPregimei = vpfull;
		decl idx1 = idx0+sizeof(m_aBreakParIdx[i]);
		vPregimei[m_aBreakParIdx[i]] += vpbreakall[idx0:idx1-1];
		idx0 = idx1;

		decl tau, vC2obs, vC2lat, vbetatil, mC1obs, mC1lat, vdelta, mOm1g2Ch;			// obtain the parameters first
		[tau, vC2obs, vC2lat, vbetatil, mC1obs, mC1lat, vdelta, mOm1g2Ch] = this.GetModelPar(vPregimei);
		decl tau2, arg2, mOm1g2inv, mXi1inv;
		[mOm1g2inv, mXi1inv, tau2, arg2] = this.GetAuxPars(mOm1g2Ch, vbetatil, vdelta, tau);
		apars[i] = {tau, vC2obs, vC2lat, vbetatil, mC1obs, mC1lat, vdelta, mOm1g2Ch, mOm1g2inv, mXi1inv, tau2, arg2}; 
	}	
	return apars;
}
CKSVAR::GetModelPar(const vpall)
/** Computes VAR coefficients and variance parameters for the likelihood
	vpall: full (unrestricted) parameter vector
	returns array[8] of {tau (double), vC2obs[m_cX][1],	vC2lat[m_cY2lags+m_fdiffY2][1],
						vbetatilde[m_cY1][1], mC1obs[m_cX][m_cY1], mC1lat[m_cY2lags][m_cY1],
						vdelta[m_cY1][1, mOm1g2Ch[m_cY1][m_cY1]}
**/
{
	decl k = m_cY1+1;				// # of variables	
	decl cpsi2 = 1+m_cX+!m_fCSVAR*m_fLatLags*(m_cY2lags+m_fdiffY2); 	// # of pars in Y2 equation, i.e., count of vpsi2 = tau, vC2obs (of sizec(m_mX)), vC2lat (of size m_cY2lags+m_fdiffY2) if m_fLatLags==TRUE (uncless censoring is imposed)
 	decl vpsi2 = vpall[:cpsi2-1];		
	decl vpsi1 = m_cY1 ? vpall[cpsi2:] : <>;
	decl tau = m_fOlsen ? 1/vpsi2[0] : vpsi2[0];					// RF error st dev (sqrt(omega22))
	if(m_fOlsen)
		vpsi2[1:] *= tau;
	
	decl idx1 = 1, idx2 = idx1+m_cX-m_cY2lags;
	decl vC21 = idx2>idx1 ? vpsi2[idx1:idx2-1] : <>;		// RF coefficients on observable regressors X1 = (exog,Y1lags) in eq. for Y2*
	idx1 = idx2; idx2 = idx1+m_cY2lags;
	decl vC2obs, vC2lat;
	decl vC22bar = idx2>idx1 ? vpsi2[idx1:idx2-1] : <>;		// RF coefficients on observable regressors Y2lags in eq. for Y2*
	if(!m_fCSVAR){
		vC2lat = m_fLatLags ? vpsi2[idx2:] : <>;			// RF coefficients on latent regressors (lagged Y2*) in eq. for Y2*
		vC2obs = vC21 | (vC22bar-(m_fLatLags ? (m_fdiffY2 ?	// RF coefficients on observable regressors m_mX in eq. for Y2*
		dropr(vC2lat,sizeof(vC2lat)-1) : vC2lat) : 0));		// if fdiffY2=TRUE, last element of vC2lat contains coefficient on lagged level and should be removed
	}
	else{
		vC2lat = vC22bar | (m_fdiffY2 ? 1 : <>);			// if censored VAR imposed C2* = C22bar (levels) or C22bar | 1 (diff)
		vC2obs = vC21 | 0*vC22bar;							// .. and C22 = 0
	}
	idx1 = 0; 
	decl vbetabar, vbetatil, dxi = m_dzeta*m_dlambda;
	if(!m_fCSVAR){
		idx2 = m_cY1;
		decl arg = idx2>idx1 ? vpsi1[idx1:idx2-1] : <>;		// betatilde or betabar coefficients in RF
		if(!m_fUseBetabar) vbetatil = arg;	else vbetabar = arg;
	}
	else{
		idx2 = 0;
		vbetatil = zeros(m_cY1,1);							// betabar unidentified in this case, so UseBetabar irrelevant
	}
	idx1 = idx2; idx2 = idx1+m_cY1*m_cX;
	decl mC11, mC12bar, mC1lat, mC1obs;
	if(idx2==idx1)											// if there are no Xs
		mC1lat = mC1obs = <>;		// set all to empty (note, if X=<>, then so is X* in Y1 eqs, by assumption)
	else{
		mC1obs = shape(vpsi1[idx1:idx2-1],m_cX,m_cY1);			// RF coeffs on observable regressors X in Y1 eqs
		mC1lat = <>;
		if(m_fLatLags){
			mC11 = mC1obs[:m_cX-m_cY2lags-1][];					// RF coeffs on observable regressors X1 in Y1 eqs
			mC12bar = m_cY2lags ? mC1obs[m_cX-m_cY2lags:][] : <>;	// RF coeffs on observable regressors X2 in Y1 eqs for uncensored obs
			if(!m_fCSVAR){
				idx1 = idx2; idx2 = idx1+m_cY2lags*(m_cY1);
				mC1lat = m_cY2lags ? shape(vpsi1[idx1:idx2-1],m_cY2lags,m_cY1) : <>;	// RF coeffs on latent regressors X2* in Y1 eqs
			}
			else
				mC1lat = m_cY2lags ? shape(mC12bar,m_cY2lags,m_cY1) : <>;		// RF coeffs on latent regressors X2* in Y1 eqs = C12bar under CSVAR
			mC1obs = mC11 | (mC12bar-mC1lat);					// RF coeffs on observable regressors in Y1 eqs
		}
	}
	idx1 = idx2; idx2 = idx1+m_cY1;
	decl vdelta = idx2>idx1 ? vpsi1[idx1:idx2-1] : <>;		// regression coefficients of u_2 (error in Y2star) in Y1 eqs
	decl mOm1g2Ch = <>; 									// choleski factor of Covariance matrix of u_1 given u_2
	if(m_cY1){
		mOm1g2Ch = unvech(vpsi1[idx2:]);
		mOm1g2Ch = setupper(mOm1g2Ch,0);						// set all elements above diagonal to zero
	}
	if(!m_fUseBetabar)
		return {tau, vC2obs, vC2lat, vbetatil, mC1obs, mC1lat, vdelta, mOm1g2Ch};
	//esle
	decl mOmega = GetOmega(tau, vdelta, mOm1g2Ch);			// variance matrix of reduced form errors for all Y
	vbetatil = (1-dxi) * invert(unit(m_cY1)-dxi*vbetabar*GetGamma(vbetabar, mOmega))*vbetabar;
	return {tau, vC2obs, vC2lat, vbetatil, mC1obs, mC1lat, vdelta, mOm1g2Ch}; 
}  
CKSVAR::SelectIRF(const vIRFs, const mscale)
/**	if multiple IRFs, select the one that is closest to satisfying the sign restrictions
	@param vIRFs [m_iSol*m_cY*(m_cHorizon+1)] vector of vectorized IRFs for the m_iSol different solutions, as returned by funcIRF
	@param mscale [m_iSol*m_cY][m_cHorizon+1] matrix of sign restrictions set by SetIRFsigns(.) 
	@return [m_cY][m_cHorizon+1] matrix of selected IRFs
**/
{
	decl mIRFs = shape(vIRFs,m_iSol*(m_cY1+1), m_cHorizon+1);		// obtain all IRFs and reshape in cY timex cHorizon + 1 matrix (or multiple such matrices concatenated vertically if serveral solutions)
	if(m_iSol > 1){
		decl amIRFs = new array[m_iSol];						// collect each IRF as separate matrix
		decl vsatisfied = new matrix[m_iSol][1];				// check if any of the solutions satisfies sign restrictions
		for(decl isol = 0; isol < m_iSol; ++isol){
			amIRFs[isol] = mIRFs[isol*(m_cY1+1):(isol+1)*(m_cY1+1)-1][];
			vsatisfied[isol] = amIRFs[isol].*mscale >= 0;		// true if all sign restrictions hold
		}
		if(sumc(vsatisfied) == 1)								// if only one solution satisfies restrictions
			mIRFs = amIRFs[vecindex(vsatisfied)];				// return the one that satisfies restrictions
		else if(sumc(vsatisfied) == 0){							// if no solution satisfies restrictions
			for(decl isol = 0; isol < m_iSol; ++isol){
				decl mtest = amIRFs[isol].*mscale;
				mtest = mtest .< 0 .? mtest .: 0;				// only keep the ones that violate restrictions
				vsatisfied[isol] = sumsqrc(vec(mtest));			// return (one-sided) Eucledian norm
			}
			mIRFs = amIRFs[mincindex(vsatisfied)];				// return the one with the smallest Eucledian distance from positive
		}
		else{
			decl vidx = vecindex(vsatisfied);					// collect indices of solutions that satisfy restrictions
			for(decl isol = 0; isol < sizerc(vidx); ++isol){
				decl mtest = amIRFs[vidx[isol]].*mscale;
				vsatisfied[isol] = sumsqrc(vec(mtest));			// return Eucledian norm
			}
			mIRFs = amIRFs[maxcindex(vsatisfied)];				// return the one farthest from the boundary (WARNING: might induce discontinuity during optimization, but not clear if anything else would be continuous)
		}	
	}
	return mIRFs;
}	
CKSVAR::cfunc_gt0(const avF, const vP)
/**	inequality constraints (sign restrictions), must be greater than zero, see MaxSQP for syntax
	@param avF: address of vector, function output
	@param vP: (unrestricted) parameter values
	@return TRUE if successful
**/
{
	if(!m_mSigns)
		return;													// no sign restrictions imposed, function will not be used
	decl mscale = (m_dImpulse > 0 ? 1 : -1) * m_mSigns;		
	decl mIRFs, ret;
	if(m_mrannIRForigin){
		decl mX, mXst, mExog;
		[mX, mXst] = this.SimData(m_mrannIRForigin);
		decl aorigin = this.GetIRForigin();							// get this to restore later
 		mExog = this.GetGroup(X_VAR)[0:min(m_cHorizon,m_cT-1)][];	// set exogenous as they appear in the observed sample
		if(m_cHorizon > m_cT-1)
			mExog |= ones(m_cHorizon-m_cT+1,1)*mExog[m_cT-1][];
		SetIRForigin(/*Origin array*/{mX[0][],mXst[0][], mExog});	// vary the IRF origin over the simulated sample
		ret = funcIRF(&mIRFs, vP);									// evaluate the desired IRFs, store return value
		if(!ret)   													// if funcIRF fails, stop
			return;
		mIRFs = SelectIRF(mIRFs,mscale);
		mIRFs = mIRFs .* mscale;									// change sign if necessary, so they are all positive under sign restrictions
		decl mIRF0 = mIRFs;
		for(decl t = 1; t < sizeof(mX); ++t){
			if(t > m_cT-1)
				mExog = ones(m_cHorizon+1,1)*this.GetGroup(X_VAR)[m_cT-1][];	// keep fixed at last observed value
			else{
				mExog = this.GetGroup(X_VAR)[t:min(t+m_cHorizon,m_cT-1)][];	
				if(m_cHorizon > m_cT-1-t)
					mExog |= ones(m_cHorizon-m_cT+1+t,1)*this.GetGroup(X_VAR)[m_cT-1][];
			}
			SetIRForigin(/*Origin array*/{mX[t][], mXst[t][], mExog});	// vary the IRF origin over the simulated sample
			ret *= funcIRF(&mIRFs, vP);									// evaluate the desired IRFs, store return value 
			if(!ret)   													// if funcIRF fails, stop
				return;
			mIRFs = SelectIRF(mIRFs,mscale);
			mIRFs = mIRFs .* mscale;									// change sign if necessary, so they are all positive under sign restrictions
			mIRFs = mIRFs .< mIRF0 .? mIRFs .: mIRF0;					// return the minimum of this and previous value
			mIRF0 = mIRFs;												// reset last period's IRF
		}
		SetIRForigin(aorigin);										// restore previous selection
	}
	else{
		ret = funcIRF(&mIRFs, vP);									// evaluate the desired IRFs, store return value 
		if(!ret)   													// if funcIRF fails, stop
			return;
		mIRFs = SelectIRF(mIRFs,mscale);
		mIRFs = mIRFs .* mscale;									// change sign if necessary, so they are all positive under sign restrictions
	}
	avF[0] = selectifr(vec(mIRFs),vec(mscale));	 					// select only those IRFs with nonzero value after scaling (i.e., drop all those with unrestricted sign)
	return ret;
}
CKSVAR::Covar()
/**	Modelbase function that computes the estimates of the covariance matrix of the parameter estimates
	the result is stored in the m_mCovar data member
**/
{
	decl mHessian;
	Num2Derivative(this.LogLik, this.GetFreePar(), &mHessian);
	m_mCovar = invertgen(-mHessian,30);						// use generalized inverse for the case that Hessian is singular, e.g., when you do OLS and the model is over-parametrized
}
CKSVAR::SetSimAnn(const fSimAnn, const aStep)
/**	Set Simulated Annealing parameters if needed
	@param fSimAnn, flag, if TRUE, use MaxSA in DoEstimation
	@param aStep, array[6] of MaxSA parameters: {iNS, iNT, dRt, vM, vC, dT}
**/
{
  	m_fSimAnn = fSimAnn;
	if(m_fSimAnn){
		decl iNS = aStep[0];
		decl iNT = aStep[1];
		decl dRt = aStep[2];
		decl vM = aStep[3];
		decl vC = aStep[4];
  		MaxSAControlStep(iNS, iNT, dRt, vM, vC);
		m_dTemp = aStep[5];
	}
}

CKSVAR::DoEstimation(vP)
/**	Modelbase function that is needed to estimate the model
	In: vP: vector of free (unrestricted) parameters
**/
{
	SetNumMCreps(m_cR);
	decl ires;
	if(!m_fbounded)
		ires = m_fSimAnn ? MaxSA(this.LogLik, &vP, &m_dLogLik, &m_dTemp) :	MaxBFGS(this.LogLik, &vP, &m_dLogLik, 0, TRUE);
	else{
		decl vlo, vhi;
		vlo = m_vLo[m_vIdxFreePar];
		vhi =  m_vHi[m_vIdxFreePar];
//		ires = MaxSQP(this.LogLik, &vP, &m_dLogLik, 0, /*fNumerical*/TRUE, /*Ineq Restr*/this.cfunc_gt0,/*Equality Restr*/0, vlo, vhi);
		ires = MaxSQPF(this.LogLik, &vP, &m_dLogLik, 0, /*fNumerical*/TRUE, /*Ineq Restr*/this.cfunc_gt0,/*Equality Restr*/0, vlo, vhi);
	}	
	SetResult(ires);
    return vP;
}
CKSVAR::GetParNames()
{
 	// set names
	decl k = columns(m_mY);				// # of variables

	decl asNames2 = new array[1+m_cX+!m_fCSVAR*m_fLatLags*(m_cY2lags+m_fdiffY2)];		// names of parameters in Y2 equation = tau, observed regressors, and m_cY2lags latent lags (+1 if Y2 in first difference)
	asNames2[0] = m_fOlsen ? "$1/\\tau$" : "$\\tau$";
	decl asx0, asY1lags, asY2lags;
	this.GetGroupNames(X_VAR,&asx0);
	this.GetGroupLagNames(Y1_VAR,1, this.GetMaxGroupLag(Y1_VAR), &asY1lags);
	this.GetGroupLagNames(Y2_VAR,1, this.GetMaxGroupLag(Y2_VAR), &asY2lags);
	decl asx = asx0~asY1lags~asY2lags;			// all observed exogenous regressors
	for(decl j = 0; j < m_cX; ++j)
		asNames2[j+1] = sprint("Eq.", k, " ", asx[j]);	// append equation number
	decl aslatent = {};
	if(m_fLatLags && !m_fCSVAR){
		if(m_fdiffY2)								// need an extra name
			this.GetGroupLagNames(Y2_LAGLEV, 1, 1, &aslatent);
		aslatent = asY2lags~aslatent;
		for(decl j = 0; j < sizeof(aslatent); ++j)
			asNames2[j+m_cX+1] = sprint("Eq.", k, " l", aslatent[j]);		// names of coefficients on latent variables, indicated with prefix "l"
	}

	decl asbeta = m_fUseBetabar ? "\\bar{\\beta}" : "\\tilde{\\beta}";
	m_asbeta = new array[k-1];				// names of beta parameters in Y1 equation(s) 
	for(decl j = 0; j < sizeof(m_asbeta); ++j)
		m_asbeta[j] = k==2 ? sprint("$",asbeta,"$") : sprint("$", asbeta, "_", j+1,"$");		// index refers to equation
	decl asCobs = new array[(k-1)*m_cX];		// coefficients on X in Y1 eqs
	for(decl ieq = 0; ieq < k-1; ++ieq){
		for(decl j = 0; j < m_cX; ++j){
			asCobs[ieq*m_cX+j] = sprint("Eq.", ieq+1, " ", asx[j]);
		}
	}
	decl asClat = {};
	if(m_fLatLags && !m_fCSVAR){
		asClat = new array[(k-1)*m_cY2lags];	// coefficients on latent lags in Y1 eqs 
		for(decl ieq = 0; ieq < k-1; ++ieq){
			for(decl j = 0; j < m_cY2lags; ++j)
				asClat[ieq*m_cY2lags+j] = sprint("Eq.", ieq+1, " l", asY2lags[j]);
		}
	}
	decl asdelta = new array[k-1];				// delta (k-1) parameters 
	for(decl j = 0; j < sizeof(asdelta); ++j)
		asdelta[j] = k==2 ? "$\\delta$" : sprint("$\\delta_", j+1,"$");		// index refers to equation
	decl asChol = new array[(k-1)*k/2];	// choleski of Omega1g2 (k-1)*k/2 
	decl aidx = new array[k];			// array of indices
	aidx[0] = 0;
	for(decl j = 1; j < sizeof(aidx); ++j)
		aidx[j] = aidx[j-1] + k-j;
	for(decl j = 0; j < sizeof(asChol); ++j){
		for(decl i = 0; i < sizeof(aidx)-1; ++i)
			if(j >= aidx[i] && j < aidx[i+1])
				asChol[j] = sprint("Ch_", j-aidx[i]+i+1, i+1);
	}
	if(!m_fCSVAR){							// if you're not imposing CSVAR, keep track of indices of latent and observed Y2 lags for TestObservedEqualLatentLags() and TestPureCensoredVAR()
		if(m_fLatLags)
			m_vIdxLatLags = (1+m_cX+range(0,m_cY2lags+m_fdiffY2-1)')
						| (m_cY1 && sizeof(asClat) ? (sizeof(asNames2~m_asbeta~asCobs)
						+ range(0,sizeof(asClat)-1)') : <>); 			// record position of coefficients on latent lags
		decl cX1 = m_cX-m_cY2lags;										// # of observed X's and lags of Y1
		m_vIdxObsLags = 1 + cX1 + range(0,m_cY2lags-1)';
		if(m_cY1){
			for(decl j = 0; j < m_cY1; ++j)
				m_vIdxObsLags |= sizeof(asNames2~m_asbeta)+j*m_cX+cX1+range(0, m_cY2lags-1)';
			m_vIdxbetatilde = sizeof(asNames2) + range(0,m_cY1-1);
		}
	}
	decl asnames = m_fCSVAR ? asNames2 ~ asCobs ~ asdelta ~ asChol :			// if you impose censoring, you don't need betatilde and c*
			asNames2~m_asbeta~asCobs~asClat~asdelta~asChol;
	for(decl i = 0; i < m_cBreaks; ++i)
		for(decl j = 0; j < sizeof(m_aBreakParIdx[i]); ++j)
			asnames |= sprint(asnames[m_aBreakParIdx[i][j]],"_",GetDateByIndex(this.GetSelStart()+m_vBreakDates[i]));

	return asnames;
}
CKSVAR::Reparameterize(const avPall, const vPfree)
/**	Map free parameters to full parameter vector
	@param avPall address, 	out: vector [m_cPar]
	@param vPfree vector of free parameters (not necessarily m_cParFree)
	@returns 1 or 0	**/
{
	decl vp, ret = 1;	  
//	if(m_frestricted && (m_iRestrict == RES_LIN || m_iRestrict == RES_GEN))		// applies only to restrictions for testing
//	{	
//		if(m_iRestrict == RES_LIN)		/** 1. check if any linear restrictions **/
//			vp = sizerc(vPfree) ?  m_mS * vPfree + m_vS : m_vS;
//		else if(m_iRestrict == RES_GEN)
//		{	decl vpr;
//			vp = m_vFixedPar;
//			vp[selectr(m_vIdxUnres,m_vIdxFreePar)] = vPfree;
//			ret = GeneralRestrictions(&vpr, vp[m_vIdxUnres]);	// allows for FixPar to be combined with restrictions
//			if(sizerc(matrix(vpr)) != m_cRestrict)
//			{	if(m_fPrintError) println("*** Error in setting GeneralRestrictions.");
////				println(vpr, vPfree);
//				return 0;
//			}
//			vp[m_vIdxRestr] = vpr;
//		}
//	}
//	else					/** check if any parameters are fixed	**/
	{	if(sizerc(matrix(vPfree)) != sizerc(m_vIdxFreePar)) // check that dimensions are correct
		{	println("*** Error in Reparameterize. sizerc(vPfree) = ", sizerc(vPfree),
					" # of free parameters = ", sizerc(m_vIdxFreePar));
			return 0;
		}
		vp = m_vFixedPar;			// full parameter vector
		vp[m_vIdxFreePar] = vPfree;
	}
	avPall[0] = vp;
	return ret;
}
CKSVAR::Fixbeta(const vbeta)
/**	Fixes the value of beta in the model, removes it from FreePar list
	@param vbeta, vector of values of beta
	@returns 1 or 0	**/
{
	decl vidx = strfind(this.GetParNames(), m_asbeta);			// locate the betas within the parameter vector
	if(any(vidx .== -1))										// if any of the indices is negative, skip it
		return;
	for(decl i = 0; i < sizerc(vidx);++i)
		FixPar(vidx[i], vbeta[i]);
}
CKSVAR::Freebeta()
/**	Frees beta in the model, i.e., adds it to list of free parameters
	@param vbeta, vector of values of beta
	@returns 1 or 0	**/
{
	decl vidx = strfind(this.GetParNames(), m_asbeta);			// locate the betas within the parameter vector
	for(decl i = 0; i < sizerc(vidx);++i)
		FreePar(vidx[i]);
}
CKSVAR::GetOLS()
/**	 Returns OLS estimates	of reduced-form VAR coefficients (ignoring the ZLB) **/
{
    if (m_iModelStatus < MS_DATA)
        InitData();
    if (m_iModelStatus < MS_DATA)                              /* need data */
        return FALSE;
	decl mC,xxinv;											// coefficients
	olsc(m_mY, m_mX, &mC,&xxinv);
	decl Omega = variance(m_mY-m_mX*mC);				// OLS estimate	of RF error variance
	decl u = m_mY-m_mX*mC;
	decl mV = variance(u)**xxinv;
	return {mC, Omega};
}
CKSVAR::SimLik(const vPar, const fSimulate)
/** Returns m_cR draws of {Y2*_t} from its conditional distribution given Y1_t, Y2_t and lags, IS weights and joint density of Y,Y2*
	@param vPar: [m_cPar][1] vector of model parameters
	@param fSimulate: flag 
	returns array {vLik, mY2star, mweights} where vLik [m_cT][1], prediction error density of Y_t|t-1 per observation
			and, if fSimulate == TRUE, mY2star, mweights are [m_cT][m_cR] matrices of draws of Y2star and IS weights so that, e.g., Y2*_{t|t} = sumr(mY2star.*mweights), Y2*_{t|T} = sumr(mY2star.*mweights[T-1][])
			otherwise mY2star, mweights are empty matrices
**/
{
	decl vP = m_cBreakPar ? vPar[:(m_cPar-m_cBreakPar-1)] : vPar;					// get parameters before the first break (default regime)
	decl tau, vC2obs, vC2lat, vbetatil, mC1obs, mC1lat, vdelta, mOm1g2Ch;			// obtain the parameters first
	[tau, vC2obs, vC2lat, vbetatil, mC1obs, mC1lat, vdelta, mOm1g2Ch] = this.GetModelPar(vP);
	decl tau2, arg2, mOm1g2inv, mXi1inv;
	[mOm1g2inv, mXi1inv, tau2, arg2] = GetAuxPars(mOm1g2Ch, vbetatil, vdelta, tau);
	if(mOm1g2inv == 0 || mXi1inv == 0)			// if either inversion failed return nothing
		return {.NaN, <>, <>};
	decl apars = new array[m_cBreaks+1];
	decl atau, avC2obs, avC2lat, avbetatil, amC1obs, amC1lat, avdelta,
	amOm1g2Ch, atau2, aarg2, amOm1g2inv, amXi1inv;							// array to hold parameter values in different regimes when there are breaks
	atau = avC2obs = avC2lat = avbetatil = amC1obs = amC1lat = avdelta = amOm1g2Ch =
	atau2 = aarg2 = amOm1g2inv = amXi1inv = apars;					// array to hold parameter values in different regimes when there are breaks
	if(m_cBreaks)
		apars = this.GetBreakPar(vPar);			// get parameter values at different regimes
	[atau[0], avC2obs[0], avC2lat[0], avbetatil[0], amC1obs[0], amC1lat[0], 
		avdelta[0], amOm1g2Ch[0], amOm1g2inv[0], amXi1inv[0], atau2[0], aarg2[0]] =
	{tau, vC2obs, vC2lat, vbetatil, mC1obs, mC1lat, 
			vdelta, mOm1g2Ch, mOm1g2inv, mXi1inv, tau2, arg2};		  // first regime is the benchmark
	for(decl i = 1; i < m_cBreaks+1; ++i)
		[atau[i], avC2obs[i], avC2lat[i], avbetatil[i], amC1obs[i], amC1lat[i], 
		avdelta[i], amOm1g2Ch[i], amOm1g2inv[i], amXi1inv[i], atau2[i], aarg2[i]] = apars[i-1];

	if(!fSimulate){							// if you only need to compute exact likelihood without latent lags, avoid importance sampling step for much faster execution	(faster than setting R=1, because we avoid the for loop over T)
		decl mmu1, vmu2, mmu1b, vmu2b;
		mmu1 = mmu1b = new matrix[m_cT][m_cY1];
		vmu2 = vmu2b = new matrix[m_cT][1];
		for(decl i = 0; i < m_cBreaks + 1; ++i){
			decl idxr = m_aRegimeDates[i];			// for short
			mmu1[idxr][] = m_cY1 ? (m_mX[idxr][]*amC1obs[i]) : <>;					// mean of Y1
			vmu2[idxr][] = m_mX[idxr][]*avC2obs[i];									// mean of Y2*
			mmu1b[idxr][] = m_cY1 ? mmu1[idxr][] + (m_vb[idxr]-vmu2[idxr][])*avbetatil[i]' : <>;	// mean of Y1 when Y2 = b
			vmu2b[idxr][] = vmu2[idxr][] + (m_cY1 ? (m_mY1[idxr][]-mmu1b[idxr][])*aarg2[i]' : 0);	// mean of Y2 | Y1 and lags (when Y2=b)
		}
		decl vlik1 = 1;											// contribution of Y1 to the likelihood
		if(m_cY1){
			decl mmu1g2 = new matrix[m_cT][m_cY1];				// mean of Y1 | Y2 when Y2>b
			decl marg1gb = new matrix[m_cT][m_cY1];				// contribution of Y1 to lik when Y2>b
			decl marg1gbSt = new matrix[m_cT][m_cY1];			// standardized varg1gb
			for(decl i = 0; i < m_cBreaks + 1; ++i){
				decl idxr = m_aRegimeDates[i];			// for short
				mmu1g2[idxr][] = mmu1[idxr][] + (m_mY2[idxr][]-vmu2[idxr][])*avdelta[i]';			
				marg1gb[idxr][] = m_mY1[idxr][] - mmu1g2[idxr][];					
				marg1gbSt[idxr][] = marg1gb[idxr][]*amOm1g2inv[i];
			}
//			println("(m_mY1 - mmu1g2)", (m_mY1 - mmu1g2)');
//			println("mOm1g2inv", mOm1g2inv);
//			println(mOm1g2inv~amOm1g2inv);
			decl varg1gb = sumr(marg1gb .* marg1gbSt);			// kernel of contribution of Y1 to lik when Y2>b, i.e., (Y1-mu1g2)'Om1g2^(-1)*(Y1-mu1g2)		
			// now move to case Y2 = b
			decl marg1b = new matrix[m_cT][m_cY1];				// contribution of Y1 to lik when Y2=b
			decl marg1bSt = new matrix[m_cT][m_cY1];			// standardized varg1b
			for(decl i = 0; i < m_cBreaks + 1; ++i){
				decl idxr = m_aRegimeDates[i];			// for short
				mmu1[idxr][] += (m_vb[idxr][]-vmu2[idxr][])*avbetatil[i]';	// mean of Y1 when Y2 = b
				marg1b[idxr][] = (m_mY1[idxr][] - mmu1[idxr][]);							
				marg1bSt[idxr][] = marg1b[idxr][]*amXi1inv[i];
			}
			decl varg1b = sumr(marg1b .* marg1bSt);			// kernel of contribution of Y1 to lik when Y2=b, i.e., (Y1-mu1)'Xi^(-1)*(Y1-mu1)
			// finally, put together the likelihood function
//			vlik1 = m_mY2 .> m_vb .? sqrt(determinant(mOm1g2inv))*exp(-varg1gb/2)
//											.: sqrt(determinant(mXi1inv))*exp(-varg1b/2);	// contribution of Y1 to the likelihood
//			println("varg1gb", varg1gb', "varg1b", varg1b');
			decl sign;
			vlik1 = new matrix[m_cT][1];
			for(decl i = 0; i < m_cBreaks + 1; ++i){
				decl idxr = m_aRegimeDates[i];			// for short
				vlik1[idxr] = m_mY2[idxr] .> m_vb[idxr] .? exp(logdet(amOm1g2inv[i],&sign)/2)*exp(-varg1gb[idxr]/2)
											.: exp(logdet(amXi1inv[i],&sign)/2)*exp(-varg1b[idxr]/2);		// contribution of Y1 to the likelihood
			}
			if(sign != 1)
				return {.NaN, <>, <>};				// if the determinant is zero, negative or unreliable (0,-1,-2 or 2), return nothing
		}
//		decl varg2gb = new matrix[m_cT][1];
		decl vlik2 = new matrix[m_cT][1];
		for(decl i = 0; i < m_cBreaks + 1; ++i){
			decl idxr = m_aRegimeDates[i];			// for short
//			varg2gb[idxr] = (m_mY2[idxr] - vmu2[idxr])/(i ? atau[i-1] : tau);			// kernel of contribution of Y2 to the likelihood when Y2>b	 (gb stands for greater than b)
			decl varg2gb = (m_mY2[idxr] - vmu2[idxr])/atau[i];			// kernel of contribution of Y2 to the likelihood when Y2>b	 (gb stands for greater than b)
			// now move to case Y2 = b
			decl vmu2i;
			if(m_cY1)
				vmu2i = vmu2[idxr]+(m_mY1-mmu1)[idxr][]*aarg2[i]';				// mean of Y_2* given Y1 when Y2* < b, stacked over all R MC draws, [k-1][R]			 
			decl varg2b = (m_vb[idxr] - vmu2i)/atau2[i];			// kernel of contribution of Y2 to the likelihood when Y2=b, i.e., Pr(Y2*<b|Y1,lags)	

			vlik2[idxr] = m_mY2[idxr] .> m_vb[idxr] .? densn(varg2gb)/atau[i] .: probn(varg2b);		// contribution of Y2 to the likelihood
		}
		return {vlik1.*vlik2, <>, <>};
	}
	// else
	decl mY2star = new matrix[m_cT][m_cR];				// to store simulated m_cR draws of Y2* for the censored observations
	decl mw = ones(m_cT, m_cR);							// IS weights, initialized at 1

	decl vS = new matrix[m_cT][1];						// predictive density f(Y_t|Y_{t-1})
	decl vlaglevst, vlaglev, mY1meanObs, vY2stmeanObs;
	decl mXstar = new matrix[m_cY2lags][m_cR];			// latent regressors (lags of Y2*) Xt*, [m_cY2lags][R]
	if(m_fFAPF && m_fLatLags){							// if FAPF, initialize X* 
		for(decl j = 0; j < m_cY2lags; ++j)				// not executed if m_cY2lags = 0
			mXstar[j][] = m_vY2stinit[][j]*ones(1,m_cR);		
	}
	for(decl t = 0; t< m_cT; ++t){
		//if there are breaks, set time-varying parameters
		if(m_cBreaks){
			decl i = sumc(m_vBreakDates .<= t);
			[tau, vC2obs, vC2lat, vbetatil, mC1obs, mC1lat,	vdelta,
					mOm1g2Ch, mOm1g2inv, mXi1inv, tau2, arg2] =
			{atau[i], avC2obs[i], avC2lat[i], avbetatil[i], amC1obs[i], amC1lat[i], 
			avdelta[i], amOm1g2Ch[i], amOm1g2inv[i], amXi1inv[i], atau2[i], aarg2[i]};				
		}
		decl b = m_vb[t];
		vlaglevst = m_fdiffY2 ? (t ? mY2star[t-1][] : m_vY2stinit[][m_cY2lags]*ones(1,m_cR)) : <>;
		vlaglev = m_fdiffY2 ? vlaglevst-m_vY2laglev[t] : <>; 

		decl mY1meanObs = m_cY1 ? (m_mX[t][]*mC1obs)'ones(1,m_cR) : 0;		// to hold mean of Y1 coming from observed regressors, stacked over all R MC draws, [m_cY1][R]	  (if m_cY1=0, set to zero)
		decl vY2stmeanObs = ones(1,m_cR)*(m_mX[t][]*vC2obs);				// to hold mean of Y2* coming from observed regressors, stacked over all R MC draws, [1][R]

		decl mY1mean = mY1meanObs;											// to hold mean of Y1, stacked over all R MC draws, [m_cY1][R]	  (if m_cY1=0, set to zero)
		decl vY2stmean = vY2stmeanObs;										// to hold mean of Y2*, stacked over all R MC draws, [1][R]
		if(m_fLatLags){
			if(!m_fFAPF){
				if(!m_fdiffY2)
					for(decl j = 0; j < m_cY2lags; ++j)				// not executed if m_cY2lags = 0
						mXstar[j][] = t>j ? mY2star[t-j-1][] : m_vY2stinit[][j-t]*ones(1,m_cR);		
				else
					for(decl j = 0; j < m_cY2lags; ++j)				// not executed if m_cY2lags = 0
						mXstar[j][] = t>j ? mY2star[t-j-1][] - (t-j>1 ? mY2star[t-j-2][] : m_vY2stinit[][m_cY2lags])
								: m_vY2stinit[][j-t]*ones(1,m_cR);		
			}		
			if(m_cY1)
				mY1mean += mC1lat'mXstar;
			vY2stmean += vC2lat'(mXstar | vlaglev);
		}
		decl mY1meanb = m_cY1 ? mY1mean + vbetatil.*(b-vY2stmean) : 0;		// mean of Y1 when Y2 = b, stacked over all R MC draws, [m_cY1][R]
		decl vY2stmeanb = vY2stmean + (m_cY1 ?
							arg2*(m_mY1[t][]'ones(1,m_cR)-mY1meanb) : 0);	// mean of Y2* | Y1 and lags (when Y2=b)	(assuming kappa=1, which is innocuous, since it drops out in the relevant calculations)
		decl vlik1_t=1;													// contribution of Y1 to the likelihood (only compute if m_cY1>0)
		if(!m_fFAPF){
			decl veta = m_mY2[t] > b ? .NaN : tau2*quann(m_mxi[t][].*probn((b-vY2stmeanb)/tau2));	 // draw shocks from distribution of Y2* | Y1 and lags
			mY2star[t][] = (m_mY2[t] > b ? m_mY2[t]*ones(1,m_cR) : m_dkappa*(vY2stmeanb+veta)+(1-m_dkappa)*b)
						+ (m_fdiffY2 ? m_vY2laglev[t] : 0);							// add lagged level if Y2 is in first difference
		}
		if(m_cY1){
			decl mY1gY2mean = mY1mean + vdelta*(m_mY2[t]-vY2stmean);	// mean of Y1 | Y2 when Y2>b
			decl varg1gb = (m_mY1[t][]'ones(1,m_cR) - mY1gY2mean);		// kernel of contribution of Y1 to lik when Y2>b..
			varg1gb = sumc(varg1gb .* (mOm1g2inv*varg1gb));				// ..., i.e., (Y1-E(Y1|Y2))'Om1g2^(-1)*(Y1-E(Y1|Y2)), stacked in [1][R] vector		
			// now move to case Y2 = b
			decl varg1b = (m_mY1[t][]'ones(1,m_cR) - mY1meanb);			// kernel of contribution of Y1 to lik when Y2=b..
			varg1b = sumc(varg1b .* (mXi1inv*varg1b));					// ... i.e., (Y1-E(Y1))'Xi^(-1)*(Y1-E(Y1)), stacked in [1][R] vector		

			// finally, put together the likelihood function for all MC replications [1][R]
//			vlik1_t = m_mY2[t] > b ? sqrt(determinant(mOm1g2inv))*exp(-varg1gb/2)
//										: sqrt(determinant(mXi1inv))*exp(-varg1b/2);	// contribution of Y1 to the likelihood
			decl sign;
			vlik1_t = m_mY2[t] > b ? exp(logdet(mOm1g2inv,&sign)/2)*exp(-varg1gb/2)
										: exp(logdet(mXi1inv,&sign)/2)*exp(-varg1b/2);	// contribution of Y1 to the likelihood
//			if(fabs(sign)>1 || sign == 0){		// if any of the two is "unreliable", print the result -- ONLY FOR DEBUGGING, REMOVE IN FINAL VERSION! ************
//				if(m_mY2[t] > b)
//					println("mOm1g2inv", mOm1g2inv, "determinant(mOm1g2inv)", exp(logdet(mOm1g2inv,&sign))~determinant(mOm1g2inv), "sign = ", sign,
//					"%c", {"min lik1", "max lik1"}, min(vlik1_t)~max(vlik1_t), "%c", {"m_cZLB", "m_cPar", "1+m_fLatLags*cY2lat"}, m_cZLB~m_cPar~(1+m_fLatLags*(m_cY2lags+m_fdiffY2)));
//				else
//					println("mXi1inv", mXi1inv, "determinant(mXi1inv)", exp(logdet(mXi1inv,&sign))~determinant(mXi1inv), "sign = ", sign,
//					"%c", {"min lik1", "max lik1"}, min(vlik1_t)~max(vlik1_t), "%c", {"m_cZLB", "m_cPar", "1+m_fLatLags*cY2lat"}, m_cZLB~m_cPar~(1+m_fLatLags*(m_cY2lags+m_fdiffY2)));
//			}
			if(sign != 1)
				return {.NaN, <>, <>};				// if the determinant is zero, negative or unreliable (0,-1,-2 or 2), return nothing

		}
		decl varg2gb = (m_mY2[t] - vY2stmean)/tau;					// kernel of contribution of Y2 to the likelihood when Y2>b	 (gb stands for greater than b)
		decl varg2b = (b - vY2stmeanb)/tau2;						// kernel of contribution of Y2 to the likelihood when Y2=b, i.e., Pr(Y2*<b|Y1,lags)	
	
		decl vlik2_t = m_mY2[t] > b ? densn(varg2gb)/tau : probn(varg2b);				// contribution of Y2 to the likelihood

		decl vlik_t = vlik1_t.*vlik2_t;
		if(!m_fFAPF){
			vS[t] = meanr(vlik_t.*mw[t][]);
			if(t<m_cT-1)
				mw[t+1][] = vlik_t .* mw[t][]/vS[t];
		}
		else{
			mw[t][] = vlik_t;											// w_{t-1|t}^j, j = 1:m_cR
			decl vprob = vlik_t ./ sumr(vlik_t);						// pi_{t-1|t}^j, j = 1:m_cR
			decl vc = ranmultinomial(m_cR, vprob);
			decl vidx = new matrix[m_cR][1];
			decl j0 = 0;
			for(decl j = 0; j < max(vc); ++j){
				decl arg = vecindex(vc.>j);
				decl carg = sizeof(arg);
				vidx[j0:j0+carg-1] = arg;
				j0 += carg;
			}
			// FAPF
			mXstar = mXstar[][vidx];									// resample the X* to get Xtilde*
			if(m_fdiffY2){
				vlaglevst = vlaglevst[][vidx];							// resample Y2*_1 to get Y2*tilde_1
				vlaglev = vlaglevst - m_vY2laglev[t];					// same for the lagged level if present
			}
	
			if(m_cY1)
				mY1mean = mY1meanObs + mC1lat'mXstar;					// recompute the means using Xtilde*
			vY2stmean = vY2stmeanObs + vC2lat'(mXstar | vlaglev);
	
			// then get mutilde_1 (mean of Y1 | Y2=b and Xstartilde rather than Xstar)
			mY1mean = m_cY1 ? (m_mX[t][]*mC1obs)'ones(1,m_cR) : 0;		// to hold mean of Y1, stacked over all R MC draws, [m_cY1][R]	  (if m_cY1=0, set to zero)
			mY1meanb = m_cY1 ? mY1mean + vbetatil.*(b-vY2stmean) : 0;		// mean of Y1 when Y2 = b, stacked over all R MC draws, [m_cY1][R]
			if(m_cY1)
				mY1mean += mC1lat'mXstar;
			mY1meanb = m_cY1 ? mY1mean + vbetatil.*(b-vY2stmean) : 0;		// mean of Y1 when Y2 = b, stacked over all R MC draws, [m_cY1][R]
			vY2stmeanb = vY2stmean + (m_cY1 ?
								arg2*(m_mY1[t][]'ones(1,m_cR)-mY1meanb) : 0);		// mean of Y2* | Y1 and lags (when Y2=b)	(assuming kappa=1, which is innocuous, since it drops out in the relevant calculations)
			decl veta = m_mY2[t] > b ? .NaN : tau2*quann(m_mxi[t][].*probn((b-vY2stmeanb)/tau2));	 // draw shocks from distribution of Y2* | Y1 and lags
			mY2star[t][] = (m_mY2[t] > b ? m_mY2[t]*ones(1,m_cR) : m_dkappa*(vY2stmeanb+veta)+(1-m_dkappa)*b)
						+ (m_fdiffY2 ? m_vY2laglev[t] : 0);							// add lagged level if Y2 is in first difference
			for(decl j = m_cY2lags-1; j > 0; --j)							 // reset Xstar
				mXstar[j][] = mXstar[j-1][];
			if(m_cY2lags)
				mXstar[0][] = m_fdiffY2 ? mY2star[t][]-vlaglevst : mY2star[t][];			// first element of X* is either Y2* or DY2* depending on m_fdiffY2
			
		}
		if(t==m_cOrigin) 												// if you require draws for computation of IRF at origin 0<=m_cOrigin<=T-1...
			m_mXstar = mXstar;											// .. store them in m_mXstar to be used in funcIRF
	}
	return {m_fFAPF ? meanr(mw) : vS, mY2star, mw};
}
CKSVAR::GetWeights(const fPlot)
/** Plots the weights and the effective sample size of the particle filter algorithm.
	@param fPlot, boolean, plot if TRUE
	returns effective sample size
**/
{
	decl mw = SimLik(m_vPar, TRUE)[2];
	mw ./= sumr(mw);										// convert to probabilities
	decl ess = 1 ./ sumsqrr(mw);							// effective sample size for each t
	if(fPlot){
		decl vidx = round(0~range(1,10)*m_cT/10-1);				// plot weights at every tenth of sample.
		decl aidx = new array[sizerc(vidx)];
		for(decl j = 0; j < sizerc(aidx); ++j)
			aidx[j] = sprint("t=", vidx[j]+1);
		SetDrawWindow("Weights");
		DrawTMatrix(0, mw[vidx][], aidx); 
		println("Smallest effective sample size: ", min(ess));
		DrawMatrix(1, ess', "ESS", 1, 1);
		ShowDrawWindow();
	}
	return ess;
}
CKSVAR::TestExcludeLatentLags()
/** Performs Wald test of exclusion of all latent lags from the model
	Prints test result
	returns array {Wald statistic, p-value}
**/
{
	if(!m_fLatLags){											// if there are no latent lags in the model, return nothing
		println("TestExcludeLatentLags(): No latent lags to exclude.");
		return;
	}//else
	decl vSel = sumr(unit(m_cPar)[][m_vIdxLatLags]);			// vector with 1 in position of coefficients on latent lags to be excluded
	return TestRestrictions(vSel);
}
CKSVAR::TestObservedEqualLatentLags()				
/** Performs Wald test of equality of coefficients on observed and latent lags of Y2 (necessary but not sufficient condition for purely censored VAR)
	Prints test result
	returns array {Wald statistic, p-value}
**/
{
	if(!m_fLatLags){											// if there are no latent lags in the model, return nothing
		println("TestObservedEqualLatentLags(): No latent lags to test.");
		return;
	}//else
	this.GetParNames();					// call this function in order to make sure m_vIdxLatLags has been set, because this is only set in GetParNames(), and the latter is not called when m_fPrint==FALSE, so it would produce a runtime error below
	decl vIdxLatLags = m_fdiffY2 ? dropr(m_vIdxLatLags, m_cY2lags)	// if contrained variable is in first difference, drop the coef in the lagged latent level in the Y2 eqn, because it is not constrained to be zero
				 : m_vIdxLatLags;									// otherwise, use them all
	decl mR = unit(m_cPar)[vIdxLatLags][]-unit(m_cPar)[m_vIdxObsLags][];	 // restrictions of the form mR*b=0
	decl vr = zeros(sizer(mR),1);
	if(m_fdiffY2){
		mR |= unit(m_cPar)[m_vIdxLatLags[m_cY2lags]][];					// m_fdiffY2 == TRUE, also impose restriction that coef on lagged latent lever = 1.
		vr |= 1;
	}
	decl fPrint = m_fPrint;
	this.SetPrint(FALSE); 											// disable printing to avoid printing R and r which could be big
	decl atest = TestRestrictions(mR, vr);							// array that holds {test stat, pval}
	this.SetPrint(fPrint);
	if(fPrint)
		PrintTestVal(atest[0], sizer(mR), 0, "Wald test of equality of coefficients on observed and latent lags");
	return atest;
}
CKSVAR::TestPureCensoredVAR()				
/** Performs Wald test of pure CSVAR (no kink). This entails betatilde = 0 and equality of coefficients on observed and latent lags of Y2. 
	Prints test result
	returns array {Wald statistic, p-value}
**/
{
	if(m_cY2lags && !m_fLatLags){					// if there are observed Y2 lags and latent lags in the model, the hypothesis cannot be true (unless the observed coefs on Y2lags = 0, which is an additional restriction)
		println("TestPureCSVAR(): No latent lags, CSVAR hypothesis cannot hold.");
		return;
	}//else

	decl mR = unit(m_cPar)[m_vIdxbetatilde][];	 // restrictions of the form mR*b=0
	decl vr = zeros(sizer(mR),1);
	if(m_fLatLags){
		decl vIdxLatLags = m_fdiffY2 ? dropr(m_vIdxLatLags, m_cY2lags)	// if contrained variable is in first difference, drop the coef in the lagged latent level in the Y2 eqn, because it is not constrained to be zero
					 : m_vIdxLatLags;									// otherwise, use them all
		mR |= unit(m_cPar)[vIdxLatLags][]-unit(m_cPar)[m_vIdxObsLags][]; 
		vr |= zeros(sizerc(vIdxLatLags),1);
		if(m_fdiffY2){
			mR |= unit(m_cPar)[m_vIdxLatLags[m_cY2lags]][];				// m_fdiffY2 == TRUE, also impose restriction that coef on lagged latent lever = 1.
			vr |= 1;
		}
	}
	decl fPrint = m_fPrint;
	this.SetPrint(FALSE); 											// disable printing to avoid printing R and r which could be big
	decl atest = TestRestrictions(mR, vr);							// array that holds {test stat, pval}
	this.SetPrint(fPrint);
	PrintTestVal(atest[0], sizer(mR), 0, "Wald test of purely censored VAR");
	return atest;
}
CKSVAR::TestAgainstUnrestr()
/** Perform LR test against unrestricted model stored in m_oUnrestricted. Should be called after SetUnrestrictedModel, otherwise it does nothing
**/
{
	if(m_iModelStatus < MS_ESTIMATED){				// you need to do this after Estimate()
		println("TestAgainstUnrestr(): Estimate restricted model first.");
		return;
	}
	if(isnan(m_oUnrestricted)){						// if unrestricted model object hasn't been stored, nothing to test against
		println("TestAgainstUnrestr(): Set unrestricted model object first using SetUnrestrictedModel(<object>).");
		return FALSE;
	}
	if(m_oUnrestricted.GetModelStatus()<MS_ESTIMATED){
		if(m_fPrint)
			println("\n**** UNRESTRICTED ESTIMATION *****\n");
		m_oUnrestricted.Estimate();					// estimate unretricted model, unless it was already done
	}
	decl dLR = m_oUnrestricted.GetLogLik();			// obtain unrestricted log likelihood
	decl df = m_oUnrestricted.GetFreeParCount();	// count # of free parameters in unrestricted model
	dLR = 2*(dLR - this.GetLogLik());				// compute LR statistic for the estimated model against the unrestricted model
	df -= this.GetFreeParCount();					// compute # of restrictions (degrees of freedom for asymptotic pvalue)
	if(dLR < 0 || df <= 0){							// stop and print warning in case models appear nonnested
		println("TestAgainstUnrestr(): Unrestricted does not nest restricted model, or estimation failed to attain global maximum likelihood.",
		"\nMax Unrestr Log Lik: ", double(m_oUnrestricted.GetLogLik()), ", Max Restr Log Lik: ", double(this.GetLogLik()), ", df = ", df);
		return FALSE;
	}
		
	decl bpval = .NaN;							// bootstrap p-value (if desired)
	
	if(m_fPrint)										 
		PrintTestVal(dLR, df, 0, "LR test of restrictions");
	if(m_cBootRep){								// do bootstrap test if m_cBootRep>0 (this is set using SetBootstrap)
		if(!m_fBootstrap || !m_fTestUnrestr){	// if you have not already computed the bootstrap sample under null and alternative models,
			m_fTestUnrestr = TRUE;				// (re-)do the bootstrap, making sure you compute it for both null and alternative model
			Bootstrap();
		}
		m_mBootData[0][find(m_asBootNames,"LR")] = dLR;			// replace missing value with just computed LR statistic from actual data (remaining rows contain bootstrap replications calculated within Bootstrap())
		bpval = double(1-meanc(deleter(m_mBootData[][find(m_asBootNames,"LR")] .< dLR)));		// tail probability of observed dLR in boostrap distribution of LR
		if(m_fPrint)
			println("Bootstrap p-value using ", min(m_cBootRep,sizeof(deleter(m_mBootData[][find(m_asBootNames,"LR")]))-1), " bootstrap replications: ", bpval);
		if(m_fStoreBoot)						// save the file again to include the LR statistic on actual data just computed
			savemat(m_sStoreBootFile, m_mBootData, m_asBootNames);

	}
	return {dLR, df, tailchi(dLR,df), bpval};
}
CKSVAR::GetBetabar(const vbetatil, const xi, const mOmega)
/**	Finds betabar by solving the equation vbetatil = (1-xi)*(I-xi*betabar*gammabar)^(-1)*betabar
	@param vbetatil [m_cY1][1] vector of coefficients of kink in reduced form
	@param xi, double
	@param mOmega, [m_cY][m_cY], reduced form error variance
	@return	mbetabar [m_cY1][L], L<=m_cY solutions that satisfy coherency condition, empty if L=0
**/
{
	decl k = sizeof(mOmega)-1;
	if(k>2){
		println("CKSVAR::GetBetabar: cases with k > 3 not supported yet.");
		return <>;
	}
	decl om11 = mOmega[:k-1][:k-1], om12 = mOmega[:k-1][k], om22 = mOmega[k][k];
	if(xi == 0 || isfeq(vbetatil,0))
//		return GetGamma(vbetatil,mOmega)*vbetatil < 1 ? vbetatil : .NaN;
		return vbetatil;				// disable coherency check

	decl flold = oxwarning(0); 			// disable warnings
	decl om11inv = invertsym(om11);
    oxwarning(flold);      		// reset the previous warning flags.
	if(om11inv == 0)			// if inverse fails
		return <>;				// return nothing
	decl ident = unit(sizeof(om11));
	decl btil = -om11inv*( ((om12'om11inv*om12-om22)*ident - om12*om12'om11inv) * vbetatil*xi - (1-xi)*om12);
	decl Atil = vbetatil*om12'om11inv+(xi*om12'om11inv*vbetatil+1-xi)*ident;
	decl vA = new matrix[4][1], vz, L;			// coefficients of cubic, roots, z = btil'vbetabar
	decl vw;									// holds solutions for bperp'vbetabar
	decl vinfinityrange = range(-10,10);		// an arbitrary sequence of number to "span" the real line
	decl mbetabar, mgammabeta;									// spans space orthogonal to btil
	if(isfeq(btil,0)){							// if zero (can happen), solution is unique
		flold = oxwarning(0); 			// disable warnings
		decl arg = invert(Atil);
		oxwarning(flold);     	 		// reset the previous warning flags.
		if(arg == 0)					// if inverse fails
			return <>;					// return nothing (no solution exists, except in the nongeneric case betatil==0, in which there are infinite solutions)
		mbetabar = arg*vbetatil;
		flold = oxwarning(0); 			// disable warnings
		arg = invert(om11-om12*mbetabar');
		oxwarning(flold);     	 		// reset the previous warning flags.
		if(arg == 0)					// if inverse fails
			return <>;					// return nothing (no solution exists)
		mgammabeta = (om12-om22*mbetabar)'arg*mbetabar;
	}
	else if(k==1){								// if m_cY1 = 1, simple quadratic
		vA = vbetatil|-Atil|btil|0;
		cubic(vA, &vz, &L);					  	// this gives the solutions for z = btil'vbetabar, in this case quadratic
		if(L==0)					// no real solutions in quadratic
			return <>;				// implies no overall solution
		mbetabar = vz;
		decl dfold = fuzziness(0);		// keep track of current fuzziness, to temporarily reset it
		fuzziness(1e-9);
		mbetabar = selectifc(mbetabar, !isdotfeq(om12'om11inv*mbetabar,1));				// drop solution with singularity of om12'om11inv*mbetabar-1
		mgammabeta = new matrix[1][columns(mbetabar)];
		fuzziness(dfold);			    // restore previous fuzziness
		mgammabeta = new matrix[1][columns(mbetabar)];
		for(decl i = 0; i < sizec(mgammabeta); ++i)
			mgammabeta[i] = (om12-om22*mbetabar[][i])'invert(om11-om12*mbetabar[][i]')*mbetabar[][i];
	}
	else{
		decl bperp = btil[1]|-btil[0];			// easy because it's 2x2, in general use SVD
		bperp /= sqrt(bperp'bperp);					// standardize to simplify the rest
		decl n = btil'btil;
		decl a11 = btil'Atil*btil/n;
		decl a12 = btil'Atil*bperp;
		decl a21 = bperp'Atil*btil/n;
		decl a22 = bperp'Atil*bperp;
		decl c1 = btil'vbetatil;
		decl c2 = bperp'vbetatil;
		if(isfeq(a12,0)){		// solve quadratic
			vA[0] = c1;
			vA[1] = -a11;
			vA[2] = 1;
			vA[3] = 0;
			cubic(vA, &vz, &L);		  	// this gives the solutions for z = btil'vbetabar, in this case quadratic
			if(L==0)					// no real solutions in quadratic
				return <>;				// implies no overall solution
			else if(L==1){				// one double real root
				if(!isfeq(a22,vz))
					vw = (c2-a21*vz)/(a22-vz);	// unique solution
				else if(isfeq(c2,a21/n*vz))
					vw = vinfinityrange;			// w is completely indeterminate, infinite solutions, return a banch of values
				else				
					return <>;			// no solution
				mbetabar = btil/(btil'btil)*vz + bperp*vw;		// obtain generic solutions for betabar given vz and vw	(1 or "infinite")
			}
			else{
				if(!any(isdotfeq(vz-a22,0))){							// if none of the solutions of z are equal to a22,
					vw = (c2-a21*vz) ./ (a22-vz);		// then exactly two solutions
					mbetabar = btil/(btil'btil)*vz + bperp*vw;		// obtain two generic solutions for betabar given vz and vw
				}
				else{
					mbetabar = <>;				// initialize dimensions
					for(decl iz = 0; iz < L; ++iz){
						if(!isfeq(a22,vz[iz])){
							vw =  (c2-a21*vz[iz]) / (a22-vz[iz]);
							mbetabar ~= btil/(btil'btil)*vz[iz] + bperp*vw;
						}
						else if(isfeq(c2,a21/n*vz[iz])){
							vw = vinfinityrange;			// w is completely indeterminate, infinite solutions, return a banch of values
							mbetabar ~= btil/(btil'btil)*vz[iz] + bperp*vw;
						}
						//else, no solution, there has to be at least one, so mbetabar won't be empty
					}
				}
			}
		}					
		else{
			vA[0] = (a22*c1-a12*c2);
			vA[1] = (a12*a21-a11*a22)-c1;
			vA[2] = a11+a22;
			vA[3] = -1;
			
			cubic(vA, &vz, &L);		  // this gives the solutions for z = btil'vbetabar
			vw = (c1-a11*vz+sqr(vz))/a12;
			mbetabar = btil/(btil'btil)*vz + bperp*vw;		// obtain 1 up to 3 generic solutions for betabar given vz and vw
		}
		decl dfold = fuzziness(0);		// keep track of current fuzziness, to temporarily reset it
		fuzziness(1e-9);
		mbetabar = selectifc(mbetabar, !isdotfeq(om12'om11inv*mbetabar,1));				// drop solution with singularity of om12'om11inv*mbetabar-1
		mgammabeta = new matrix[1][columns(mbetabar)];
		fuzziness(dfold);
		for(decl i = 0; i < sizec(mgammabeta); ++i)
			mgammabeta[i] = (om12-om22*mbetabar[][i])'invert(om11-om12*mbetabar[][i]')*mbetabar[][i];
	}
	decl vcoh = ((1-mgammabeta).*(1-xi*mgammabeta)) .> 0;

	return selectifc(mbetabar, vcoh);			// return only those solutions that satisfy coherency
}

CKSVAR::GetGamma(const vbetabar, const mOmega)
/**	Returns vgammabar [1][m_cY1] coefficients on Y1 in Y2* equation
	@param vbetabar [m_cY1][1] vector of sum of contemporaneous coefficients of Y2* and Y2 on Y1.
	@param mOmega [m_cY][m_cY] reduced form variance matrix
**/
{
	decl k = sizeof(mOmega)-1;
	return (mOmega[:k-1][k]'-mOmega[k][k]*vbetabar')*invert(mOmega[:k-1][:k-1]-mOmega[:k-1][k]*vbetabar');
}
CKSVAR::SetIRFsigns(const mSigns)
/** Sets the sign of IRFs for use in the sign restrictions
	@param mSigns: matrix [m_cY][1+m_cHorizon] of 1s (positive), -1s (negative) and 0s (unrestricted) signs corresponding to each IR_i,h
	Note: will be somewhat inefficient if there are many zeros, since we will be computing IRs that don't contribute to the sign restrictions, so try to avoid as many zeros as posisble in mSigns, e.g., by keeping cHorizon to a minimum
**/
{
	decl arg = mSigns .> 0 .? 1 .: mSigns .< 0 .? -1 .: 0;		// make sure it's always 1s, -1s and 0s
	decl r = rows(arg), c = columns(arg);
	if(r != m_cY1+1 ||  c != m_cHorizon+1){
		println("SetIRFsigns: input dimension incorrect.");
		return;
	}
	format(1000);
	println("Sign restrictions (horizon in columns starting at 0): 1 means positive, -1 negative, 0 unrestricted",
			"%r", m_asY1~m_asY2, mSigns);
	println("Total # of sign restrictions: ", int(sumc(sumr(fabs(mSigns)))));
	m_mSigns = mSigns;			// matrix of 1 (positive), -1 (negative) and 0 (unrestricted)
}
CKSVAR::SetRannIRForigin(const cT)
/** Sets the random numbers to use when computing origin of IRFs for sign restrictions
	@param cT: integer>0: sample size to draw
**/
{
	if(!isint(cT) || !(cT>0)){
		eprint("CKSVAR::SetIRFrannOrigin: argument must be positive integer.");
		return;
	}
	m_mrannIRForigin = rann(cT, m_cY1+1);			// matrix of random numbers to use in SimData()
}
CKSVAR::funcIRF(const avIRF, const vP)
/** Computes IRF to structural shock to Y2*, in form needed for NumJacobian. It expects Y2 in levels.
	@param avR, address of vector of desired IRFs
	@param vP, vector of unrestricted (free) parameters
	returns TRUE if function evaluation successful
**/
{
	avIRF[0] = <>;										// initialize at null, as default value when function evaluation fails
	decl vpall;
	if(!Reparameterize(&vpall, vP))						// obtain full parameter vector
		return 0;
	decl tau, vC2obs, vC2lat, vbetatil, mC1obs, mC1lat, vdelta, mOm1g2Ch;			// obtain the parameters first
	[tau, vC2obs, vC2lat, vbetatil, mC1obs, mC1lat, vdelta, mOm1g2Ch] = this.GetModelPar(vpall);
	decl vgamma, vbetabar, mOm11bar_ch, mgammabar, mOmega;										
	if(m_cY1){
		decl mOm1g2 = mOm1g2Ch*mOm1g2Ch';					// variance of u_1|u_2
		decl mOm11 = mOm1g2 + tau^2*vdelta*vdelta';	// unconditional variance matrix for Y1
		mOmega = mOm11 ~ vdelta*tau^2
				| vdelta'tau^2 ~ tau^2;			// variance matrix of reduced form errors for all Y
		if(m_fCholeski){
			decl vidx = sortcindex(vec(m_vorder));
			decl idxY2 = vidx[m_cY1];					// locate the position of Y2_VAR in m_vorder -- will determine which slice of the Choleski corresponds to beta and gamma
			decl vidxY1 = vidx[:m_cY1-1];				// index of Y1_VAR
			mOmega = mOmega[m_vorder][m_vorder];		// re-order according to the desired Choleski ordering
			decl chol = choleski(mOmega);
			decl stdevs = diagonalize(chol); 			// st deviations of structural errors in Choleski factorization
			decl mat = stdevs*invert(chol);				// matrix of recursive regression coefficients
			vbetabar = -mat[vidxY1][idxY2];
			vgamma = -mat[idxY2][vidxY1];
		}
		else{
			if(m_fUseBetabar)
				m_mbetabar = vpall[m_vIdxbetatilde];
			else
				m_mbetabar = this.GetBetabar(vbetatil,m_dlambda*m_dzeta,mOmega);
			m_iSol = columns(matrix(m_mbetabar));			// record the number of solutions
			if(!m_iSol)									// if there's no solution, return nothing
				return 0;	
		}
	}
	else
		vgamma = vbetabar = 0;							// set to zero if empty, to avoid conditional statements where these appear below

	decl crep = max(m_crepIRF,1);						// # of MC replications to compute IRF. If 0, set to 1 to execute for loop below
	decl k = m_cY1+1;									// dimension of VAR
	decl p = max(this.GetMaxGroupLag(Y1_VAR),
			this.GetMaxGroupLag(Y2_VAR)+m_fdiffY2);		// VAR order, assumes all Y variables have same number of lags
	decl T = m_cHorizon+1;								// # of time periods to simulate data
	decl cX0 = columns(this.GetGroup(X_VAR));			// # of exogenous and deterministic regressors

	if(m_fdiffY2){										// if TRUE recover the coefficients on constrained variable in levels
		decl cX1 = m_cX-m_cY2lags;						// # of regressors excluding lagged Y2 on RHS
		decl mC = mC1obs ~ vC2obs;						// collect all coefficints in [k][m_cY] matrix
		decl mCY1 = cX1 ? mC[:cX1-1][] : <>;			// coefficients up to Y1 lags (if there are any)
		decl mCY2 = p>m_fdiffY2 ? mC[cX1:][] : (zeros(1,k-1)~1);			// coefficients on Y2 lags
		if(p==1){
			mCY2 = zeros(1,k-1)~1;
			mC1lat = zeros(1,k-1);						// the coefficients on the latent lags in Y1 eqns must be zero
			//vC2lat includes only coef on lagged level in this case, so does not need to change
		}
		else if (k==1){
			mCY2 = (1+mCY2[0])|reversec(-diff0(0|reversec(mCY2))[1:]);
			vC2lat = m_fLatLags ? (vC2lat[p-1] + vC2lat[0])|reversec(-diff0(0|reversec(vC2lat[:p-2]))[1:]) : 0;
		}
		else{
			mCY2 = (mCY2[0][:k-2] ~ (1+mCY2[0][k-1]) )
				|reversec(-diff0(zeros(1,k)|reversec(mCY2))[1:][]);
			vC2lat = m_fLatLags ? (vC2lat[p-1] + vC2lat[0])|reversec(-diff0(0|reversec(vC2lat[:p-2]))[1:]) : 0;
			mC1lat = (m_cY1 && m_fLatLags) ? (mC1lat[0][] | reversec(-diff0(zeros(1,k-1)|reversec(mC1lat))[1:][])) : 0;
		}
		decl mClev = mCY1 | mCY2;					// collect coefficients again, but where Y2 is in levels
		if(m_cY1)
			mC1obs = mClev[][:k-2];					// update coefficients on lags in Y1 equations
		vC2obs = mClev[][k-1];
	}
	if(!m_fLatLags)
		vC2lat = mC1lat = 0;

	decl b = m_vb[(m_cOrigin < 0 || m_cOrigin > m_cT-1) ? m_cT-1 : m_cOrigin];

	if(m_fLatLags && !m_fSetOriginIRF){					// if necessary, compute initial values of latent lags by IS smoother
		decl p2 = sizerc(vC2lat);
		decl idx = sizerc(m_vObsOriginIRF)-p2;			// index of beginning of lagged Y2s in m_vObsOrigin
		decl vX2 = m_vObsOriginIRF[idx:];				// lagged Y2s
		decl vLatLags = new matrix[1][p2];
		if(any(vX2 .<= b)){								// if any of the lagged Y2 is on the ZLB, compute smoothed estimate of latent Y2*
			decl vlik, mY2draws, mw;
			[vlik, mY2draws, mw] = this.SimLik(vpall, TRUE);
			decl weights = m_fFAPF ? 1 : mw[m_cT-1][];	// if FAPF, using particles from filtering density, otherwise if SIS, use the smoothing weights
			if(!m_fdiffY2)
				vLatLags = meanr(m_mXstar.*weights)';
			else{
				vLatLags[0][0] = meanr(mY2draws[m_cOrigin][].*weights);
				for(decl ip = 1; ip < p2; ++ip)
					vLatLags[0][ip] = vLatLags[0][ip-1] - meanr(m_mXstar[ip-1][].*weights);				// compute Y2*[t-ip] = Y2*[t-ip+1] - \Delta Y2*[t-ip+1]
			}
		}
		m_vLatOriginIRF = vX2 .> b .? vX2 .: vLatLags;		  // use actual values for Y2*>b and smooth estimates otherwise
	}

	decl mIRFs = new matrix[m_iSol*k][T];
	
	decl fdraw = m_crepIRF>0;								// if cMCrep>0, you will need to draw from the distribution of the other shocks, otherwise, set those shocks to zero
	for(decl j = 0; j < m_iSol; ++j){
		vbetabar = m_mbetabar[][j];
		vgamma = GetGamma(vbetabar,mOmega);
		decl mOm11bar = (unit(m_cY1)~(-vbetabar))*mOmega*(unit(m_cY1)|-vbetabar');
		mOm11bar_ch = choleski(mOm11bar);
		decl mInv = invert(unit(k-1)-vbetabar*vgamma);		// to be used inside parallel loop
	
		decl mX = new matrix[crep][sizerc(m_vObsOriginIRF)];				// initialize observed lags when shock = dImpulse (force row vector)
		decl mXst = m_fLatLags ? new matrix[crep][sizerc(m_vLatOriginIRF)] : 0;	// initialize latent lags (lags of Y2*) Xt*, only necessary if m_cLatLags == TRUE
		decl mX0 = mX, mXst0 = mXst;					// corresponding regressors when shock=0
		decl amY = new array[p], amY0=amY;				// holds each of the lags Y
		for(decl i = 0; i < p; ++i)
			amY[i] = amY0[i] = new matrix[crep][k];		// all draws of Y for lag t-i-1
		decl my2star = new matrix[crep][T];				// all draws of Y2*
		decl my2star0 = my2star;						// under no impulse
		decl mYt = new matrix[crep][k], mY0t = mYt;		// draws of Y at t
		
		decl id1 = crep*(k-1);
		decl meps1bar = m_cY1 && fdraw ? shape(m_mrannIRF[:id1-1], crep, k-1)*mOm11bar_ch' : zeros(1,k-1); 		//new matrix[crep][k-1];
		decl id2 = id1 + crep*(T-1);
		decl mu2 = T>1 ? (fdraw ? tau*shape(m_mrannIRF[id1:id2-1], crep, T-1) : zeros(1,T-1)) : 0; 			//new matrix[crep][T-1];
				
		decl u2 = (m_dImpulse + meps1bar*vgamma')/(1-vgamma*vbetabar);	// reduced form error in Y2* equation in period 0 when shock = dImpulse
		decl u20 = (meps1bar*vgamma')/(1-vgamma*vbetabar);				// reduced form error in Y2* equation in period 0 when shock = 0
		for(decl t = 0; t < T; ++t){
			decl mu1g2t;
			if(t){
				id1 = id2;
				id2 += crep*(k-1);
				mu1g2t = m_cY1 && fdraw ? shape(m_mrannIRF[id1:id2-1], crep, k-1)*mOm1g2Ch' : zeros(1,k-1);
			}
			for(decl ik = 0; ik < k; ++ik){
				for(decl j = 0; j < p; ++j){
					mX[][ik*p+j] = t>j ? amY[j][][ik] : ones(crep,1)*vec(m_vObsOriginIRF[ik*p+j-t])';
					mX0[][ik*p+j] = t>j ? amY0[j][][ik] : ones(crep,1)*vec(m_vObsOriginIRF[ik*p+j-t])';
				}
			}
			if(m_fLatLags){
				for(decl j = 0; j < p; ++j){					// not executed if m_cY2lags = 0
					mXst[][j] = t>j ? my2star[][t-j-1] : ones(crep,1)*vec(m_vLatOriginIRF[j-t])';		
					mXst0[][j] = t>j ? my2star0[][t-j-1] : ones(crep,1)*vec(m_vLatOriginIRF[j-t])';		
				}
			}
			decl arg  = (cX0 ? m_mExogIRF[t][]*vC2obs[:cX0-1] : 0) + (sizeof(vC2obs)>cX0 ?  mX*vC2obs[cX0:] : 0) + mXst*vC2lat	+ (t ? mu2[][t-1] : u2);
			decl arg0 = (cX0 ? m_mExogIRF[t][]*vC2obs[:cX0-1] : 0) + (sizeof(vC2obs)>cX0 ? mX0*vC2obs[cX0:] : 0) + mXst0*vC2lat+ (t ? mu2[][t-1] : u20);
			my2star[][t] = (arg .> b) .? arg .: m_dkappa*arg + (1-m_dkappa)*b; 			// keep in mind the kink in the latent equation
			my2star0[][t] = (arg0 .> b) .? arg0 .: m_dkappa*arg0 + (1-m_dkappa)*b; 		// keep in mind the kink in the latent equation
			mYt[][k-1] = my2star[][t] .> b .? my2star[][t] .: b;
			mY0t[][k-1] = my2star0[][t] .> b .? my2star0[][t] .: b;
			if(k>1){
				mYt[][:k-2] = (cX0 ? m_mExogIRF[t][]*mC1obs[:cX0-1][] : 0) + (sizeof(mC1obs)>cX0 ?  mX*mC1obs[cX0:][] : 0) + mXst*mC1lat
						- (arg .< b) .* (arg - b) .* vbetatil'
						+ (t ? (mu1g2t+mu2[][t-1]*vdelta') : (meps1bar+m_dImpulse*vbetabar')*mInv');		// simulated Y1 when shock = dImpulse
				mY0t[][:k-2] = (cX0 ? m_mExogIRF[t][]*mC1obs[:cX0-1][] : 0) + (sizeof(mC1obs)>cX0 ? mX0*mC1obs[cX0:][] : 0) + mXst0*mC1lat
						- (arg0 .< b) .* (arg0 - b) * vbetatil'
						+ (t ? (mu1g2t+mu2[][t-1]*vdelta') : meps1bar*mInv');						// simulated Y1 when shock = 0
			}
			mIRFs[j*k:((j+1)*k-1)][t] = meanc(mYt-mY0t)';
			for(decl j = p-1; j >= 0; --j){
				amY[j] = j ? amY[j-1] : mYt;
				amY0[j] = j ? amY0[j-1] : mY0t;
			}
		}
	}
/*		//analytical formula for impact effect
    	decl coeffu2 = vgamma*mInv*vbetabar*(1-vgamma*vbetabar); // the coefficient in front of u2t in the expression for omegabarsq
    	decl omegabarsq =	1/((1-vgamma*vbetabar)^2)*(vgamma*mOm11*vgamma' + coeffu2*(tau^2)*(coeffu2)' - 2*vgamma*vdelta*tau^2*coeffu2');
    	decl omegabar = sqrt(omegabarsq);
    	decl normarg1 = (b-arg-dImpulse/(1-vgamma*vbetabar))/omegabar;
    	decl normarg2 = (b-arg)/omegabar;
    	decl diff_prob = probn(normarg1)-probn(normarg2);
    	decl diff_dens = densn(normarg1)-densn(normarg2);
    	//calculate the g0 functions			
    	mg1IRF[][j] = mInv*vbetabar*dImpulse - diff_prob*vbetatil*(arg - b) - probn(normarg1)*vbetatil*dImpulse/(1 - vgamma*vbetabar) + omegabar*diff_dens*vbetatil;
        mg2IRF[j] = dImpulse/(1 - vgamma*vbetabar)*(1 - probn(normarg1)) - diff_prob*(arg-b) + omegabar*diff_dens;
*/
	avIRF[0] = vec(mIRFs);
	return !isnan(avIRF[0]);

}
CKSVAR::SetCholeskiIRF(const fCholeski, const vorder)
/** Sets flag for Choleski and order 
	@param fCholeski, TRUE or FALSE
	@param vorder, vector [1][m_cY] with the order of the selected variables in the choleski factorization
**/
{
	m_fCholeski = fCholeski;			// if TRUE, report IRF using Choleski factorization
	if(fCholeski == TRUE)
		m_vorder = vorder;				// store Wold causal ordering
}
CKSVAR::GetBootIRFfromStored(const mBootPar, const dcoverage)
/** Computes Bootstrap confidence bands for IRF at desired coverage using stored data
	@param mBootPar, matrix [m_cBootRep+1][m_cFreePar]
	@param dcoverage, double, in (0,1).
	returns {mIRF, mIRFlo, mIRFup} array of [m_cY][m_cHorizon+1] matrices of IRF and dcoverage error bands
**/
{

	decl amIRF = new array[m_iSol*(m_cY1+1)];			// array to hold results per variable, repeating for each solution if more than one
	for(decl i = 0; i < sizeof(amIRF); ++i)
		amIRF[i] = new matrix[sizeof(mBootPar)][m_cHorizon+1];	// bootstrap reps in rows for IRF of var i at horizon in columns. 
	decl arg;
	for(decl j = 0; j < sizeof(mBootPar); ++j){
		decl vParBoot = mBootPar[j][]';					// extract parameter estimates from jth bootstrap replication
		funcIRF(&arg, vParBoot);
		arg = shape(arg,m_iSol*(m_cY1+1), m_cHorizon+1);
		for(decl i = 0; i < sizeof(amIRF); ++i)
			amIRF[i][j][] = arg[i][];
	}		
	decl mIRF, mIRFlo, mIRFup;
	mIRF=mIRFlo=mIRFup=new matrix[sizeof(amIRF)][m_cHorizon+1];
	for(decl i = 0; i < sizeof(amIRF); ++i){
		mIRF[i][] = amIRF[i][0][];									// first entry contains estimated IRF for var i
		decl vquants = quantilec(amIRF[i]-mIRF[i][],(1-dcoverage)/2~(1-(1-dcoverage)/2));		// alpha/2 and 1-alpha/2 quantiles of demeaned bootstraps of IRFs
		mIRFup[i][] = mIRF[i][]-vquants[0][];
		mIRFlo[i][] = mIRF[i][]-vquants[1][];
	}
	return {mIRF, mIRFlo, mIRFup};
}
CKSVAR::GetIRFset(const vrange, ...)
/** Computes identified set of IRFs over the required range of lambda*zeta
	@param 	vector, range of values to compute IDset
	optional: vector of values of lambda for each IRF
	optional: matrix of values of betabar for each IRF
	returns array[m_cY] of [cSet][cHorizon+1] matrices of IRFs, where cSet is the number of IRFs computed over all elements of vrange
**/
{
	decl k = m_cY1+1;
	decl amIRF = new array[k];
	decl vlambda = <>;								// to hold the lambda corresponding to each IRF, implicitly giving the admissible range of lambda
	decl ambetabar = {};								// to hold the betabar corresponding to each lambda
	decl avlam = <>, aambetabar = <>;					// address to return vlambda and mbetabar if required -- initialize to empty
	if(sizeof(va_arglist()))
		avlam = va_arglist()[0];
	if(sizeof(va_arglist())==2)
		aambetabar = va_arglist()[1];
	for(decl j = 0; j < sizeof(amIRF); ++j)
		amIRF[j] = <>;								// initialize at empty matrices, since you don't know dimension beforehand

	decl lambda = m_dlambda, zeta = m_dzeta;		// store original values to restore at the end
	this.Setzeta(1);								// set zeta=1 (no change in variance during ZLB regime)
	for(decl i = 0; i < sizerc(vrange); ++i){
		ranseed(111113);
		this.Setlambda(vrange[i]);
		decl mIRFnl = this.GetIRF(/*coverage*/0,/*Asy*/FALSE,/*Boot*/FALSE);
		if(mIRFnl==<>)
			continue;
		for(decl j = 0; j < sizeof(amIRF); ++j)
			amIRF[j] |= mIRFnl[j+range(0,m_iSol-1)*k][];						// collect IRFs for variable j
		if(sizeof(avlam))
			vlambda |= vrange[i]*ones(m_iSol,1);
		if(sizeof(aambetabar))						
			ambetabar ~= {m_mbetabar};				// optionally store betabar too
		
	}
	
	this.Setzeta(zeta);								// restore to original parameters
	this.Setlambda(lambda);							// restore to original parameters
	if(sizeof(avlam))
		avlam[0] = vlambda;							// return the admissible lambdas
	if(sizeof(aambetabar))						
		aambetabar[0] = ambetabar;					// return the admissible betabars associated with each lambda
	return amIRF;
}
CKSVAR::SetIRF(const cHorizon, const dImpulse, const cMCreps, const origin)
/** Stores settings of IRF computation to class data members
	@param cHorizon, count # of periods to compute IRF for
	@param dImpulse, double, magnitude of impulse on Y2* shock from zero (epsilon2)
	@param cMCreps, count, to take expectations over future shocks. If set to zero, all other shocks are set to zero
	@param origin, either: integer, point in sample to use as origin for IRF calculation. If t<0 or t>m_cT use latest observation in the sample
					or: array[3] with vector of initial values of Observed, Latent and Exogenous regressors
	returns TRUE if successful
**/
{
////////////////////////////////////////////////////////
//	set member variables, to be used inside funcIRF(...)
	m_cHorizon = cHorizon;
	m_dImpulse = dImpulse;
	m_fSetOriginIRF = FALSE;
	if(isarray(origin) && sizeof(origin)==3){
		m_fSetOriginIRF = TRUE;
		this.SetIRForigin(origin);
//		m_vObsOriginIRF = vec(origin[0])';			// force row vector
//		decl idx = sizerc(origin[0])-m_cY2lags;		// index of beginning of lagged Y2s in origin[0]
//		m_vLatOriginIRF = m_vObsOriginIRF[][idx:] .> m_db .? m_vObsOriginIRF[][idx:]
//						.: vec(origin[1])';		// force Y2*=Y2 when Y2>b
//		m_mExogIRF = origin[2];
	}
	else if(isint(origin)){
		decl t = origin;
		m_cOrigin = t>m_cT-1 || t<0 ? m_cT-1 : t;			// origin of IRF computation (indexing starts at 0)
		m_vObsOriginIRF = this.GetGroupLag(Y1_VAR,1,this.GetMaxGroupLag(Y1_VAR))[m_cOrigin][]		// lags of Y1, starting at t, t-1, t-p+
						~ (m_fdiffY2 ? this.GetGroupLag(Y2_LAGLEV,0,this.GetMaxGroupLag(Y2_LAGLEV))	  
						: this.GetGroupLag(Y2_VAR,1,this.GetMaxGroupLag(Y2_VAR)))[m_cOrigin][];	// lags of Y2
		m_mExogIRF = ones(cHorizon+1,1)*this.GetGroup(X_VAR)[m_cOrigin][];							// exogenous regressors treated as fixed throughout the horizon (not good for trends and seasonals)
		m_vLatOriginIRF = <>;						// leave empty, to be computed by smoothing inside funcIRF, at given parameter values
	}
	else{
		println("IRF origin hasn't been set correctly.");
		return;
	}
	m_crepIRF = cMCreps;							// # of MC reps for IRF computation
	if(cMCreps)
		m_mrannIRF = rann(cMCreps*(m_cHorizon+1)*m_cY1+cMCreps*m_cHorizon,1);		// all random draws for IRF	

	return !any(isnan(m_cHorizon~m_dImpulse~m_fSetOriginIRF~m_crepIRF)~isnan(m_vObsOriginIRF)~isnan(m_mExogIRF)~isnan(m_mrannIRF));
}
CKSVAR::SetIRForigin(const aOrigin)
/** Sets the origin for computation of IRFs, i.e., m_vObsOriginIRF, m_vLatOriginIRF, m_mExogIRF
	@param aOrigin, array[3] with vector of initial values of Observed, Latent and Exogenous regressors
**/
{
	m_vObsOriginIRF = vec(aOrigin[0])';				// force row vector
	decl idx = sizerc(aOrigin[0])-m_cY2lags;		// index of beginning of lagged Y2s in origin[0]
	m_vLatOriginIRF = m_vObsOriginIRF[][idx:] .> meanc(m_vb) .? m_vObsOriginIRF[][idx:]			// evaluate at average value of lower bound if time-varying, for lack of a better alternative
						.: vec(aOrigin[1])';		// force Y2*=Y2 when Y2>b
	m_mExogIRF = aOrigin[2];
}
CKSVAR::GetIRForigin()
/** @return array {m_vObsOriginIRF, m_vLatOriginIRF, m_mExogIRF}
**/
{
	return {m_vObsOriginIRF, m_vLatOriginIRF, m_mExogIRF};
}

CKSVAR::GetIRF(const dcoverage, const fAsy, const fBoot)
/** Computes Impulse Response function to structural shock to Y2*
	@param dcoverage, double, in (0,1).
	@param fAsy, if TRUE, return asymptotic error bands
	@param fBoot, if TRUE, return bootstrap error bands
	returns if fAsy=fBoot=FALSE, or if m_iSol > 1, matrix [m_cY*m_iSol][cHorizon+1] of IRFs (m_iSol is number of solutions of betabar from given RF parameters)
			if any fAsy,fBoot true AND m_iSol =1, array[1+2*fAsy+2*fBoot] of [m_cY][cHorizon+1] matrices of IRFs, and lower/upper error bands 
**/
{
	decl mIRFs;
	funcIRF(&mIRFs, this.GetFreePar());
	mIRFs = shape(mIRFs,m_iSol*(m_cY1+1), m_cHorizon+1);
	if(m_iSol > 1 && m_fPrint)
		println("GetIRF: ", m_iSol, " solutions for betabar at this parameterization.");
	if(m_iSol > 1 || (!fAsy && !fBoot) || isfeq(dcoverage,0))	// if more than one solutions, you cannot guarantee the ordering of solutions is continuous wrt the parameters, or across boostraps, so numerical differentiation or boostrap quantiles may be off 
		return mIRFs;											// or if you don't want error bands, or if error band coverage is zero, just return the point estimates
	//else
	decl aresult = new array[1+2*fAsy+2*fBoot];
	aresult[0] = mIRFs;
	if(fAsy){
		decl vJac;
		decl ires = NumJacobian(funcIRF, this.GetFreePar(), &vJac);
		decl mIRFvar = vJac*this.GetCovar()*vJac';			// variance matrix of IRFs by the delta method
		decl msterr = ires ? shape(sqrt(diagonal(mIRFvar)), m_cY1+1, m_cHorizon+1) : constant(.NaN, m_cY1+1, m_cHorizon+1);	// to hold matrix of asymptotic standard errors
		decl dcrit = tailn(1-(1-min(dcoverage,1))/2);		// compute critical value for asymptotically normal error bands
		aresult[1] = mIRFs - dcrit*msterr;					// lower bound
		aresult[2] = mIRFs + dcrit*msterr;					// upper bound
	}
	if(fBoot && m_fBootstrap){
		decl mBootPar = m_mBootData[][m_fTestUnrestr+range(1,this.GetFreeParCount())];		// retrieve bootstrapped parameter estimates
		aresult[2*fAsy+<1:2>] = GetBootIRFfromStored(mBootPar, dcoverage)[1:2];
	}
	return aresult;
}
CKSVAR::funcLinearIRF(const avR, const vP)
/** Computes pseudo linear IRF to structural shock to Y2*, in form needed for NumJacobian
	@param avR, address of vector of desired IRFs
	@param vP, vector of unrestricted (free) parameters
	returns TRUE if function evaluation successful
**/
{
	decl vpall;
	if(!Reparameterize(&vpall, vP))									 // obtain full parameter vector
		return 0;
	decl tau, vC2obs, vC2lat, vbetatil, mC1obs, mC1lat, vdelta, mOm1g2Ch;			// obtain the parameters first
	[tau, vC2obs, vC2lat, vbetatil, mC1obs, mC1lat, vdelta, mOm1g2Ch] = this.GetModelPar(vpall);
	decl mOm11 = mOm1g2Ch*mOm1g2Ch' + tau^2*vdelta*vdelta';	// unconditional variance matrix for Y1
	decl mOmega = mOm11 ~ vdelta*tau^2
				| vdelta'tau^2 ~ tau^2;			// variance matrix of reduced form errors for all Y
	decl vgamma = this.GetGamma(vbetatil,mOmega);
	decl vR0 = (vbetatil | 1)/(m_cY1 ? (1-vgamma*vbetatil) : 1);
	decl mComp = this.GetCompanionMatrix(vpall);
	decl k = sizer(vR0);			// dimension of VAR
	decl p = sizer(mComp)/k;		// order of VAR
	decl mIRF = new matrix[k][m_cHorizon+1];
	decl vR0tilde = vR0 | zeros(sizer(mComp)-k,1);
	decl mSel = unit(k)~zeros(k,k*(p-1));		// select the first k elements of companion vector
	for(decl i = 0; i <= m_cHorizon; ++i)
		mIRF[][i]= mSel*mComp^i*vR0tilde; 		 // compute i-th irf
	avR[0] = m_dImpulse*vec(mIRF);
	return !isnan(mIRF);

}
CKSVAR::GetLinearIRF(const cHorizon, const amStErr)
/** Computes pseudo linear IRF to structural shock to Y2* (above the ZLB)
	@param cHorizon, count # of periods to compute IRF for
	@param amStErr, address of matrix of standard errors
	returns mIRF [m_cY][cHorizon+1] matrix of IRFs
**/
{
	m_cHorizon = cHorizon;
	decl vIRF;
	funcLinearIRF(&vIRF, this.GetFreePar());
	if(amStErr){							// if nonzero, compute standard errors for each IRF
		decl vJac;
		decl ires = NumJacobian(funcLinearIRF, this.GetFreePar(), &vJac);
		decl mIRFvar = vJac*this.GetCovar()*vJac';			// variance matrix of IRFs by the delta method
		amStErr[0] = ires ? shape(sqrt(diagonal(mIRFvar)), m_cY1+1, cHorizon+1) : constant(.NaN, m_cY1+1, cHorizon+1);
	}
	return shape(vIRF, m_cY1+1, cHorizon+1);
}
CKSVAR::GetCompanion(mC)
/** computes companion matrix by re-arranging coefficients on VAR lags
	@param mC coefficient on lags
	@returns [m_cY*maxlag] square matrix, where maxlag is the highest lag in the model
**/
{
	decl p = max(this.GetMaxGroupLag(Y1_VAR),
			this.GetMaxGroupLag(Y2_VAR)+m_fdiffY2);		// VAR order, assumes all Y variables have same number of lags
	decl cX0 = columns(this.GetGroup(X_VAR));			// # of exogenous and deterministic regressors
	mC = sizeof(mC)>cX0 ? mC[cX0:][] : mC;				// drop deterministics, exogenous stuff
	decl k = columns(mC);								// # of variables in the VAR
	decl cX1 = m_cX-cX0-m_cY2lags;						// total # of Y1 lagged regressors on RHS
	decl mCY1 = cX1 ? mC[:cX1-1][] : <>;				// coefficients on Y1 lags (if there are any)
	decl mCY2 = p>m_fdiffY2 ? mC[cX1:][] : (zeros(1,k-1)~1);			// coefficients on Y2 lags
	if(m_fdiffY2){										// if TRUE recover the coefficients on constrained variable in levels
		if(p==1)
			mCY2 = zeros(1,k-1)~1;
		else if (k==1)
			mCY2 = 1+mCY2[0]|reversec(-diff0(0|reversec(mCY2))[1:]);
		else
			mCY2 = (mCY2[0][:k-2] ~ (1+mCY2[0][k-1]) )
				|reversec(-diff0(zeros(1,k)|reversec(mCY2))[1:][]);
	}
	decl mClags = mCY1 | mCY2;
	decl mComp = new matrix[k*p][k*p];			// Companion matrix of VAR lag polynomial
	for(decl i = 0; i < k; ++i){				// to compute it, you need to rearrange C
		for(decl j = 0; j < k*p; ++j){
			decl s = imod(j,k)*p+idiv(j,k);
			mComp[i][j] = mClags[s][i];
		}
	}
	if(p>1)
		mComp[k:][] = unit((p-1)*k)~zeros((p-1)*k,k);
	return mComp;
}
CKSVAR::GetCompanionMatrix(const vPar)
/** computes companion matrix by re-arranging coefficients on VAR lags
	@param vPar full parameter vector
	@returns [m_cY*maxlag] square matrix, where maxlag is the highest lag in the model
**/
{
	decl p = max(this.GetMaxGroupLag(Y1_VAR),
			this.GetMaxGroupLag(Y2_VAR)+m_fdiffY2);		// VAR order, assumes all Y variables have same number of lags
	if(p == 0)
		return zeros(m_cY1+1, m_cY1+1);					// if no lags, the companion matrix is trivial
	//else
	decl tau, vC2obs, vC2lat, vbetatil, mC1obs, mC1lat, vdelta, mOm1g2Ch;			// obtain the parameters first
	[tau, vC2obs, vC2lat, vbetatil, mC1obs, mC1lat, vdelta, mOm1g2Ch] = this.GetModelPar(vPar);
	vC2obs[sizeof(vC2obs)-m_cY2lags:][] = vC2lat;
	mC1obs[sizeof(vC2obs)-m_cY2lags:][] = mC1lat;
	decl mC = mC1obs~vC2obs;							// only interested in coefficients on observed lags
	return GetCompanion(mC);
	
}
CKSVAR::GetCompanionRoots(const vPar)
/** computes roots of companion matrices in both regimes
	@param vPar full parameter vector
	@returns array[2] of matrices of eigenvalues
**/
{
	decl p = max(this.GetMaxGroupLag(Y1_VAR),
			this.GetMaxGroupLag(Y2_VAR)+m_fdiffY2);		// VAR order, assumes all Y variables have same number of lags
	if(p == 0)
		return zeros(m_cY1+1, m_cY1+1);					// if no lags, the companion matrix is trivial
	//else
	decl tau, vC2obs, vC2lat, vbetatil, mC1obs, mC1lat, vdelta, mOm1g2Ch;			// obtain the parameters first
	[tau, vC2obs, vC2lat, vbetatil, mC1obs, mC1lat, vdelta, mOm1g2Ch] = this.GetModelPar(vPar);
	decl eig0, eig1;
	// first compute roots of companion matrix of non-ZLB regime
	if(m_fLatLags){
		vC2obs[sizeof(vC2obs)-m_cY2lags:][] += vC2lat;	// get C2bar by adding coef on latent lags on to coef on observed Y2 lags
		mC1obs[sizeof(vC2obs)-m_cY2lags:][] += mC1lat;	// get C1bar by adding coef on latent lags on to coef on observed Y2 lags
	}
	decl mC = mC1obs~vC2obs;							// only interested in coefficients on observed lags
	decl mComp = GetCompanion(mC);
	eigen(mComp, &eig0);
	// next, compute roots of companion matrix at ZLB regime
	if(m_fLatLags){
		vC2obs[sizeof(vC2obs)-m_cY2lags:][] = vC2lat;
		mC1obs[sizeof(vC2obs)-m_cY2lags:][] = mC1lat;
	}
	else{
		vC2obs[sizeof(vC2obs)-m_cY2lags:][] *= 0;
		mC1obs[sizeof(vC2obs)-m_cY2lags:][] *= 0;
	}	
	mC = mC1obs~vC2obs;									// now replace coefs of observed with latent lags (or zeros if no latent lags present)
	mComp = GetCompanion(mC);
	eigen(mComp, &eig1);
	return {eig0, eig1};
}

CKSVAR::distfunc(const avF, const vx)
{
	decl vDF;				// distribution function
	vDF = meanr((m_mY2draws .< vx) .*m_mw);			// compute DF by IS
	avF[0] = vDF - m_dquant;
	return !isnan(vDF);
}
CKSVAR::GetY2star(const quantiles)
/**	Computes smoothed espectation of Y2* and quantiles.
	@param quantiles vector of requested quantiles
	@return matrix [m_cT][1+sizerc(quantiles)]
**/
{
	decl vlik, mY2draws, mw;
	[vlik, mY2draws, mw] = this.SimLik(m_vPar, TRUE);
	decl vsmoothed = meanr(mY2draws.*mw[m_cT-1][]);			// E(Y2*_t|Y_T)
	decl mquants = <>;
//	MaxControl(1000, 10);
	for(decl i = 0; i < sizerc(quantiles); ++i){
		if(quantiles[i] >= 1 || quantiles[i] <= 0)
			continue;
		m_dquant = quantiles[i];
		decl vquants = quantiler(mY2draws, m_dquant);		//start with naive quantiles, not correct
		for(decl t=0; t < m_cT; ++t){
			if(mY2draws[t][] > m_vb[t])						// if unconstrained, you don't need to do anything
				continue;
			m_mY2draws = mY2draws[t][];
			m_mw = mw[m_cT-1][];
			decl quant = vquants[t], DF;
			this.distfunc(&DF, quant);
			decl ires = SolveNLE(this.distfunc, &quant);
			vquants[t] = quant;
		}
		mquants ~= vquants;
	}
	return vsmoothed ~ mquants;
}																 
CKSVAR::Output()
{
    if (!OutputHeader(m_fCSVAR ? "CSVAR" : GetPackageName()))     // returns FALSE if no estimation
        return;

    OutputPar();
    OutputLogLik();
	println("ZLB: ", !isfeq(variance(m_vb),0) ? "variable" : double(m_vb[0]), ", no of observations at ZLB: ", m_cZLB, " (", "%.1d", 100*m_cZLB/m_cT, "%)");
	if(m_fLatLags){
		println("Likelihood evaluated using ", m_fFAPF ? "fully-adapted particle filter." : "sequential importance sampling.");
		if(m_fSimAnn)
			println("Maximization using Simulated Annealing");
		println("no. of MC repl. for simul. lik. ", m_cR);
		println("minimum effective MC size ", min(this.GetWeights(FALSE)));
	}
	if(m_fSimAnn)
		println("MaxSA final temperature ", m_dTemp);
}
CKSVAR::CBRT(const Z)
/**	To be used in cubic(..) Not my code
**/
{
	decl sign = (Z > 0) ? 1 : (Z < 0) ? -1 : 0;
	return fabs(pow(fabs(Z),1/3)) * sign;
}
CKSVAR::cubic(const vA, const avX, const aL)
/**in: vA vector of four coefficients
	out: vX	vector with 1 to 3 real roots
		L = number of real roots
**/
{
	decl THIRD = 1/3;
	decl U = new matrix[3],W, P, Q, DIS, PHI;
	decl i;

	//define cubic root as statement function
	// In C, the function is defined outside of the cubic fct

	// ====determine the degree of the polynomial ====

	if (!isfeq(vA[3], 0.0))
	{
		//cubic problem
		W = vA[2]/vA[3]*THIRD;
		P = pow((vA[1]/vA[3]*THIRD - pow(W,2)),3);
		Q = -.5*(2.0*pow(W,3)-(vA[1]*W-vA[0])/vA[3] );
		DIS = pow(Q,2)+P;
		if ( !(DIS > 0.0 || isfeq(DIS,0)) )	 	// "isfuzzy negative"
		{
			//three real solutions!
			//Confine the argument of ACOS to the interval [-1;1]!
			PHI = acos(min(1.0,max(-1.0,Q/sqrt(-P))));
			P=2.0*pow((-P),(5.e-1*THIRD));
			for (i=0;i<3;i++)	U[i] = P*cos((PHI+2*(i)*M_PI)*THIRD)-W;
			avX[0] = min(U[0], min(U[1], U[2]))
					~ max(min(U[0], U[1]),max( min(U[0], U[2]), min(U[1], U[2])))
					~ max(U[0], max(U[1], U[2]));
			aL[0] = 3;
		}
		else
		{
			// only one real solution!
			DIS = isfeq(DIS,0) ? 0 : sqrt(DIS);
			avX[0] = matrix(CBRT(Q+DIS)+CBRT(Q-DIS)-W);
			aL[0]=1;
		}
	}
	else if (!isfeq(vA[2],0.0))
	{
		// quadratic problem
		P = 0.5*vA[1]/vA[2];
		DIS = pow(P,2)-vA[0]/vA[2];
		if (DIS > 0.0)
		{
			// 2 real solutions
			avX[0] = -P - sqrt(DIS)
					~ -P + sqrt(DIS);
			aL[0]=2;
		}
		else if(isfeq(DIS,0))
		{
			avX[0] = matrix(-P);
			aL[0]=1;
			return 0;		// avoid newton-raphson
		}
		else
		{
			// no real solution
			aL[0]=0;
		}
	}
	else if (!isfeq(vA[1], 0.0))
	{
		//linear equation
		avX[0] =matrix(vA[0]/vA[1]);
		aL[0]=1;
	}
	else
	{
		//no equation
		aL[0]=0;
	}
 /*
  *     ==== perform one step of a newton iteration in order to minimize
  *          round-off errors ====
  */
	for (i=0;i<aL[0];i++)
		avX[0][i] = avX[0][i] - (vA[0]+avX[0][i]*(vA[1]+avX[0][i]*(vA[2]+avX[0][i]*vA[3])))
			/(vA[1]+avX[0][i]*(2.0*vA[2]+avX[0][i]*3.0*vA[3]));

	return 0;
}
