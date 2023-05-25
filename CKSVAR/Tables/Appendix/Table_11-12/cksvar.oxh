#include <modelbase.oxh>

class CKSVAR : Modelbase
{
	decl m_fOlsen;								// Olsen's (1978) parametrization of the Tobit to make likelihood concave. Default = FALSE
	decl m_fbounded;							// flag to indicate whether to do bounded estimation
	decl m_vLo, m_vHi;							// lower and upper bound for each parameter
	decl m_mY1, m_mY2;							// data
	decl m_cY1;									// # of Y1 eqns (could be zero, in which case you get Tobit)
	decl m_cY2lags;								// columns of m_mY2lags
	decl m_fdiffY2;								// flag, TRUE if Y2 appears in first difference of original censored variable
	decl m_vY2laglev;							// if m_fdiffY2 = TRUE, you need to stare the first lag of Y2 in levels, to adjust the bound
	decl m_vY2stinit;							// initial values of latent process
	decl m_mxi;									// uniform random draws for IS
	decl m_mY2star;								// [T][R] matrix holding draws of Y2* from the its dist conditional on Y1, Y2 and exogenous regressors, computed within LogLik
	decl m_db, m_vb;							// lower bound
	decl m_cR;									// sample size, # of MC replications for Simulated Likelihood, IRFs, Y2*
	decl m_dkappa;								// kink in Y2* eqn -- not identified seperately from the coefficients on the lags. 
	decl m_fLatLags;							// flags, if TRUE, the model contains latent lags and importance sampling is used to calculate likelihood
	decl m_fUseOLS;								// if TRUE (default), compute OLS estimates of VAR and use them as starting values, setting all other coefficients on latent variables and betatilde to zero
	decl m_vIdxLatLags, m_vIdxObsLags;			// vector of indices of coefficients on latent regressors Y2*_j and observed Y2lags
	decl m_vIdxbetatilde;						// vector of indices of coefficients betatilde
	decl m_asbeta; 								// names for contemporaneous impact coefficients of Y2 on Y1
	decl m_aslambda; 							// name for weight on shadow rate in PKSVAR
	decl m_cHorizon; 							// # of periods for computing IRF
	decl m_cOrigin; 							// date for computing IRF
	decl m_mXstar;								// X* particles at IRF origin (to be used to filter out X* if necessary)
	decl m_dImpulse, m_vObsOriginIRF, m_vX2origin, m_vLatOriginIRF, m_mExogIRF, m_crepIRF;		// variables needed to compute (nonlinear) IRF
	decl m_mrannIRF;							// random numbers for IRF computation
	decl m_mrannIRForigin;						// random numbers for origin of IRFs to be used when imposing sign restrictions
	decl m_fSetOriginIRF;						// TRUE if initial values for IRF are set via the GetIRF function
	decl m_vorder; 								// vector [1][m_cY] holding the index of the selected Y variables in the Choleski factorization. Set in GetCholeskiIRF( _ )
	decl m_fCholeski;							// flag, if TRUE, use Choleski identification for IRFs
	decl m_dquant, m_mY2draws, m_mw;						// variables needed to compute smoothed quantiles of Y2*
	decl m_asY1, m_asY2, m_asX;					// array of strings, to hold variable names
	decl m_cZLB;								// # of observations on the ZLB
	decl m_fFAPF;								// flag, if TRUE, use Fully Adapted Particle filter (default is FALSE)
	decl m_fSimAnn;								// flag, if TRUE, use Simulated Annealing (MaxSA)
	decl m_dTemp;								// double, temperature of simulated annealing algorithm
	decl m_fCSVAR;								// flag, if TRUE, impose CSVAR restrictions
	decl m_fPKSVAR;								// flag, if TRUE, impose PKSVAR restrictions
	decl m_dlambda;								// double, coef on shadow rate, efficacy of uncov policy
	decl m_dzeta;								// double, ratio of vol of eps_2 in ZLB vs non-ZLB regime
	decl m_iSol;								// number of solutions for betabar as a function of betatilde, xi=lambda*zeta, and Omega
	decl m_cBootRep;							// number of bootstrap replications
	decl m_fInitTrue;							// if TRUE, use "true" values to initialize bootstrap MLE, otherwise, use OLS
	decl m_fStoreBoot;							// if TRUE, store bootstrap replications, default FALSE
	decl m_sStoreBootFile;						// file name to store bootstrap repls if desired, default "bootstrap.in7", can be altered
	decl m_fBootstrap;							// flag, if TRUE, bootstrap has been performed, results stored in m_mBootData
	decl m_mBootData;							// to hold bootstrap data, parameter estimates, LogLik, LR against CKSVAR (if desired)
	decl m_mIntermBoot;							// variable to hold intermediate parallel bootstrap results, to be used in InsideLoop and ParBoot
	decl m_fStoreBootCovar;						// flag, if TRUE compute and store Covar()
	decl m_asBootNames;							// to hold bootstrap data names, "LogLik", "LR", FreeParNames
	decl m_dBootSeed;							// holds seed for the bootstrap random numbers (optional), default = 0;
	decl m_oUnrestricted;						// holds unrestricted model object
	decl m_fTestUnrestr;						// flag, if TRUE, compute bootstrap simulations for m_oUnrestricted model 
	decl m_cParallel;							// Size of blocks of bootstraps to put in parallel for loop. 1 means no parallelization
	decl m_fSavePerSet;							// if TRUE, and if m_fStoreBoot==TRUE, save bootstrap data in file for every set of bootstrap replications
	decl m_fUseBetabar;							// if TRUE, use the parameterization betabar instead of betatilde in GetModelPar and estimation -- can be used to impose sign restrictions for given value of xi = lambda*zeta, default = FALSE
	decl m_mSigns;								// matrix of 1,-1 and 0 for signs of IRF (0 means unrestricted)
	decl m_vBreakDates, m_aBreakParIdx;			// vector of indices of break dates, array of vectors of indices of parameters that change at each candidate break date 
	decl m_aRegimeDates;						// array of vectors of indices of dates for each regime defined by each break date
	decl m_cBreakPar;							// # of parameters that are allowed to change
	decl m_cBreaks;								// # of breaks (sizeof m_vBreakDates and m_aBreakParIdx)
	decl m_mbetabar;							// to store values of betabar generated within funcIRF when m_fUseBetabar==FALSE;
	
	SetBreaks(const vBreakDates, const aBreakParIdx);	// sets the dates when breaks are allowed and the coefficients that change at those dates
	GetBreakPar(const vpall);					// Returns set of full parameter vectors one for each regime
	SetUseBetabar(const fUseBetabar); 			// sets flag m_fUseBetabar
	GetUseBetabar();				 			// returns flag m_fUseBetabar
	SetOlsen(const fOlsen);						// sets flag m_fOlsen
	GetOlsen();									// returns m_fOlsen
	SetParBounds(const viPar, const vLo, const vHi);   // sets bounds on individual or group of parameters
	GetParBounds();									   // gets parameter bounds
	ClearParBounds();								   // clears parameter bounds
	Estimate();									// Prepare and estimate the model -- exactly as in Modelbase with eprint replaced by print (so it can be disabled during bootstrap)
	SetUnrestrictedModel(const oUnrestricted);	// sets m_oUnrestricted
	DeleteUnrestrictedModel();					// deletes m_oUnrestricted
	GetZLBcount();								// return # of observations at the ZLB
	ClearBootstrap();							// sets m_fBootstrap to FALSE and clears bootstrap data
	SetBootstrap(const cBootRep, const fInitTrue,	// sets # of bootstrap replications, whether to initialize bootstrap MLE at true (original data MLE)
			 const cParallel, const fSavePerSet);	// .. and settings for parallel bootstrap computations
	GetBootstrapData(const fSave);				//	returns the bootstrap data
	SetBootstrapData(const arg);				//	sets the bootstrap data computed and stored in a previous run
	SetBootstrapSeed(const dBootSeed);			// sets a seed for the random numbers used for the bootstrap
	SetStoreBoot(const fStoreBoot, ...);		// sets flag to determine whether to store boostrap estimates, optional ... whether to store covariance
	SetStoreBootFile(const sStoreBootFileName);	// sets file name to store boostrap estimates
	Fixbeta(const vbeta);						// Fixes the value of beta in the model, removes it from list of free parameters
	Freebeta();									// Frees beta in the model, i.e., adds it to list of free parameters
	Setlambda(const dlambda);					// Sets the value of lambda in the model
	Setzeta(const dlambda);						// Sets ratio of vols of eps2 in ZLB vs nonZLB regime, if 1 (default), no change across regimes
	static GetBetabar(const vbetatilde, const xi, const mOmega);		// solves m_cY-th order polynomial equation to find betabar given xi=lambda*zeta
	static GetGamma(const vbetabar, const mOmega);		// return coefficient gammabar
	static GetOmega(const tau, const vdelta, const mOm1g2Ch); // Returns reduced form variance matrix
	virtual GetPackageName();					// Modelbase function
	GetParNames();								// Sets and returns parameter names
	InitData();									// Modelbase function, checks data
	virtual InitPar();							// Modelbase function, initializes parameters for estimation
	SetStartPar(const vPar);					// Modelbase function: allows to set starting values
	SetDiffY2(const fdiffY2);
	SetFAPF(const fFAPF);						// sets flag m_fFAPF
	SetLatentLags(const flatlags);
	SetNumMCreps(const cR);
	Setkappa(const kappa);
	virtual Reparameterize(const avPall, const vPfree);	// Map free parameters to full parameter vector
	CKSVAR();
	LogLik(const vp, const adFunc, const avScore, const amHessian);
	Covar();
	SetSimAnn(const fSimAnn, const aStep);		// if TRUE, estimate using simulated annealing
	SetCSVAR(const fCSVAR);						// sets whether to impose CSVAR restriction (TRUE) or not. Default is FALSE.
	virtual DoEstimation(vPar);
	virtual GetModelPar(const vp);				// returns VAR coefficients and variance parameters for the likelihood
	GetOLS();									// estimate reduced-form VAR coefficients by OLS
	SimLik(const vPar, const fSimulate);		// computes joint density of observables, draws of Y2*, and Importance Sampling weights
	GetWeights(const fPlot);					// plots the weights and effective sample size of the particle filter
	GetAuxPars(const mOm1g2Ch, const vbetatil,
					const vdelta, const tau);	// computes parameters Omega_{1|2}^{-1} (mOm1g2inv), tau2, Xi1inv, arg2, to be used in SimLik
	distfunc(const avF, const x);				// distribution function of Y2*_t conditional on Y_{1:T} (smoothing distribution)
	GetY2star(const quantiles);					// returns smoothed estimates of Y2* and requested quantiles (for error bands)
	TestExcludeLatentLags();					// Wald test for exclusion of latent lags
	TestObservedEqualLatentLags();				// Wald test of equality of coefficients on observed and latent lags of Y2 (necessary but not sufficient condition for purely censored VAR)
	TestPureCensoredVAR();						// Wald test of pure CSVAR (no kink). This entails betatilde = 0 and equality of coefficients on observed and latent lags of Y2. 
	TestAgainstUnrestr();						// Perform LR test against unrestricted model stored in m_oUnrestricted. Should be called after SetUnrestrictedModel, otherwise it does nothing
	SetCholeskiIRF(const fCholeski, const vorder);	// sets flag and order for Choleski identification of IRFs
	SetIRF(const cHorizon, const dImpulse, const cMCreps,	const t);	// store settings to class data members for computation of Impulse Response function to structural shock to Y2* at origin given by data at time t
	GetIRF(const dcover, const fAsy, const fBoot);		// Computes Impulse Response function to structural shock to Y2* at origin given by data at time t
	GetIRFset(const vrange, ...);				// Computes identified set of IRFs over the required range of lambda*zeta
	SetIRForigin(const aOrigin);				// Sets the origin for computation of IRFs, i.e., m_vObsOriginIRF, m_vLatOriginIRF, m_mExogIRF
	GetIRForigin();								// returns array {m_vObsOriginIRF, m_vLatOriginIRF, m_mExogIRF}
	GetBootIRFfromStored(const mBootPar, const dcoverage);	//Computes Bootstrap confidence bands for IRF at desired coverage using stored data
	SelectIRF(const vIRFs, const mscale);		// if multiple IRFs, select the one that is closest to satisfying the sign restrictions
	SetIRFsigns(const mSigns); 					// sets the matrix of signs of IRFs -- matrix of 1s, -1s and 0s (no restrictions)
	SetRannIRForigin(const cT);					// sets the random numbers to use when computing origin of IRFs for sign restrictions
	funcIRF(const avIRF, const vP);				// computes true (nonlinear) IRF, to be used with NumJacobian to obtain st error bands	
	funcLinearIRF(const avIRF, const vP);		// computes pseudo linear IRF, to be used with NumJacobian to obtain st error bands	
	GetLinearIRF(const cHorizon, const amStErr);// obtains pseudo linear IRF using the identified parameters from the kinked SVAR model ignoring the ZLB constraint
	GetCompanion(mC);							// intermediate function used in GetCompanionMatrix
	GetCompanionMatrix(const vPar);				// computes and returns companion VAR matrix, to be used in GetLinearIRF
	GetCompanionRoots(const vPar);				// computes roots of companion matrices in both regimes
	SetBound(const b);							// (Re)sets the lower bound on Y2
	Output();
	SimData(const mrann);						// simulates data from the model, to be used in the bootstrap, and for setting sign restrictions on the IRFs
	Bootstrap();								// performs parametric bootstrap
	ParBoot(const csize, const mEstimData);		// performs parallel for loop for bootstraps -- auxiliary function used twice in Bootstrap()
	InsideLoop(const ib);						// contents of for loop in ParBoot
	virtual cfunc_gt0(const avF, const vP);		// inequality constraints (e.g., sign restrictions), must be greater than zero

	static CBRT(const Z);
	static cubic(const vA, const avX, const aL);
	
public:
	enum { Y1_VAR=50, Y2_VAR, Y2_LAGLEV, Y2_BOUNDS};

};
