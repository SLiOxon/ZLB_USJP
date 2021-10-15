#ifndef CKSVAR_UI_INCLUDED
#define CKSVAR_UI_INCLUDED

#include "cksvar.oxh"

/*-------------------------- CKSVAR_UI : cksvar ---------------------------*/
class CKSVAR_UI : CKSVAR
{
    CKSVAR_UI();															   //initialise
	AddOxCall(const sFunc, ...args);										   //store actions

	// OxPack related
	virtual Output(); 															  //set display of output
	virtual Buffering(const bBufferOn);											 //set buffer for output display
	SendMenu(const sMenu);														//set up dropdown menu for test
	SendSpecials();																 //constants etc. in formulate dialog
	SendVarStatus();															//variable status in formulate dialog
	virtual DoFormulateDlg(const iLagMode, ...);							   //dialog for formulate
	virtual DoSettingsDlg();												   //dialog for model settings
	virtual DoOption(const sOpt, const val);								   //execute the options(blank)
	virtual DoOptionsDlg(const aMoreOptions);								  //dialog for options...
	virtual DoEstimateDlg(const iFirstMethod, const cMethods, const sMethods, //estimation dialog
		const bForcAllowed, const bRecAllowed, const bMaxDlgAllowed);
	TestEstimate();														//Estimate using DoEstimation in CKSVAR
	DoBootstrapTest(m_fKV, m_fCV, m_fEX, m_fboottest);									    //Do Bootstrap test against KSVAR and CSVAR
	virtual DoAutoOutput();											   //call Do AutoGraphics
	virtual DoAutoGraphics();											//***
	virtual DoTestGraphicsDlg();									   //dialog for test/graphic analysis
	virtual DoTestTestDlg();										   //dialog for test/test... and do test
	DoStore(const sOpt, const fQuery=TRUE, const sRename="");		   //store residuals
	virtual DoTestStoreDlg();										  //dialog for test/store residuals etc.
	virtual DoTestFurtherOutputDlg();									//dialog for further output
	virtual ReceiveMenuChoice(const sMenu);								//receive the choice from the test dropdown menu and open corresponding dialogs
	ReceiveModel();													  //receive model specifications and sample
	virtual GetModelSettings();										  //get model setting for model history
	virtual SetModelSettings(const aValues);						  //recall settings from the previous model
	LoadOptions();													  //load options for model setting
	SaveOptions();													  //save options for model setting
	virtual GetOxCode();											  //generate ox code

	virtual DoEstimation(vP);										 //	display and execute of estimation results
	CompanionRoots();												 // compute companion roots
	FAPFLikelihood();
	GetResiduals();
//	SimulateY();
	ImpulseResponses(const irfchoice, const setiden, const setirf);					//draw IRFs
	GetIRFbounds(const dcoverage, const ximin, const ximax, const bootband);
	funcIRF0(const avIRF, const vP);
	fcritIM(const avF, const vX);
	fJacoIM(const avF, const vX);
	
	decl m_dterm1, m_dterm2;												// constants in equation of Imbens and Manski critical value
	decl m_asActions, m_asActionsArgs;									  //to store actions and setting to generate ox code
	decl m_sY2;															 //latent variable
	decl m_cLags, m_xLags;											     //endogenous and exogenous lags
	decl m_dOutlierFactor;												//outlier factor in further output
	decl m_iNS, m_iNT, m_dRT, m_vM, m_vC, m_dT;							//simulated annealing parameters
	decl m_bsR, m_bsT, m_bsP, m_bsST;									//bootstrap parameters
	decl model0;
	decl mBootPar;
};
/*------------------------ END CKSVAR_UI : cksvar -------------------------*/

#endif /* CKSVAR_UI_INCLUDED */
