#include <oxstd.oxh>
#include <oxdraw.oxh>
#include <oxprob.oxh>
#import <maximize>
#import <packages/maxsa/maxsa>
#import <solvenle>

//#import <packages/CKSVAR/CKSVAR>
#import "cksvar"
#include "CKSVAR_UI.oxh"


///////////////////////////////////////////////////////////////////////
// CKSVAR_UI : CKSVAR
CKSVAR_UI::CKSVAR_UI()
{
    CKSVAR::CKSVAR();          // intialize base class

	println("\n---- ", GetPackageName(), " ", GetPackageVersion(),
        " session started at ", time(), " on ", date(), " ----\n");

	SetDrawWindow("CKSVAR");

	m_sY2 = "";
	m_cLags = 2;
	m_xLags = 2;
	m_dOutlierFactor = 3;
	m_iNS=20;
	m_iNT=5;
	m_dRT=0.4;
	m_vM=1;
	m_vC=2;
	m_dT=2;
	m_bsR=99;
	m_bsT=1;
	m_bsP=15;
	m_bsST=1;

	m_asActions = m_asActionsArgs = {};
}
CKSVAR_UI::AddOxCall(const sFunc, ...args)
{
	m_asActions ~= sFunc;
	decl sargs = "(", sep = "", val;
	foreach (val in args)
	{
		sargs ~= sprint(sep, "%v", val);
		sep = ", ";
	}
	sargs ~= ");";
	m_asActionsArgs ~= sargs;

	"OxPackAddOxCode"("\tmodel." ~ sFunc ~ sargs ~ "\n");
}

///////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////
// OxPack related
CKSVAR_UI::Output()
{
	Buffering(11);
	Buffering(1);

	CKSVAR::Output();

	Buffering(0);
	Buffering(10);
}

CKSVAR_UI::Buffering(const iBufferOn)
{
	if (iBufferOn == 11)
		"OxPackSetMarker"(TRUE);
	else if (iBufferOn == 1)
		"OxPackBufferOn"();
	else
		"OxPackBufferOff"();
}

CKSVAR_UI::SendSpecials()
{
	  return {"Constant", "Trend", "CSeasonal" };
}

CKSVAR_UI::SendVarStatus()
{
	return
        {{ "&Y endogenous", 			'Y', STATUS_GROUP + STATUS_MULTIVARIATE, Y_VAR},
         { "&X weakly exogenous", 		'X', STATUS_GROUP + STATUS_DEFAULT, X_VAR},
		 { "&b lower bound", 			'b', STATUS_GROUP + STATUS_ONEONLY, Y2_BOUNDS}
//         { "&U short-run variable", 	'U', STATUS_GROUP + STATUS_SPECIAL, U_VAR}
//         { "&Select By",  'S', STATUS_GROUP + STATUS_ONEONLY, SEL_VAR}
		};
}

CKSVAR_UI::SendMenu(const sMenu)
{
	if (sMenu == "ModelClass")
	{
		return 
			{{ "CKSVAR", 	"modelclass0", 1}
			 };
	}
	
	else if (sMenu == "Test")
	{	return 
			{{ "&Graphic Analysis...", 		"OP_TEST_GRAPHICS"},
//			 { "&Recursive Graphics...", 	"OP_TEST_GRAPHREC"},
//			 0,
//			 { "&Dynamic Analysis...", 		"OP_TEST_DYNAMICS"},
//			 { "&Impulse Responses...", 	"Impulses"},
			 0,
			 { "&Tests...", 				"Tests"},
			 { "Further &Output...", 		"Further Output"},
			 0,
			 { "Store Residuals etc. in D&atabase...", "OP_TEST_STORE"}
			 };
	}
}

CKSVAR_UI::DoFormulateDlg(const iLagMode, ...)
{
	return Modelbase::DoFormulateDlg(-1);
}
CKSVAR_UI::DoSettingsDlg()
{
	decl asoptions, avalues, adlg, asy = {};
	decl cy = "OxPackGetData"("GetGroupLagCount", Y_VAR, 0, 0);
	decl ay = "OxPackGetData"("SelGroup", Y_VAR, 0, 0);
	decl freq = "OxPackGetData"("DbSample")[0];

	// get the names
	for (decl i = 0; i < sizeof(ay); i += 3)
		asy ~= ay[i];
	
	adlg =
	{	
		{ "Lower bound", CTL_GROUP, 1},
		{ "Constrained variable Y2*", 		CTL_SELECT, max(find(asy, m_sY2), 0), asy ~ "- None", 0, "slower" },
		{ "Value of lower bound", 			CTL_DOUBLE,    m_db, "m_db" },
		{ "Lag lengths", CTL_GROUP, 1},
		{ "Latent variable has lags", 		CTL_CHECK, m_fLatLags, "m_fLatLags" },
		{ "Endogenous variables",           CTL_INT,   m_cLags,    "m_cLags"},
		{ "Exogenous variables",            CTL_INT,   m_xLags,    "m_xLags"},
//		{ "Deterministic terms", CTL_GROUP, 1 },
//		{ "Constant and trend", 				CTL_SELECT, m_iDetTrend, m_asDetTrend, 1/*base 1*/, "m_iDetTrend" },
//		{ "Seasonals", 							CTL_SELECT, freq <= 1 ? -1 : GetSeasonals(), "None|Centred seasonals|Seasonal dummies", 0, "m_iSeasType" },
		{ "Options", CTL_GROUP, 0},
		{ "Number of MC replications", 		CTL_VARIABLE, &m_cR, 		"m_cR" },
		{ "Kink in Y2* equation",	 		CTL_VARIABLE, &m_dkappa, 	"m_dkappa" },
		{ "CSVAR model",                    CTL_CHECK, m_fCSVAR, "	m_fCSVAR" }		
	};

	if ("OxPackDialog"("Model Settings", adlg, &asoptions, &avalues))
	{
		m_asActions = m_asActionsArgs = {};

		m_sY2 = asy[avalues[0]];
		m_db = avalues[1];
		m_fLatLags = avalues[2];
		m_cLags = avalues[3];
		m_xLags = avalues[4];
		m_fCSVAR = avalues[7];
		return TRUE;
	}
	return FALSE;
}	
CKSVAR_UI::DoOption(const sOpt, const val)
{
	if (sOpt == "FAPF")
	  m_fFAPF=val[];
	if (sOpt == "m_iNS")
	  m_iNS=val[];
	if (sOpt == "m_iNT")
	  m_iNT=val[];
	if (sOpt == "m_dRT")
	  m_dRT=val[];
	if (sOpt == "m_vM")
	  m_vM=val[];
	if (sOpt == "m_vC")
	  m_vC=val[];
	if (sOpt == "m_dT")
	  m_dT=val[];
	if (sOpt == "m_bsR")
	  m_bsR=val[];
	if (sOpt == "m_bsT")
	  m_bsT=val[];
	if (sOpt == "m_bsP")
	  m_bsP=val[];
	if (sOpt == "m_bsST")
	  m_bsST=val[];
	else
       return 0;
	return 1;
}
CKSVAR_UI::DoOptionsDlg(const aMoreOptions)
{
	decl adlg =
		{	{ "Further options", CTL_GROUP, 1 },
			{ "Evaluate Likelihood by Fully Adapted Particle Filter", CTL_CHECK, m_fFAPF, "FAPF" },
			{ "Settings for Simulated Annealing", CTL_GROUP, 0 },
			{ "Number of cycles", CTL_INT, m_iNS, "m_iNS" },
			{ "Number of iterations before temperature reduction", CTL_INT, m_iNT, "m_iNT" },
			{ "Temperature reduction factor", CTL_DOUBLE, m_dRT, "m_dRT" },
			{ "Step length vector used in initial step", CTL_DOUBLE, m_vM,  "m_vM"  },
			{ "Step length adjustment", CTL_DOUBLE, m_vC,  "m_vC"  },
			{ "Initial temperature", CTL_DOUBLE, m_dT,  "m_dT"  },
			{ "Settings for Bootstrap", CTL_GROUP, 0 },
			{ "Bootstrap Replications", CTL_INT, m_bsR, "m_bsR" },
			{ "use true values to initialize bootstrap MLE", CTL_CHECK, m_bsT,  "m_bsT"  },
			{ "number of iterations in each parallel for loop", CTL_INT, m_bsP,  "m_bsP"  },
			{ "store intermediate results", CTL_CHECK, m_bsST,  "m_bsST"  }
			
		};
	return Modelbase::DoOptionsDlg(adlg);
}
CKSVAR_UI::DoEstimateDlg(const iFirstMethod, const cMethods, const sMethods,
		const bForcAllowed, const bRecAllowed, const bMaxDlgAllowed)
{
	decl ret_val = FALSE;
	
	if (Modelbase::DoEstimateDlg(0, 2, "BFGS|Simulated Annealing", FALSE, FALSE, FALSE))
	{
	    SetLatentLags(m_fLatLags);
		if(m_iMethod == 1)
		   MaxSAControl(1e6, 1),
		   SetSimAnn(TRUE, {m_iNS, m_iNT, m_dRT, m_vM, m_vC, m_dT});
		else
		   SetSimAnn(FALSE, {m_iNS, m_iNT, m_dRT, m_vM, m_vC, m_dT});
		ret_val=TRUE;
	}
	return ret_val;
}

CKSVAR_UI::DoAutoOutput()
{
//	CKSVAR::DoAutoOutput();
	DoAutoGraphics();
}
CKSVAR_UI::DoAutoGraphics()
{
//	SetDrawWindow("CKSVAR Model");
//	GraphicAnalysis(0, TRUE, FALSE, 2, 0, 1);
//	::DrawAdjust(ADJ_AREAMATRIX, 3, m_cY);
//	ShowDrawWindow();
//
//	SetDrawWindow("CKSVAR CVec");
//	GraphicAnalysis(0, 0, 0, 0, 2, 0);
//	ShowDrawWindow();
}
CKSVAR_UI::DoTestGraphicsDlg()
{
	decl asoptions, avalues;
	decl t0 = GetSelStart();		// index over entire db of first obs in estimation sample
	decl iYear0 = ObsYear(t0);
	decl iquarter0 = ObsPeriod(t0);
	decl t1 = GetSelEnd();		// index over entire db of first obs in estimation sample
	decl iYear1 = ObsYear(t1);
	decl iquarter1 = ObsPeriod(t1);
	decl ifre = GetFrequency();
	decl asalem = iYear0~iquarter0~iYear1~iquarter1~ifre;
	decl adlg =
	{	{ "Impulse responses"},
		{ "OLS, ignoring lower bound", 		CTL_CHECK, 1, "OLS" },
		{ "Recursive CKSVAR",  		CTL_CHECK, 1, "CKSVAR" },
		{ "Recursive CSVAR", CTL_CHECK, m_fCSVAR ? 0 : -1, "Chol" },
		{ "CKSVAR set identification",      CTL_CHECK, (m_cY1 < 3) && m_fLatLags ? 0 : -1, "setID"},
		{ "...with sign restrictions",      CTL_CHECK, (m_cY1 < 3) && m_fLatLags ? 0 : -1, "signres"},
		{ "Settings for IRF set identification", CTL_GROUP, 0},
		{ "Lowest admissible value of xi", 		CTL_DOUBLE, 0.0, "lambdalo" },
		{ "Highest admissible value of xi", 		CTL_DOUBLE, 1.0,   "lambdahi" },
		{ "Increment step", 				CTL_DOUBLE, 0.01,   "lambdastep" },
		{ "Number of horizons to impose sign restrictions",	CTL_INT, 1, "signhorz" },
		{ "Impose sign restrictions over period...",   CTL_CHECK, 0, "periodres"},
		{ "",	CTL_SAMPLERANGE, t0, t1, asalem, "restfrom" },
		{ "Sign restrictions on...",             CTL_SELECT, 0, "only Y2|all variables", "sr"},
		{ "Settings for IRFs", CTL_GROUP, 1},
		{ "Number of replications", 		CTL_INT, 1000, "reps" },
		{ "Horizon", 						CTL_INT, 24,   "horz" },
		{ "IRF origin", 	CTL_SAMPLEINDEX, t1, asalem,  "iorigin" },
		{ "Confidence bands",             CTL_SELECT, 0, "None|Asymptotic|Bootstrapped", "bt"},
		{ "Confidence band width(0 ~ 1)",             CTL_DOUBLE, 0.9,   "cfi"},
		{ "Impulse size",             CTL_DOUBLE, -0.25,   "dimpulse"}
	};
		
	if ("OxPackDialog"("Graphic Analysis", adlg, &asoptions, &avalues))
	{
		decl irfchoice = avalues[0:4];
		decl setiden = avalues[5:11];
		decl setirf = avalues[12:17];
		ImpulseResponses(irfchoice, setiden, setirf);
		AddOxCall("ImpulseResponses", irfchoice, setiden, setirf);
		
		return TRUE;
	}
	return FALSE;
}
CKSVAR_UI::DoTestTestDlg()
{
	decl asoptions, avalues;
	decl adlg =
	{	{ "Wald Tests" },
		{ "Exclude latent lags", 		CTL_CHECK, m_fLatLags ? 0 : -1, "exlat" },
		{ "Observed equal latent lags", CTL_CHECK, m_fLatLags ? 0 : -1, "eqlat" },
		{ "Pure censored VAR", 			CTL_CHECK, m_fLatLags ? 0 : -1, "cens" },
		{ "Likelihood Ratio Tests" },
		{ "KSVAR against CKSVAR", 		CTL_CHECK, m_fLatLags && !m_fCSVAR ? 0 : -1, "aksvar" },
		{ "CSVAR against CKSVAR", 		CTL_CHECK, m_fLatLags && !m_fCSVAR ? 0 : -1, "acsvar" },
		{ "Exclude Y2 lags from Y1 equations", 		CTL_CHECK, m_cY1 ? 0 : -1, "aexcl" },
		{ "Bootstrap for Test", CTL_CHECK, 0, "boottest"}
	};
		
	if ("OxPackDialog"("Tests", adlg, &asoptions, &avalues))
	{
		Buffering(TRUE);
		if (avalues[0] == 1)
		{
			TestExcludeLatentLags();
			AddOxCall("TestExcludeLatentLags");
		}
		if (avalues[1] == 1)
		{
			TestObservedEqualLatentLags();
			AddOxCall("TestObservedEqualLatentLags");
		}
		if (avalues[2] == 1)
		{
			TestPureCensoredVAR();
			AddOxCall("TestPureCensoredVAR");
		}
		DoBootstrapTest(avalues[3]==1, avalues[4]==1, avalues[5]==1, avalues[6]==1);

		Buffering(FALSE);

		return TRUE;
	}
	return FALSE;
}
CKSVAR_UI::DoTestFurtherOutputDlg()
{
	decl asoptions, avalues;
	decl adlg =
	{	{ "Further Output" },
		{ "Roots of companion matrix", 		CTL_CHECK, 0, "roots" },
		{ "Re-evaluate the likelihood by FAPF", CTL_CHECK, 0, "evaFAPF" }
	};
		
	if ("OxPackDialog"("Further Output", adlg, &asoptions, &avalues))
	{
		Buffering(TRUE);
		if (avalues[0] == 1)
		{
			CompanionRoots();
			AddOxCall("CompanionRoots");
		}
		if (avalues[1] == 1)
		{
			FAPFLikelihood();
			AddOxCall("FAPFLikelihood");
		}		
		Buffering(FALSE);

		return TRUE;
	}
	return FALSE;
}
CKSVAR_UI::DoStore(const sOpt, const fQuery, const sRename)
{
	decl i, t1 = m_iT1est, t2 = m_iT1est + m_cT - 1;
	decl brename = isstring(sRename) && sRename != "";
	decl asy = m_asY1 ~ m_asY2;

	if (sOpt == "residuals")
	{
		decl mres = GetResiduals();
		for (i = 0; i < m_cY1 + 1; ++i)
		{
			decl s = brename ? sprint(sRename, i + 1) : "V" ~ asy[i];
			"OxPackStore"(mres[][i], t1, t2, s, fQuery);
		}
	}
	else if (sOpt == "fitted")
	{
		decl my = m_mY, mres = GetResiduals();
		for (i = 0; i < m_cY1 + 1; ++i)
		{
			decl s = brename ? sprint(sRename, i + 1) : "F" ~ asy[i];
			"OxPackStore"(my[][i] - mres[][i], t1, t2, s, fQuery);
		}
	}
	else if (sOpt == "bootstrap")
	{
	    decl time = timer();
		decl sbootfile = sprint("bootstrap_", !m_fLatLags ? "k" : m_fCSVAR ? "c" : "ck", "svar_p_",m_cLags,"_R_",m_cR,".in7");
		SetBootstrap(m_bsR, m_bsT, m_bsP, m_bsST);
		SetStoreBoot(TRUE);
		SetStoreBootFile(sbootfile);
		SetBootstrapSeed(13);
		Bootstrap();
		println("\nTime to complete Bootstraps ", timespan(time),"\n");
	}
	else
		return 0;
	return 1;
}
CKSVAR_UI::DoTestStoreDlg()
{
	decl asoptions, avalues;
	if ("OxPackDialog"("Store in Database",
			{	{ "Store in database" },
				{ "Residuals", CTL_CHECK, 0, "residuals" },
				{ "Fitted values", CTL_CHECK, 0, "fitted" },
				{ "Bootstrap replications", CTL_CHECK, 0, "bootstrap"}
			},
			&asoptions, &avalues)
		)
	{
		decl i, s;
		foreach (s in asoptions[i])
			if (avalues[i] > 0)
				DoStore(s);
		return TRUE;
	}
	return FALSE;
}

CKSVAR_UI::ReceiveModel()
{
	// get the selection from OxPack in terms of Y_VAR,X_VAR
	// translate Y_VAR to Y1_VAR,Y2_VAR
	decl ay = "OxPackGetData"("SelGroup", Y_VAR);
	decl ax = "OxPackGetData"("SelGroup", X_VAR);
	decl ab = "OxPackGetData"("SelGroup", Y2_BOUNDS);

	for (decl i = 0; i < sizeof(ay); i += 3)
	{
		if (ay[i] == m_sY2)
			Select(Y2_VAR, {ay[i], 0, m_cLags});
		else
			Select(Y1_VAR, {ay[i], 0, m_cLags});
	}

	for (decl j = 0; j < sizeof(ax); j += 3)
	{
		if (ax[j] == "Constant")
			Select(X_VAR, {ax[j], 0, 0});
		else
			Select(X_VAR, {ax[j], 0, m_xLags});
	}

	Select(Y2_BOUNDS, {ab, 0, 0});

	// default to the full sample
	SetSelSampleMode(SAM_ALLVALID);
	SetSelSample(-1, 0, -1, 0);

	// lag length in CKSVAR not specified through selection, reduce by lag length of Y
	m_iT1sel += m_cLags - max(GetMaxGroupLag(Y1_VAR), GetMaxGroupLag(Y2_VAR));

	SetLatentLags(m_fLatLags);
}

CKSVAR_UI::ReceiveMenuChoice(const sMenu)
{
	if (sMenu == "OP_TEST_GRAPHICS")
	{
		DoTestGraphicsDlg();
	}
	else if (sMenu == "OP_TEST_DYNAMICS")
	{
		println("\nDynamic analysis not implemented yet");
	}
	else if (sMenu == "OP_TEST_SUMMARY")
	{
		AddOxCall("TestSummary");
//		CKSVAR::TestSummary();
	}
	else if (sMenu == "Further Output")
	{
		DoTestFurtherOutputDlg();
	}
	else if (sMenu == "Tests")
	{
		DoTestTestDlg();
	}
	else if (sMenu == "OP_TEST_STORE")
	{
		DoTestStoreDlg();
	}
	else
	{	// allow base class to process unhandled cases
		return Modelbase::ReceiveMenuChoice(sMenu);
	}
	return TRUE;
}
CKSVAR_UI::LoadOptions()
{
	Modelbase::LoadOptions();
	
	m_sY2		= "OxPackReadProfileString"(0,	"m_sY2", 		m_sY2		);
	m_db		= "OxPackReadProfileDouble"(0,	"m_db", 		m_db		);
	m_cLags		= "OxPackReadProfileInt"(0,		"m_cLags", 		m_cLags		);
	m_xLags		= "OxPackReadProfileInt"(0,		"m_xLags", 		m_xLags		);
	m_fLatLags	= "OxPackReadProfileInt"(0,		"m_fLatLags", 	m_fLatLags	);
	m_cR		= "OxPackReadProfileInt"(0,		"m_cR", 		m_cR		);
	m_dkappa	= "OxPackReadProfileDouble"(0,	"m_dkappa", 	m_dkappa	);
    m_fCSVAR    = "OxPackReadProfileInt"(0,	    "m_fCSVAR", 	m_fCSVAR	);
	m_fFAPF     = "OxPackReadProfileInt"(0,	    "m_fFAPF", 	    m_fFAPF	    );
	m_iNS       = "OxPackReadProfileInt"(0,	    "m_iNS", 	    m_iNS	    );
	m_iNT       = "OxPackReadProfileInt"(0,	    "m_iNT", 	    m_iNT	    );
	m_dRT   	= "OxPackReadProfileDouble"(0,	"m_dRT",     	m_dRT   	);
	m_vM    	= "OxPackReadProfileDouble"(0,	"m_vM",     	m_vM    	);
	m_vC    	= "OxPackReadProfileDouble"(0,	"m_vC", 	    m_vC    	);
	m_dT    	= "OxPackReadProfileDouble"(0,	"m_dT", 	    m_dT    	);
	m_bsR       = "OxPackReadProfileInt"(0,	    "m_bsR", 	    m_bsR	    );
	m_bsT       = "OxPackReadProfileInt"(0,	    "m_bsT", 	    m_bsT	    );
	m_bsP       = "OxPackReadProfileInt"(0,	    "m_bsP", 	    m_bsP	    );
	m_bsST      = "OxPackReadProfileInt"(0,	    "m_bsST", 	    m_bsST	    );
}
CKSVAR_UI::SaveOptions()
{
	"OxPackWriteProfileString"(0,	"m_sY2", 		m_sY2		);
	"OxPackWriteProfileDouble"(0,	"m_db", 		m_db		);
	"OxPackWriteProfileInt"(0,		"m_cLags", 		m_cLags		);
	"OxPackWriteProfileInt"(0,		"m_xLags", 		m_xLags		);
	"OxPackWriteProfileInt"(0,		"m_fLatLags", 	m_fLatLags	);
	"OxPackWriteProfileInt"(0,		"m_cR", 		m_cR		);
	"OxPackWriteProfileDouble"(0,	"m_dkappa", 	m_dkappa	);
	"OxPackWriteProfileInt"(0,		"m_fCSVAR", 	m_fCSVAR	);
	"OxPackWriteProfileInt"(0,		"m_fFAPF", 	    m_fFAPF 	);
	"OxPackWriteProfileInt"(0,		"m_iNS", 	    m_iNS    	);
	"OxPackWriteProfileInt"(0,		"m_iNT", 	    m_iNT   	);
	"OxPackWriteProfileDouble"(0,	"m_dRT", 	    m_dRT   	);
	"OxPackWriteProfileDouble"(0,	"m_vM", 	    m_vM    	);
	"OxPackWriteProfileDouble"(0,	"m_vC", 	    m_vC    	);
	"OxPackWriteProfileDouble"(0,	"m_dT", 	    m_dT    	);
	"OxPackWriteProfileInt"(0,		"m_bsR", 	    m_bsR   	);
	"OxPackWriteProfileInt"(0,		"m_bsT", 	    m_bsT   	);
	"OxPackWriteProfileInt"(0,		"m_bsP", 	    m_bsP   	);
	"OxPackWriteProfileInt"(0,		"m_bsST", 	    m_bsST   	);

	Modelbase::SaveOptions();
}
CKSVAR_UI::GetModelSettings()
{
	decl aval = Modelbase::GetModelSettings();

	return aval ~
	{	{	m_cLags,		"m_cLags"		}
	};
}
CKSVAR_UI::SetModelSettings(const aValues)
{
	if (!isarray(aValues))
		return;

	Modelbase::SetModelSettings(aValues);

	for (decl i = 0; i < sizeof(aValues); ++i)
	{
		switch_single (aValues[i][1])
		{	case "m_cLags":			m_cLags 		= aValues[i][0];
		}
	}
}
CKSVAR_UI::GetOxCode()
{
	decl sobj = "model", sbatch = "";
	decl sclass = "CKSVAR";

	sbatch += GetOxDecl(sobj, sclass);
	sbatch += GetOxDatabase(sobj, sclass);
	sbatch += GetOxModel(sobj, sclass);
	sbatch += GetOxModelSettings(sobj, sclass);
	sbatch += GetOxEstimate(sobj, sclass);

	decl s, i;
	foreach (s in m_asActions[i])
		sbatch += "\tmodel." + s + m_asActionsArgs[i] + "\n";

	return "//#import <packages/CKSVAR/cksvar>\n"	+ sbatch;
}
///////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////
CKSVAR_UI::DoEstimation(vP)
{
	decl time = timer();

	//	Report OLS estimates ignoring the ZLB -- for comparison
	decl mc, omega;								// red form VAR coefs and error variance matrix
	[mc, omega] = GetOLS();						// obtain OLS coefficients and variance matrix
	println("\n\nSimple OLS reduced form estimates, ignoring the lower bound",
		"%r", m_asX, "%c", m_asY1 ~ m_asY2, mc);

    println("Number observations at the ZLB: ", double(sumc(m_mY2 .<= m_vb)));
		
	print("\nStarting estimation . . .");

	SetStartPar("OLS");
	
	decl retval = CKSVAR::DoEstimation(vP);

	println(" done in: ", timespan(time), "\n");
	print("Starting covariance evaluation . . .");
	
	return retval;
}

CKSVAR_UI::TestEstimate()
{//Do estimation using DoEstimation in CKSVAR rather than here.
    decl  vpstart, vpfree, estout;
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
    estout = CKSVAR::DoEstimation(vpstart); // do the estimation
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

CKSVAR_UI::DoBootstrapTest(m_fKV, m_fCV, m_fEX, m_fboottest)
{
    model0 = clone(this);
	SetUnrestrictedModel(model0);
	decl Latlags_init = m_fLatLags;
	if(!Latlags_init)
			println("The unrestricted model is KSVAR!");
	if(m_fCSVAR)
			println("The unrestricted model is CSVAR!");
	m_cBootRep = 0;
//    SetBootstrap(m_bsR, m_bsT, m_bsP, m_bsST);

	if(m_fKV){
	    SetLatentLags(FALSE);
		decl time = timer();
    	println("\n*****   Started estimation of KSVAR   *****");
		TestEstimate();
		if(m_fboottest){
			decl sbootfile = sprint("bootstrap_ksvar_p_",m_cLags,"_R_",m_cR,".in7");
			decl fbootfile = fopen(sbootfile, "v");
			SetBootstrap(m_bsR, m_bsT, m_bsP, m_bsST);
			if(!isfile(fbootfile)){				
				SetStoreBoot(TRUE);
				SetStoreBootFile(sbootfile);
				SetBootstrapSeed(13);
				ClearBootstrap();
			}
			else{
				m_mBootData = loadmat(sbootfile, &m_asBootNames);
				m_fBootstrap = TRUE;
				m_fTestUnrestr = TRUE;
			}
			TestAgainstUnrestr();
			m_fBootstrap = FALSE;
			m_fTestUnrestr = FALSE;
			m_cBootRep = 0;
			println("\nTime to complete Bootstraps for KSVAR: ", timespan(time),"\n");
		}
		else{
			TestAgainstUnrestr();
		}
		SetLatentLags(TRUE);
		InitData();
		InitPar();
	}

	if(m_fCV){
	    SetLatentLags(TRUE);
		if(m_cLags==0)
	    	SetLatentLags(FALSE);			// if VAR(0), policy rate cannot be in first difference and model does not have latent lags
		SetCSVAR(TRUE);
		decl time = timer();
    	println("\n*****   Started estimation of CSVAR   *****");
		TestEstimate();
		if(m_fboottest){
			decl sbootfile = sprint("bootstrap_csvar_p_",m_cLags,"_R_",m_cR,".in7");
			decl fbootfile = fopen(sbootfile, "v");
			SetBootstrap(m_bsR, m_bsT, m_bsP, m_bsST);
			if(!isfile(fbootfile)){
				SetStoreBoot(TRUE);
				SetStoreBootFile(sbootfile);
				SetBootstrapSeed(13);
				ClearBootstrap();
			}
			else{
				m_mBootData = loadmat(sbootfile, &m_asBootNames);
				m_fBootstrap = TRUE;
				m_fTestUnrestr = TRUE;		
			}
			TestAgainstUnrestr();
			m_fBootstrap = FALSE;
			m_fTestUnrestr = FALSE;
			m_cBootRep = 0;
			println("\nTime to complete Bootstraps for CSVAR: ", timespan(time),"\n");
		}
		else{
			TestAgainstUnrestr();
		}
		SetLatentLags(TRUE);
		SetCSVAR(FALSE);
		InitData();
		InitPar();
	}

	if(m_fEX){
		decl cY1 = columns(GetGroupLag(CKSVAR::Y1_VAR, 0,0));
		decl p = columns(GetGroupLag(CKSVAR::Y2_VAR, 1, GetMaxGroupLag(CKSVAR::Y2_VAR)));
		Fixbeta(zeros(cY1,1));
		if(p){
			decl vIdxExcludeLags = <>;
			decl asY2lags;
			GetGroupLagNames(CKSVAR::Y2_VAR, 1, p, &asY2lags);
			for(decl i = 0; i < cY1; ++i){
				for(decl j = 0; j < sizeof(asY2lags); ++j)
					vIdxExcludeLags |= find(GetParNames(), sprint("Eq.", i+1, " ", asY2lags[j]));
				if(Latlags_init){
					for(decl k = 0; k < sizeof(asY2lags); ++k)
						vIdxExcludeLags |= find(GetParNames(), sprint("Eq.", i+1, " ", "l", asY2lags[k]));
				}
			}
			FixPar(vIdxExcludeLags, zeros(sizerc(vIdxExcludeLags), 1));
		}
		decl time = timer();
    	println("\n*****   Started estimation of the restricted model   *****");
		TestEstimate();
		if(m_fboottest){
			decl sbootfile = sprint("bootstrap_zlbirr_p_",m_cLags,"_R_",m_cR,".in7");
			decl fbootfile = fopen(sbootfile, "v");
			SetBootstrap(m_bsR, m_bsT, m_bsP, m_bsST);
			if(!isfile(fbootfile)){
				SetStoreBoot(TRUE);
				SetStoreBootFile(sbootfile);
				SetBootstrapSeed(13);
				ClearBootstrap();
			}
			else{
				m_mBootData = loadmat(sbootfile, &m_asBootNames);
				m_fBootstrap = TRUE;
				m_fTestUnrestr = TRUE;
			}
			TestAgainstUnrestr();
			m_fBootstrap = FALSE;
			m_fTestUnrestr = FALSE;
			m_cBootRep = 0;
			println("\nTime to complete Bootstraps for the restricted model ", timespan(time),"\n");
		}
		else{
			TestAgainstUnrestr();
		}
		SetLatentLags(TRUE);	
		InitData();
		InitPar();
	}
}

CKSVAR_UI::CompanionRoots()
{
	decl mcomp = GetCompanionRoots(GetPar());					// obtains companion matrix based on MLE of VAR parameters
	println("Roots of companion matrices in constrained and unconstrained regimes", sortr(cabs(mcomp[0])), sortr(cabs(mcomp[1])));	 // prints abs val of eigenvalues in ascending order
}

CKSVAR_UI::FAPFLikelihood()
{
	SetFAPF(TRUE);
	decl dlikFAPF;
	LogLik(GetPar(), &dlikFAPF, 0, 0);
	decl essFAPF = min(GetWeights(FALSE));
	SetFAPF(FALSE);
	println("likelihood using FAPF ", dlikFAPF);
	println("Minimum effective MC size ", essFAPF);
}

CKSVAR_UI::GetResiduals()
{
	decl tau, vC2obs, vC2lat, vbetatil, mC1obs, mC1lat, vdelta, mOm1g2Ch;			// obtain the parameters first
	[tau, vC2obs, vC2lat, vbetatil, mC1obs, mC1lat, vdelta, mOm1g2Ch] = GetModelPar(GetPar());
	decl m_cY1 = sizec(GetGroupLag(CKSVAR::Y1_VAR, 0, 0));
	decl T = sizer(GetGroupLag(CKSVAR::Y1_VAR, 0, 0));												// # of periods to simulate	 sizer(model.GetGroupLag(CKSVAR::Y1_VAR, 0, 0));
	
	decl mY1lags = GetGroupLag(CKSVAR::Y1_VAR, 1, GetMaxGroupLag(CKSVAR::Y1_VAR));
	decl mY2lags = GetGroupLag(CKSVAR::Y2_VAR, 1, GetMaxGroupLag(CKSVAR::Y2_VAR));
	decl mY1 = GetGroupLag(CKSVAR::Y1_VAR, 0, 0);
	decl mY2 = GetGroupLag(CKSVAR::Y2_VAR, 0, 0);

	decl p = sizec(mY2lags);

	decl vb = m_vb[0:p-1]|m_vb;
	for(decl idx=0; idx<sizer(mY2lags); ++idx){
		for(decl iyx=0; iyx<sizec(mY2lags); ++iyx)
			mY2lags[idx][iyx] = mY2lags[idx][iyx] < vb[idx+p-iyx-1] ? vb[idx+p-iyx-1] : mY2lags[idx][iyx];
	}
	decl mX = GetGroup(CKSVAR::X_VAR)~mY1lags~mY2lags;
	decl mY2st = GetY2star(<>);
	decl mX_star = <>;
	for(decl ip = 0; ip < p; ++ ip){
		mX_star = mX_star~lag(mY2st,ip+1);
	}
	decl varg1b = mY1 - mX*mC1obs - mX_star*mC1lat;
	decl varg2b = mY2st - mX*vC2obs - mX_star*vC2lat;
	decl vzlb = vecindex(mY2 .< vb[p:]);
	varg1b[vzlb][] = mY1[vzlb][] - mX[vzlb][]*mC1obs - mX_star[vzlb][]*mC1lat + (mY2st[vzlb] - vb[vzlb+p])*vbetatil';	
	decl mu = varg1b~varg2b;

	return mu;
	
}

//CKSVAR_UI::SimulateY()
//{
//	decl tau, vC2obs, vC2lat, vbetatil, mC1obs, mC1lat, vdelta, mOm1g2Ch;			// obtain the parameters first
//	[tau, vC2obs, vC2lat, vbetatil, mC1obs, mC1lat, vdelta, mOm1g2Ch] = GetModelPar(GetPar());
//	decl p = sizec(mY2lags);
//	decl mu = GetResiduals()[p:][];
//	decl k = sizec(mu);
//	decl T = sizer(GetGroupLag(CKSVAR::Y1_VAR, 0, 0));
//	decl mY2lags = GetGroupLag(CKSVAR::Y2_VAR, 1, GetMaxGroupLag(CKSVAR::Y2_VAR));
//	decl vb = m_vb[0:p-1]|m_vb;
//	decl mY = new matrix[T+p][k];						// simulated sample must include p initial values to maintain size T. Simulations conditional on observed initial values.
//	decl vy2star = mY[][k-1];
//	decl Y0 = (GetGroupLag(CKSVAR::Y1_VAR, 1, GetMaxGroupLag(CKSVAR::Y1_VAR))
//			~ GetGroupLag(CKSVAR::Y2_VAR, 1, GetMaxGroupLag(CKSVAR::Y2_VAR)))[p][];
//	decl mX0 = GetGroup(CKSVAR::X_VAR);					// exogenous (typically deterministic) regressors, to be treated as fixed in the simulation, not valid unless strictly exogenous
//	decl cX0 = columns(GetGroup(CKSVAR::X_VAR));							// # of exogenous and deterministic regressors
//	decl mX = new matrix[T][columns(Y0)];				// initialize observed lags
//	decl mXst = m_fLatLags ? new matrix[T][p] : 0;	// initialize latent lags (lags of Y2*) Xt*, only necessary if m_cLatLags == TRUE
//	for(decl t = -p; t < T - p - 1; ++t){
//		if(t<0){									// first, store initial values in first p rows of mY
//			mY[t+p][] = Y0[][range(0,k-1)*p-t-1];	// recall that Y0 has Y1_-1,...,Y1_-p, Y2_-1,...,Y2_-p, ..., Yk_-1,...,Yk_-p
//			vy2star[t+p] = Y0[][(k-1)*p-t-1]; 	// and similarly for Y2star
//			continue;
//		}
//		decl it = t;
//		for(decl ik = 0; ik < k; ++ik){
//			for(decl j = 0; j < p; ++j)
//				mX[it][ik*p+j] = mY[t+p-j-1][ik];
//		}
//		if(m_fLatLags){
//			for(decl j = 0; j < p; ++j)					// not executed if m_cY2lags = 0
//				mXst[it][j] = vy2star[t+p-j-1];		
//		}
//		vy2star[t+p] = (cX0 ? mX0[t][]*vC2obs[:cX0-1] : 0) + (sizeof(vC2obs)>cX0 ? mX[it][]*vC2obs[cX0:] : 0) + ( (vC2lat != <>) ? mXst[it][]*vC2lat : 0) + mu[t][k-1];
//		mY[t+p][k-1] = max(vy2star[t+p],vb[t+p]);
//		if(k>1)
//			mY[t+p][:k-2] = (cX0 ? mX0[t][]*mC1obs[:cX0-1][] : 0) + (sizeof(mC1obs)>cX0 ?  mX[it][]*mC1obs[cX0:][] : 0) + ( (mC1lat != <>) ? mXst[it][]*mC1lat : 0)
//					- min(vy2star[t+p] - vb[t+p],0) * vbetatil' +mu[t][:k-2];
//	}
//	return mY;
//}

CKSVAR_UI::ImpulseResponses(const irfchoice, const setiden, const setirf)
{
	decl t0 = GetSelStart();		// index over entire db of first obs in estimation sample
	decl iYear0 = ObsYear(t0);
	decl iquarter0 = ObsPeriod(t0);
	decl iT_irf = setirf[2];		// origin of IRF expressed as observation index
	decl sirf_origin = sprint(ObsYear(iT_irf+t0), "q", ObsPeriod(iT_irf+t0));		// date at which to compute IRF
    SetIRF(/*horizon*/setirf[1], /*dImpulse*/setirf[5], /*#MC reps*/setirf[0], /*Origin t*/iT_irf);
	
	decl mIRF, mIRFlo, mIRFup, amIRF, amIRFlo, amIRFup, amIRFset;
	decl aamIRF = {};
	decl vxi, vxirange, vxirestrange;	    				// hold identified range of xi
	decl asY, asY1, asY2;
	GetGroupLagNames(CKSVAR::Y1_VAR, 0, 0, &asY1);
	GetGroupLagNames(CKSVAR::Y2_VAR, 0, 0, &asY2);
	asY = asY1|asY2;

	if((irfchoice[1] || irfchoice[3] || irfchoice[4]) && setirf[3] == 2)
	{
		println("Getting bootstrapped data...");
		decl sbootfile = sprint("bootstrap_", !m_fLatLags ? "k" : m_fCSVAR ? "c" : "ck", "svar_p_",m_cLags,"_R_",m_cR,".in7");
		decl fbootfile = fopen(sbootfile, "v");

		if(isfile(fbootfile)){
			m_mBootData = loadmat(sbootfile, &m_asBootNames);
			m_fBootstrap = TRUE;
		}
		else{
			decl time = timer();
			SetBootstrap(m_bsR, m_bsT, m_bsP, m_bsST);
			SetStoreBoot(TRUE);
			SetStoreBootFile(sbootfile);
			SetBootstrapSeed(13);
			Bootstrap();
			println("\nTime to complete Bootstraps ", timespan(time),"\n");
		}	
	}

	decl impose0 = int(setiden[5][0]);
	decl impose1 = int(setiden[5][1]);
	decl SRper = range(0,setiden[3]-1);					// range of periods over which you impose sign restrictions

	
	if(irfchoice[4] == 1 && setiden[4]){
		if(impose0 > iT_irf){
			println("The IRF origin is chosen before the starting point for imposing sign restrictions! \nPlease choose the correct origin or the period for imposing the sign restrictions.");
			return;
		}
//		get the narrowest identified IRF set and the corresponding IRFs 
		vxi = range(setiden[0],setiden[1],setiden[2]);
		decl loopend = int(setiden[5][1]-setiden[5][0]+1);
		vxirestrange = vxi;
		for(decl iori = 0; iori < loopend; ++iori){		    
			SetIRF(/*horizon*/setirf[1], /*dImpulse*/setirf[5], /*#MC reps*/setirf[0], /*Origin t*/(impose0+iori));
			SetPrint(FALSE);
			amIRF = GetIRFset(vxirestrange, &vxirange);			// computes IRFs over range of values for xi
			if(iori == 0)
				println("Identified xi without sign restrictions: [", min(vxirange), ", ", max(vxirange), "]");
			decl msignrestri = 0;
			if(!setiden[6] || (m_cY1 == 0)){
				msignrestri = (!isdotfeq(amIRF[m_cY1][][SRper],0) .&& amIRF[m_cY1][][SRper] .> 0);
			}
			else{
				for(decl ik = 0; ik < m_cY1; ++ik)
					msignrestri = msignrestri .|| (!isdotfeq(amIRF[ik][][SRper],0) .&& amIRF[ik][][SRper] .< 0);
				msignrestri = msignrestri .|| (!isdotfeq(amIRF[m_cY1][][SRper],0) .&& amIRF[m_cY1][][SRper] .> 0);
			}
			vxirestrange = deleteifr(vxirange,msignrestri);
			vxirestrange = range(min(vxirestrange),max(vxirestrange),setiden[2]);		// do this to remove duplicates
						
			println(ObsYear(impose0+t0+iori), "q", ObsPeriod(impose0+t0+iori),
			": Identified set for lambda with sign restriction: [", min(vxirestrange), ", ", max(vxirestrange),"]");
		}
		vxi = range(min(vxirestrange),max(vxirestrange),setiden[2]);			// make the grid finer for US to fill the space
		println("The narrowest set of xi : [", min(vxirestrange), ", ", max(vxirestrange), "]");

		println("Now computing IRFs at all dates over the above range of xi...");

		decl msignrestr = 0;
		// first collect all IRFs and sign-restrictions across all dates
		for(decl iori = 0; iori < loopend; ++iori){		    
			SetIRF(/*horizon*/setirf[1], /*dImpulse*/setirf[5], /*#MC reps*/setirf[0], /*Origin t*/impose0+iori);
			amIRF = GetIRFset(vxi, &vxirange);			// computes IRFs over range of values for xi 
 			aamIRF ~= {amIRF};									// store IRFs to avoid repeating later
			// evaluate sign restrictions
			decl msignrestri = 0;
			if(!setiden[6] || (m_cY1 == 0)){
				msignrestri = (!isdotfeq(amIRF[m_cY1][][SRper],0) .&& amIRF[m_cY1][][SRper] .> 0);
			}
			else{
				for(decl ik = 0; ik < m_cY1; ++ik)
					msignrestri = msignrestri .|| (!isdotfeq(amIRF[ik][][SRper],0) .&& amIRF[ik][][SRper] .< 0);
				msignrestri = msignrestri .|| (!isdotfeq(amIRF[m_cY1][][SRper],0) .&& amIRF[m_cY1][][SRper] .> 0);
			}
			msignrestr = msignrestri .|| msignrestr;  			// combine restrictions from previous cases
		}
		// then remove IRFs that fail the joint sign restrictions
		for(decl iori = 0; iori < loopend; ++iori)	    
		    for(decl j = 0; j < m_cY1+1; ++j){					
		    	aamIRF[iori][j] = deleteifr(aamIRF[iori][j], msignrestr);
			}
	}
 
   if (irfchoice[3] == 1 || irfchoice[4] == 1)
   {
		SetIRF(/*horizon*/setirf[1], /*dImpulse*/setirf[5], /*#MC reps*/setirf[0], /*Origin t*/iT_irf);
		Setlambda(0);
		mIRF = GetIRF(0, FALSE, FALSE);
		decl vlambda = range(setiden[0],setiden[1],setiden[2]);					// admissible values of efficacy of unconv policy 
//		vlambda = range(0,0.45,.01)~range(0.451, 0.49, 0.001)~range(0.4901,0.506, 0.0001)~range(0.507, 1, 0.001);		// (adjusted range to make sure graph doesn't take too much memory)				
		SetPrint(FALSE);										// disable printing to avoid printing number of solutions in intermediate steps
		amIRFset = GetIRFset(vlambda, &vlambda);			// computes IRFs over range of values for xi
		SetPrint(TRUE);
		println("Identified set for xi without sign restrictions: [", min(vlambda), ", ", max(vlambda),"]");
		SetDrawWindow(sprint("SetID IRF CKSVAR(",m_cLags,")"));
		decl mlimits = new matrix[m_cY1+1][2];
		for(decl j = 0; j < m_cY1+1; ++j){
			mlimits[j][] = min(amIRFset[j])~max(amIRFset[j]);					// store limits for next plot
			DrawTitle(j, sprint("Monetary policy shock to ", asY[j], " in ", sirf_origin));
			DrawMatrix(j, mIRF[j][], "$\\lambda=0$", 0, 1, 2, 3);								 // IRF with unconventional policy weight lambda=0
			DrawMatrix(j, amIRFset[j], "", 0, 1, 0, 2*ones(sizeof(amIRFset[j]),1));					 // IRFs at all other values of lambda
//			DrawMatrix(j,  minc(amIRF[j])|maxc(amIRF[j]), "", 0, 1, 0, <2;2>);				 // use this line instead if you just want to plot the boundaries of the ID set
			DrawAxis(j, TRUE, .NaN, 0, setirf[1], 2, 2, 0, 0);
			DrawLegend(j, 400,0, j);
			DrawAdjust(ADJ_AREA_Y, j, mlimits[j][0] > 0 ? mlimits[j][0]/1.1 : mlimits[j][0]*1.1, mlimits[j][1] > 0 ? mlimits[j][1]*1.1 : mlimits[j][1]/1.1);
		}
		ShowDrawWindow();
	if (irfchoice[4] == 1)
	{
	// plot identified set when you impose the restirction that the effect of MPS on interest rate is nonnegative
		if(setiden[4]){
			decl mlimits = new matrix[m_cY1+1][2];
			if(setirf[3]){
				[amIRFlo, amIRFup] = GetIRFbounds(setirf[4], min(vxi), max(vxi), setirf[3]-1);
				if(m_cY1)
					amIRFlo[:m_cY1-1][SRper] = amIRFlo[:m_cY1-1][SRper] .> 0 .? amIRFlo[:m_cY1-1][SRper] .: 0;			// impose sign restrictions on C.I. as in Granziera et al step 2 on p.1096
				amIRFup[m_cY1][SRper] = amIRFup[m_cY1][SRper] .< 0 .? amIRFup[m_cY1][SRper] .: 0;				// impose sign restrictions on C.I. as in Granziera et al step 2 on p.1096
			}
			SetDrawWindow(sprint("SetID IRF CKSVAR(",m_cLags,") with sign restriction"));
			for(decl j = 0; j < m_cY1+1; ++j){
//				DrawTitle(2*j+1, "... with sign restriction");
				DrawTitle(j, sprint("Monetary policy shock to ", asY[j], " in ", sirf_origin));
				mlimits[j][] = min(setirf[3] ? amIRFlo[j][] : aamIRF[iT_irf-impose0][j])~max(setirf[3] ? amIRFup[j][] : aamIRF[iT_irf-impose0][j]);							// store limits for next plot
				mlimits[j][0] = mlimits[j][0].>0 .? mlimits[j][0]/1.1 .: mlimits[j][0]*1.1;		// set range for plots
				mlimits[j][1] = mlimits[j][1].>0 .? mlimits[j][1]*1.1 .: mlimits[j][1]/1.1;		// set range for plots
				DrawMatrix(j, aamIRF[iT_irf-impose0][j], "", 0, 1, 0, 2*ones(sizeof(aamIRF[iT_irf-impose0][j]),1));					 // IRFs at all other values of xi
				if(setirf[3]){
					DrawMatrix(j, amIRFlo[j][], sprint(100*setirf[4],"\% error bands"), 0, 1, 0, 5);							
					DrawMatrix(j, amIRFup[j][], "", 0, 1, 0, 5);							 
				}
				DrawAxis(j, TRUE, .NaN, 0, setirf[1], 2, 2, 0, 0);
				DrawLegend(j, 400,0, j);
				DrawAdjust(ADJ_AREA_Y, j, mlimits[j][0], mlimits[j][1]);
			}
			ShowDrawWindow();
		}
		else{	
			decl msignrestri = 0;
			if(!setiden[6] || (m_cY1 == 0)){
				msignrestri = (!isdotfeq(amIRFset[m_cY1][][SRper],0) .&& amIRFset[m_cY1][][SRper] .> 0);
			}
			else{
				for(decl ik = 0; ik < m_cY1; ++ik)
					msignrestri = msignrestri .|| (!isdotfeq(amIRFset[ik][][SRper],0) .&& amIRFset[ik][][SRper] .< 0);
				msignrestri = msignrestri .|| (!isdotfeq(amIRFset[m_cY1][][SRper],0) .&& amIRFset[m_cY1][][SRper] .> 0);
			}
			vlambda = deleteifr(vlambda,msignrestri);
			vlambda = range(min(vlambda),max(vlambda),setiden[2]);
			println("Identified set for lambda with sign restriction: [", min(vlambda), ", ", max(vlambda),"]");
			decl mlimits = new matrix[m_cY1+1][2];
			if(setirf[3]){
				[amIRFlo, amIRFup] = GetIRFbounds(setirf[4], min(vlambda), max(vlambda), setirf[3]-1);
				if(m_cY1)
					amIRFlo[:m_cY1-1][SRper] = amIRFlo[:m_cY1-1][SRper] .> 0 .? amIRFlo[:m_cY1-1][SRper] .: 0;			// impose sign restrictions on C.I. as in Granziera et al step 2 on p.1096
				amIRFup[m_cY1][SRper] = amIRFup[m_cY1][SRper] .< 0 .? amIRFup[m_cY1][SRper] .: 0;				// impose sign restrictions on C.I. as in Granziera et al step 2 on p.1096
			}
			for(decl j = 0; j < m_cY1+1; ++j)
				amIRFset[j] = deleteifr(amIRFset[j], msignrestri);			// delete IRFs with negative impact effect of policy shock on interest rate
			SetDrawWindow(sprint("SetID IRF CKSVAR(",m_cLags,") with sign restriction"));
			for(decl j = 0; j < m_cY1+1; ++j){
//				DrawTitle(j, sprint("Monetary policy shock to ", {"Inflation", "Unemployment", "Fed Funds rate"}[j], " in ", sirf_origin));
//				DrawTitle(2*j+1, "... with sign restriction");
				DrawTitle(j, sprint("Monetary policy shock to ", asY[j], " in ", sirf_origin));
				mlimits[j][] = min(setirf[3] ? amIRFlo[j][] : amIRFset[j])~max(setirf[3] ? amIRFup[j][] : amIRFset[j]);							// store limits for next plot
				mlimits[j][0] = mlimits[j][0].>0 .? mlimits[j][0]/1.1 .: mlimits[j][0]*1.1;		// set range for plots
				mlimits[j][1] = mlimits[j][1].>0 .? mlimits[j][1]*1.1 .: mlimits[j][1]/1.1;		// set range for plots
				DrawMatrix(j, amIRFset[j], "", 0, 1, 0, 2*ones(sizeof(amIRFset[j]),1));				// IRFs at all other values of lambda
//				DrawMatrix(2*j+1,  minc(amIRF[j])|maxc(amIRF[j]), "", 0, 1, 0, <2;2>);				// use this line instead if you just want to plot the boundaries of the ID set
				if(setirf[3]){
					DrawMatrix(j, amIRFlo[j][], sprint(100*setirf[4],"\% error bands"), 0, 1, 0, 5);							
					DrawMatrix(j, amIRFup[j][], "", 0, 1, 0, 5);	
				}
				DrawAxis(j, TRUE, .NaN, 0, setirf[1], 2, 2, 0, 0);
				DrawLegend(j, 400,0, TRUE);
				DrawAdjust(ADJ_AREA_Y, j, mlimits[j][0], mlimits[j][1]);
	 		}
		}
		ShowDrawWindow();
	}
   }

   SetPrint(TRUE);

   SetIRF(/*horizon*/setirf[1], /*dImpulse*/setirf[5], /*#MC reps*/setirf[0], /*Origin t*/iT_irf);
   
   if (irfchoice[1] == 1)
   {
	if(setirf[3] == 0)
		mIRF = GetIRF(0, FALSE, FALSE);
	if(setirf[3] == 1)
		[mIRF,mIRFlo,mIRFup] = GetIRF(setirf[4], TRUE, FALSE);
	if(setirf[3] == 2)
		[mIRF,mIRFlo,mIRFup] = GetIRF(setirf[4], FALSE, TRUE);
	SetDrawWindow(sprint("IRFs CKSVAR(", m_cLags, ")"));
	for(decl j = 0; j < m_cY1+1; ++j){
		DrawTitle(j, sprint("Monetary policy shock to ", asY[j], " in ", sirf_origin));
		DrawMatrix(j, mIRF[j][], "CKSVAR", 0, 1, 0, 2);						 // draw IRF at point estimates
		if(setirf[3] == 1){
			DrawMatrix(j, mIRFlo[j][], sprint(100*setirf[4],"\% asymptotic error bands"), 0, 1, 0, 5);							 // draw lower bound of bootstrap CI
			DrawMatrix(j, mIRFup[j][], "", 0, 1, 0, 5);							 // draw upper bound of bootstrap CI
		}
		if(setirf[3] == 2){
			DrawMatrix(j, mIRFlo[j][], sprint(100*setirf[4],"\% bootstrap error bands"), 0, 1, 0, 5);							 // draw lower bound of bootstrap CI
			DrawMatrix(j, mIRFup[j][], "", 0, 1, 0, 5);							 // draw upper bound of bootstrap CI
		}		
		DrawAxis(j, TRUE, .NaN, 0, setirf[1], 2, 2, 0, 0);
		DrawLegend(j, 400,0, j);
	 }
//	   ShowDrawWindow();
   }

   if (irfchoice[2] == 1)
   {// add choleski IRFs
//	model.Fixbeta(zeros(sizeof(asY1),1));				// impose recursive ID (policy rate ordered last)
	SetCSVAR(TRUE);								// impose CSVAR
	Estimate();									// reestimate the model
	decl mIRFrec = GetIRF(/*setirf[4]*/0,/*Asy error bands*/FALSE,/*Bootstrap error bands*/FALSE);		// return only point estimates of recursive ID, for comparison
	// compute IRFs from naive OLS estimation
	SetCSVAR(FALSE);								// remove all latent variables
    for(decl j = 0; j < m_cY1+1; ++j){
		DrawMatrix(j, mIRFrec[j][], "Recursive, CSVAR", 0, 1, 0, 3);					 // draw recursive IRF 
		DrawAxis(j, TRUE, .NaN, 0, setirf[1], 2, 2, 0, 0);
		DrawAdjust(ADJ_LEGEND, 0, 2, -1, -1, 1);		
	}
//	ShowDrawWindow();
   }

   if (irfchoice[0] == 1)
   {
	SetCSVAR(FALSE);								// remove all latent variables
	SetLatentLags(FALSE);							// disable latent lags, to estimate by OLS
	InitPar();									// computes OLS estimates
	decl mIRFols = GetLinearIRF(setirf[1],/*&<st. errors>*/0);		// compute IRFs with choleski factorizations
	for(decl j = 0; j < m_cY1+1; ++j){
		DrawMatrix(j, mIRFols[j][], "Recursive, OLS", 0, 1, 2, 4);						 // draw pseudo-IRF using OLS
		DrawAxis(j, TRUE, .NaN, 0, setirf[1], 2, 2, 0, 0);
		DrawAdjust(ADJ_LEGEND, 0, 2, -1, -1, 1);		
	}
	SetLatentLags(TRUE);
//	ShowDrawWindow();
   }
   ShowDrawWindow();

   SetPrint(TRUE);
   
}
CKSVAR_UI::fcritIM(const avF, const vX)
/** equation that determines the critical value in Imbens and Manski (2004)
**/
{
	avF[0] = probn(vX+m_dterm1) - probn(-vX) - m_dterm2;
	return TRUE;
}
CKSVAR_UI::fJacoIM(const avF, const vX)
/** Jacobian of equation for Imbens and Manski critical value
**/
{
	avF[0] = densn(vX+m_dterm1) + densn(vX);
	return TRUE;
}
CKSVAR_UI::funcIRF0(const avIRF, const vP)
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
CKSVAR_UI::GetIRFbounds(const dcoverage, const ximin, const ximax, const bootband)
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
		if(!m_fBootstrap || !bootband){			// if you don't have bootstrapped samples, use delta method (prior estimation suggests Delta is accurate for Japan, i.e., ~ bootstrap, but not for US -- too big s.e.s)
			decl vJac;
			decl ires = NumJacobian(funcIRF0, this.GetFreePar(), &vJac);
			decl mIRFvar = vJac*this.GetCovar()*vJac';			// variance matrix of IRFs by the delta method
			msterr = ires ? shape(sqrt(diagonal(mIRFvar)), m_cY1+1, m_cHorizon+1) : constant(.NaN, m_cY1+1, m_cHorizon+1);	// to hold matrix of asymptotic standard errors
		}
		if(m_fBootstrap && bootband){
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
///////////////////////////////////////////////////////////////////////
