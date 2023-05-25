% Testing the effectiveness of unconventional monetary policy in Japan and the United States
% Daisuke Ikeda, Shangshang Li, Sophocles Mavroeidis, and Francesco Zanetti
% American Economic Journal: Macroeconomics

%-----------------------------------------------------------------------------------------------
% Define variables
%-----------------------------------------------------------------------------------------------
var yU,cuU,crU,hU,huU,hrU,wU,iU,iTU,RbLU,piU,bLU,bLuU,bLrU,dbLU,bU,FpuU,FprU,KpuU,KprU,zetaU;  % 21 endogenous variables.
var zbU,zaU;                                                                                   % 2 shocks
var ZU;                                                                                        % 1 variable related to forward guidance
var psi_p, zetab_p, betar_p, kap_p;                                                            % parameter pinnded down by the target values
varexo e_b, e_a, e_i;                                                                          % 3 exogenous disturbances

parameters betau_p sigma_p nu_p PbL_p b_p alpha_p omegau_p xip_p lambda_p phipi_p phiy_p iotap_p pi_p;
parameters rhob_p rhoa_p rhoL_p rhozeta_p rhoi_p sigb_p siga_p sigi_p;
parameters RL_R_p cc_p du_p gam_p;
parameters rhoZ_p alphaZ_p;
parameters eff_p;

@#define ms = 0
% 0: lambda=0 -- ineffective long-term bond purchases (elasticity=0)
% 1: lamdba=1 -- perfectly effective long-term bond purchases (elasticitiy>0) <In the case of no forward guidance (alphaZ_p=0)>
% 2: lambda=eff_p
% 3: between (calibrated values set in paramfile_model.m) 

%@# include "params.mod"

%------------------------------------------------------------------------------------------------
% Model
%------------------------------------------------------------------------------------------------
model; 
  
  % Eq.1: U-Household, short-term bonds
  # betauUp1 = betau_p*exp(zbU(+1));
  %1 = betauUp1*(cuU(+1)/cuU)^(-sigma_p)*(1+(iTU-alphaZ_p*ZU))/piU(+1);
  1 = betauUp1*(cuU(+1)/cuU)^(-sigma_p)*(1+iU)/piU(+1);

  % Eq.2: U-Household, long-term bonds
  1+zetaU = betauUp1*(cuU(+1)/cuU)^(-sigma_p)*RbLU(+1)/(piU(+1))*(RbLU-kap_p)/(RbLU(+1)-kap_p);
  
  % Eq.3: U-Household, labor supply
  wU = psi_p*(huU)^(nu_p)*(cuU)^(sigma_p);
  
  % Eq.4: R-Household, long-term bonds
  # betarUp1 = betar_p*exp(zbU(+1));
  1 = betarUp1*(crU(+1)/crU)^(-sigma_p)*RbLU(+1)/(piU(+1))*(RbLU-kap_p)/(RbLU(+1)-kap_p);
  
  % Eq.5: R-Household, labor supply
  wU = psi_p*(hrU)^(nu_p)*(crU)^(sigma_p);

  % Eq.6: R-Household budget constraint
  # PLU = 1/(RbLU-kap_p);
  # profU = omegau_p*zetaU*PLU*bLuU + yU-wU*yU^(1/alpha_p);
  %crU + PLU*bLrU = PLU*RbLU/(piU)*(-dbLU(-1)) + wU*hrU + profU - ((1+(iTU(-1)-alphaZ_p*ZU(-1)))*bU(-1)/(piU) - (bU + PLU*bLU));
  crU + PLU*bLrU = PLU*RbLU/(piU)*(-dbLU(-1)) + wU*hrU + profU - ((1+iU(-1))*bU(-1)/(piU) - (bU + PLU*bLU));

  % Eq.7: Differece between bLU and bLrU
  dbLU = bLU - bLrU;

  % Eq.8: Price setting, piIU
  # omegar_p = 1-omegau_p;
  # PipU = (piU(-1))^(iotap_p)*(pi_p)^(1-iotap_p)/(piU);
  ((1-xip_p*PipU^(1/(1-lambda_p)))/(1-xip_p))^(1-lambda_p) = (lambda_p/alpha_p*(omegau_p*KpuU+omegar_p*KprU)/(omegau_p*FpuU+omegar_p*FprU))^((1-lambda_p)*alpha_p/(alpha_p-lambda_p));

  % Eq.9: Price setting, FpuU
  # PipUp1 = ((piU)^(iotap_p)*(pi_p)^(1-iotap_p))/piU(+1);
  # beta_p = (1-omegau_p)*betar_p+omegau_p*betau_p;
  FpuU = (cuU)^(-sigma_p)*yU + beta_p*xip_p*PipUp1^(1/(1-lambda_p))*FpuU(+1);

  % Eq.10: Price setting, FprU
  FprU = (crU)^(-sigma_p)*yU + beta_p*xip_p*PipUp1^(1/(1-lambda_p))*FprU(+1);
 
  % Eq.11: Price setting, KpuU
  KpuU = (cuU)^(-sigma_p)*(yU/exp(zaU))^(1/alpha_p)*wU + beta_p*xip_p*PipUp1^(lambda_p/((1-lambda_p)*alpha_p))*KpuU(+1);

  % Eq.12: Price setting, KprU
  KprU = (crU)^(-sigma_p)*(yU/exp(zaU))^(1/alpha_p)*wU + beta_p*xip_p*PipUp1^(lambda_p/((1-lambda_p)*alpha_p))*KprU(+1);
  
  % Eq.13: Firms, production
  yU = exp(zaU)*hU^(alpha_p);
  
  % Eq.14: Long-term bonds market clearing
  bLU = omegau_p*bLuU + (1-omegau_p)*bLrU;

  % Eq.15: Labor market clearing
  hU = omegau_p*huU + (1-omegau_p)*hrU;
  
  % Eq.16: Goods market clearing
  yU = omegau_p*cuU + (1-omegau_p)*crU;

  % Eq.17: Shadow rate
  # iSS = pi_p/betau_p-1;
  %iTU-iSS = rhoi_p*(iTU(-1)-alphaZ_p*ZU(-1)-iSS)+(1-rhoi_p)*(phipi_p*log(piU/pi_p)+phiy_p*log(yU/steady_state(yU))) + e_i;
   @# if ms == 0
      iTU-iSS = rhoi_p*(iU(-1)-iSS)+(1-rhoi_p)*(phipi_p*log(piU/pi_p)+phiy_p*log(yU/steady_state(yU))) + e_i;
  @# endif
  @# if ms == 1
      iTU-iSS = rhoi_p*(iTU(-1)-alphaZ_p*ZU(-1)-iSS)+(1-rhoi_p)*(phipi_p*log(piU/pi_p)+phiy_p*log(yU/steady_state(yU))) + e_i;
  @# endif
  @# if ms == 2
      iTU-iSS = rhoi_p*((1-eff_p)*(iU(-1)-iSS)+eff_p*(iTU(-1)-alphaZ_p*ZU(-1)-iSS))+(1-rhoi_p)*(phipi_p*log(piU/pi_p)+phiy_p*log(yU/steady_state(yU))) + e_i;
  @# endif

  % Eq.18: monetary policy rule
  iU = iTU - alphaZ_p*ZU;

  % Eq.19: Supply of short-term bonds
  bU = b_p;

  % Eq.20: Supply of long-term bonds
  bLU = steady_state(bLU);
  
  % Eq.21: Costs of trading long-term bonds
  %zetaU = zetab_p*(PLU*bLU/bU)^(rhozeta_p);  
  zetaU = zetab_p*(bLU/steady_state(bLU))^(rhozeta_p);

  % Eq.22: Forwarg guidance related term
  ZU = rhoZ_p*ZU(-1) + (iU - iTU);
  %ZU = 0;

  % Shocks
  %------------------------------------------------------------------------
  % Eq.20: Technology shock
  zaU = rhoa_p*zaU(-1) + e_a;

  % Eq.21: Preference shock
  zbU = rhob_p*zbU(-1) + e_b;

  % Paramaters pinned down by target values
  psi_p   = steady_state(psi_p);
  zetab_p = steady_state(zetab_p);
  betar_p = steady_state(betar_p);
  kap_p = steady_state(kap_p);

end;

%-------------------------------------------------------------------------
% Shocks (e_a, e_b, e_i)
%-------------------------------------------------------------------------
shocks;
  var e_b;   stderr 1;    % preference shock
  var e_a;   stderr 1;    % TFP shock
  var e_i;   stderr 1;    % monetary policy shock
end;


steady;


stoch_simul(order=1,nocorr,nomoments,irf=30,print,nograph);
