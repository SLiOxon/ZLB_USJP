betau_p = (1.01)^(-1/4); %(1.0025)^(-1/4); % preference discount factor, unrestricted households
sigma_p = 2;               % relative risk aversion
nu_p    = 2;               % inverse of labor supply elasticity 
PbL_p = 1*4;               % long-term bonds (relative to output)
b_p = 0.25*4;              % short-term bonds (relative to output) 
alpha_p = 1;               % curvature of production function
omegau_p = 1e-7;           % the share of U-households
xip_p = 0.75;              % Calvo parameter
lambda_p = 1.15;           % markup
iotap_p = 0;               % indexation
phipi_p = 1.5;             % Taylor rule, inflation
phiy_p = 0.5;              % Taylor rule, output
pi_p = 1 + 0/400;          % inflation 

rhoL_p = 0.9;              % persistence, long-term bonds supply
rhoa_p = 0.9;              % persistence, technology shock
rhob_p = 0.9;              % persistence, preference shock
rhozeta_p = 6;             % elasticity of costs of trading long-term bonds
rhoi_p = 0.7;              % interest rate smoothing

RL_R_p = 1 + 1/400;        % RL/R
cc_p = 1;                  % cc_p = cuU/crU in steady state     
du_p = 10*4;               % duration of long-term bonds    

gam_p = 20;                 % -1 annual percent interest rate = gam_p/4 percent increase in asset purchase. (gam_p = 20)

rhoZ_p = 0;                % forward guidance persistence 
alphaZ_p = 0;              % forward guidance strength (0 = no forward guidance)

eff_p = 0.5;              % Effectiveness of UMP (= lambdas)

sigb_p = 0.25/400; 
siga_p = 0.25/100;
sigi_p = 0.25/400;

warning off

