%[ys,params,check] = NK_baseline_steadystate(ys,exo,M_,options_)
function [ys,params,check]=model_steadystate(ys,exo,M_,options_)

paramfile_model

%NumberOfParameters = M_.param_nbr;
%for ii = 1:NumberOfParameters
%  paramname = M_.param_names{ii};
%  eval([ paramname ' = M_.params(' int2str(ii) ');']);
%end

NumberOfParameters = M_.param_nbr;
for icount = 1:NumberOfParameters
    M_.params(icount) = eval([M_.param_names{icount}]);
    %eval(['M_.params(' int2str(icount) ')  =  M_.param_names{' int2str(icount) '};'])
end

check=0;

% THIS BLOCK IS MODEL SPECIFIC.
% Here the user has to define the steady state.

zbU = 0;                     % preference shock
zaU = 0;                     % TFP growth rate
yU = 1;                      % output, normalization
hU = 1;                      % hours worked
piU = pi_p;                  % inflation
iU = piU/betau_p-1;          % short-term bonds yield
iTU = iU;
RbLU = RL_R_p*(1+iU);        % long-term bonds yield
zetaU = RL_R_p-1;            % premium of long-term bonds
betar_p = piU/RbLU;          % (*1) prefence discount factor of R-household
kap_p = (du_p-1)*RbLU/du_p;  % (*2) duration parameter
PLU = 1/(RbLU-kap_p);        % price of long-term bonds
bU = b_p;                    % short-term bonds
bLU = PbL_p/PLU;             % long-term bonds
zetab_p = zetaU;             % (*3) premium parameter

wU = alpha_p/lambda_p;                                   % real wage
hrU = (omegau_p*cc_p^(-sigma_p/nu_p)+(1-omegau_p))^(-1); % hours worked, R-household
huU = cc_p^(-sigma_p/nu_p)*hrU;                          % hours worked, U-household
crU = (omegau_p*cc_p+(1-omegau_p))^(-1);                 % consumption, R-household
cuU = cc_p*crU;                                          % consumption, U-household
psi_p = wU/((hrU)^(nu_p)*(crU)^(sigma_p));               % (*4) parameter of the disutility of hours worked
bLrU = (crU-wU*hrU-(1-wU)+(RbLU/piU-1)*PLU*bLU + ((1+iU)/piU-1)*bU - zetaU*PLU*bLU)/((RbLU/piU-1-zetaU*(1-omegau_p))*PLU); % long-term bonds, R-household
bLuU = (bLU-(1-omegau_p)*bLrU)/omegau_p;                 % long-term bonds, U-household

beta_p = (1-omegau_p)*betar_p+omegau_p*betau_p;
FpuU = (cuU)^(-sigma_p)*yU/(1-beta_p*xip_p);
FprU = (crU)^(-sigma_p)*yU/(1-beta_p*xip_p);
KpuU = (cuU)^(-sigma_p)*yU^(1/alpha_p)*wU/(1-beta_p*xip_p);
KprU = (crU)^(-sigma_p)*yU^(1/alpha_p)*wU/(1-beta_p*xip_p);

dbLU = bLU - bLrU;
ZU = 0;

%
% END OF THE MODEL SPECIFIC BLOCK.


% ---- DO NOT CHANGE THIS PART --------------------------------------------
% Here we define the steady state values of the endogenous variables of the model.
%
params=NaN(NumberOfParameters,1);
for iter = 1:length(M_.params) %update parameters set in the file
  eval([ 'params(' num2str(iter) ') = ' M_.param_names{iter} ';' ])
end

NumberOfEndogenousVariables = M_.orig_endo_nbr; %auxiliary variables are set automatically
for ii = 1:NumberOfEndogenousVariables
  varname = M_.endo_names{ii};
  eval(['ys(' int2str(ii) ') = ' varname ';']);
end
%
% ---- END OF THE SECOND MODEL INDEPENDENT BLOCK --------------------------

