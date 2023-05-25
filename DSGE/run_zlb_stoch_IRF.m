% Testing the effectiveness of unconventional monetary policy in Japan and the United States
% Daisuke Ikeda, Shangshang Li, Sophocles Mavroeidis, and Francesco Zanetti
% American Economic Journal: Macroeconomics

% This code computes impulse responses using OccBin (occbin_20140630).
% The saved impulse responses are used in run_IRF_prop2.
% Set the effectiveness of UMP "lambdas_p = 0, 0.5, or 1" in line 12 below.

clear

%% Set variables and draw shocks
lambdas_p = 1;         % = 0, 0.5, or 1
% ************************************************************************
% 1. In model.mod and model_zlb.modl, "ms" has to be set accordingly. 
%    For example, if lambdas_p=0, "ms" should be set as "ms=0" in line 17.
% 2. Run run_zlb.m first to save the solution using OccBin.
% ************************************************************************

T = 30;                % size of simulated sample size
N = 1000;              % number of simulations
E = zeros(3,T,N);

rng(1)                 
E(:,2:end,:) = randn(3,T-1,N); % shocks
E0 = E;              
IRF = zeros(T,N,5);   % IRF with e_j > 0; [iU piU yU iTU ZU]
IRF0 = zeros(T,N,5);  % IRF with e_j = 0

%% Compute IRFs using OccBin
global M_ opt
opt = 1;               % 0: solve the model every time of the iteration
                       % 1: use the saved solution

modnam = 'model';
modnamstar = 'model_zlb';

constraint = 'iTU-alphaZ_p*ZU<-(pi_p/betau_p-1)';
constraint_relax ='iTU-alphaZ_p*ZU>-(pi_p/betau_p-1)';

% Pick innovation for IRFs
irfshock =char('e_b','e_a','e_i');      % label for innovation for IRFs

maxiter = 100;
tol0 = 1e-8;
paramfile_model

h = timebar('Loop counter','Progress');
tic;
for n=1:N
    
    % Solve nonlinear model
    Et = [E(2,:,n);E(1,:,n);E(3,:,n)]; 
    Et0 = [E0(2,:,n);E0(1,:,n);E0(3,:,n)];
    Et(1,1) = 5;   % demand shock
    Et0(1,1) = 5;  % demand shock
    Et(3,2) = -1;  % monetary policy shock
    Et0(3,2) = 0;  % no monetary policy shock

    SHOCKS = Et'*diag([sigb_p siga_p sigi_p]); 
    SHOCKS0 = Et0'*diag([sigb_p siga_p sigi_p]);
    shockssequence = [SHOCKS;zeros(30,3)];
    shockssequence0 = [SHOCKS0;zeros(30,3)];
    nperiods = T+30;
    
    % Solve model, generate model IRFs
    [zdatalinear,zdatapiecewise,zdatass,oobase_,Mbase_] = ...
        solve_one_constraint(modnam,modnamstar,...
        constraint, constraint_relax,...
        shockssequence,irfshock,nperiods);
    
     [zdatalinear0,zdatapiecewise0,zdatass0,oobase0_,Mbase0_] = ...
        solve_one_constraint(modnam,modnamstar,...
        constraint, constraint_relax,...
        shockssequence0,irfshock,nperiods);
    
    % unpack the IRFs
    for i=1:M_.endo_nbr
        eval([deblank(M_.endo_names{i}),'_p=zdatapiecewise(:,i);']);
        eval([deblank(M_.endo_names{i}),'_ss=zdatass(i);']);
        eval([deblank(M_.endo_names{i}),'_p0=zdatapiecewise0(:,i);']);
    end
    
    
    % Save data [iU piU yU isU]
    IRF(:,n,1) = iU_p(1:T).*400;        % gross nominal short-term rate (APR diff)
    IRF(:,n,2) = piU_p(1:T).*400;       % gross inflation rate (APR diff)
    IRF(:,n,3) = yU_p(1:T)./yU_ss.*100; % output (%)
    IRF(:,n,4) = iTU_p(1:T).*400;       % Taylor rate
    IRF(:,n,5) = ZU_p(1:T).*400;        % FG term
    IRF0(:,n,1) = iU_p0(1:T).*400;
    IRF0(:,n,2) = piU_p0(1:T).*400;
    IRF0(:,n,3) = yU_p0(1:T)./yU_ss.*100;  
    IRF0(:,n,4) = iTU_p0(1:T).*400;       
    IRF0(:,n,5) = ZU_p0(1:T).*400;        
    
    timebar(h,n/N)
end
close(h)
toc


% IRFs
IRFm = zeros(T,3);
for j=1:3
    IRFm(:,j) = mean(IRF(:,:,j),2) - mean(IRF0(:,:,j),2);
end

% Save IRFs
IRFm_ob = IRFm;
IRF_init = IRF(1,1,:);
if lambdas_p == 0
    save IRFlambda0 IRFm_ob IRF_init 
elseif lambdas_p == 1
    save IRFlambda1 IRFm_ob IRF_init
elseif lambdas_p == 0.5
    save IRFlambda05 IRFm_ob IRF_init
end



