% Testing the effectiveness of unconventional monetary policy in Japan and the United States
% Daisuke Ikeda, Shangshang Li, Sophocles Mavroeidis, and Francesco Zanetti
% American Economic Journal: Macroeconomics

close all
clear

% ************* Set path *******************
dir = 'occbin_20140630\toolkit_files';
path(dir,path);
% ******************************************

lambdas_p = 1;         % = 0, 0.5, or 1
% For saving results used in run_zlb_prop2.m (IRF under Proposition 2),
% "ms" in model.mod and model_zlb.mod should be set as "ms=0." 
% For saving results used in run_zlb_stoch_IRF.m (IRF under OccBin), "ms"
% should be set in accordance with lambdas_p. For example, if lambdas_p=1,
% "ms=1."

global M_ opt

% option
opt = 0;               % 0: solve the model every time of the iteration

modnam = 'model';
modnamstar = 'model_zlb';

constraint = 'iTU -alphaZ_p*ZU<-(pi_p/betau_p-1)';
constraint_relax ='iTU -alphaZ_p*ZU>-(pi_p/betau_p-1)';

% Pick innovation for IRFs
irfshock =char('e_b','e_a','e_i');      % label for innovation for IRFs

maxiter = 100;
tol0 = 1e-8;


% Solve nonlinear model
paramfile_model  % read parameter values

SD = diag([sigb_p siga_p sigi_p]);
SHOCKS = [ zeros(5,3)
   8 0 0
   0 0 0
  zeros(20,3) ]*SD;

shockssequence = SHOCKS;
nperiods = 30;


% Solve model, generate model IRFs
[zdatalinear, zdatapiecewise, zdatass, oobase_, Mbase_] = ...
  solve_one_constraint(modnam,modnamstar,...
  constraint, constraint_relax,...
  shockssequence,irfshock,nperiods);


% unpack the IRFs                          
for i=1:M_.endo_nbr
  eval([deblank(M_.endo_names{i}),'_u=zdatalinear(:,i);']);
  eval([deblank(M_.endo_names{i}),'_p=zdatapiecewise(:,i);']);
  eval([deblank(M_.endo_names{i}),'_ss=zdatass(i);']);
end

%% Check Condition 1 (ABCs (and Ds) of understanding VARs, AER) and derive matrices in Proposition 2
As = oobase_.dr.ghx;
Bs = oobase_.dr.ghu;
[ny,nk] = size(As);
if nk == 7 % if ms = 0 in model.mod and model_zlb.mod
    % remove [dbLU bU] from the set of state variables because they are irrelevant for [isU yU piU]
    kindex = [oobase_.dr.state_var(1) oobase_.dr.state_var(4:end)];
    As = [As(:,1) As(:,4:end)];
    nk = nk-2;
    Cs = zeros(nk,ny);
    Ds = zeros(3,ny); Ds2 = zeros(3,ny);
    for j=1:nk
        Cs(j,oobase_.dr.inv_order_var(kindex(j)))=1;
    end
    [~,i_iU] = select('iU',Mbase_.endo_names);
    [~,i_yU] = select('yU',Mbase_.endo_names);
    [~,i_piU] = select('piU',Mbase_.endo_names);
    [~,i_RbLU] = select('RbLU',Mbase_.endo_names);
    yindex = [i_yU; i_piU; i_iU];
    yindex2 = [i_yU; i_piU; i_RbLU];
    for j=1:3
        Ds(j,oobase_.dr.inv_order_var(yindex(j)))=1;
        Ds2(j,oobase_.dr.inv_order_var(yindex2(j)))=1;
    end
    A = Cs*As;
    B = Cs*Bs;
    C = Ds*As; C2 = Ds2*As;
    D = Ds*Bs; D2 = Ds2*Bs;
    E = A - B/D*C;  E2 = A - B/D2*C2;
    [~,F] = eig(E); [~,F2] = eig(E2);
    eigE = diag(F); eigE2 = diag(F2);
    max_eig = max(abs(eigE));
    max_eig2 = max(abs(eigE2));
    if max_eig > 1
        disp('Condition A is not satisfied (explosive eigenvalues) with [iU yU piU]')
    else
        disp('Condition A is satisfied (|eigenvalues|<1) with [iU yU piU]')
    end
    if max_eig2 > 1
        disp('Condition A is not satisfied (explosive eigenvalues) with [RbLU yU piU]')
    else
        disp('Condition A is satisfied (|eigenvalues|<1) with [RbLU yU piU]')
    end
    
    %Proposition 1 and 2
    % endogenous variables Y = [y pi is]
    % state variables X = [is za zb]
    % shock eps = [monetary supply demand]
    C = [C(:,1:2) C(:,5)];
    D = [D(:,3) D(:,2) D(:,1)];
    A = [C(3,:);0 rhoa_p 0; 0 0 rhob_p];
    B = [D(3,:); 0 1 0; 0 0 1];
    E = A - B/D*C;
    eigE = diag(F);
    max_eig = max(abs(eigE));
    if max_eig > 1
        disp('Condition A is not satisfied (explosive eigenvalues) with [iU yU piU]')
    end
    G = C*B/D;
    
    % Long-term yield; E(RbL(t+1)) = EC(1)*is(t) + EC(2)*za(t) + EC(3)*zb(t)
    RbL.EC = [C2(3,1)*C(3,1) (C2(3,1)*C(3,2)+C2(3,2))*rhoa_p (C2(3,1)*C(3,3)+C2(3,5))*rhob_p];
    RbL.Rho = diag([rhoa_p rhob_p]);
    RbL.w = kap_p_ss/RbLU_ss;
    RbL.ss = RbLU_ss;
    
    % Proposition 2 (see Appendix A5)
    betar_p = betau_p/(RL_R_p);
    beta_p = (1-omegau_p)*betar_p+omegau_p*betau_p;
    kap_p = (1-beta_p*xip_p)*(1-xip_p)/xip_p*(sigma_p+1/(1/nu_p));
    chib_p = rhob_p/sigma_p;
    chia_p = (1-beta_p*xip_p)*(1-xip_p)/xip_p*(1+(1/nu_p))/(1/nu_p);
    H1 = [1 0;-kap_p 1];
    H2 = [-1/sigma_p+C(1,1)+C(2,1)/sigma_p;beta_p*C(2,1)];
    H4 = [rhoa_p*D(1,2)+rhoa_p*D(2,2)/sigma_p, rhob_p*D(1,3)+rhob_p*D(2,3)/sigma_p-chib_p;
        beta_p*rhoa_p*D(2,2)-chia_p, beta_p*rhob_p*D(2,3)];
    H3 = H4*diag([rhoa_p rhob_p]);
    
    J = [(beta_p*G(2,1)+kap_p)/chia_p (beta_p*G(2,2)-1)/chia_p beta_p*G(2,3)/chia_p;
        (G(1,1)+G(2,1)/sigma_p-1)/chib_p (G(1,2)+G(2,2)/sigma_p)/chib_p (G(1,3)+G(2,3)/sigma_p-1/sigma_p)/chib_p];
    
    % Y1 = [output inflation]', Y2=[interest rate], Y=[Y2' Y1']', Y2*=[shadow rate]
    % Equations (16) and (17) in Mavroeidis (2021)
    A11 = H4\H1;
    A12 = -H4\H2.*(1-lambdas_p);
    A12s = -H4\H2.*lambdas_p;
    BB = H4\H3*J;
    B11 = [BB(:,1:2) BB(:,3).*(1-lambdas_p)];
    B11s = BB(:,3).*lambdas_p;
    A22s = 1;
    A22 = 0;
    A21 = -(1-rhoi_p).*[phiy_p phipi_p]./(1+iU_ss);
    B21 = [0 0 rhoi_p*(1-lambdas_p)];
    B21s = rhoi_p*lambdas_p;
    b = -iU_ss/(1+iU_ss);
    
    save var_matrices_common A11 H4 H2 BB A22s A22 A21 rhoi_p b iU_ss piU_ss yU_ss RbL
    
    if lambdas_p == 0
        save var_matrices_lambdas0 lambdas_p A11 A12 A12s B11 B11s A22s A22 A21 B21 B21s b iU_ss piU_ss yU_ss
    elseif lambdas_p == 0.5
        save var_matrices_lambdas05 lambdas_p A11 A12 A12s B11 B11s A22s A22 A21 B21 B21s b iU_ss piU_ss yU_ss
    elseif lambdas_p == 1
        save var_matrices_lambdas1 lambdas_p A11 A12 A12s B11 B11s A22s A22 A21 B21 B21s b iU_ss piU_ss yU_ss
    end    
end