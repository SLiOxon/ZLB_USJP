%% set inputs for solution below
%  The program produces responses to the shocks selected below under
%  irfshock. Paths for all the endogenous variables in the model selected
%  are produced. VariableName_difference holds the piece-wise linear
%  solution for
%  VariableName.  VariableName_uncdifference holds the linear solution for
%  VariableName.

clear


global M_ oo_
global betap rhop sigmap phip r_lowerbar

% modname below chooses model
% directory. But simple param choices are made from paramfile in current
% directory.
modnam = 'model_simple';
modnamstar = 'model_simple_lb';


compute_euler_residuals = 0;



% express the occasionally binding constraint
% in linearized form
% one can use any combination of endogenous variables and parameters
% declared in the the dynare .mod files
% constraint1 defines the first constraint
% if constraint1 is true, solution switches to model2
% but if constraint_relax is true, solution reverts to model1


constraint = 'r<r_lowerbar';
constraint_relax ='rnot>r_lowerbar';

% Pick innovation for IRFs
irfshock =char('u');      % label for innovation for IRFs
% needs to be an exogenous variable in the
% dynare .mod files

maxiter = 10;


% Option=1: impulse responses
% Option=2: random simulation
% Option=3: comparison with nonlinear model when available

option=1;

%%%%%%%%%%%%%%%% Inputs stop here %%%%%%%%%%%%%%%%%%%%%
%%

if option==1
    nper=1;
    
    shockssequence = [0
        0.05*ones(nper,1)/nper
        zeros(7,1)
        -0.05*ones(nper,1)/nper
        zeros(7,1)
        ];         % scale factor for simulations
    nperiods = size(shockssequence,1);            %length of IRFs
    
    
end

if option==2
    nperiods = 25;
    randn('seed',3);
    shockssequence = 1*randn(nperiods,1)*0.03 ;
end


% Solve model, generate model IRFs
[zdata zdataconcatenated zdatass oobase_ Mbase_ oostar_ Mstar_] = ...
    solve_one_constraint_temp1(modnam,modnamstar,...
    constraint, constraint_relax,...
    shockssequence,irfshock,nperiods,maxiter);


% unpack the IRFs
for i=1:M_.endo_nbr
    eval([deblank(M_.endo_names(i,:)),'_uncdifference=zdata(:,i);']);
    eval([deblank(M_.endo_names(i,:)),'_difference=zdataconcatenated(:,i);']);
end


%% Modify to plot IRFs and decision rules

% modify to plot IRFs
lbss=0;

titlelist = char('q','r');
percent = 'Percent dev. from steady state';
level = 'Level';
ylabels = char(percent,level);

%i_difference = log(exp(k_difference+log(kss))-(1-DELTAK)*exp(klag_difference+log(kss)))-log(iss);
%i_uncdifference = log(exp(k_uncdifference+log(kss))-(1-DELTAK)*exp(klag_uncdifference+log(kss)))-log(iss);


figtitle = '';
line1=100*[q_difference,r_difference];
line2=100*[q_uncdifference,r_uncdifference];

legendlist = char('Piecewise linear','Linear','Wrong Solution');


% create line3 -- mashup of 

[zdata_wrong zdataconcatenated_wrong zdatass oobase_ Mbase_ oostar_ Mstar_] = ...
    solve_one_constraint_wrong(modnam,modnamstar,...
    constraint, constraint_relax,...
    shockssequence,irfshock,nperiods,maxiter);

% unpack the IRFs
for i=1:M_.endo_nbr
    eval([deblank(M_.endo_names(i,:)),'_uncwrong=zdata_wrong(:,i);']);
    eval([deblank(M_.endo_names(i,:)),'_wrong=zdataconcatenated_wrong(:,i);']);
end
line3=100*[q_wrong,r_wrong];

makechart9(titlelist,legendlist,figtitle,-1000,ylabels,line1,line2,line3);


