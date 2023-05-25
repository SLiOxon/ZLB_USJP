%% set inputs for solution below
%  The program produces responses to the shocks selected below under
%  irfshock. Paths for all the endogenous variables in the model selected
%  are produced. VariableName_difference holds the piece-wise linear solution for
%  VariableName.  VariableName_uncdifference holds the linear solution for
%  VariableName.

global M_ oo_

% modname below chooses model
% directory. But simple param choices are made from paramfile_borrcon in current
% directory.
modnam = 'borrcon';
modnamstar = 'borrconnotbinding';

% Enter parameter values and save them in PARAM_EXTRA.mat file.
% This choices overwrite any choice that are made in the paramfile_borrcon file.
R = 1.05;
BETA = 0.945;
RHO   = 0.9;
SIGMA = 0.05;
M = 1;
GAMMAC = 1;
save PARAM_EXTRA R BETA RHO SIGMA GAMMAC


% see notes 1 and 2 in README. Note that 
% lb_ss is the steady-state value of the multiplier lb (see borrcon_steadystate.m).
%    lb_ss is calculated automatically as additional auxiliary parameter
%    around line 49 of solve_one_constraint
constraint = 'lb<-lb_ss';
constraint_relax ='b>M*y';

% Pick innovation for IRFs
irfshock =char('u');      % label for innovation for IRFs

shockssequence(10,1) = 0.03;
shockssequence(30,1) = -0.03;

nperiods = 60;

    
[zdatalinear zdatapiecewise zdatass oobase_ Mbase_  ] = ...
  solve_one_constraint(modnam,modnamstar,...
  constraint, constraint_relax,...
  shockssequence,irfshock,nperiods);


% unpack the IRFs
for i=1:Mbase_.endo_nbr
  eval([deblank(Mbase_.endo_names(i,:)),'_l=zdatalinear(:,i);']);
  eval([deblank(Mbase_.endo_names(i,:)),'_p=zdatapiecewise(:,i);']);
  eval([deblank(Mbase_.endo_names(i,:)),'_ss=zdatass(i);']);
end
   

titlelist = char('c (consumption)','b (borrowing)','y (income)','lb (multiplier)');
percent = 'Percent';
level = 'Level';
ylabels = char(percent,percent,percent,level);
figtitle = 'Simulated variables';
legendlist = cellstr(char('Piecewise Linear','Linear'));

line1=100*[c_p/c_ss,b_p/b_ss,y_p/y_ss,(lb_p+lb_ss)/100];
line2=100*[c_l/c_ss,b_l/b_ss,y_l/y_ss,(lb_l+lb_ss)/100 ];

makechart(titlelist,legendlist,figtitle,ylabels,line1,line2);

 
