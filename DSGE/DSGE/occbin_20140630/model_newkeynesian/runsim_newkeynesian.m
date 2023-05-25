%% set inputs for solution below
%  The program produces responses to the shocks selected below under
%  irfshock. Paths for all the endogenous variables in the model selected
%  are produced. VariableName_difference holds the piece-wise linear solution for
%  VariableName.  VariableName_uncdifference holds the linear solution for
%  VariableName.


clear


global M_ oo_

modnam = 'fv';
modnamstar = 'fv_zlb';



constraint = 'r<RZLB-INFSS/BETA';
constraint_relax ='r>RZLB-INFSS/BETA';

% Pick innovation for IRFs
irfshock =char('eps_a','eps_u');      % label for innovation for IRFs

maxiter = 10;
tol0 = 1e-8;






% Solve nonlinear model

SHOCKS = [ zeros(5,2)
   0  0.024
  zeros(20,2) ] ;


shockssequence = SHOCKS;
nperiods = 30;



% Solve model, generate model IRFs
[zdatalinear zdatapiecewise zdatass oobase_ Mbase_  ] = ...
  solve_one_constraint(modnam,modnamstar,...
  constraint, constraint_relax,...
  shockssequence,irfshock,nperiods);


% unpack the IRFs                          
for i=1:M_.endo_nbr
  eval([deblank(M_.endo_names(i,:)),'_u=zdatalinear(:,i);']);
  eval([deblank(M_.endo_names(i,:)),'_p=zdatapiecewise(:,i);']);
  eval([deblank(M_.endo_names(i,:)),'_ss=zdatass(i);']);
end


%% Modify to plot IRFs and decision rules
paramfile_fv

figure(1)

subplot(3,2,1)
plot(BETA*(ut_p+ut_ss),'k')
title('Discount Factor')
ylabel('Level')

subplot(3,2,2)
plot(400*(r_p+r_ss-1),'ok')
title('Interest Rate')
ylabel('Annualized Level, PPt')

subplot(3,2,3)
plot(100*y_p/y_ss,'k')
title('Output')
ylabel('% from s.s.')

subplot(3,2,4)
plot(v_p+v_ss,'k')
title('Price Dispersion')
ylabel('Level')

subplot(3,2,5)
plot(400*(p_p+p_ss-1),'k')
hold on
title('Inflation')
ylabel('Annualized Level, PPt')





