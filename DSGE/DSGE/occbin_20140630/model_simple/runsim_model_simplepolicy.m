%% Reproduces figure 1 in the paper

clear

global M_ oo_
global betap rhop sigmap phip r_lowerbar rho_pu;

% modname below chooses model
% directory. But simple param choices are made from paramfile in current
% directory.
modnam = 'model_simple';
modnamstar = 'model_simple_lb';




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


maxiter = 10;
nperiods = 50;



ishock=1;

for x=[-0.2:0.02:0.2]
    
        
    [zdata_wrong zdataconcatenated_wrong zdatass oobase_ Mbase_ oostar_ Mstar_] = ...
        solve_one_constraint_wrong(modnam,modnamstar,...
        constraint, constraint_relax,...
        x,irfshock,nperiods,maxiter);
    
    [zdata zdataconcatenated zdatass oobase_ Mbase_ oostar_ Mstar_] = ...
        solve_one_constraint_temp1(modnam,modnamstar,...
        constraint, constraint_relax,...
        x,irfshock,nperiods,maxiter);
    
    
    
    
    % unpack the IRFs
    for i=1:M_.endo_nbr
        eval([deblank(M_.endo_names(i,:)),'_w=zdataconcatenated_wrong(:,i);']);
        eval([deblank(M_.endo_names(i,:)),'_l=zdata(:,i);']);
        eval([deblank(M_.endo_names(i,:)),'_p=zdataconcatenated(:,i);']);
    end
    line3=100*[q_w,r_w];
    
    qpol_w(ishock)=q_w(1);
    rpol_w(ishock)=r_w(1);
    qpol_l(ishock)=q_l(1);
    rpol_l(ishock)=r_l(1);
    qpol_p(ishock)=q_p(1);
    rpol_p(ishock)=r_p(1);
    upol(ishock)=x;
    
    ishock=ishock+1;
    
end

r_ss = 1/betap-1;

        
close all
figure
subplot(2,3,1:3)
plot(upol,100*(rpol_p+r_ss),'k-','LineWidth',2); hold on
plot(upol,100*(rpol_l+r_ss),'r--','LineWidth',2); hold on
plot(upol,100*(rpol_w+r_ss),'b-.','LineWidth',2)
title('Interest Rate - Level, %')
legend('Piecewise Linear','Linear','Naive Piecewise Linear','Location','SouthEast')
ylim([-1.5 3.5])

subplot(2,3,4:5)
breakpoint=find(diff(qpol_w)<0);

plot(upol,100*qpol_p,'k-','LineWidth',2); hold on
plot(upol,100*qpol_l,'r--','LineWidth',2); hold on
plot(upol(1:breakpoint),100*qpol_w(1:breakpoint),'b-.','LineWidth',2); hold on
plot(upol(breakpoint+1:end),100*qpol_w(breakpoint+1:end),'b-.','LineWidth',2)
ylim([-30 70])
title('Asset Price - Percent Dev. from Steady State')
xlabel('Shock')

hold all
subplot(2,3,6)
plot(upol(1:end-1),(qpol_p(2:end)-qpol_p(1:end-1))./(upol(2:end)-upol(1:end-1)),'k-','LineWidth',2); hold on
plot(upol(1:end-1),(qpol_l(2:end)-qpol_l(1:end-1))./(upol(2:end)-upol(1:end-1)),'r--','LineWidth',2)
plot(upol(1:breakpoint-1),(qpol_w(2:breakpoint)-qpol_w(1:breakpoint-1))./(upol(2:breakpoint)-upol(1:breakpoint-1)),'b-.','LineWidth',2); hold on
plot(upol(breakpoint+1:end-1),(qpol_w(breakpoint+2:end)-qpol_w(breakpoint+1:end-1))./(upol(breakpoint+2:end)-upol(breakpoint+1:end-1)),'b-.','LineWidth',2)
title('Slope of Asset Price')

xlabel('Shock')

delete *static* *auxiliary* *dynamic* *.log *.asv *results*
