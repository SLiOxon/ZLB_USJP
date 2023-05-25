var q r rnot pu;    
 
varexo u;

parameters betap, rhop, rho_pu, sigmap, phip, r_lowerbar;

betap = 0.99;
rhop = 0.5;
sigmap = 5;

phip = 0.2;
r_lowerbar = -(1/betap-1);
rho_pu = 0.5;

model;

q = betap*(1-rhop)*q(1)+rhop*q(-1)-sigmap*r+pu;

pu = rho_pu*pu(-1)+u;

rnot = phip*q;

r = r_lowerbar;

end;

initval;
pu = 0;
q = -sigmap*phip*r_lowerbar/(1-betap*(1-rhop)-rhop);
rnot = phip*q;
r = r_lowerbar;
end;

steady;


shocks;
var u; stderr 1;
end;


stoch_simul(order=1,nocorr,nomoments,irf=0);
