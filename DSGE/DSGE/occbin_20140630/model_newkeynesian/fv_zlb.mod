
// Endogenous variables
var at c ics p ps r ut v y ;

varexo eps_a eps_u ;

parameters BETA TETA GYSS EPSILON PHI TAYLOR_P TAYLOR_Y TAYLOR_R INFSS RZLB RHO_A RHO_U;



//////////////////
model;

#XSS = EPSILON/(EPSILON-1);
#pss = INFSS;
#psss =  ((1-TETA*pss^(EPSILON-1))/(1-TETA))^(1/(1-EPSILON)) ;
#vss = (1-TETA)*psss^(-EPSILON)/(1-TETA*pss^EPSILON);
#rss = pss/BETA;
#icsss = psss/(1-TETA*BETA*pss^(EPSILON-1))/(1-GYSS);
#PSI = icsss*(1-TETA*BETA*pss^EPSILON)/(XSS*vss^PHI) ;
#YSS = (icsss*(1-TETA*BETA*pss^EPSILON)/(XSS*PSI*vss^PHI))^(1/(1+PHI)) ;

c^(-1) - BETA*ut*r*(c(+1)^(-1)/(p(+1)));

ics - ( XSS*PSI*v^PHI*((y/at)^(1+PHI)) + TETA*BETA*ut*((p(+1)^EPSILON*ics(+1))) );

ics/ps - ( y/c+TETA*BETA*ut*((p(+1)^(EPSILON-1)/ps(+1)*ics(+1))) );

r - RZLB ;

1 - (TETA*p^(EPSILON-1)+(1-TETA)*ps^(1-EPSILON));

v - (TETA*p^EPSILON*v(-1)+(1-TETA)*ps^(-EPSILON));

y - ( c+GYSS*y )  ;

log(ut) = RHO_U*log(ut(-1))+eps_u ;

log(at) = RHO_A*log(at(-1))+eps_a ;

end;



