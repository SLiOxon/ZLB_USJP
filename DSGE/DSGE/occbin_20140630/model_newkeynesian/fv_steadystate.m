
function [ys,check]=fv_steadystate(junk,ys);

global M_

paramfile_fv

nparams = size(M_.param_names,1);
for icount = 1:nparams
    eval(['M_.params(icount) = ',M_.param_names(icount,:),';'])
end


check=0;

pss = INFSS;
psss =  ((1-TETA*pss^(EPSILON-1))/(1-TETA))^(1/(1-EPSILON)) ;
vss = (1-TETA)*psss^(-EPSILON)/(1-TETA*pss^EPSILON);
rss = pss/BETA;
icsss = psss/(1-TETA*BETA*pss^(EPSILON-1))/(1-GYSS);

XSS = EPSILON/(EPSILON-1);
PSI = icsss*(1-TETA*BETA*pss^EPSILON)/(XSS*vss^PHI) ;

YSS = (icsss*(1-TETA*BETA*pss^EPSILON)/(XSS*PSI*vss^PHI))^(1/(1+PHI)) ;
yss = YSS ;
css = (1-GYSS)*yss;


atss=1;
utss=1;

ys = [ atss
css  
icsss
pss
psss
rss
utss 
vss
yss ] ;

save ys css icsss pss psss rss vss yss

