function [zdata zdataconcatenated zdatass oobase_ Mbase_ oostar_ Mstar_ ys_ ] = ...
    solve_one_constraint_wrong(modnam,modnamstar,...
    constraint, constraint_relax,...
    shockssequence,irfshock,nperiods,maxiter,init_orig)

global M_ oo_

errlist = [];

% solve model
eval(['dynare ',modnam,' noclearall'])
oobase_ = oo_;
Mbase_ = M_;
setss

eval(['dynare ',modnamstar,' noclearall '])
oostar_ = oo_;
Mstar_ = M_;





% [decrulea,decruleb]=get_pq(oobase_.dr);
% [decrulea_star,decruleb_star]=get_pq(oostar_.dr);


if isfield(Mbase_,'nfwrd')
    % the latest Dynare distributions have moved nstatic and nfwrd
    [decrulea,decruleb]=get_pq(oobase_.dr,Mbase_.nstatic,Mbase_.nfwrd);
    [decrulea_star,decruleb_star]=get_pq(oostar_.dr,Mstar_.nstatic,Mstar_.nfwrd);
else
    [decrulea,decruleb]=get_pq(oobase_.dr,oobase_.dr.nstatic,oobase_.dr.nfwrd);
    [decrulea_star,decruleb_star]=get_pq(oostar_.dr,oostar_.dr.nstatic,oostar_.dr.nfwrd);
end




endog_ = M_.endo_names;
exog_ =  M_.exo_names;

steady_state = oobase_.dr.ys;
steady_state_star = oostar_.dr.ys;
zdatass = steady_state;

% processes the constrain so as to uppend a suffix to each
% endogenous variables
constraint_difference = process_constraint(constraint,'_difference',Mbase_.endo_names,0);

constraint_relax_difference = process_constraint(constraint_relax,'_difference',Mbase_.endo_names,0);


nvars = Mbase_.endo_nbr;

nshocks = size(shockssequence,1);
if ~exist('init_orig')
    init = zeros(nvars,1);
    init_orig = init;
else
    init = init_orig;
end
zdataconcatenated = zeros(nperiods,nvars);
wishlist = endog_;
nwishes = size(wishlist,1);

violvecbool = zeros(nperiods,1);
for ishock = 1:nshocks
    
    changes=1;
    iter = 0;
    
    
    while (changes & iter<maxiter)
        iter = iter +1;
        
        
        [regime regimestart]=map_regime(violvecbool);
        
        [zdata] = mkdata_mixture(nperiods, decrulea, decruleb, ...
                                 decrulea_star, decruleb_star, ...
                                 endog_, exog_, endog_, irfshock, ...
                                 shockssequence(ishock,:),violvecbool, ...
                                 steady_state,steady_state_star,init);

        for i=1:nwishes
            eval([deblank(wishlist(i,:)),'_difference=zdata(:,i);']);
        end
        
        
        
        newviolvecbool = eval(constraint_difference);
        relaxconstraint = eval(constraint_relax_difference);
        
        
        
        % check if changes
        if (max(newviolvecbool-violvecbool>0)) | sum(relaxconstraint(find(violvecbool==1))>0)
            changes = 1;
        else
            changes = 0;
        end
        
        
        violvecbool = (violvecbool|newviolvecbool)-(relaxconstraint & violvecbool);
        
        
    end
    
    init = zdata(1,:);
    zdataconcatenated(ishock,:)=init;
    init= init';
    
    % reset violvecbool for next period -- consistent with expecting no
    % additional shocks
    violvecbool=[violvecbool(2:end);0];
    
end

if changes ==1
    display('Did not converge -- increase maxiter')
end

zdataconcatenated(ishock+1:end,:)=zdata(2:nperiods-ishock+1,:);

zdata = mkdata(max(nperiods,size(shockssequence,1)),decrulea,decruleb,endog_,exog_,wishlist,irfshock,shockssequence,init_orig);

