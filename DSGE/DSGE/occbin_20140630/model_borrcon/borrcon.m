%
% Status : main Dynare file 
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

tic;
global M_ oo_ options_ ys0_ ex0_ estimation_info
options_ = [];
M_.fname = 'borrcon';
%
% Some global variables initialization
%
global_initialization;
diary off;
M_.exo_names = 'u';
M_.exo_names_tex = 'u';
M_.endo_names = 'b';
M_.endo_names_tex = 'b';
M_.endo_names = char(M_.endo_names, 'c');
M_.endo_names_tex = char(M_.endo_names_tex, 'c');
M_.endo_names = char(M_.endo_names, 'ec');
M_.endo_names_tex = char(M_.endo_names_tex, 'ec');
M_.endo_names = char(M_.endo_names, 'lb');
M_.endo_names_tex = char(M_.endo_names_tex, 'lb');
M_.endo_names = char(M_.endo_names, 'y');
M_.endo_names_tex = char(M_.endo_names_tex, 'y');
M_.param_names = 'RHO';
M_.param_names_tex = 'RHO';
M_.param_names = char(M_.param_names, 'BETA');
M_.param_names_tex = char(M_.param_names_tex, 'BETA');
M_.param_names = char(M_.param_names, 'M');
M_.param_names_tex = char(M_.param_names_tex, 'M');
M_.param_names = char(M_.param_names, 'R');
M_.param_names_tex = char(M_.param_names_tex, 'R');
M_.param_names = char(M_.param_names, 'SIGMA');
M_.param_names_tex = char(M_.param_names_tex, 'SIGMA');
M_.param_names = char(M_.param_names, 'GAMMAC');
M_.param_names_tex = char(M_.param_names_tex, 'GAMMAC');
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 5;
M_.param_nbr = 6;
M_.orig_endo_nbr = 5;
M_.aux_vars = [];
M_.Sigma_e = zeros(1, 1);
M_.H = 0;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
erase_compiled_function('borrcon_static');
erase_compiled_function('borrcon_dynamic');
M_.lead_lag_incidence = [
 1 3 0;
 0 4 8;
 0 5 0;
 0 6 0;
 2 7 0;]';
M_.nstatic = 2;
M_.nfwrd   = 1;
M_.npred   = 2;
M_.nboth   = 0;
M_.equations_tags = {
};
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(5, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(6, 1);
M_.NNZDerivatives = zeros(3, 1);
M_.NNZDerivatives(1) = 14;
M_.NNZDerivatives(2) = -1;
M_.NNZDerivatives(3) = -1;
%
% SHOCKS instructions
%
make_ex_;
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (M_.params(5))^2;
M_.sigma_e_is_diagonal = 1;
steady;
options_.irf = 0;
options_.nocorr = 1;
options_.nomoments = 1;
options_.order = 1;
var_list_=[];
info = stoch_simul(var_list_);
save('borrcon_results.mat', 'oo_', 'M_', 'options_');


disp(['Total computing time : ' dynsec2hms(toc) ]);
