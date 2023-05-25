%
% Status : main Dynare file 
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

tic;
global M_ oo_ options_ ys0_ ex0_ estimation_info
options_ = [];
M_.fname = 'model_simple_lb';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('model_simple_lb.log');
M_.exo_names = 'u';
M_.exo_names_tex = 'u';
M_.exo_names_long = 'u';
M_.endo_names = 'q';
M_.endo_names_tex = 'q';
M_.endo_names_long = 'q';
M_.endo_names = char(M_.endo_names, 'r');
M_.endo_names_tex = char(M_.endo_names_tex, 'r');
M_.endo_names_long = char(M_.endo_names_long, 'r');
M_.endo_names = char(M_.endo_names, 'rnot');
M_.endo_names_tex = char(M_.endo_names_tex, 'rnot');
M_.endo_names_long = char(M_.endo_names_long, 'rnot');
M_.endo_names = char(M_.endo_names, 'pu');
M_.endo_names_tex = char(M_.endo_names_tex, 'pu');
M_.endo_names_long = char(M_.endo_names_long, 'pu');
M_.param_names = 'betap';
M_.param_names_tex = 'betap';
M_.param_names_long = 'betap';
M_.param_names = char(M_.param_names, 'rhop');
M_.param_names_tex = char(M_.param_names_tex, 'rhop');
M_.param_names_long = char(M_.param_names_long, 'rhop');
M_.param_names = char(M_.param_names, 'rho_pu');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_pu');
M_.param_names_long = char(M_.param_names_long, 'rho_pu');
M_.param_names = char(M_.param_names, 'sigmap');
M_.param_names_tex = char(M_.param_names_tex, 'sigmap');
M_.param_names_long = char(M_.param_names_long, 'sigmap');
M_.param_names = char(M_.param_names, 'phip');
M_.param_names_tex = char(M_.param_names_tex, 'phip');
M_.param_names_long = char(M_.param_names_long, 'phip');
M_.param_names = char(M_.param_names, 'r_lowerbar');
M_.param_names_tex = char(M_.param_names_tex, 'r\_lowerbar');
M_.param_names_long = char(M_.param_names_long, 'r_lowerbar');
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 4;
M_.param_nbr = 6;
M_.orig_endo_nbr = 4;
M_.aux_vars = [];
M_.Sigma_e = zeros(1, 1);
M_.Correlation_matrix = eye(1, 1);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
erase_compiled_function('model_simple_lb_static');
erase_compiled_function('model_simple_lb_dynamic');
M_.lead_lag_incidence = [
 1 3 7;
 0 4 0;
 0 5 0;
 2 6 0;]';
M_.nstatic = 2;
M_.nfwrd   = 0;
M_.npred   = 1;
M_.nboth   = 1;
M_.nsfwrd   = 1;
M_.nspred   = 2;
M_.ndynamic   = 2;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(4, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(6, 1);
M_.NNZDerivatives = zeros(3, 1);
M_.NNZDerivatives(1) = 11;
M_.NNZDerivatives(2) = -1;
M_.NNZDerivatives(3) = -1;
M_.params( 1 ) = 0.99;
betap = M_.params( 1 );
M_.params( 2 ) = 0.5;
rhop = M_.params( 2 );
M_.params( 4 ) = 5;
sigmap = M_.params( 4 );
M_.params( 5 ) = 0.2;
phip = M_.params( 5 );
M_.params( 6 ) = (-(1/M_.params(1)-1));
r_lowerbar = M_.params( 6 );
M_.params( 3 ) = 0.5;
rho_pu = M_.params( 3 );
%
% INITVAL instructions
%
options_.initval_file = 0;
oo_.steady_state( 4 ) = 0;
oo_.steady_state( 1 ) = (-M_.params(4))*M_.params(5)*M_.params(6)/(1-M_.params(1)*(1-M_.params(2))-M_.params(2));
oo_.steady_state( 3 ) = M_.params(5)*oo_.steady_state(1);
oo_.steady_state( 2 ) = M_.params(6);
if M_.exo_nbr > 0;
	oo_.exo_simul = [ones(M_.maximum_lag,1)*oo_.exo_steady_state'];
end;
if M_.exo_det_nbr > 0;
	oo_.exo_det_simul = [ones(M_.maximum_lag,1)*oo_.exo_det_steady_state'];
end;
steady;
%
% SHOCKS instructions
%
make_ex_;
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (1)^2;
M_.sigma_e_is_diagonal = 1;
options_.irf = 0;
options_.nocorr = 1;
options_.nomoments = 1;
options_.order = 1;
var_list_=[];
info = stoch_simul(var_list_);
save('model_simple_lb_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('model_simple_lb_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('model_simple_lb_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('model_simple_lb_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('model_simple_lb_results.mat', 'estimation_info', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
