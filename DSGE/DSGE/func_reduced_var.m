function [C, Cs, G, Ct, Cts, ct, Gt, iU_ss, piU_ss, yU_ss, b, RbL] = func_reduced_var(lambdas_p)

% Input: lambdas_p  -- the parameter of the effectiveness of UMP 
% Output: reduced from VAR matrices
       
load var_matrices_common A11 H4 H2 BB A22s A22 A21 rhoi_p b iU_ss piU_ss yU_ss RbL  % generated from run_zlb.m
A12 = -H4\H2.*(1-lambdas_p);
A12s = -H4\H2.*lambdas_p;
B11 = [BB(:,1:2) BB(:,3).*(1-lambdas_p)];
B11s = BB(:,3).*lambdas_p;
B21 = [0 0 rhoi_p*(1-lambdas_p)];
B21s = rhoi_p*lambdas_p;
Ab = [A11 A12+A12s; A21 A22+A22s];
As = [A11 A12s; A21 A22s];
B = [B11;B21];
Bs = [B11s; B21s];
G = inv(Ab);
C = G*B;
Cs = G*Bs;
Gt = inv(As);
Ct = Gt*B;
Cts = Gt*Bs;
ct = - Gt*[A12;A22].*b;
    
end

