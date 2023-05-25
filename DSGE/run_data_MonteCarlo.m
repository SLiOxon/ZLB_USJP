% Testing the effectiveness of unconventional monetary policy in Japan and the United States
% Daisuke Ikeda, Shangshang Li, Sophocles Mavroeidis, and Francesco Zanetti
% American Economic Journal: Macroeconomics

% This code generates data the theoretical model using Monte Carlo (MC) simulations
% Data to be generated and saved
% Data_MC_lambdas**.mat for ** = 70, 75, 80, 85, 90, 95, and 99
% where **/100 = lambdas_p: the effectiveness of UMP
% The mat file includes
%    Data:  Time * [output gap, inflation, interest rate, shadow rate, long rate] * N (N = number of samples)
%    ELB :  N * 1 -- the frequency of ELB periods for each sample   

clear
close all

%% System of equations for a different value of lambdas
Lambdas100 = [70:5:95 99];
Lambdas_p = Lambdas100./100;   % effectiveness of UMP
J = length(Lambdas100);
rootname = 'var_reduced_lambdas';  
extension = '.mat';
for j = 1:J
    % Variables: Ys = [Y1' Y2s']' = [y pi is]'; X = Y(t-1); Y = [Y1' Y2']' = [y pi i]'; Xs = Y2s(t-1) = is;
    %            eps = [supply, demand, monetary]';
    % Model: Ys = C*X + Cs*Xs + G*eps if D=0 (non-ZLB)
    %           = Ct*X + Cts*Xs + ct + Gt*eps (ZLB)
    lambdas_p = Lambdas_p(j);
    [C, Cs, G, Ct, Cts, ct, Gt, iU_ss, piU_ss, yU_ss, b, RbL] = func_reduced_var(lambdas_p);
    filename = [rootname, num2str(Lambdas100(j)), extension]; % name of the file to be saved
    save(filename, 'lambdas_p', 'C', 'Cs', 'G', 'Ct', 'Cts', 'ct', 'Gt', 'RbL');
end

%% Simulate the model and generate data for each value of lambdas
N = 100;                      % number of simulations
T1 = (2019-1960+1)*4-3;       % sample period 1960q1 - 2019q1 (US)
T0 = 100;                     % data up to period T0 are discarded
T = T0 + T1;                  % total length of sample period
Data = zeros(T1,5,N);         % T1 * [y, pi, i, is, RbL] * N
ELB = zeros(N,1);             % frequency of ELB periods
ELBj = zeros(J,3);            % ELB frequency [mean min max]
Ysd = zeros(J,5);             % mean SD of [y, pi, i, is, RbL]
b1 = (b+iU_ss)*400*(1+iU_ss); % ELB level
rng(1)
SD = diag([0.5/100, 0.5/100, 0.25/400]);          % SD of shocks [supply, demand, monetary]
rootname2 = 'Data_MC_lambdas';
for j = 1:J
     filename = [rootname, num2str(Lambdas100(j)), extension];
     load(filename);
     lambdas_p = Lambdas_p(j);
     for n = 1:N
         rng(n)
         E0 = randn(3,T);  % random disturbances
         En = SD*E0;        % 3 * T matrix
         Datan = func_simul_prop2(C,Cs,G,Ct,Cts,ct,Gt,En,T,iU_ss,yU_ss,b,RbL,lambdas_p);
         Data(:,:,n) = Datan(T0+1:end,:);
         ELB(n) = length(find(Data(:,3,n)==b1))/T1;
     end
     filename = [rootname2, num2str(Lambdas100(j)), extension];
     save(filename, 'Data', 'ELB', 'lambdas_p');
     ELBj(j,:) = [mean(ELB), min(ELB), max(ELB)];
     Ysd(j,:) = [mean(std(Data(:,1,:))), mean(std(Data(:,2,:))), mean(std(Data(:,3,:))), mean(std(Data(:,4,:))), mean(std(Data(:,5,:)))];
end