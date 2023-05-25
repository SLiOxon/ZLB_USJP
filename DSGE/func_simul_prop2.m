function Datan = func_simul_prop2(C,Cs,G,Ct,Cts,ct,Gt,En,T,iU_ss,yU_ss,b,RbL,lambdas_p)

% Output: Datan = [y, pi, i, is, RL] (output gap, inflation, short-rate, long-rate)

Datan = zeros(T,5);
X = zeros(3,1);       % initial condition for [y pi i]
Xs = 0;               % initial condition for iT
Zm = zeros(2,1);      % initial condition for [za zb]

for t=1:T   
    % Ys=[y, pi, i]';
    Ys = C*X + Cs*Xs + G*En(:,t);
    if Ys(3) < b
        Ys = Ct*X + Cts*Xs + ct + Gt*En(:,t);
    end
    X = [Ys(1);Ys(2);max(Ys(3),b)];
    Xs = Ys(3);
    Datan(t,3) = (X(3)+iU_ss)*400*(1+iU_ss);
    Datan(t,4) = (Xs+iU_ss)*400*(1+iU_ss);
    Datan(t,1) = X(1)/yU_ss*100;
    Datan(t,2) = X(2)*400; 
    
    % Long-term yield
    Z = RbL.Rho*Zm + En(1:2,t);
    ERbL = RbL.EC*[Ys(3) Z']';
    if Ys(3) >= b
        RbLU = (1-RbL.w)*Ys(3) + RbL.w*ERbL;
    else
        RbLU = (1-RbL.w)*(lambdas_p*Ys(3) + (1-lambdas_p)*b) + RbL.w*ERbL;
    end
    Datan(t,5) = ((RbLU + RbL.ss)*RbL.ss-1)*400;
    Zm = Z;
end

end

