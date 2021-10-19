function [rank,S_new] = getOptimalRank(S,Waf,Rs_full,std_W)
%GETOPTIMALRANK Summary of this function goes here
%   Estimate the optimal rank for post uncertainty propagation

D=Func_q_inv(S);
[U,Sigma,V] = svd(D,'econ');
[N,F] = size(D);
N = N/3;
r = 1;
R = full(Rs_full);
while(r < max(F)-1)
    Sigma_tmp = Sigma;
    Sigma_tmp(r+1:end,r+1:end)=0;
    D=U*Sigma_tmp*V';
    S_new=Func_q(D);
    Z1 = abs(Waf - R*S_new);
    Z1 = reshape(Z1,1,2*F*N);
    rate =  sum(Z1 < 1.959*std_W) / size(Z1,2);
    if(rate > 0.95)
        break;
    end
    r = r + 1;
end
rank = r;
end

