function [Var_W_post] = getVarNRSfM(S3_original,r,var_W)


[T,n] = size(S3_original); T = T/3;
S_sharp_tmp = Func_q_inv(S3_original);
rank = r;
[X,Y,U,Sigma,V]= svdXY_rank(S_sharp_tmp,rank);
% tmpX = X*inv(X'*X)*X';
% tmpY = Y*inv(Y'*Y)*Y';
tmpX = zeros(size(X,1),1);
tmpY = zeros(size(Y,1),1);
inv_XX = inv(X'*X);
for i = 1 : size(X,1)
    tmpX(i) = X(i,:)*inv_XX*X(i,:)';
end
inv_YY = inv(Y'*Y);
for i = 1 : size(Y,1)
    tmpY(i) = Y(i,:)*inv_YY*Y(i,:)';
end

Var_W_post = zeros(3*n,T);
for i = 1 : 3*n
    for j = 1 : T
%         Var_W_post(i,j) = 1.5*var_W*(tmpX(i,i) + tmpY(j,j));
        Var_W_post(i,j) = 1.5*var_W*(tmpX(i) + tmpY(j));
    end
end
Var_W_post = Func_q(Var_W_post);

end


function [X,Y,U,Sigma,V]= svdXY_rank( S ,rank)

[U,Sigma,V] = svd(S,'econ');

%   Shrink U Sigma V by rank
tmp = diag(Sigma);
U = U(:,1:rank);
Sigma = diag(tmp(1:rank));
V = V(:,1:rank);

[m,n] = size(Sigma);
Sigma_sqrt = sqrt(Sigma);
if(m >= n)
    X = U*Sigma_sqrt;
    Y = Sigma_sqrt(1:n,1:n)*V';
else
    X = U*Sigma_sqrt(1:m,1:m);
    Y = Sigma_sqrt*V';
end
Y = Y';
end
