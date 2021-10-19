clc;clear all;close all;

%   ===========Parameter control============ %
std_W = 0.07;
i_data = 4;
%   ===========Parameter control============ %


var_W = std_W^2;
ind_x = [1 3 5 30 40];  %   For drawing Q-Q map
ind_y = [1 4 8 30 10];  %   For drawing Q-Q map
result_coverage = zeros(1,2);  % mean + std
result_swtest = zeros(1,1);  % mean + std

load(['Main_batchDataset_',num2str(100*std_W),'_',num2str(i_data),'.mat']);

element_xy = zeros(batch_iter,1);
%   Generate post statistics
S3_original_mean = zeros(size(S0));
for iter = 1 : batch_iter
    S3_original_mean = S3_original_mean + S3_original{1,iter};
end
S3_original_mean = S3_original_mean / batch_iter;

[T,n] = size(S0); T = T/3;
rate_original = zeros(3*n,T);
Var_Wxx = zeros(3*n,T);
for iter = 1 : batch_iter
    S_sharp_tmp = Func_q_inv(S3_original{1,iter});
    rank = r{1,iter};
    [X,Y,U,Sigma,V]= svdXY_rank(S_sharp_tmp,rank);

    tmpX = X*inv(X'*X)*X';
    tmpY = Y*inv(Y'*Y)*Y';
    
    S3_original_subspace = S3_original{1,iter};
    err3D_original_post = S3_original_subspace - S3_original_mean;
    err3D_original_post = Func_q_inv(err3D_original_post);
    for i = 1 : 3*n
        for j = 1 : T
            Var_W_post_ij = 1.5*var_W*(tmpX(i,i) + tmpY(j,j));
            %                 Var_Wxx(i,j) = sqrt(Var_W_post_ij);
            if(abs(err3D_original_post(i,j))<1.96*sqrt(Var_W_post_ij))
                rate_original(i,j) = rate_original(i,j) + 1;
            end
            if(i == ind_x(i_data)&& j == ind_y(i_data))
                element_xy(iter) = S3_original{1,iter}(i,j);
            end
        end
    end
    
end
rate_original = rate_original / batch_iter;

mean_rate_original = mean(mean(rate_original));
std_rate_original = std(rate_original(:));
result_coverage(1,:) = [mean_rate_original std_rate_original];

figure(i_data);
qqplot(element_xy);grid on;
set(gcf,'color','w');
[H,result_swtest(1)] = swtest(element_xy);
disp(['--------------------- Dataset' num2str(i_data) '  ---------------------']);
disp('Ratio of error within the bound: "mean", "std"');
disp(['                   ' num2str(result_coverage(1)) '   ' num2str(result_coverage(2))]);
disp('SW test result:');
disp(num2str(1));





function [X,Y,U,Sigma,V]= svdXY_threshold( S ,threshold_rank)

[U,Sigma,V] = svd(S);

%   Shrink U Sigma V by rank
tmp = diag(Sigma);
rank = sum(tmp > threshold_rank);
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

function [X,Y,U,Sigma,V]= svdXY_rank( S ,rank)

[U,Sigma,V] = svd(S);

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

function [S_subspace]= project( S ,rank)

[U,Sigma,V] = svd(S);

%   Shrink U Sigma V by rank
tmp = diag(Sigma);
U = U(:,1:rank);
Sigma = diag(tmp(1:rank));
V = V(:,1:rank);


S_subspace = U*Sigma*V';

end