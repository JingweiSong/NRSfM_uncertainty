
clear all;close all;

%   ===========Parameter control============ %
i_data = 4;     %   1-5 Data set ID
std_W = 0.07;
batch_iter = 4; %   Number of batches in experiment
%   ===========Parameter control============ %

var_W = std_W^2;




strings = { 'drink', 'pickup', 'stretch','yoga','dance'};


%   Load data
load([ './data/' strings{i_data} '.mat' ]);
S0 = S;
[T,n] = size(W); T = T / 2;          % number of frames and points
C = size(list,1);
M = zeros(C,n);
for i = 1 : size(list,1)
    M(i,list(i,1)) = 1;
    M(i,list(i,2)) = -1;
end


% Estimate camera matrices using PTA's Euclidean upgrade method:
Rs_full = sparse(2*T, 3*T);
for f = 1:T
    f2 = 2*f-[1 0]; f3 = 3*f-[2 1 0];
    Rs_full(f2,f3) = Rs(f2,:);
end
disp('The Monte Carlo test requires heavy computation, please be patient.');
disp('Consider reducing the number of iteration: "batch_iter".');
parfor i = 1 : batch_iter
    %   Add Gaussian noise
    Waf = W + std_W*randn(size(W));
    
    % %   My method
    S3_original{i} = shape_recovery_fpca_s_sharp(Waf,full(Rs_full),pinv(full(Rs_full))*Waf,10);
    [r{i},S3_original{i}] = getOptimalRank(S3_original{i},Waf,Rs_full,std_W);
    err3D_original{i} = S0 - S3_original{i};
    disp(['Finish model: ',num2str(i_data),' in iteration: ',num2str(i)]);
end
save(['Main_batchDataset_',num2str(100*std_W),'_',num2str(i_data),'.mat'],'S3_original','S0','K','batch_iter','r');





