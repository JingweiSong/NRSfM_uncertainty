
clear all;close all;

% =========================  Parameter setting ===================%
std_W = 0.01;
num_subtraj = 4;
overlap_percent = 0.2;      %   Percentage of overlap of subtrajectories
ind_dataset = 4;            % 1-5 sparse points | 6-9 dense points
% =================================================================%

var_W = std_W^2;


strings = { 'drink', 'pickup', 'stretch','yoga','dance',...
    'SYNTHETIC_FACES\Seq3','SYNTHETIC_FACES\Seq4'};

%   Load data
switch (ind_dataset > 5)
    case 1  %   Dense
        [W,S,Rs,Rs_full] = loadData(strings{ind_dataset});
    case 0  %   Sparse
        load([ './data/' strings{ind_dataset} '.mat' ]);
end


S0 = S;
[T,n] = size(W); T = T / 2;          % number of frames and points



% Estimate camera matrices using PTA's Euclidean upgrade method:
Rs_full = sparse(2*T, 3*T);
for f = 1:T
    f2 = 2*f-[1 0]; f3 = 3*f-[2 1 0];
    Rs_full(f2,f3) = Rs(f2,:);
end

%   Add noise
Waf = W + std_W*randn(size(W));

%   Segmenting the trajectory into submaps
tmp1 = round(size(Waf,1)/(2*num_subtraj));     %   Frames within patch
tmp2 = round(overlap_percent*tmp1/2);   %   Half of the overlap
ind = zeros(num_subtraj,2);
for i = 1 : num_subtraj
    ind(i,:) = [1+(i-1)*tmp1-tmp2 i*tmp1+tmp2];
end
ind(1,1) = 1;
ind(num_subtraj,2) = size(Waf,1)/2;
tic
parfor i = 1 : num_subtraj
    % -----------------------------------------------------------------------------
    tol = 1e-6;
    maxind_dataset = 100;
    
    Waf_sub = Waf(2*ind(i,1)-1:2*ind(i,2),:);
    Rs_sub  = Rs_full(2*ind(i,1)-1:2*ind(i,2),3*ind(i,1)-2:3*ind(i,2));
    
    S3_sub{i} = shape_recovery_fpca_s_sharp(Waf_sub,full(Rs_sub),pinv(full(Rs_sub))*Waf_sub,10);
    [r{i},S3_sub{i}] = getOptimalRank(S3_sub{i},Waf_sub,Rs_sub,std_W);
    [Var_W_post{i}] = getVarNRSfM(S3_sub{i},r{i},var_W);
end
S3_original = zeros(size(W,1)*1.5,size(W,2));
S3_original(3*ind(1,1)-2:3*ind(1,2),:) = S3_sub{1};
for i = 2 : num_subtraj
    num_overlap = ind(i-1,2) - ind(i,1) + 1;
    W1 = (1./Var_W_post{i-1}(end-3*num_overlap+1:end,:)) ./ (1./Var_W_post{i-1}(end-3*num_overlap+1:end,:) + 1./Var_W_post{i}(1:3*num_overlap,:));
    W2 = (1./Var_W_post{i  }(1:3*num_overlap,:)        ) ./ (1./Var_W_post{i-1}(end-3*num_overlap+1:end,:) + 1./Var_W_post{i}(1:3*num_overlap,:));
    S3_original(3*ind(i,1)-2:3*ind(i-1,2),:) = S3_sub{i-1}(end-3*num_overlap+1:end,:) .* W1 + ...
        S3_sub{i  }(1:3*num_overlap,:) .* W2;
    S3_original(3*ind(i-1,2)+1:3*ind(i,2),:) = S3_sub{i  }(3*num_overlap+1:end,:);
end
toc

S0_tmp = [];
S3_original_tmp = [];
Rs_tmp = [];
for i = 2 : num_subtraj
    S0_tmp = [S0_tmp;S0(3*ind(i,1)-2:3*ind(i-1,2),:)];
    S3_original_tmp = [S3_original_tmp;S3_original(3*ind(i,1)-2:3*ind(i-1,2),:)];
    Rs_tmp = [Rs_tmp;Rs(2*ind(i,1)-1:2*ind(i-1,2),:)];
end
[err3D_original,S1,S2] = getAccuracy( S0_tmp,S3_original_tmp,  Rs_tmp,size(Rs_tmp,1)/2,n );
disp(['[Overlap]Result of weighted submap based NRSfM is: ' num2str(err3D_original) '   Number of submap: ' num2str(num_subtraj) '   Percentage of overlap: ' num2str(overlap_percent)]);
[err3D_original,S1,S2] = getAccuracy( S0,S3_original,  Rs_tmp,size(Rs,1)/2,n );
disp(['[Entire]Result of weighted submap based NRSfM is: ' num2str(err3D_original) '   Number of submap: ' num2str(num_subtraj) '   Percentage of overlap: ' num2str(overlap_percent)]);

S3_average = zeros(size(W,1)*1.5,size(W,2));
S3_average(3*ind(1,1)-2:3*ind(1,2),:) = S3_sub{1};
for i = 2 : num_subtraj
    num_overlap = ind(i-1,2) - ind(i,1) + 1;
    W1 = (1./Var_W_post{i-1}(end-3*num_overlap+1:end,:)) ./ (1./Var_W_post{i-1}(end-3*num_overlap+1:end,:) + 1./Var_W_post{i}(1:3*num_overlap,:));
    W2 = (1./Var_W_post{i  }(1:3*num_overlap,:)        ) ./ (1./Var_W_post{i-1}(end-3*num_overlap+1:end,:) + 1./Var_W_post{i}(1:3*num_overlap,:));
    W1(:) = 0.5;
    W2(:) = 0.5;
    S3_average(3*ind(i,1)-2:3*ind(i-1,2),:) = S3_sub{i-1}(end-3*num_overlap+1:end,:) .* W1 + ...
        S3_sub{i  }(1:3*num_overlap,:) .* W2;
    S3_average(3*ind(i-1,2)+1:3*ind(i,2),:) = S3_sub{i  }(3*num_overlap+1:end,:);
end


S0_tmp = [];
S3_average_tmp = [];
Rs_tmp = [];
for i = 2 : num_subtraj
    S0_tmp = [S0_tmp;S0(3*ind(i,1)-2:3*ind(i-1,2),:)];
    S3_average_tmp = [S3_average_tmp;S3_average(3*ind(i,1)-2:3*ind(i-1,2),:)];
    Rs_tmp = [Rs_tmp;Rs(2*ind(i,1)-1:2*ind(i-1,2),:)];
end
[err3D_average,S1,S2] = getAccuracy( S0_tmp,S3_average_tmp,  Rs_tmp,size(Rs_tmp,1)/2,n );
disp(['[Overlap]Result of averaging submap based NRSfM is: ' num2str(err3D_average) '   Number of submap: ' num2str(num_subtraj) '   Percentage of overlap: ' num2str(overlap_percent)]);
[err3D_average,S1,S2] = getAccuracy( S0,S3_average,  Rs,size(Rs,1)/2,n );
disp(['[Entire]Result of averaging submap based NRSfM is: ' num2str(err3D_average) '   Number of submap: ' num2str(num_subtraj) '   Percentage of overlap: ' num2str(overlap_percent)]);


function [W,S0,Rs,Rs_full] = loadData(filename)

if(filename=="SYNTHETIC_FACES\Seq3" || filename=="SYNTHETIC_FACES\Seq4")
    filepath = ['E:\dataset\NRSfM\nnrsfm_datasets\' filename '\MATLAB_original_files\matlab.mat'];
    load(filepath);
    S0 = S;
    Rs_full = R_GT;
    ind = 3:3:size(Rs_full,1);
    Rs_full(ind,:)=[];
    Rs = zeros(size(Rs_full,1),3);
    for i = 1 : size(Rs_full,1)/2
        Rs(2*i-1:2*i,:) = Rs_full(2*i-1:2*i,3*i-2:3*i);
    end
end
end

