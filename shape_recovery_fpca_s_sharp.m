%%

function S = shape_recovery_fpca_s_sharp(W,R,S0,K)

[F P] = size(W);
F = F/2;  %Frames

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R_Recover = zeros(2*F,3);

%obtain the rotation matrix 2Fx3
for f=1:F
    R_Recover(2*f-1:2*f,:) = R(2*f-1:2*f,3*f-2:3*f);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu_0 = 4;                  %Initial continuation parameter mu
eta_mu = 0.25;             %0.25
mu_th = 1e-10;             %threshold on mu
tau = 0.2;                 %gradient step size
espilon = 1e-8;            %Threshold

S_k = S0;                  %Initial shape estimation 3F X P
%------------------------------------------------------------------------
%Transform from 3F*P matrix S_k to F*3P matrix S_sharp_k [X Y Z]
S_sharp_k = [S_k(1:3:end,:) S_k(2:3:end,:) S_k(3:3:end,:)];
%--------------------------------------------------------------------------
outer_iteration = 1; %boolean variable to control the outer loop

for l=1:20
    if(~outer_iteration)
        break;
    end
%     disp('Outer loop count: out of 20 iterations')
    
    %decrease mu in each iteration
    mu = max(mu_th,mu_0*eta_mu^l); 
    
    if (mu == mu_th)
        outer_iteration = 0;
    end
    
    Energy_Residual = [];
    non_converged = 1;
    it = 1;
    
    %Inner loop with gradient descent and singular value shrinkage
    while(non_converged)
        %Step1: Compute gradient g(S^k_sharp)
        g_S = R'*(R*S_k - W);    %gradient for S_k (3FxP)
        
        %transform g_S to g_S_sharp (Fx3P)        
        g_S_sharp = [g_S(1:3:end,:) g_S(2:3:end,:) g_S(3:3:end,:)];
        
        %record the last S_sharp
        S_sharp_k_1 = S_sharp_k;
        
        %Gradient descent
        Y_k_sharp = S_sharp_k - tau*g_S_sharp;
        
        %Step2: Singular value shrinkage operation
        [U_Y D_Y V_Y] = svd(Y_k_sharp);
        
        for i = 1:min(size(D_Y,1),size(D_Y,2))
            D_Y(i,i) = D_Y(i,i) - tau*mu;
            
            if(D_Y(i,i) < 0)
                D_Y(i,i) = 0;
            end
        end
        
        %Update S_k_sharp
        S_sharp_k = U_Y*D_Y*V_Y';
        
        %Transform S_k_sharp F x 3P to S_k 3F x P        
        S_k(1:3:end,:) = S_sharp_k(:,1:P);
        S_k(2:3:end,:) = S_sharp_k(:,P+1:2*P);
        S_k(3:3:end,:) = S_sharp_k(:,2*P+1:3*P);
                
        %Evaluate the residual
        Energy_Residual(it) = norm(S_sharp_k - S_sharp_k_1,'fro')/max(1,norm(S_sharp_k_1,'fro'));
%         Energy_Residual(it)
        if(Energy_Residual(it) < espilon)
            non_converged = 0;
        end
                
        it = it + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project to K-shape basis 
[U D V] = svd(S_sharp_k);
D(K+1:end,K+1:end) = 0;
S_sharp = U*D*V';

S = zeros(3*F,P);
S(1:3:end,:) = S_sharp(:,1:P);
S(2:3:end,:) = S_sharp(:,P+1:2*P);
S(3:3:end,:) = S_sharp(:,2*P+1:3*P);
