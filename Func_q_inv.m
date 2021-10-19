function [ S ] = Func_q_inv( D )
%===========================================%
%   3D shape recovery based 2D tracks
%   This is to convert D(N * 3F) to S(3N * F)
%===========================================%
[F,N] = size(D);
F = F / 3;

S = zeros(3*N,F);
for i = 1 : F
    for j = 1 : N
        S(j,i)      = D(3*(i-1)+1,j);
        S(N+j,i)    = D(3*(i-1)+2,j);
        S(2*N+j,i)  = D(3*i,j);
    end
end


end

