function [ D ] = Func_q( S )
%===========================================%
%   3D shape recovery based 2D tracks
%   This is to convert S(3N * F) to D(3F)
%===========================================%
[N,F] = size(S);
N = N / 3;

D = zeros(3*F,N);
for i = 1 : F
    for j = 1 : N
        x = S(j,i);
        y = S(N+j,i);
        z = S(2*N+j,i);
        D(3*(i-1)+1:3*i,j) = [x y z]';
    end
end


end

