function [ err3D,S0R,S3R ] = getAccuracy( S0, S3, Rs,T,n )
%GETACCURACY Summary of this function goes here
%   Detailed explanation goes here

% Align original and recovered 3D shapes
[S0R,S3R] = alignStruct( S0, S3, Rs ,1); 
       
% Compute normalized, average 3D error
errS = zeros(T,n);                                  % 3D reconstruction errors
for t = 1:T, t3 = 3*t-[2 1 0];
    errS(t,:) = sqrt(sum( (S0R(t3,:)-S3R(t3,:)).^2 )); % 3D Euclidean distance
end
s = mean( std(S0R,1,2) );                           % "scale" (avg row std.dev.)
err3D = mean(errS(:)) / s;


end

