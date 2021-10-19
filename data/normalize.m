clc;clear all;close all;
filename = 'yoga';
load(filename);

scale = max(max(abs(S)));
S = S/scale;
W = W/scale;

save(filename);
disp(['Finish scaling: ' filename]);