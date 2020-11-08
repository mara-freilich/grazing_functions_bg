function [mld,x1,x2,tm,ml_min,ml_max] = mldmodel(t)

% function to compute the mixed layer depth used in Freilich et al (in
% review) Biogeosciences. This function is used to approximate the average
% mixed layer depth observed in BGC Argo floats by Mignot et al 2018 Nature
% Communications.

% Input:
% t : time in days
% Output:
% mld: mixed layer depth (meters)
% x1 and x2: functions of time (unitless)
% tm: day of year (days)
% ml_min: minimum mixed layer depth (meters)
% ml_max: maximum mixed layer depth (meters)

%% user adjustable parameters
ml_min = 25; ml_max = 600;
N = 80; 
M = 80;
sl = 25;
%% compute mixed layer depth
tm=mod(t+365-N,365);
x1=pi*(tm/M.*(2-tm/M)).*(tm<(M+0.0001));
x2=(pi+pi/(365-sl)^3*(tm-sl).^3).*(tm>M);
mld=(ml_min+ml_max)/2+(ml_max-ml_min)/2*cos(x1+x2);
