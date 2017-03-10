%
% run a single large species on three resources
%

baseparameters

%Fishing rate
param.F = 0.0 * [0 1 10]';
%Length of run (yrs)
param.tEnd = 500;
%Result
result = poem(param);
%Plot
plotPoem(param, result)