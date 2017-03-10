%
% run a single large species on three resources
%

baseparameters

param.F = 0.0*[0 0 0 0 0 0 0.1 1]';
param.tEnd = 200;
result = poem(param);
plotPoem(param, result)