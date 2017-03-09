%
% run a single large species on three resources
%

baseparameters

param.F = .0*[0 1 10]';
param.tEnd = 500;
result = poem(param);
plotPoem(param, result)