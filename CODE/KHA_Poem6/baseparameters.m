

param.tEnd = 100;
%
% Resources:
%
param.ixR = [1 2 3];
param.wc(param.ixR) = [0.00001 0.001 0.1];
param.r = 10*[10 1 0.1];
param.K = 1*[1 1 1];
%
% Species:
%
param.nSpecies = 1;
% large pelagics
param.ix1(1) = 4;
param.ix2(1) = 6;

param.ixFish = param.ix1(1):param.ix2(end);

% Weight classes
param.beta = 500;
for i = 1:param.nSpecies
    param.nClasses(i) = param.ix2(i)-param.ix1(i)+1;
    ix = param.ix1(i):param.ix2(i);
    param.w(ix) = 0.001*param.beta.^(0:(param.ix2(i)-param.ix1(i)));
end

% Init conditions
param.B0 = [10 0 0];
%
% Interactions
%
param.theta = [0,   0,   0,   0, 0, 0; 
               0,   0,   0,   0, 0, 0;
               0,   0,   0,   0, 0, 0;
               1,   0,   0,   0, 0, 0;
               0,   1,   0,   1, 0, 0;
               0,   0,   1,   0, 1, 0];
%param.theta = [0,   0,   0; 
%               1,   0,   0;
%               0.25,   1,   0;
%               ];
           
%
% Physiology:
%
param.wc(param.ixFish) = param.w(param.ixFish)*sqrt(param.beta); % central sizes
param.wu(param.ixFish) = param.w(param.ixFish)*param.beta; % Upper sizes
param.fc = 0.2;
param.h = 20;
param.epsAssim = 0.7;
param.q = 0.8;
param.n = 0.75;
param.mort0 = [10 0.5 0.1]';%1*param.wc(param.ixFish).^(-0.25)';
param.kappa = [1 1 0.5];
param.F = 0*param.ixFish';
param.gamma = 70;
param.eRepro = [.001];

param.z = (param.w./param.wu)';
param.Cmax = param.h*param.wc.^param.n;
param.V = param.gamma*param.wc.^param.q;
