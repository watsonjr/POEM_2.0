param.tEnd = 100;
%
% Resources:
%
param.ixR = [1 2 3 4 5];
param.wc(param.ixR) = [0.00001 0.001 0.00001 0.001 0.1];
param.r =    [10 2     5   1  0.2];
param.K = 10*[1  1    0.5 0.5 0.5];
%
% Species:
%
param.nSpecies = 3;
% Indices for small pelagics
param.ix1(1) = 6;
param.ix2(1) = 7;
% Indices for large pelagics:
param.ix1(2) = 8;
param.ix2(2) = 10;
% Indices for large demersals:
param.ix1(3) = 11;
param.ix2(3) = 13;

param.ixFish = param.ix1(1):param.ix2(end);  % Index for all fish

% Weight classes
param.beta = 500; % Pred-prey mass ratio
for i = 1:param.nSpecies
    param.nClasses(i) = param.ix2(i)-param.ix1(i)+1;
    ix = param.ix1(i):param.ix2(i);
    param.w(ix) = 0.001*param.beta.^(0:(param.ix2(i)-param.ix1(i)));
end
param.wc(param.ixFish) = param.w(param.ixFish)*sqrt(param.beta); % central sizes
param.wu(param.ixFish) = param.w(param.ixFish)*param.beta; % Upper sizes

% Initial conditions 
param.B0 = [   6.1773    0.1055    0.2001    0.0112    0.4827    0.1688    0.0045 0.7033];
%
% Interactions
%
J = 0.5;  % Juvenile foraging reduction
D = 0.5;  % Demersal foraging on pelagics
param.theta = ...% Z   Benthos    SPel   LPelag    Ldemersal   
              [0, 0,   0, 0, 0,   0, 0,  0, 0, 0,  0, 0, 0; 
               0, 0,   0, 0, 0,   0, 0,  0, 0, 0,  0, 0, 0;
               % Benthos
               0, 0,   0, 0, 0,   0, 0,  0, 0, 0,  0, 0, 0;
               0, 0,   0, 0, 0,   0, 0,  0, 0, 0,  0, 0, 0;
               0, 0,   0, 0, 0,   0, 0,  0, 0, 0,  0, 0, 0;
               % Small pelagics:
               1, 0,   0, 0, 0,   0, 0,  0, 0, 0,  0, 0, 0;
               0, 1,   0, 0, 0,   1, 0,  1, 0, 0,  1, 0, 0;
               % Large pelagics:
               1, 0,   0, 0, 0,   0, 0,  0, 0, 0,  0, 0, 0;
               0, .5,  0 ,0, 0,   J, 0,  J, 0, 0,  J, 0, 0;
               0, 0,   0, 0, 0,   0, 1,  0, 1, 0,  0, D, 0;
               % Large demersals:
               1, 0,   0, 0, 0,   0, 0,  0, 0, 0,  0, 0, 0;
               0, 0,   0.25,1,0,  0, 0,  0 ,0, 0,  0 ,0, 0;
               0, 0,   0, 0, 1,   0,D*J, 0, D*J,0, 0, 1, 0];
           
%
% Physiology:
%
param.fc = 0.2;
param.h = 20;
param.epsAssim = 0.7;
param.q = 0.8;
param.n = 0.75;
param.mort0 = [0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1]';%1*param.wc(param.ixFish).^(-0.25)';
param.F = 0*param.ixFish';  % Fishing mortality
param.gamma = 70;
param.eRepro = [.5 .01 .01];
param.kappa = [1 0.5 1 1 0.5 1 1 0.5];

param.z = (param.w./param.wu)';
param.Cmax = param.h*param.wc.^param.n;  % Max consumption rate
param.V = param.gamma*param.wc.^param.q; % Clearance rate
