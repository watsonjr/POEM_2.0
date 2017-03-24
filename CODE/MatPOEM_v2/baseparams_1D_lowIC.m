%Length of run (years)
param.tEnd = 100;
%Time step (years)
param.dt = 1/365;
%
% Environment:
%
param.Ddep = 200.0;
%
% Resources:
%
%Indices in feeding matrix
param.ixR = [1 2 3 4 5];
%Mean weights (only used for plotting)
param.wc(param.ixR) = [0.00001 0.001 0.00001 0.001 0.1];
%Transfer efficiency of detritus to benthic prey
param.be = 0.05;
%Growth rates (yr-1)
param.r = [10    3     10      3        1];
%Carrying capacities (gWW/m2)
param.K = [10   10     0.5    0.5      0.5];
%
% Species:
%
% Number of types
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
% Index for all fish
param.ixFish = param.ix1(1):param.ix2(end);  

% Weight classes
param.beta = 500; % Pred-prey mass ratio
for i = 1:param.nSpecies
    param.nClasses(i) = param.ix2(i)-param.ix1(i)+1;
    ix = param.ix1(i):param.ix2(i);
    param.w(ix) = 0.001*param.beta.^(0:(param.ix2(i)-param.ix1(i)));
end
%Central sizes
param.wc(param.ixFish) = param.w(param.ixFish)*sqrt(param.beta); 
%Upper sizes
param.wu(param.ixFish) = param.w(param.ixFish)*param.beta; 

% Initial conditions (from some run Ken did)
param.B0 = [1e-2,1e-3,1e-2,1e-3,1e-4,1e-2,1e-3,1e-4];
%
% Interactions
%
J = 0.5;   % Juvenile foraging reduction
D = 0.5;   % Demersal foraging on pelagics
Sm = 0.25; % 2 size classes down
param.theta = ...% Z   Benthos    F      Pelag     Demersal
              ...%MZ LZ   SB MB LB   SF MF  SP MP LP  SD MD LD
              ...%1  2    3  4  5    6  7   8  9  10  11 12 13 
                 [0, 0,   0, 0, 0,   0, 0,  0, 0, 0,  0, 0, 0; %1 MZ
                  0, 0,   0, 0, 0,   0, 0,  0, 0, 0,  0, 0, 0; %2 LZ
                  % Benthos
                  0, 0,   0, 0, 0,   0, 0,  0, 0, 0,  0, 0, 0; %3 SB
                  0, 0,   0, 0, 0,   0, 0,  0, 0, 0,  0, 0, 0; %4 MB
                  0, 0,   0, 0, 0,   0, 0,  0, 0, 0,  0, 0, 0; %5 LB
                  % Forage:
                  1, 0,   0, 0, 0,   0, 0,  0, 0, 0,  0, 0, 0; %6 SF
                 Sm, 1,   0, 0, 0,   1, 0,  1, 0, 0,  1, 0, 0; %7 MF
                  % Pelag:
                  1, 0,   0, 0, 0,   0, 0,  0, 0, 0,  0, 0, 0; %8 SP
               J*Sm, J,   0 ,0, 0,   J, 0,  J, 0, 0,  J, 0, 0; %9 MP
                  0, 0,   0, 0, 0,   0, 1,  0, 1, 0,  0, 0, 0; %10 LP
                  % Demers:
                  1, 0,   0, 0, 0,   0, 0,  0, 0, 0,  0, 0, 0; %11
                  0, 0,   Sm,1, 0,   0, 0,  0, 0, 0,  0 ,0, 0; %12
                  0, 0,   0, 0, 1,   0, D,  0, D, 0,  0, 1, 0];%13
param.D = D;
param.J = J;
param.Sm = Sm;

%
% Physiology:
%
%fcrit = fraction Cmax used for metab
param.fc = 0.2;
%Cmax coeff
param.h = 20;
%Assimilation effic
param.epsAssim = 0.7;
%Exponent on clearance/search
param.q = 0.8;
%Exponent on Cmax
param.n = 0.75;
%Background mortality (d^-1)
param.nmrt = [0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1];
param.mort0 = param.dt * [0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1]; %1*param.wc(param.ixFish).^(-0.25)';
%Fishing mortality (d^-1)
param.F = 0*param.ixFish; 
%Clearance/search coeff
param.gamma = 70;
%Repro effic
param.RE = 0.1;
param.eRepro = param.RE*[1 1 1];
%Fraction energy to growth
param.kappa = [1 0.5    1 1 0.5     1 1 0.5];
%Ratio initial:final size
param.z = (param.w./param.wu);
%Max consumption rate (g/g/d)
param.Cmax = param.dt * param.h * param.wc .^ param.n;  
%Clearance rate (m2/d)
param.V = param.dt * param.gamma * param.wc .^ param.q; 
