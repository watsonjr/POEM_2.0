%
% Interactions
%
J = 0.5;   % Juvenile foraging reduction
%D = 0.5;   % Demersal foraging on pelagics
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