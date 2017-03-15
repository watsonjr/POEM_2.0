function result = poem(param, result)

%
% Init
%
if (nargin==1)
    % Default initial conditions:
    R0 = 0.1*param.K;
    B0 = param.B0;
    y0 = [R0 B0];
else
    % Use "results" argument for initial conditions:
    y0 = result.y(end,:);
end
 
%
% Run:
%
[t,y] = ode45(@poem_deriv, [0 param.tEnd], y0, odeset('NonNegative',1:length(y0)), param);
%
% Construct output:
%
%Diff eq result
result.y = y;
%Resources
result.R = y(:,param.ixR);
%Fish biomass
result.B = y(:,param.ixFish);
%Time
result.t = t;

%Fishing yield
result.Y = result.B .* (ones(length(t),1)*param.F');
%Consumption
[result.f, result.mortpred, result.Eavail] = calcEncounter(y(end,:)', param);
%Mortality
result.mort = result.mortpred(param.ixFish)' + param.mort0 + param.F;
%Growth
[result.v  result.nu] = calcNu(result.Eavail, result.mort, param);

