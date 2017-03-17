function result = poem_1D_jhp(param,ENVR,t,result)

%
% Init
%
if (nargin==3)
    % Default initial conditions:
    R0 = 0.1*param.K;
    B0 = param.B0;
    y = [R0 B0];
else
    % Use "results" argument for initial conditions:
    y = result.y(end,:);
end
 
%
% Run:
%
[y, f, mortpred, Eavail] = poem_deriv_1D_jhp(y,param,ENVR);
%
% Construct output:
%
%Diff eq result
result.y(t,:) = y;
%Resources
result.R(t,:) = y(:,param.ixR);
%Fish biomass
result.B(t,:) = y(:,param.ixFish);
%Time
result.t(t,:) = t;

%Fishing yield
result.Y(t,:) = result.B(t,:) .* param.F;
%Consumption
%[result.f(t,:), result.mortpred(t,:), result.Eavail(t,:)] = calcEncounter(y, param);
result.f(t,:) = f; 
result.mortpred(t,:) = mortpred;
result.Eavail(t,:) = Eavail;
%Mortality
result.mort(t,:) = result.mortpred(t,param.ixFish) + param.mort0 + param.F;
%Growth
[result.v(t,:),  result.nu(t,:)] = calcNu(result.Eavail(t,:), result.mort(t,:), param);

