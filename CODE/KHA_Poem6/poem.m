function result = poem(param, result)

%
% Init
%
if (nargin==1)
    % Default initial conditions:
    R0 = param.K;
    B0 = param.B0;
    y0 = [R0 B0];
else
    % Use end results for initial conditions:
    y0 = result.y(end,:);
end
    
%
% Run:
%
[t,y] = ode45(@poem_deriv, [0 param.tEnd], y0, odeset('NonNegative',1:length(y0)), param);

result.y = y;
result.R = y(:,param.ixR);
result.B = y(:,param.ixFish);
result.t = t;

result.Y = result.B .* (ones(length(t),1)*param.F');

[result.f, result.mortpred, result.Eavail] = calcEncounter(y(end,:)', param);
result.mort = result.mortpred(param.ixFish)' + param.mort0 + param.F;
[result.v  result.nu] = calcNu(result.Eavail, result.mort, param);

