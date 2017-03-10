function [f, mort, Eavail] = calcEncounter(y, param)
%
% Consumption
%
for i = 1:length(y)
    Enc(i,:) = param.V(i) .* param.theta(i,:) .* y';
end

Encspecies = sum(Enc');
f = Encspecies ./ (param.Cmax + Encspecies);
f(isnan(f)) = 0;
Eavail = param.Cmax .* (param.epsAssim * f - param.fc);
%
% Mortality:
%
for i = 1:length(y)
    mort(i) = sum( Enc(:,i) .* (y ./ param.wc') .* (1-f)' )/(eps + y(i));
end
