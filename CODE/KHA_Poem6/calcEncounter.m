function [f, mort, Eavail] = calcEncounter(y, param)
%
% Consumption
%
%Single prey encounters
for i = 1:length(y)
    Enc(i,:) = param.V(i) .* param.theta(i,:) .* y';
end
%Total encounter
Encspecies = sum(Enc');
%Consumption
f = Encspecies ./ (param.Cmax + Encspecies);
f(isnan(f)) = 0;
%Metabolism
Eavail = param.Cmax .* (param.epsAssim * f - param.fc);
%
% Mortality:
%
%Predation mortality = (biom enc by all preds * feeding level) / (prey biom)
for i = 1:length(y)
    mort(i) = sum(Enc(:,i) .* (y ./ param.wc') .* (1-f)') / (eps + y(i));
end
