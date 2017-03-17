function [f, mort, Eavail] = calcEncounter(y, param)
%
% Encounter, Consumption, & Predation mortality
%

% Single prey encounters
for i = 1:length(y)
    Enc(i,:) = param.V(i) .* param.theta(i,:) .* y;
end

% Total encounter
Encspecies = sum(Enc');

% Consumption
f = Encspecies ./ (param.Cmax + Encspecies);
f(isnan(f)) = 0;
%Con = param.Cmax .* f;

% Metabolism
%Eavail =    param.epsAssim*Con - param.fc*param.Cmax;
Eavail = param.Cmax .* (param.epsAssim * f - param.fc);

% Mortality:
%Predation mortality = (biom enc by all preds * feeding level) / (prey biom)
for i = 1:length(y)
    mort(i) = sum(Enc(:,i)' .* (y ./ param.wc) .* (1-f)) / (eps + y(i));
end

%%May need to add code to prevent eating more than HP losses
conM = sum((Enc(:,1)' .* (y ./ param.wc) .* (1-f)));
if (conM > ENVR.dZm)
    frac = (Enc(:,1)' .* (y ./ param.wc) .* (1-f)) / conM;
    conMz = frac * ENVR.dZm; 
    %Where do I put this in and recalc energy avail?
    encMz = conMz ./ (y ./ param.wc);
    Enc(:,1) = encMz;
    % Total encounter
    Encspecies = sum(Enc');
    % Consumption
    f = Encspecies ./ (param.Cmax + Encspecies);
    f(isnan(f)) = 0;
    % Metabolism
    Eavail = param.Cmax .* (param.epsAssim * f - param.fc);
    % Mortality
    mort(1) = sum(Enc(:,1)' .* (y ./ param.wc) .* (1-f)) / (eps + y(1));
end
conL = sum((Enc(:,2)' .* (y ./ param.wc) .* (1-f)));
if (conL > ENVR.dZl)
    frac = (Enc(:,2)' .* (y ./ param.wc) .* (1-f)) / conL;
    conLz = frac * ENVR.dZl; 
    %Where do I put this in and recalc energy avail?
    encLz = conLz ./ (y ./ param.wc);
    Enc(:,1) = encLz;
    % Total encounter
    Encspecies = sum(Enc');
    % Consumption
    f = Encspecies ./ (param.Cmax + Encspecies);
    f(isnan(f)) = 0;
    % Metabolism
    Eavail = param.Cmax .* (param.epsAssim * f - param.fc);
    % Mortality
    mort(2) = sum(Enc(:,2)' .* (y ./ param.wc) .* (1-f)) / (eps + y(2));
end
