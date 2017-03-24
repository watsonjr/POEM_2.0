function [dydt, f, mortpred, Eavail] = poem_deriv_1D_jhp2_1B(y,param,ENVR)
% Non-neg
y(y<0) = 0;

% Resources
R = y(param.ixR);
R(1) = ENVR.Zm;
R(2) = ENVR.Zl;
det = ENVR.det;
%Logistic growth
dRdt(3) = param.be*det .* (1 - R(3)./param.K(3));
%Chemostat
%dRdt(3) = param.be*det .*(param.K(3)-R(3));
R(3) = R(3) + dRdt(3);
y(param.ixR) = R;
y(y<0) = 0;

% Fish biomass
B = y(param.ixFish);

% Feeding 
%Pelagic-demersal coupling
if (ENVR.H > param.Ddep)
    D = 0;
    set_thetaD_1B
end
%Encounter, Consumption, & Predation mortality
[f, mortpred, Eavail] = calcEncounter_jhp(y,param,ENVR);

% Fish:
ixFish = param.ixFish;

% Total mortality:
mort = mortpred(ixFish) + param.mort0 + param.F;

% Flux out of the size group:
v = (Eavail(param.ixFish) ./ param.wc(param.ixFish));
vplus = max(0,v);
gamma = (param.kappa .* vplus - mort) ./ ...
    (1 - param.z(param.ixFish) .^ (1 - mort ./ (param.kappa .* vplus)));
Fout = gamma .* B;

% Reproduction
Repro = (1-param.kappa).*vplus.*B;

% Flux into the size group:
for i = 1:param.nSpecies
    ix = (param.ix1(i):param.ix2(i)) - length(R);
    ixPrev = [ix(end) ix(1:(end-1))];
    Fin(ix) = Fout(ixPrev);
    % Reproduction = RE * (would-be growth + repro)
    Fin(ix(1)) = param.eRepro(i) * (Fin(ix(1)) + Repro(ix(end)));
end

% Mass balance
%Fish
dBdt = Fin - Fout + (v - mort).*B - Repro;
B = B + dBdt;
%Resource
dRdt = -mortpred(param.ixR).*R;
R = R + dRdt;

dydt = [R, B];

