%%%  Encounter rates for demersal
function enc = sub_enc_dem(Tp,Tb,wgt,prey,tpel,tprey,pref)
    % Tp: pelagic temp
    % Tb: bottom temp
    % wgt: ind weight of size class
    % pred: pred biomass density,
    % prey: prey biomass density,
    % A: predator search rate,
    % tpel: time spent in pelagic,
    % tprey: time spent in area with that prey item.
    % pref: preference for prey item
    temp = (Tp.*tpel) + (Tb.*(1.0-tpel));
    % Specific clearance rates from Blanchard et al. 2008 (m3/g/day)
    A = (exp(0.063*(temp-15.0)) .* 64 .* wgt^(-0.25)) ./365.0;
    %Encounter per predator, mult by biomass later
    enc = prey*A*tprey*pref;
end
