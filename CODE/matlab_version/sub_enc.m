%%%  Encounter rates
function enc = sub_enc(Tp,Tb,wgt,prey,tpel,tprey,pref)
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
    % Specific clearance rates from Kiorboe & Hirst (m3/g/day)
    % divide by 100m to put in m2/g/day b/c zoop is integrated over 100m depth
%     A = (exp(0.063.*(temp-15.0)) .* 10^(3.24) .* wgt^(-0.24)) .* (24e-3/9); %./ 100.0;
    % Hartvig et al (m3/g/day) gamma = 0.8e4 ref to 10C
    %A = (exp(0.063*(temp-10.0)) .* 0.8e4 .* wgt^(0.8-1.0)) ./365.0;
    % mizer et al (m3/g/day) gamma = 2.9e3? ref to 10C
    A = (exp(0.063*(temp-10.0)) .* 2.9e3 .* wgt^(0.8-1.0)) ./365.0;
    % J&C15 et al (m3/g/day) gamma = 2.9e3? ref to 10C
    %A = (exp(0.063*(temp-10.0)) .* 1.37e4 .* wgt^(0.9-1.0)) ./365.0;
    %Encounter per predator, mult by biomass later
    enc = prey.*A.*tprey.*pref;
end
