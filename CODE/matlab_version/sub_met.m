%%% Metabolism
function met = sub_met(Tp,Tb,tdif,wgt)
    %Tp: pelagic temp
    %Tb: bottom temp
    %tdif: frac pelagic time
    %wgt: ind weight of size class
    %fcrit: feeding level to meet resting respiration rate
    %cmax: max consumption rate
    %U: swimming speed
    
    global fcrit
    
    temp = (Tp.*tdif) + (Tb.*(1.0-tdif));
    
    %Cmax
    % Specific ingestion rate from Hartvig et al (g/g/day) ref to 15C
%     cmax = (exp(0.063.*(temp-15.0)) .* 60.0 .* wgt^(-0.25)) ./365.0;
    
    % Specific ingestion rate from Hartvig et al (g/g/day) ref to 10C
    cmax = (exp(0.063*(temp-10.0)) .* 85.0 .* wgt^(-0.25)) ./365.0;
    % Specific ingestion rate from mizer (g/g/day) ref to 10C
    %cmax = (exp(0.063*(temp-10.0)) .* 40.0 .* wgt^(-0.25)) ./365.0;
    % Specific ingestion rate from J&C15 (g/g/day) ref to 10C
    %cmax = (exp(0.063*(temp-10.0)) .* 25.0 .* wgt^(-0.33)) ./365.0;
    
    %Metabolism
    bas = fcrit .* cmax;
    met = bas;
    
    
    % Hartvig et al (g/g/day) k=10 ref to 10C
    %met = (exp(0.063*(temp-10.0)) .* 10 .* wgt^(-0.25)) ./365.0;
    % mizer et al (g/g/day) k=4.8 ref to 10C
    %met = (exp(0.063*(temp-10.0)) .* 4.8 .* wgt^(-0.25)) ./365.0;
    % J&C15 et al (g/g/day) k=2.5 ref to 10C
    %met = (exp(0.063*(temp-10.0)) .* 2.5 .* wgt^(-0.25)) ./365.0;

end
