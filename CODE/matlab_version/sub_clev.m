%%% Consumption/Cmax
function clev = sub_clev(con,Tp,Tb,tdif,wgt)
    % calculates consumption rate of first element of enc
    %Cmax
    temp = (Tp.*tdif) + (Tb.*(1.0-tdif));
    % Specific ingestion rate from Hartvig et al (g/g/day) ref to 15C
    cmax = (exp(0.063.*(temp-15.0)) .* 60.0 .* wgt^(-0.25)) ./365.0;
    % Hartvig et al (g/g/day) ref to 10C
    %cmax = (exp(0.063*(temp-10.0)) .* 85.0 .* wgt^(-0.25)) ./365.0;
    % mizer (g/g/day) ref to 10C
    %cmax = (exp(0.063*(temp-10.0)) .* 40.0 .* wgt^(-0.25)) ./365.0;
    % J&C15 (g/g/day) ref to 10C
    %cmax = (exp(0.063*(temp-10.0)) .* 25.0 .* wgt^(-0.33)) ./365.0;
    % Specific ingestion rate from Kiorboe & Hirst (g/g/day)
    %cmax = (exp(0.063*(temp-15.0)) .* 10^(0.4) .* wgt^(-0.51)) .* 24e-3;
    %clev
    clev = con./cmax;
end
