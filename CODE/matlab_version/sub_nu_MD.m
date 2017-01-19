%%% ENERGY AVAILABLE FOR GROWTH NU for medium demersal
function [nu, prod] = sub_nu_MD(I,B,met)
    % convert to biomass specific ingestion
    %nu = ((I/B)*Lambda) - met
    %nu = 0.5*I
    % Already in biomass specific ingestion
    nu = (I*0.6) - met;
    prod = nu .* B;
end
