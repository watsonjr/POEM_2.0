%%% ENERGY AVAILABLE FOR GROWTH NU for large demersal
function [nu, prod] = sub_nu_LD(B,met,c1,c2,c3,c4)

    global Lambda
    
    % convert to biomass specific ingestion
    %nu = ((I/B)*Lambda) - met
    %nu = 0.5*I
    % Already in biomass specific ingestion
    Ip = c1+c2+c3;
    Ib = c4;
    nu = (Ip*Lambda + Ib*0.6) - met;
    prod = nu .* B;
end
