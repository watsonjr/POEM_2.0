%%% ENERGY AVAILABLE FOR SOMATIC GROWTH
function gamma = sub_gamma(K,Z,nu,d,B,nmrt,Frate,selec)
    % d = predation loss
    % nmort = natural mortality rate
    % Frate = fishing mortality rate
    % selec = if harvested
    
    % convert predation mortality to biomass specific rate
    if (selec == 1)
        D = (d./B) + nmrt + Frate;
    else
        D = (d./B) + nmrt;
    end
    kap=K;
    gg = ((kap.*nu) - D) ./ (1 - (Z.^(1 - (D ./ (kap.*nu)))));
    gneg = (gg < 0);
    gamma = min(gg,nu);
    gamma(gneg) = 0.0;
end
