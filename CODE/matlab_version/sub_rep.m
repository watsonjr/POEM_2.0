%%% BIOMASS MADE FROM REPRODUCTION
function [nu, rep, egg] = sub_rep(nu,K,S,egg)
    %nu: energy for growth or spawning
    %K: proportion allocated to growth
    %S: fraction of pop spawning at that time
    %egg: energy stored for later repro
    
    % NOTE: Still never going to accumulate biomass as muscle tissue
    % If it is spawning season, it gets spawned
    % If it is not spawning season, it gets stored as repro energy "egg"
    
    if K<1.0
        if nu > 0.0
            rho = (1.0-K) .* nu;  %energy available for from eating
        else
            rho = 0.0;
        end
        if S>0.0
            %rep = eggs from energy now + eggs from stored energy
            rep = S.*(rho+egg);         %fraction of pop reproducing now
            egg = (1.0-S).*(rho+egg);   %rest gets stored for later
        else
            rep = 0.0;
            egg = egg + rho;
        end
    else
        rep = 0.0;
        egg = 0.0;
    end
    if nu > 0.0
        %nu is now split into used for repro (nu) and stored (egg)
        nu = rep;  
    end
end
