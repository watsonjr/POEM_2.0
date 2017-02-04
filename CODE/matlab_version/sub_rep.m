%%% BIOMASS MADE FROM REPRODUCTION
function [nu, rep, egg] = sub_rep(nu,K,S,egg)
%nu: energy for growth or spawning
%K: proportion allocated to growth
%S: fraction of pop spawning at that time
%egg: energy stored for later repro

global NX

% NOTE: Still never going to accumulate biomass as muscle tissue
% If it is spawning season, it gets spawned
% If it is not spawning season, it gets stored as repro energy "egg"

    if K<1.0
        rho = zeros(NX,1);
        id = (nu > 0.0);
        rho(id) = (1.0-K) .* nu(id);  %energy available for from eating
        
        id2 = (S>0.0);
        %rep = eggs from energy now + eggs from stored energy
        rep(id2) = S(id2).*(rho(id2)+egg(id2));         %fraction of pop reproducing now
        egg(id2) = (1.0-S(id2)).*(rho(id2)+egg(id2));   %rest gets stored for later
        
        id3 = (S<=0.0);
        rep(id3) = 0.0;
        egg(id3) = egg(id3) + rho(id3);
        
        id4 = (nu > 0.0);
            %nu is now split into used for repro (nu) and stored (egg)
        nu(id4) = rep(id4);
    else
        rep = zeros(NX,1);
        egg = zeros(NX,1);
    end

end
