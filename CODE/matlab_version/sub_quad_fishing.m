%%%% Fishing
function [bio, caught, fmort] = sub_quad_fishing(bio,F,selec)
    %bio = fish biomass
    %F = fishing rate per day
    %selec = fishery selectivity 
    %NOTE: selec only 1 or 0 now, but could update code (here & gamma calc) so it is a fraction
        
    % Quadratic fishing mortality
    caught = bio.^2 .* selec .* F;
    fmort = caught ./ bio;
    bio = bio - caught;
    

end
