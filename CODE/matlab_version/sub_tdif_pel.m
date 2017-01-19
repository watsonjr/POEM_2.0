%%% Fraction of time spent in pelagic (for piscivore)
function tdif = sub_tdif_pel(Z,bio1,bio2,biod)
    % bio1, bio2: pelagic prey
    % biod: demersal prey
    biop = bio1+bio2;
    if Z < PI_be_cutoff
        tdif = biop ./ (biop+biod);
    else
        tdif = 1.0;
    end
end
