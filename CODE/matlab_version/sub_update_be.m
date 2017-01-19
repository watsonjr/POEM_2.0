%%% Update biomass
function bio_out = sub_update_be(bio_in,con,bio)
    %bio_in = benthic biomass
    %con = biomass specific consumption rate by MD and LD
    %bio = biomass of MD and LD
    die = con.*bio;
    bio_out = bio_in - sum(die);
end
