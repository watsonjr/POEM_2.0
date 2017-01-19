%%% Forward Euler checks
function bio = sub_check(bio)
    ID = (bio < 0);
    bio(ID) = eps();
    %ID2 = find(isnan(bio))
    %bio(ID2) = eps()
end
