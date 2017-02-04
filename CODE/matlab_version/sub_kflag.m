%%% SPAWNING FLAG
function [S, dd] = sub_kflag(S,dd,dthresh,dtot)
    global SP DAYS
    dur=59;
    id = (dd >= dthresh);
    %Change spawning flag
    if ((dtot+dur) <= DAYS)
        S(id,dtot:(dtot+dur)) = SP(id,:);
    else
        dleft = DAYS - dtot + 1;
        S(id,dtot:DAYS) = SP(id,1:dleft);
    end
    %Reset cumulative deg days
    dd(id) = 0.0;
    
end
