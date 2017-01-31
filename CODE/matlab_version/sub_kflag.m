%%% SPAWNING FLAG
function [S, dd] = sub_kflag(S,dd,dthresh,dtot)
    global SP DAYS
    if (dd >= dthresh)
        dur=59;
        %Change spawning flag
        if ((dtot+dur) <= DAYS)
            S(:,dtot:(dtot+dur)) = SP;
        else
            dleft = DAYS - dtot + 1;
            S(:,dtot:DAYS) = SP(1:dleft);
        end
        %Reset cumulative deg days
        dd = 0.0;
    end
end
