%%% Offline coupling
function [out_1, out_2, zf] = sub_offline_zl(enc_1,enc_2,bio_1,bio_2,dZ)
    % ADD FLAG FOR COUNTING HOW MANY TIMES THIS HAPPENS
    % offline switch
    con_1 = enc_1 .* bio_1;
    con_2 = enc_2 .* bio_2;
    if ((con_1 + con_2) > dZ)
        frac1 = con_1 ./ (con_1 + con_2);
        frac2 = con_2 ./ (con_1 + con_2);
        out_1 = (frac1 .* dZ) ./ bio_1;
        out_2 = (frac2 .* dZ) ./ bio_2;
    else
        out_1 = enc_1;
        out_2 = enc_2;
    end
    zf = (out_1*bio_1 + out_2*bio_2) ./ dZ;
end
