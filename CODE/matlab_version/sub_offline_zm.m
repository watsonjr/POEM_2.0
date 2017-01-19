%%% Offline coupling
function [out_1, out_2, out_3, out_4, out_5, zf] = sub_offline_zm(enc_1,enc_2,enc_3,enc_4,enc_5,bio_1,bio_2,bio_3,bio_4,bio_5,dZ)
    % ADD FLAG FOR COUNTING HOW MANY TIMES THIS HAPPENS
    % offline switch
    con_1 = enc_1 .* bio_1;
    con_2 = enc_2 .* bio_2;
    con_3 = enc_3 .* bio_3;
    con_4 = enc_4 .* bio_4;
    con_5 = enc_5 .* bio_5;
    if ((con_1 + con_2 + con_3 + con_4 + con_5) > dZ)
        frac1 = con_1 ./ (con_1 + con_2 + con_3 + con_4 + con_5);
        frac2 = con_2 ./ (con_1 + con_2 + con_3 + con_4 + con_5);
        frac3 = con_3 ./ (con_1 + con_2 + con_3 + con_4 + con_5);
        frac4 = con_4 ./ (con_1 + con_2 + con_3 + con_4 + con_5);
        frac5 = con_5 ./ (con_1 + con_2 + con_3 + con_4 + con_5);
        out_1 = (frac1 .* dZ) ./ bio_1;
        out_2 = (frac2 .* dZ) ./ bio_2;
        out_3 = (frac3 .* dZ) ./ bio_3;
        out_4 = (frac4 .* dZ) ./ bio_4;
        out_5 = (frac5 .* dZ) ./ bio_5;
    else
        out_1 = enc_1;
        out_2 = enc_2;
        out_3 = enc_3;
        out_4 = enc_4;
        out_5 = enc_5;
    end
    zf = (out_1*bio_1 + out_2*bio_2 + out_3*bio_3 + out_4*bio_4 + out_5*bio_5) ./ dZ;
end
