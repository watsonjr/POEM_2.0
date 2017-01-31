%%% Offline coupling
function [out_1, out_2, bf] = sub_offline_bent(enc_1,enc_2,bio_1,bio_2,B)
    con_1 = enc_1 .* bio_1;
    con_2 = enc_2 .* bio_2;
    if ((con_1 + con_2) > B)
        frac1 = con_1 ./ (con_1 + con_2);
        frac2 = con_2 ./ (con_1 + con_2);
        out_1 = (frac1 .* B) ./ bio_1;
        out_2 = (frac2 .* B) ./ bio_2;
    else
        out_1 = enc_1;
        out_2 = enc_2;
    end
    bf = (out_1.*bio_1 + out_2.*bio_2) ./ B; 
end
