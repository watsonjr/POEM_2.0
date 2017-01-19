%%% DEGREE DAYS
function dd = sub_degday(dd,Tp,Tb,tdif,Tref,S,dtot)
    %if (S==0.0) %Don't accumulate temp while spawning, DD represents recovery after
    %if (sum(S(1:dtot)) < dtot) %Only spawn once per year
    if (sum(S(1:dtot)) > 0.0) %Only spawn once per year
        dd = 0.0;
    else
        Tavg = (Tp.*tdif) + (Tb.*(1.0-tdif));
        dd = dd + max((Tavg-Tref),0.0);
    end
end
