%%% Temp-dep natural mortality
function nmort = sub_nmort(Tp,Tb,tpel,wgt)
    %Tp: pelagic temp
    %Tb: bottom temp
    %tpel: frac pelagic time
    
    global MORT Nat_mrt NX
    
    if (MORT==0) % None
        nmort = 0.0;
    end
    if (MORT==1) % Constant
        nmort = Nat_mrt * ones(NX,1);
    end
    if (MORT==2) % Hartvig Temperature-dependent mortality
        temp = (Tp.*tpel) + (Tb.*(1.0-tpel));
        % Hartvig
        nmort = exp(0.063*(temp-10.0)) .* 0.84 .* wgt^(-0.25) /365.0;
    end
    if (MORT==3) % mizer Temperature-dependent mortality
        temp = (Tp.*tpel) + (Tb.*(1.0-tpel));
        % mizer
        nmort = exp(0.063*(temp-10.0)) .* 3.0 .* wgt^(-0.25) /365.0;
    end
    if (MORT==4) % Jennings & Collingridge Temperature-dependent mortality
        temp = (Tp.*tpel) + (Tb.*(1.0-tpel));
        % J&C
        nmort = exp(0.063*(temp-10.0)) .* 0.5 .* wgt^(-0.33) /365.0;
    end
    if (MORT==5) % Peterson & Wrob Temperature-dependent mortality
        temp = (Tp.*tpel) + (Tb.*(1.0-tpel));
        % Peterson & Wroblewski (daily & uses dry weight)
        nmort = exp(0.063*(temp-15.0)) * 5.26e-3 * (wgt/9.0)^(-0.25); 
    end
end
