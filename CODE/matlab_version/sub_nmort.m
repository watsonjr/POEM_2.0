%%% Temp-dep natural mortality
function nmort = sub_nmort(Tp,Tb,tpel,wgt)
    %Tp: pelagic temp
    %Tb: bottom temp
    %tpel: frac pelagic time
    
    global MORT Nat_mrt M_l
    
    if (MORT==0) % None
        nmort = 0.0;
    end
    if (MORT==1) % Constant
        nmort = Nat_mrt;
    end
    if (MORT==2) % Temperature-dependent mortality
        temp = (Tp.*tpel) + (Tb.*(1.0-tpel));
        % Const with size
        %nmort = exp(0.063*(temp-15.0)) .* Nat_mrt;
        % Hartvig
        %nmort = exp(0.063*(temp-10.0)) .* 0.84 .* wgt^(-0.25) /365.0;
        % mizer
        %nmort = exp(0.063*(temp-10.0)) .* 3.0 .* wgt^(-0.25) /365.0;
        % J&C
        %nmort = exp(0.063*(temp-10.0)) .* 0.5 .* wgt^(-0.33) /365.0;
        % Intermediate
        %nmort = exp(0.063*(temp-10.0)) .* 1.5 .* wgt^(-0.33) /365.0;
    end
    if (MORT==3) % Large fishes only
        if (wgt == M_l)
            nmort = Nat_mrt;
        else
            nmort = 0.0;
        end
    end
    if (MORT==4) % Large fishes only w/ temp-dep
        if (wgt == M_l)
            temp = (Tp.*tpel) + (Tb.*(1.0-tpel));
            nmort = exp(0.063*(temp-15.0)) .* Nat_mrt;
        else
            nmort = 0.0;
        end
    end
end
