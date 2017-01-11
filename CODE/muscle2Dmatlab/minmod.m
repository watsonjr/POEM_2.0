% MINMOD corrector for AdvTVD.m
% J Guiet
% 19-12-2016
%
function delta=minmod(a,b)

delta=0;

if (a*b > 0)
    if abs(a)<abs(b)
        delta=a;
    elseif abs(a)>abs(b)
        delta=b;
    end
end

end