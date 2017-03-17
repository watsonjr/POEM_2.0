%%% Get COBALT data
function ENVR = get_COBALT(COBALT,GRD,ID,DY)

    NX = length(ID);
    % Get data
    ENVR.Tp  = COBALT.Tp(ID,DY);
    ENVR.Tb  = COBALT.Tb(ID,DY);
    ENVR.Zm  = COBALT.Zm(ID,DY);
    ENVR.Zl  = COBALT.Zl(ID,DY);
    ENVR.det = COBALT.det(ID,DY);
    ENVR.dZm = COBALT.dZm(ID,DY);
    ENVR.dZl = COBALT.dZl(ID,DY);
    %ENVR.U   = COBALT.U(ID,DY);
    %ENVR.V   = COBALT.V(ID,DY);
    ENVR.fZm = zeros(NX,1);
    ENVR.fZl = zeros(NX,1);
    ENVR.fB  = zeros(NX,1);
    ENVR.H   = GRD.Z(ID);
    ENVR.A   = GRD.AREA(ID);

    % Make sure no negatives
    ENVR.det = sub_neg(ENVR.det);
    ENVR.Zm  = sub_neg(ENVR.Zm);
    ENVR.Zl  = sub_neg(ENVR.Zl);
    ENVR.dZm = sub_neg(ENVR.dZm);
    ENVR.dZl = sub_neg(ENVR.dZl);
end
