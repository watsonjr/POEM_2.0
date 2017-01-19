%%% Get COBALT data
function ENVR = get_COBALT(COBALT,ID,DY)

    global GRD NX TrefP TrefB Dthresh

    %% Get data
    ENVR.Tp(:,1)  = COBALT.Tp(ID,DY);
    ENVR.Tb(:,1)  = COBALT.Tb(ID,DY);
    ENVR.Zm(:,1)  = COBALT.Zm(ID,DY);
    ENVR.Zl(:,1)  = COBALT.Zl(ID,DY);
    ENVR.det(:,1) = COBALT.det(ID,DY);
    ENVR.dZm(:,1) = COBALT.dZm(ID,DY);
    ENVR.dZl(:,1) = COBALT.dZl(ID,DY);
    %ENVR.U(:,1)   = COBALT.U(ID,DY);
    %ENVR.V(:,1)   = COBALT.V(ID,DY);
    ENVR.T0p(:,1) = TrefP(ID);
    ENVR.T0b(:,1) = TrefB(ID);
    ENVR.Dthresh(:,1) = Dthresh(ID);
    ENVR.fZm(:,1) = int64(zeros(NX,1));
    ENVR.fZl(:,1) = int64(zeros(NX,1));
    ENVR.fB(:,1)  = int64(zeros(NX,1));
    ENVR.H(:,1)   = GRD.Z(ID);
    ENVR.A(:,1)   = GRD.AREA(ID);
end
