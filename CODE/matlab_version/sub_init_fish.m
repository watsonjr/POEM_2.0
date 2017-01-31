%============== INITIAL CONDITIONS =============%
function [Sf,Sp,Sd,Mf,Mp,Md,Lp,Ld,BENT] = sub_init_fish(ID,phen,DAYS)

    %===== VARIABLES =====%
    %%%! Number of spatial cells
    NX = length(ID);

    %%%! Fish biomass
    X = 1.0e-5; % very small amount
    Sf.bio = ones(NX,1) * X;
    Sp.bio = ones(NX,1) * X;
    Sd.bio = ones(NX,1) * X;
    Mf.bio = ones(NX,1) * X;
    Mp.bio = ones(NX,1) * X;
    Md.bio = ones(NX,1) * X;
    Lp.bio = ones(NX,1) * X;
    Ld.bio = ones(NX,1) * X;

    %%%! Fraction of time in pelagic
    Sf.td = ones(NX,1);
    Sp.td = ones(NX,1);
    Sd.td = ones(NX,1);
    Mf.td = ones(NX,1);
    Mp.td = ones(NX,1);
    Md.td = ones(NX,1);
    Lp.td = ones(NX,1);
    Ld.td = ones(NX,1);

    nzero = {'met' 'enc_f' 'enc_p' 'enc_d' 'enc_zm' 'enc_zl' 'enc_be' ...
        'con_f' 'con_p' 'con_d' 'con_zm' 'con_zl' 'con_be' 'I' 'die' 'pred' ...
        'nmort' 'prod' 'nu' 'gamma' 'egg' 'rep' 'rec' 'DD' 'clev' 'caught'};
    for i = 1 : length(nzero)
        Sf.(nzero{i}) = zeros(NX,1);
        Sp.(nzero{i}) = zeros(NX,1);
        Sd.(nzero{i}) = zeros(NX,1);
        Mf.(nzero{i}) = zeros(NX,1);
        Mp.(nzero{i}) = zeros(NX,1);
        Md.(nzero{i}) = zeros(NX,1);
        Lp.(nzero{i}) = zeros(NX,1);
        Ld.(nzero{i}) = zeros(NX,1);
    end

    %%%! Spawning flag
    % frac spawning at any given time
    if (phen == 1)
        Sf.S = zeros(NX,DAYS);
        Sp.S = zeros(NX,DAYS);
        Sd.S = zeros(NX,DAYS);
        Mf.S = zeros(NX,DAYS);
        Mp.S = zeros(NX,DAYS);
        Md.S = zeros(NX,DAYS);
        Lp.S = zeros(NX,DAYS);
        Ld.S = zeros(NX,DAYS);
    else
        Sf.S = ones(NX,DAYS);
        Sp.S = ones(NX,DAYS);
        Sd.S = ones(NX,DAYS);
        Mf.S = ones(NX,DAYS);
        Mp.S = ones(NX,DAYS);
        Md.S = ones(NX,DAYS);
        Lp.S = ones(NX,DAYS);
        Ld.S = ones(NX,DAYS);
    end

    %%%! Detritus
    BENT.mass = ones(NX,1) * X;
    BENT.pred = zeros(NX,1);

end
