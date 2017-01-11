%%%% THE MODEL
%%% DEMOGRAPHIC CALCULATIONS
function sub_futbio(ID,DY,COBALT,ENVR,Sf,Sp,Sd,Mf,Mp,Md,Lp,Ld,BENT,dfrate,rfrac)

%%% COBALT information
get_COBALT(COBALT,ID,DY,ENVR)

% Update benthic biomass with new detritus avail at that time step
BENT.mass = BENT.mass + bent_eff.*ENVR.det;

% Pelagic-demersal coupling
%Lp: fraction of time large piscivores spends in pelagic
%Ld: fraction of time large demersals spends in pelagic
if (pdc == 0)
    Lp.td = 1.0;
    Ld.td = 0.0;
elseif (pdc == 1)
    Lp.td = 1.0;
    Ld.td = sub_tdif_dem(ENVR.H,Mf.bio,Mp.bio,Md.bio,BENT.mass);
elseif (pdc == 2)
    Lp.td = sub_tdif_pel(ENVR.H,Mf.bio,Mp.bio,Md.bio);
    Ld.td = sub_tdif_dem(ENVR.H,Mf.bio,Mp.bio,Md.bio,BENT.mass);
end

% Metabolism
Sf.met = sub_met(ENVR.Tp,ENVR.Tb,Sf.td,M_s,L_s);
Sp.met = sub_met(ENVR.Tp,ENVR.Tb,Sp.td,M_s,L_s);
Sd.met = sub_met(ENVR.Tp,ENVR.Tb,Sd.td,M_s,L_s);
Mf.met = sub_met(ENVR.Tp,ENVR.Tb,Mf.td,M_m,L_m);
Mp.met = sub_met(ENVR.Tp,ENVR.Tb,Mp.td,M_m,L_m);
Md.met = sub_met(ENVR.Tp,ENVR.Tb,Md.td,M_m,L_m);
Lp.met = sub_met(ENVR.Tp,ENVR.Tb,Lp.td,M_l,L_l);
Ld.met = sub_met(ENVR.Tp,ENVR.Tb,Ld.td,M_l,L_l);

% Encounter rates
%sub_enc(Tp,Tb,wgt,pred,prey,td,tprey,pref)
Sf.enc_zm = sub_enc_pel(ENVR.Tp,ENVR.Tb,M_s,ENVR.Zm,Sf.td,Sf.td,1);
Sp.enc_zm = sub_enc_pel(ENVR.Tp,ENVR.Tb,M_s,ENVR.Zm,Sp.td,Sf.td,1);
Sd.enc_zm = sub_enc_pel(ENVR.Tp,ENVR.Tb,M_s,ENVR.Zm,Sd.td,Sf.td,1);

Mf.enc_zm = sub_enc_pel(ENVR.Tp,ENVR.Tb,M_m,ENVR.Zm,Mf.td,Mf.td,MF_phi_MZ);
Mf.enc_zl = sub_enc_pel(ENVR.Tp,ENVR.Tb,M_m,ENVR.Zl,Mf.td,Mf.td,MF_phi_LZ);
Mf.enc_f  = sub_enc_pel(ENVR.Tp,ENVR.Tb,M_m,Sf.bio,Mf.td,Mf.td,MF_phi_S);
Mf.enc_p  = sub_enc_pel(ENVR.Tp,ENVR.Tb,M_m,Sp.bio,Mf.td,Mf.td,MF_phi_S);
Mf.enc_d  = sub_enc_pel(ENVR.Tp,ENVR.Tb,M_m,Sd.bio,Mf.td,Mf.td,MF_phi_S);

Mp.enc_zm = sub_enc_pel(ENVR.Tp,ENVR.Tb,M_m,ENVR.Zm,Mp.td,Mp.td,MP_phi_MZ);
Mp.enc_zl = sub_enc_pel(ENVR.Tp,ENVR.Tb,M_m,ENVR.Zl,Mp.td,Mp.td,MP_phi_LZ);
Mp.enc_f  = sub_enc_pel(ENVR.Tp,ENVR.Tb,M_m,Sf.bio,Mp.td,Mp.td,MP_phi_S);
Mp.enc_p  = sub_enc_pel(ENVR.Tp,ENVR.Tb,M_m,Sp.bio,Mp.td,Mp.td,MP_phi_S);
Mp.enc_d  = sub_enc_pel(ENVR.Tp,ENVR.Tb,M_m,Sd.bio,Mp.td,Mp.td,MP_phi_S);

Md.enc_be = sub_enc_dem(ENVR.Tp,ENVR.Tb,M_m,BENT.mass,Md.td,1-Md.td,MD_phi_BE);

Lp.enc_f  = sub_enc_pel(ENVR.Tp,ENVR.Tb,M_l,Mf.bio,Lp.td,Lp.td,LP_phi_MF);
Lp.enc_p  = sub_enc_pel(ENVR.Tp,ENVR.Tb,M_l,Mp.bio,Lp.td,Lp.td,LP_phi_MP);
Lp.enc_d  = sub_enc_pel(ENVR.Tp,ENVR.Tb,M_l,Md.bio,Lp.td,1-Lp.td,LP_phi_MD);

Ld.enc_f  = sub_enc_dem(ENVR.Tp,ENVR.Tb,M_l,Mf.bio,Ld.td,Ld.td,LD_phi_MF);
Ld.enc_p  = sub_enc_dem(ENVR.Tp,ENVR.Tb,M_l,Mp.bio,Ld.td,Ld.td,LD_phi_MP);
Ld.enc_d  = sub_enc_dem(ENVR.Tp,ENVR.Tb,M_l,Md.bio,Ld.td,1-Ld.td,LD_phi_MD);
Ld.enc_be = sub_enc_dem(ENVR.Tp,ENVR.Tb,M_l,BENT.mass,Ld.td,1-Ld.td,LD_phi_BE);

% Consumption rates
Sf.con_zm = sub_cons(ENVR.Tp,ENVR.Tb,Sf.td,M_s,Sf.enc_zm);
Sp.con_zm = sub_cons(ENVR.Tp,ENVR.Tb,Sp.td,M_s,Sp.enc_zm);
Sd.con_zm = sub_cons(ENVR.Tp,ENVR.Tb,Sd.td,M_s,Sd.enc_zm);

Mf.con_zm = sub_cons(ENVR.Tp,ENVR.Tb,Mf.td,M_m,[Mf.enc_zm,Mf.enc_zl,Mf.enc_f,Mf.enc_p,Mf.enc_d]);
Mf.con_zl = sub_cons(ENVR.Tp,ENVR.Tb,Mf.td,M_m,[Mf.enc_zl,Mf.enc_zm,Mf.enc_f,Mf.enc_p,Mf.enc_d]);
Mf.con_f  = sub_cons(ENVR.Tp,ENVR.Tb,Mf.td,M_m,[Mf.enc_f,Mf.enc_zm,Mf.enc_zl,Mf.enc_p,Mf.enc_d]);
Mf.con_p  = sub_cons(ENVR.Tp,ENVR.Tb,Mf.td,M_m,[Mf.enc_p,Mf.enc_zm,Mf.enc_zl,Mf.enc_f,Mf.enc_d]);
Mf.con_d  = sub_cons(ENVR.Tp,ENVR.Tb,Mf.td,M_m,[Mf.enc_d,Mf.enc_zm,Mf.enc_zl,Mf.enc_f,Mf.enc_p]);

Mp.con_zm = sub_cons(ENVR.Tp,ENVR.Tb,Mp.td,M_m,[Mp.enc_zm,Mp.enc_zl,Mp.enc_f,Mp.enc_p,Mp.enc_d]);
Mp.con_zl = sub_cons(ENVR.Tp,ENVR.Tb,Mp.td,M_m,[Mp.enc_zl,Mp.enc_zm,Mp.enc_f,Mp.enc_p,Mp.enc_d]);
Mp.con_f  = sub_cons(ENVR.Tp,ENVR.Tb,Mp.td,M_m,[Mp.enc_f,Mp.enc_zm,Mp.enc_zl,Mp.enc_p,Mp.enc_d]);
Mp.con_p  = sub_cons(ENVR.Tp,ENVR.Tb,Mp.td,M_m,[Mp.enc_p,Mp.enc_zm,Mp.enc_zl,Mp.enc_f,Mp.enc_d]);
Mp.con_d  = sub_cons(ENVR.Tp,ENVR.Tb,Mp.td,M_m,[Mp.enc_d,Mp.enc_zm,Mp.enc_zl,Mp.enc_f,Mp.enc_p]);

Md.con_be = sub_cons(ENVR.Tp,ENVR.Tb,Md.td,M_m,Md.enc_be);

Lp.con_f  = sub_cons(ENVR.Tp,ENVR.Tb,Lp.td,M_l,[Lp.enc_f,Lp.enc_p,Lp.enc_d]);
Lp.con_p  = sub_cons(ENVR.Tp,ENVR.Tb,Lp.td,M_l,[Lp.enc_p,Lp.enc_f,Lp.enc_d]);
Lp.con_d  = sub_cons(ENVR.Tp,ENVR.Tb,Lp.td,M_l,[Lp.enc_d,Lp.enc_p,Lp.enc_f]);

Ld.con_f  = sub_cons(ENVR.Tp,ENVR.Tb,Ld.td,M_l,[Ld.enc_f,Ld.enc_p,Ld.enc_d,Ld.enc_be]);
Ld.con_p  = sub_cons(ENVR.Tp,ENVR.Tb,Ld.td,M_l,[Ld.enc_p,Ld.enc_f,Ld.enc_d,Ld.enc_be]);
Ld.con_d  = sub_cons(ENVR.Tp,ENVR.Tb,Ld.td,M_l,[Ld.enc_d,Ld.enc_p,Ld.enc_f,Ld.enc_be]);
Ld.con_be = sub_cons(ENVR.Tp,ENVR.Tb,Ld.td,M_l,[Ld.enc_be,Ld.enc_f,Ld.enc_p,Ld.enc_d]);


% Offline coupling
%Zooplankton consumption cannot exceed amount lost to higher predation in COBALT runs
[Sf.con_zm,Sp.con_zm,Sd.con_zm,Mf.con_zm,Mp.con_zm,ENVR.fZm] = ...
    sub_offline_zm(Sf.con_zm,Sp.con_zm,Sd.con_zm,Mf.con_zm,Mp.con_zm,Sf.bio,Sp.bio,Sd.bio,Mf.bio,Mp.bio,ENVR.dZm);
[Mf.con_zl,Mp.con_zl,ENVR.fZl] = ...
    sub_offline_zl(Mf.con_zl,Mp.con_zl,Mf.bio,Mp.bio,ENVR.dZl);
%Benthic material consumption cannot exceed amount present
[Md.con_be, Ld.con_be, ENVR.fB] = ...
    sub_offline_bent(Md.con_be,Ld.con_be,Md.bio,Ld.bio,BENT.mass);

% Total consumption rates (could factor in handling times here; g m-2 d-1)
Sf.I = Sf.con_zm;
Sp.I = Sp.con_zm;
Sd.I = Sd.con_zm;
Mf.I = Mf.con_zm + Mf.con_zl + Mf.con_f + Mf.con_p + Mf.con_d;
Mp.I = Mp.con_zm + Mp.con_zl + Mp.con_f + Mp.con_p + Mp.con_d;
Md.I = Md.con_be;
Lp.I = Lp.con_f + Lp.con_p + Lp.con_d;
Ld.I = Ld.con_f + Ld.con_p + Ld.con_d + Ld.con_be;

% Consumption related to Cmax
Sf.clev = sub_clev(Sf.I,ENVR.Tp,ENVR.Tb,Sf.td,M_s);
Sp.clev = sub_clev(Sp.I,ENVR.Tp,ENVR.Tb,Sp.td,M_s);
Sd.clev = sub_clev(Sd.I,ENVR.Tp,ENVR.Tb,Sd.td,M_s);
Mf.clev = sub_clev(Mf.I,ENVR.Tp,ENVR.Tb,Mf.td,M_m);
Mp.clev = sub_clev(Mp.I,ENVR.Tp,ENVR.Tb,Mp.td,M_m);
Md.clev = sub_clev(Md.I,ENVR.Tp,ENVR.Tb,Md.td,M_m);
Lp.clev = sub_clev(Lp.I,ENVR.Tp,ENVR.Tb,Lp.td,M_l);
Ld.clev = sub_clev(Ld.I,ENVR.Tp,ENVR.Tb,Ld.td,M_l);

% Death rates (g m-2 d-1)
Sf.die = Mp.con_f.*Mp.bio + Mf.con_f.*Mf.bio;
Sp.die = Mp.con_p.*Mp.bio + Mf.con_p.*Mf.bio;
Sd.die = Mp.con_d.*Mp.bio + Mf.con_d.*Mf.bio;
Mf.die = Lp.con_f.*Lp.bio + Ld.con_f.*Ld.bio;
Mp.die = Lp.con_p.*Lp.bio + Ld.con_p.*Ld.bio;
Md.die = Lp.con_d.*Lp.bio + Ld.con_d.*Ld.bio;

% predation rates (m-2 d-1)
Sf.pred = Sf.die ./ Sf.bio;
Sp.pred = Sp.die ./ Sp.bio;
Sd.pred = Sd.die ./ Sd.bio;
Mf.pred = Mf.die ./ Mf.bio;
Mp.pred = Mp.die ./ Mp.bio;
Md.pred = Md.die ./ Md.bio;

% Natural mortality rates
Sf.nmort = sub_nmort(ENVR.Tp,ENVR.Tb,Sf.td,M_s);
Sp.nmort = sub_nmort(ENVR.Tp,ENVR.Tb,Sp.td,M_s);
Sd.nmort = sub_nmort(ENVR.Tp,ENVR.Tb,Sd.td,M_s);
Mf.nmort = sub_nmort(ENVR.Tp,ENVR.Tb,Mf.td,M_m);
Mp.nmort = sub_nmort(ENVR.Tp,ENVR.Tb,Mp.td,M_m);
Md.nmort = sub_nmort(ENVR.Tp,ENVR.Tb,Md.td,M_m);
Lp.nmort = sub_nmort(ENVR.Tp,ENVR.Tb,Lp.td,M_l);
Ld.nmort = sub_nmort(ENVR.Tp,ENVR.Tb,Ld.td,M_l);

% Degree days
Mf.DD = sub_degday(Mf.DD,ENVR.Tp,ENVR.Tb,Mf.td,ENVR.T0p,Mf.S,DY);
Lp.DD = sub_degday(Lp.DD,ENVR.Tp,ENVR.Tb,Lp.td,ENVR.T0p,Lp.S,DY);
Ld.DD = sub_degday(Ld.DD,ENVR.Tp,ENVR.Tb,1-Ld.td,ENVR.T0b,Ld.S,DY);

% Spawning flag determined from DD, dthresh
[Mf.S, Mf.DD] = sub_kflag(Mf.S,Mf.DD,ENVR.Dthresh,DY);
[Ld.S, Ld.DD] = sub_kflag(Ld.S,Ld.DD,ENVR.Dthresh,DY);
[Lp.S, Lp.DD] = sub_kflag(Lp.S,Lp.DD,ENVR.Dthresh,DY);

% Energy available for somatic growth nu
[Sf.nu, Sf.prod] = sub_nu(Sf.I,Sf.bio,Sf.met);
[Sp.nu, Sp.prod] = sub_nu(Sp.I,Sp.bio,Sp.met);
[Sd.nu, Sd.prod] = sub_nu(Sd.I,Sd.bio,Sd.met);
[Mf.nu, Mf.prod] = sub_nu(Mf.I,Mf.bio,Mf.met);
[Mp.nu, Mp.prod] = sub_nu(Mp.I,Mp.bio,Mp.met);
[Md.nu, Md.prod] = sub_nu_MD(Md.I,Md.bio,Md.met);
[Lp.nu, Lp.prod] = sub_nu(Lp.I,Lp.bio,Lp.met);
[Ld.nu, Ld.prod] = sub_nu_LD(Ld.bio,Ld.met,Ld.con_f,Ld.con_p,Ld.con_d,Ld.con_be);

% Maturation (note subscript on Kappa is larvae, juv, adult)
Sf.gamma = sub_gamma(K_l,Z_s,Sf.nu,Sf.die,Sf.bio,Sf.S(:,DY),Sf.nmort);
Sp.gamma = sub_gamma(K_l,Z_s,Sp.nu,Sp.die,Sp.bio,Sp.S(:,DY),Sp.nmort);
Sd.gamma = sub_gamma(K_l,Z_s,Sd.nu,Sd.die,Sd.bio,Sd.S(:,DY),Sd.nmort);
Mf.gamma = sub_gamma(K_a,Z_m,Mf.nu,Mf.die,Mf.bio,Mf.S(:,DY),Mf.nmort);
Mp.gamma = sub_gamma(K_j,Z_m,Mp.nu,Mp.die,Mp.bio,Mp.S(:,DY),Mp.nmort);
Md.gamma = sub_gamma(K_j,Z_m,Md.nu,Md.die,Md.bio,Md.S(:,DY),Md.nmort);
Lp.gamma = sub_gamma(K_a,Z_l,Lp.nu,Lp.die,Lp.bio,Lp.S(:,DY),Lp.nmort);
Ld.gamma = sub_gamma(K_a,Z_l,Ld.nu,Ld.die,Ld.bio,Ld.S(:,DY),Ld.nmort);

% Egg production (by med and large size classes only)
[Sf.rep,Sf.egg] = sub_rep(Sf.nu,K_l,Sf.S(:,DY),Sf.egg);
[Sp.rep,Sp.egg] = sub_rep(Sp.nu,K_l,Sp.S(:,DY),Sp.egg);
[Sd.rep,Sd.egg] = sub_rep(Sd.nu,K_l,Sd.S(:,DY),Sd.egg);
[Mf.rep,Mf.egg] = sub_rep(Mf.nu,K_a,Mf.S(:,DY),Mf.egg);
[Mp.rep,Mp.egg] = sub_rep(Mp.nu,K_j,Mp.S(:,DY),Mp.egg);
[Md.rep,Md.egg] = sub_rep(Md.nu,K_j,Md.S(:,DY),Md.egg);
[Lp.rep,Lp.egg] = sub_rep(Lp.nu,K_a,Lp.S(:,DY),Lp.egg);
[Ld.rep,Ld.egg] = sub_rep(Ld.nu,K_a,Ld.S(:,DY),Ld.egg);

% Recruitment (from smaller size class)
Sf.rec = sub_rec(Mf.rep,Mf.bio,M_s,rfrac);
Sp.rec = sub_rec(Lp.rep,Lp.bio,M_s,rfrac);
Sd.rec = sub_rec(Ld.rep,Ld.bio,M_s,rfrac);
Mf.rec = sub_rec(Sf.gamma,Sf.bio,M_m,rfrac);
Mp.rec = sub_rec(Sp.gamma,Sp.bio,M_m,rfrac);
Md.rec = sub_rec(Sd.gamma,Sd.bio,M_m,rfrac);
Lp.rec = sub_rec(Mp.gamma,Mp.bio,M_l,rfrac);
Ld.rec = sub_rec(Md.gamma,Md.bio,M_l,rfrac);

% Mass balance
BENT.mass = sub_update_be(BENT.mass,[Md.con_be,Ld.con_be],[Md.bio,Ld.bio]);

Sf.bio = sub_update_fi(Sf.bio,Sf.rec,Sf.nu,Sf.rep,Sf.gamma,Sf.die,Sf.egg,Sf.nmort);
Sp.bio = sub_update_fi(Sp.bio,Sp.rec,Sp.nu,Sp.rep,Sp.gamma,Sp.die,Sp.egg,Sp.nmort);
Sd.bio = sub_update_fi(Sd.bio,Sd.rec,Sd.nu,Sd.rep,Sd.gamma,Sd.die,Sd.egg,Sd.nmort);

Mf.bio = sub_update_fi(Mf.bio,Mf.rec,Mf.nu,Mf.rep,Mf.gamma,Mf.die,Mf.egg,Mf.nmort);
Mp.bio = sub_update_fi(Mp.bio,Mp.rec,Mp.nu,Mp.rep,Mp.gamma,Mp.die,Mp.egg,Mp.nmort);
Md.bio = sub_update_fi(Md.bio,Md.rec,Md.nu,Md.rep,Md.gamma,Md.die,Md.egg,Md.nmort);

Lp.bio = sub_update_lg(Lp.bio,Lp.rec,Lp.nu,Lp.rep,Lp.gamma,Lp.die,Lp.egg,Lp.nmort;
Ld.bio = sub_update_lg(Ld.bio,Ld.rec,Ld.nu,Ld.rep,Ld.gamma,Ld.die,Ld.egg,Ld.nmort);

% Fishing by rate
%[Mf.bio, Mf.caught] = sub_fishing_rate(Mf.bio,M_m,dfrate);
%[Mp.bio, Mp.caught] = sub_fishing_rate(Mp.bio,M_m,dfrate);
%[Md.bio, Md.caught] = sub_fishing_rate(Md.bio,M_m,dfrate);
%[Lp.bio, Lp.caught] = sub_fishing_rate(Lp.bio,M_l,dfrate);
[Ld.bio, Ld.caught] = sub_fishing_rate(Ld.bio,M_l,dfrate);

% Forward Euler checks for demographics and movement
sub_check(Sf.bio);
sub_check(Sp.bio);
sub_check(Sd.bio);
sub_check(Mf.bio);
sub_check(Mp.bio);
sub_check(Md.bio);
sub_check(Lp.bio);
sub_check(Ld.bio);

end