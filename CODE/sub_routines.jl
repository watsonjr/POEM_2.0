
#### THE MODEL
###! DEMOGRAPHIC CALCULATIONS
function sub_futbio!(ID,DY,COBALT,ENVR,Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT)

	###! COBALT information
	get_COBALT!(COBALT,ID,DY,ENVR)

	for JD = 1:NX

		#! update benthic biomass with new detritus avail at that time step
		BENT.mass[JD] = BENT.mass[JD] + bent_eff*ENVR.det[JD]

		#! fraction of time large piscivores spends in pelagic
	  #Lrg_p.td[JD] = sub_tdif_pel(GRD["Z"][ID],Med_f.bio[JD],Med_p.bio[JD],Med_d.bio[JD])
		Lrg_p.td[JD] = 1.0
		#! fraction of time large demersal spends in pelagic
	  #Lrg_d.td[JD] = sub_tdif_dem(GRD["Z"][ID],Med_f.bio[JD],Med_p.bio[JD],Med_d.bio[JD],BENT.mass[JD])
		Lrg_d.td[JD] = 0.0

		#! metabolism
		Sml_f.met[JD] = sub_met(ENVR.Tp[JD],ENVR.Tb[JD],Sml_f.td[JD],M_s,L_s)
		Sml_p.met[JD] = sub_met(ENVR.Tp[JD],ENVR.Tb[JD],Sml_p.td[JD],M_s,L_s)
		Sml_d.met[JD] = sub_met(ENVR.Tp[JD],ENVR.Tb[JD],Sml_d.td[JD],M_s,L_s)
		Med_f.met[JD] = sub_met(ENVR.Tp[JD],ENVR.Tb[JD],Med_f.td[JD],M_m,L_m)
		Med_p.met[JD] = sub_met(ENVR.Tp[JD],ENVR.Tb[JD],Med_p.td[JD],M_m,L_m)
		Med_d.met[JD] = sub_met(ENVR.Tp[JD],ENVR.Tb[JD],Med_d.td[JD],M_m,L_m)
		Lrg_p.met[JD] = sub_met(ENVR.Tp[JD],ENVR.Tb[JD],Lrg_p.td[JD],M_l,L_l)
		Lrg_d.met[JD] = sub_met(ENVR.Tp[JD],ENVR.Tb[JD],Lrg_d.td[JD],M_l,L_l)

		#! encounter rates
		#sub_enc(Tp,Tb,tdif,wgt,L,tu,pred,prey,td)
		Sml_f.enc_zm[JD] = sub_enc(ENVR.Tp[JD],ENVR.Tb[JD],M_s,L_s,Tu_s,Sml_f.bio[JD],ENVR.Zm[JD],Sml_f.td[JD],1)
		Sml_p.enc_zm[JD] = sub_enc(ENVR.Tp[JD],ENVR.Tb[JD],M_s,L_s,Tu_s,Sml_p.bio[JD],ENVR.Zm[JD],Sml_p.td[JD],1)
		Sml_d.enc_zm[JD] = sub_enc(ENVR.Tp[JD],ENVR.Tb[JD],M_s,L_s,Tu_s,Sml_d.bio[JD],ENVR.Zm[JD],Sml_d.td[JD],1)

		Med_f.enc_zl[JD] = sub_enc(ENVR.Tp[JD],ENVR.Tb[JD],M_m,L_m,Tu_m,Med_f.bio[JD],ENVR.Zl[JD],Med_f.td[JD],1)
		Med_f.enc_f[JD]  = sub_enc(ENVR.Tp[JD],ENVR.Tb[JD],M_m,L_m,Tu_m,Med_f.bio[JD],Sml_f.bio[JD],Med_f.td[JD],1)
		Med_f.enc_p[JD]  = sub_enc(ENVR.Tp[JD],ENVR.Tb[JD],M_m,L_m,Tu_m,Med_f.bio[JD],Sml_p.bio[JD],Med_f.td[JD],1)
		Med_f.enc_d[JD]  = sub_enc(ENVR.Tp[JD],ENVR.Tb[JD],M_m,L_m,Tu_m,Med_f.bio[JD],Sml_d.bio[JD],Med_f.td[JD],1)

		Med_p.enc_zl[JD] = sub_enc(ENVR.Tp[JD],ENVR.Tb[JD],M_m,L_m,Tu_m,Med_p.bio[JD],ENVR.Zl[JD],Med_p.td[JD],1)
		Med_p.enc_f[JD]  = sub_enc(ENVR.Tp[JD],ENVR.Tb[JD],M_m,L_m,Tu_m,Med_p.bio[JD],Sml_f.bio[JD],Med_p.td[JD],1)
		Med_p.enc_p[JD]  = sub_enc(ENVR.Tp[JD],ENVR.Tb[JD],M_m,L_m,Tu_m,Med_p.bio[JD],Sml_p.bio[JD],Med_p.td[JD],1)
		Med_p.enc_d[JD]  = sub_enc(ENVR.Tp[JD],ENVR.Tb[JD],M_m,L_m,Tu_m,Med_p.bio[JD],Sml_d.bio[JD],Med_p.td[JD],1)

		Med_d.enc_be[JD] = sub_enc(ENVR.Tp[JD],ENVR.Tb[JD],M_m,L_m,Tu_m,Med_d.bio[JD],BENT.mass[JD],Med_d.td[JD],1)

		Lrg_p.enc_f[JD]  = sub_enc(ENVR.Tp[JD],ENVR.Tb[JD],M_l,L_l,Tu_l,Lrg_p.bio[JD],Med_f.bio[JD],Lrg_p.td[JD],Lrg_p.td[JD])
		Lrg_p.enc_p[JD]  = sub_enc(ENVR.Tp[JD],ENVR.Tb[JD],M_l,L_l,Tu_l,Lrg_p.bio[JD],Med_p.bio[JD],Lrg_p.td[JD],Lrg_p.td[JD])
		#Lrg_p.enc_d[JD]  = sub_enc(ENVR.Tp[JD],ENVR.Tb[JD],M_l,L_l,Tu_l,Lrg_p.bio[JD],Med_d.bio[JD],Lrg_p.td[JD],1-Lrg_p.td[JD])

		#Lrg_d.enc_f[JD]  = sub_enc(ENVR.Tp[JD],ENVR.Tb[JD],M_l,L_l,Tu_l,Lrg_d.bio[JD],Med_f.bio[JD],Lrg_d.td[JD],Lrg_d.td[JD])
		#Lrg_d.enc_p[JD]  = sub_enc(ENVR.Tp[JD],ENVR.Tb[JD],M_l,L_l,Tu_l,Lrg_d.bio[JD],Med_p.bio[JD],Lrg_d.td[JD],Lrg_d.td[JD])
		Lrg_d.enc_d[JD]  = sub_enc(ENVR.Tp[JD],ENVR.Tb[JD],M_l,L_l,Tu_l,Lrg_d.bio[JD],Med_d.bio[JD],Lrg_d.td[JD],1-Lrg_d.td[JD])
		Lrg_d.enc_be[JD] = sub_enc(ENVR.Tp[JD],ENVR.Tb[JD],M_l,L_l,Tu_l,Lrg_d.bio[JD],BENT.mass[JD],Lrg_d.td[JD],1-Lrg_d.td[JD])

		#! Consumption rates
		Sml_f.con_zm[JD] = sub_cons(ENVR.Tp[JD],ENVR.Tb[JD],Sml_f.td[JD],M_s,Sml_f.enc_zm[JD])
		Sml_p.con_zm[JD] = sub_cons(ENVR.Tp[JD],ENVR.Tb[JD],Sml_p.td[JD],M_s,Sml_p.enc_zm[JD])
		Sml_d.con_zm[JD] = sub_cons(ENVR.Tp[JD],ENVR.Tb[JD],Sml_d.td[JD],M_s,Sml_d.enc_zm[JD])

		Med_f.con_zl[JD] = sub_cons(ENVR.Tp[JD],ENVR.Tb[JD],Med_f.td[JD],M_m,Med_f.enc_zl[JD])

		Med_p.con_zl[JD] = sub_cons(ENVR.Tp[JD],ENVR.Tb[JD],Med_p.td[JD],M_m,[Med_p.enc_zl[JD],Med_p.enc_f[JD],Med_p.enc_p[JD],Med_p.enc_d[JD]])
		Med_p.con_f[JD]  = sub_cons(ENVR.Tp[JD],ENVR.Tb[JD],Med_p.td[JD],M_m,[Med_p.enc_f[JD],Med_p.enc_zl[JD],Med_p.enc_p[JD],Med_p.enc_d[JD]])
		Med_p.con_p[JD]  = sub_cons(ENVR.Tp[JD],ENVR.Tb[JD],Med_p.td[JD],M_m,[Med_p.enc_p[JD],Med_p.enc_zl[JD],Med_p.enc_f[JD],Med_p.enc_d[JD]])
		Med_p.con_d[JD]  = sub_cons(ENVR.Tp[JD],ENVR.Tb[JD],Med_p.td[JD],M_m,[Med_p.enc_d[JD],Med_p.enc_zl[JD],Med_p.enc_f[JD],Med_p.enc_p[JD]])

		Med_d.con_be[JD] = sub_cons(ENVR.Tp[JD],ENVR.Tb[JD],Med_d.td[JD],M_m,Med_d.enc_be[JD])

		Lrg_p.con_f[JD]  = sub_cons(ENVR.Tp[JD],ENVR.Tb[JD],Lrg_p.td[JD],M_l,[Lrg_p.enc_f[JD],Lrg_p.enc_p[JD],Lrg_p.enc_d[JD]])
		Lrg_p.con_p[JD]  = sub_cons(ENVR.Tp[JD],ENVR.Tb[JD],Lrg_p.td[JD],M_l,[Lrg_p.enc_p[JD],Lrg_p.enc_f[JD],Lrg_p.enc_d[JD]])
		Lrg_p.con_d[JD]  = sub_cons(ENVR.Tp[JD],ENVR.Tb[JD],Lrg_p.td[JD],M_l,[Lrg_p.enc_d[JD],Lrg_p.enc_p[JD],Lrg_p.enc_f[JD]])

		Lrg_d.con_f[JD]  = sub_cons(ENVR.Tp[JD],ENVR.Tb[JD],Lrg_d.td[JD],M_l,[Lrg_d.enc_f[JD],Lrg_d.enc_p[JD],Lrg_d.enc_d[JD],Lrg_d.enc_be[JD]])
		Lrg_d.con_p[JD]  = sub_cons(ENVR.Tp[JD],ENVR.Tb[JD],Lrg_d.td[JD],M_l,[Lrg_d.enc_p[JD],Lrg_d.enc_f[JD],Lrg_d.enc_d[JD],Lrg_d.enc_be[JD]])
		Lrg_d.con_d[JD]  = sub_cons(ENVR.Tp[JD],ENVR.Tb[JD],Lrg_d.td[JD],M_l,[Lrg_d.enc_d[JD],Lrg_d.enc_p[JD],Lrg_d.enc_f[JD],Lrg_d.enc_be[JD]])
		Lrg_d.con_be[JD] = sub_cons(ENVR.Tp[JD],ENVR.Tb[JD],Lrg_d.td[JD],M_l,[Lrg_d.enc_be[JD],Lrg_d.enc_f[JD],Lrg_d.enc_p[JD],Lrg_d.enc_d[JD]])


		#! Offline coupling
		#Zooplankton consumption cannot exceed amount lost to higher predation in COBALT runs
		Sml_f.con_zm[JD], Sml_p.con_zm[JD], Sml_d.con_zm[JD], ENVR.fZm[JD] = sub_offline(Sml_f.con_zm[JD],Sml_p.con_zm[JD],Sml_d.con_zm[JD],Sml_f.bio[JD],Sml_p.bio[JD],Sml_d.bio[JD],ENVR.dZm[JD])
		Med_f.con_zl[JD], Med_p.con_zl[JD], Med_d.con_zl[JD], ENVR.fZl[JD] = sub_offline(Med_f.con_zl[JD],Med_p.con_zl[JD],Med_d.con_zl[JD],Med_f.bio[JD],Med_p.bio[JD],Med_d.bio[JD],ENVR.dZl[JD])
		#Benthic material consumption cannot exceed amount present
		Med_d.con_be[JD], Lrg_d.con_be[JD], ENVR.fB[JD] = sub_offline_bent(Med_d.con_be[JD],Lrg_d.con_be[JD],Med_d.bio[JD],Lrg_d.bio[JD],BENT.mass[JD],ENVR.det[JD])

		#! total consumption rates (could factor in handling times here; g m-2 d-1)
		Sml_f.I[JD] = Sml_f.con_zm[JD]
		Sml_p.I[JD] = Sml_p.con_zm[JD]
		Sml_d.I[JD] = Sml_d.con_zm[JD]
		Med_f.I[JD] = Med_f.con_zl[JD]
		Med_p.I[JD] = Med_p.con_zl[JD] + Med_p.con_f[JD] + Med_p.con_p[JD] + Med_p.con_d[JD]
		Med_d.I[JD] = Med_d.con_be[JD]
		Lrg_p.I[JD] = Lrg_p.con_f[JD] + Lrg_p.con_p[JD] + Lrg_p.con_d[JD]
		Lrg_d.I[JD] = Lrg_d.con_f[JD] + Lrg_d.con_p[JD] + Lrg_d.con_d[JD] + Lrg_d.con_be[JD]

		#! consumption related to Cmax
		Sml_f.clev[JD] = sub_clev(Sml_f.I[JD],ENVR.Tp[JD],ENVR.Tb[JD],Sml_f.td[JD],M_s)
		Sml_p.clev[JD] = sub_clev(Sml_p.I[JD],ENVR.Tp[JD],ENVR.Tb[JD],Sml_p.td[JD],M_s)
		Sml_d.clev[JD] = sub_clev(Sml_d.I[JD],ENVR.Tp[JD],ENVR.Tb[JD],Sml_d.td[JD],M_s)
		Med_f.clev[JD] = sub_clev(Med_f.I[JD],ENVR.Tp[JD],ENVR.Tb[JD],Med_f.td[JD],M_m)
		Med_p.clev[JD] = sub_clev(Med_p.I[JD],ENVR.Tp[JD],ENVR.Tb[JD],Med_p.td[JD],M_m)
		Med_d.clev[JD] = sub_clev(Med_d.I[JD],ENVR.Tp[JD],ENVR.Tb[JD],Med_d.td[JD],M_m)
		Lrg_p.clev[JD] = sub_clev(Lrg_p.I[JD],ENVR.Tp[JD],ENVR.Tb[JD],Lrg_p.td[JD],M_l)
		Lrg_d.clev[JD] = sub_clev(Lrg_d.I[JD],ENVR.Tp[JD],ENVR.Tb[JD],Lrg_d.td[JD],M_l)

		#! death rates (g m-2 d-1)
		Sml_f.die[JD] = Med_p.con_f[JD]*Med_p.bio[JD]
		Sml_p.die[JD] = Med_p.con_p[JD]*Med_p.bio[JD]
		Sml_d.die[JD] = Med_p.con_d[JD]*Med_p.bio[JD]
		Med_f.die[JD] = Lrg_p.con_f[JD]*Lrg_p.bio[JD] + Lrg_d.con_f[JD]*Lrg_d.bio[JD]
		Med_p.die[JD] = Lrg_p.con_p[JD]*Lrg_p.bio[JD] + Lrg_d.con_p[JD]*Lrg_d.bio[JD]
		Med_d.die[JD] = Lrg_p.con_d[JD]*Lrg_p.bio[JD] + Lrg_d.con_d[JD]*Lrg_d.bio[JD]

		#! Degree days
		Med_f.DD[JD] = sub_degday(Med_f.DD[JD],ENVR.Tp[JD],ENVR.Tb[JD],Med_f.td[JD],ENVR.T0p[JD],Med_f.S[JD,:],DY)
		Lrg_p.DD[JD] = sub_degday(Lrg_p.DD[JD],ENVR.Tp[JD],ENVR.Tb[JD],Lrg_p.td[JD],ENVR.T0p[JD],Lrg_p.S[JD,:],DY)
		#Lrg_d.DD[JD] = sub_degday(Lrg_d.DD[JD],ENVR.Tp[JD],ENVR.Tb[JD],Lrg_d.td[JD],ENVR.T0b[JD],Lrg_d.S[JD,:],DY)
		#Assume demersal spawn at same time as pelagic b/c larvae also need spring bloom
		Lrg_d.DD[JD] = sub_degday(Lrg_d.DD[JD],ENVR.Tp[JD],ENVR.Tb[JD],1-Lrg_d.td[JD],ENVR.T0b[JD],Lrg_d.S[JD,:],DY)

		#! Spawning flag determined from DD, dthresh
		Med_f.S[JD,:], Med_f.DD[JD] = sub_kflag(Med_f.S[JD,:],Med_f.DD[JD],ENVR.Dthresh[JD],DY);
		Lrg_d.S[JD,:], Lrg_d.DD[JD] = sub_kflag(Lrg_d.S[JD,:],Lrg_d.DD[JD],ENVR.Dthresh[JD],DY);
		Lrg_p.S[JD,:], Lrg_p.DD[JD] = sub_kflag(Lrg_p.S[JD,:],Lrg_p.DD[JD],ENVR.Dthresh[JD],DY);

		#! energy available for somatic growth nu
		Sml_f.nu[JD] = sub_nu(Sml_f.I[JD],Sml_f.bio[JD],Sml_f.met[JD])
		Sml_p.nu[JD] = sub_nu(Sml_p.I[JD],Sml_p.bio[JD],Sml_p.met[JD])
		Sml_d.nu[JD] = sub_nu(Sml_d.I[JD],Sml_d.bio[JD],Sml_d.met[JD])
		Med_f.nu[JD] = sub_nu(Med_f.I[JD],Med_f.bio[JD],Med_f.met[JD])
		Med_p.nu[JD] = sub_nu(Med_p.I[JD],Med_p.bio[JD],Med_p.met[JD])
		Med_d.nu[JD] = sub_nu(Med_d.I[JD],Med_d.bio[JD],Med_d.met[JD])
		Lrg_p.nu[JD] = sub_nu(Lrg_p.I[JD],Lrg_p.bio[JD],Lrg_p.met[JD])
		Lrg_d.nu[JD] = sub_nu(Lrg_d.I[JD],Lrg_d.bio[JD],Lrg_d.met[JD])

		#! maturation (note subscript on Kappa is larvae, juv, adult)
		Sml_f.gamma[JD] = sub_gamma(K_l,Z_s,Sml_f.nu[JD],Sml_f.die[JD],Sml_f.bio[JD],Sml_f.S[JD,DY])
		Sml_p.gamma[JD] = sub_gamma(K_l,Z_s,Sml_p.nu[JD],Sml_p.die[JD],Sml_p.bio[JD],Sml_p.S[JD,DY])
		Sml_d.gamma[JD] = sub_gamma(K_l,Z_s,Sml_d.nu[JD],Sml_d.die[JD],Sml_d.bio[JD],Sml_d.S[JD,DY])
		Med_f.gamma[JD] = sub_gamma(K_a,Z_m,Med_f.nu[JD],Med_f.die[JD],Med_f.bio[JD],Med_f.S[JD,DY])
		Med_p.gamma[JD] = sub_gamma(K_j,Z_m,Med_p.nu[JD],Med_p.die[JD],Med_p.bio[JD],Med_p.S[JD,DY])
		Med_d.gamma[JD] = sub_gamma(K_j,Z_m,Med_d.nu[JD],Med_d.die[JD],Med_d.bio[JD],Med_d.S[JD,DY])
		Lrg_p.gamma[JD] = sub_gamma(K_a,Z_l,Lrg_p.nu[JD],Lrg_p.die[JD],Lrg_p.bio[JD],Lrg_p.S[JD,DY])
		Lrg_d.gamma[JD] = sub_gamma(K_a,Z_l,Lrg_d.nu[JD],Lrg_d.die[JD],Lrg_d.bio[JD],Lrg_d.S[JD,DY])

		#! egg production (by med and large size classes only)
		Sml_f.rep[JD],Sml_f.egg[JD] = sub_rep(Sml_f.nu[JD],K_l,Sml_f.S[JD,DY],Sml_f.egg[JD])
		Sml_p.rep[JD],Sml_p.egg[JD] = sub_rep(Sml_p.nu[JD],K_l,Sml_p.S[JD,DY],Sml_p.egg[JD])
		Sml_d.rep[JD],Sml_d.egg[JD] = sub_rep(Sml_d.nu[JD],K_l,Sml_d.S[JD,DY],Sml_d.egg[JD])
		Med_f.rep[JD],Med_f.egg[JD] = sub_rep(Med_f.nu[JD],K_a,Med_f.S[JD,DY],Med_f.egg[JD])
		Med_p.rep[JD],Med_p.egg[JD] = sub_rep(Med_p.nu[JD],K_j,Med_p.S[JD,DY],Med_p.egg[JD])
		Med_d.rep[JD],Med_d.egg[JD] = sub_rep(Med_d.nu[JD],K_j,Med_d.S[JD,DY],Med_d.egg[JD])
		Lrg_p.rep[JD],Lrg_p.egg[JD] = sub_rep(Lrg_p.nu[JD],K_a,Lrg_p.S[JD,DY],Lrg_p.egg[JD])
		Lrg_d.rep[JD],Lrg_d.egg[JD] = sub_rep(Lrg_d.nu[JD],K_a,Lrg_d.S[JD,DY],Lrg_d.egg[JD])

		#! recruitment (from smaller size class)
		Sml_f.rec[JD] = sub_rec(Med_f.rep[JD],Med_f.bio[JD])
		Sml_p.rec[JD] = sub_rec(Lrg_p.rep[JD],Lrg_p.bio[JD])
		Sml_d.rec[JD] = sub_rec(Lrg_d.rep[JD],Lrg_d.bio[JD])
		Med_f.rec[JD] = sub_rec(Sml_f.gamma[JD],Sml_f.bio[JD])
		Med_p.rec[JD] = sub_rec(Sml_p.gamma[JD],Sml_p.bio[JD])
		Med_d.rec[JD] = sub_rec(Sml_d.gamma[JD],Sml_d.bio[JD])
		Lrg_p.rec[JD] = sub_rec(Med_p.gamma[JD],Med_p.bio[JD])
		Lrg_d.rec[JD] = sub_rec(Med_d.gamma[JD],Med_d.bio[JD])

		#! Mass balance
		BENT.mass[JD] = sub_update_be(BENT.mass[JD],[Med_d.con_be[JD],Lrg_d.con_be[JD]],[Med_d.bio[JD],Lrg_d.bio[JD]])

		Sml_f.bio[JD] = sub_update_fi(Sml_f.bio[JD],Sml_f.rec[JD],Sml_f.nu[JD],
								   Sml_f.rep[JD],Sml_f.gamma[JD],Sml_f.die[JD],Sml_f.egg[JD])
		Sml_p.bio[JD] = sub_update_fi(Sml_p.bio[JD],Sml_p.rec[JD],Sml_p.nu[JD],
								   Sml_p.rep[JD],Sml_p.gamma[JD],Sml_p.die[JD],Sml_p.egg[JD])
		Sml_d.bio[JD] = sub_update_fi(Sml_d.bio[JD],Sml_d.rec[JD],Sml_d.nu[JD],
								   Sml_d.rep[JD],Sml_d.gamma[JD],Sml_d.die[JD],Sml_d.egg[JD])

		Med_f.bio[JD] = sub_update_fi(Med_f.bio[JD],Med_f.rec[JD],Med_f.nu[JD],
								   Med_f.rep[JD],Med_f.gamma[JD],Med_f.die[JD],Med_f.egg[JD])
		Med_p.bio[JD] = sub_update_fi(Med_p.bio[JD],Med_p.rec[JD],Med_p.nu[JD],
								   Med_p.rep[JD],Med_p.gamma[JD],Med_p.die[JD],Med_p.egg[JD])
		Med_d.bio[JD] = sub_update_fi(Med_d.bio[JD],Med_d.rec[JD],Med_d.nu[JD],
								   Med_d.rep[JD],Med_d.gamma[JD],Med_d.die[JD],Med_d.egg[JD])

		Lrg_p.bio[JD] = sub_update_fi(Lrg_p.bio[JD],Lrg_p.rec[JD],Lrg_p.nu[JD],
								   Lrg_p.rep[JD],Lrg_p.gamma[JD],Lrg_p.die[JD],Lrg_p.egg[JD])
	 	Lrg_d.bio[JD] = sub_update_fi(Lrg_d.bio[JD],Lrg_d.rec[JD],Lrg_d.nu[JD],
								   Lrg_d.rep[JD],Lrg_d.gamma[JD],Lrg_d.die[JD],Lrg_d.egg[JD])

	end

	#! Fishing
	#PISC.bio,PLAN.bio,DETR.bio = sub_fishing(PISC.bio,PLAN.bio,DETR.bio,GRD_A);

	#! Forward Euler checks for demographics and movement
	sub_check!(Sml_f.bio);
	sub_check!(Sml_p.bio);
	sub_check!(Sml_d.bio);
	sub_check!(Med_f.bio);
	sub_check!(Med_p.bio);
	sub_check!(Med_d.bio);
	sub_check!(Lrg_p.bio);
	sub_check!(Lrg_d.bio);

end
