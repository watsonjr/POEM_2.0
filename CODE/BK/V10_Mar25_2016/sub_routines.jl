
#### THE MODEL
###! DEMOGRAPHIC CALCULATIONS
function sub_futbio!(ID,DY,COBALT,ENVR,PISC,PLAN,DETR,BENT)
	###! COBALT information
	get_COBALT!(COBALT,ID,DY,ENVR)

	#! temperature multiplier
	PISC.tmet = map(sub_tmet,PISC.tmet,ENVR.Tp);
	PLAN.tmet = map(sub_tmet,PLAN.tmet,ENVR.Tp);
	DETR.tmet = map(sub_tmet,DETR.tmet,ENVR.Tb);

	#! activity multiplier
	map(sub_umet_pi!,PISC.umet);
	map(sub_umet_pl!,PLAN.umet);
	map(sub_umet_de!,DETR.umet);

	#! metabolism
	map(sub_metabolism_pi!,PISC.met,PISC.tmet,PISC.umet);
	map(sub_metabolism_pl!,PLAN.met,PLAN.tmet,PLAN.umet);
	map(sub_metabolism_de!,DETR.met,DETR.tmet,DETR.umet);

	#! fraction of time piscivore spends in pelagic
	PISC.tdif = map(sub_tdif,GRD_Z[ID],PLAN.bio,DETR.bio);

	#! handling times
	map(sub_tau_pi!,PISC.tau,PISC.met);
	map(sub_tau_pl!,PLAN.tau,PLAN.met);
	map(sub_tau_de!,DETR.tau,DETR.met);

	#! encounter rates
	map(sub_enc_pipi!,PISC.enc_pi,PISC.bio,PISC.tdif);
	map(sub_enc_pipl!,PISC.enc_pl,PLAN.bio,PISC.tdif);
	map(sub_enc_pide!,PISC.enc_de,DETR.bio,PISC.tdif);
	map(sub_enc_piz!,PISC.enc_z,ENVR.Zm,ENVR.Zl,PISC.tdif);
	map(sub_enc_plz!,PLAN.enc_z,ENVR.Zm,ENVR.Zl);
	map(sub_enc_dede!,DETR.enc_de,DETR.bio);
	map(sub_enc_debe!,DETR.enc_be,BENT.bio);

	#! total biomass encountered of each group
	map(sub_ENC_pi!,PISC.ENC,PISC.enc_pi,PISC.enc_pl,PISC.enc_de,PISC.enc_z);
	map(sub_ENC_pl!,PLAN.ENC,PLAN.enc_z);
	map(sub_ENC_de!,DETR.ENC,DETR.enc_de,DETR.enc_be);

	#! reset consumption and mortality
	map(sub_reset_pisc!,PISC.I,PISC.I_z,PISC.d)
	map(sub_reset_plan!,PLAN.I,PLAN.I_z,PLAN.d)
	map(sub_reset_detr!,DETR.I,DETR.d)

	#! total biomass consumed and lost to predation SPEED UP
	map(sub_consume_pipi!,PISC.I,PISC.d,PISC.bio,PISC.enc_pi,PISC.ENC,PISC.tau);
	map(sub_consume_pipl!,PISC.I,PLAN.d,PISC.bio,PISC.enc_pl,PISC.ENC,PISC.tau);
	map(sub_consume_pide!,PISC.I,DETR.d,PISC.bio,PISC.enc_de,PISC.ENC,PISC.tau);
	map(sub_consume_dede!,DETR.I,DETR.d,DETR.bio,DETR.enc_de,DETR.ENC,DETR.tau);
	map(sub_consume_debe!,DETR.I,BENT.d,DETR.bio,DETR.enc_be,DETR.ENC,DETR.tau);

	#! zooplankton consumption
	map(sub_consume_piz!,PISC.I_z,PISC.bio,PISC.enc_z,PISC.ENC,PISC.tau);
	map(sub_consume_plz!,PLAN.I_z,PLAN.bio,PLAN.enc_z,PLAN.ENC,PLAN.tau);

	#! OFFLINE Coupling
	map(sub_offline!,PISC.I_z,PLAN.I_z,ENVR.dZm,ENVR.dZl)

	#! Add zooplankton consumption of PI and PL
	map(sub_consume_pizoo!,PISC.I,PISC.I_z);
	map(sub_consume_plzoo!,PLAN.I,PLAN.I_z);

	#! energy available for somatic growth nu
	map(sub_nu_pi!,PISC.nu,PISC.bio,PISC.I,PISC.met);
	map(sub_nu_pl!,PLAN.nu,PLAN.bio,PLAN.I,PLAN.met);
	map(sub_nu_de!,DETR.nu,DETR.bio,DETR.I,DETR.met);

	#! maturation
	map(sub_gamma_pi!,PISC.gamma,PISC.nu,PISC.bio,PISC.d);
	map(sub_gamma_pl!,PLAN.gamma,PLAN.nu,PLAN.bio,PLAN.d);
	map(sub_gamma_de!,DETR.gamma,DETR.nu,DETR.bio,DETR.d);

	#! egg production
	map(sub_rep_pi!,PISC.REP,PISC.nu,PISC.bio);
	map(sub_rep_pl!,PLAN.REP,PLAN.nu,PLAN.bio);
	map(sub_rep_de!,DETR.REP,DETR.nu,DETR.bio);

	#! total biomass somatic growth
	map(sub_grw_pi!,PISC.GRW,PISC.nu,PISC.bio);
	map(sub_grw_pl!,PLAN.GRW,PLAN.nu,PLAN.bio);
	map(sub_grw_de!,DETR.GRW,DETR.nu,DETR.bio);

	#! total biomass maturing
	map(sub_mat_pi!,PISC.MAT,PISC.gamma,PISC.bio);
	map(sub_mat_pl!,PLAN.MAT,PLAN.gamma,PLAN.bio);
	map(sub_mat_de!,DETR.MAT,DETR.gamma,DETR.bio);

	#! total biomass lost to natural mortality
	map(sub_mrt_pi!,PISC.MRT,PISC.bio);
	map(sub_mrt_pl!,PLAN.MRT,PLAN.bio);
	map(sub_mrt_de!,DETR.MRT,DETR.bio);

	#! Mass balance	
	map(sub_update_pi!,PISC.bio,PISC.REP,PISC.GRW,PISC.MAT,PISC.d,PISC.MRT);
	map(sub_update_pl!,PLAN.bio,PLAN.REP,PLAN.GRW,PLAN.MAT,PLAN.d,PLAN.MRT);
	map(sub_update_de!,DETR.bio,DETR.REP,DETR.GRW,DETR.MAT,DETR.d,DETR.MRT);
	map(sub_update_be!,BENT.bio,BENT.d,ENVR.det);

	#! Fishing
	PISC.bio,PLAN.bio,DETR.bio = sub_fishing(PISC.bio,PLAN.bio,DETR.bio,GRD_A);

	#! Forward Euler checks for demographics and movement
	map(sub_check_pi!,PISC.bio);
	map(sub_check_pl!,PLAN.bio);
	map(sub_check_de!,DETR.bio);
	map(sub_check_be!,BENT.bio);

	#! Record total biomasses
	tot_pi = sub_totbio(PISC.bio,PI_N)
	tot_pl = sub_totbio(PLAN.bio,PL_N)
	tot_de = sub_totbio(DETR.bio,DE_N)

	#! Advection / Swimming
	PISC.bio = sub_advection(PISC.bio,PISC.nu,PI_N,PI_U,ENVR.U,ENVR.V,GRD["DX"],GRD["DY"])
	PLAN.bio = sub_advection(PLAN.bio,PLAN.nu,PL_N,PL_U,ENVR.U,ENVR.V,GRD["DX"],GRD["DY"])

	#! Diffusion
	PISC.bio = sub_diffuse(PISC.bio,PI_N,PI_a,GRD["DX"],GRD["DY"])
	PLAN.bio = sub_diffuse(PLAN.bio,PL_N,PL_a,GRD["DX"],GRD["DY"])
	DETR.bio = sub_diffuse(DETR.bio,DE_N,DE_a,GRD["DX"],GRD["DY"])

	##! Movement conservation of mass
	#PISC.bio = sub_conserve(PISC.bio,PI_N,tot_pi)
	#PLAN.bio = sub_conserve(PLAN.bio,PL_N,tot_pl)
	#DETR.bio = sub_conserve(DETR.bio,DE_N,tot_de)

end



