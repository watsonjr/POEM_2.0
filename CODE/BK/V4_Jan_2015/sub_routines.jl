
###! Metabolism
function sub_metabolism(PISC,PLAN,DETR,PRM_PI,PRM_PL,PRM_DE,TEMP_p)
	@inbounds for i = 1:PRM_PL.N
		fnc_met!(PISC.met,i,PRM_PI.s[i],TEMP_p)
	end
	@inbounds for i = 1:PRM_PI.N
		fnc_met!(PLAN.met,i,PRM_PL.s[i],TEMP_p)
	end
	@inbounds for i = 1:PRM_DE.N
		## change this when changed bottom temp to account for pressure
		fnc_met!(DETR.met,i,PRM_DE.s[i],TEMP_p)
	end
	return PISC,PLAN,DETR
end

###! Encounter rate
function sub_enc(PISC,PLAN,DETR,W,ZOO,PRM_PI,PRM_PL,PRM_DE)
	@inbounds for i = 1:PRM_PI.N
		fnc_enc!(PISC.enc["pipi"],i,PISC.bio,PRM_PI.Phi_PI[:,i],PRM_PI.a[i])
		fnc_enc!(PISC.enc["pipl"],i,PLAN.bio,PRM_PI.Phi_PL[:,i],PRM_PI.a[i])
		fnc_enc!(PISC.enc["pide"],i,DETR.bio,PRM_PI.Phi_DE[:,i],PRM_PI.a[i])
		fnc_enc!(PISC.enc["piz"],i,ZOO,PRM_PI.Phi_Z[:,i],PRM_PI.a[i])
	end
	@inbounds for i = 1:PRM_PL.N
		fnc_enc!(PLAN.enc["plz"],i,ZOO,PRM_PL.Phi_Z[:,i],PRM_PL.a[i])
	end
	@inbounds for i = 1:PRM_DE.N
		fnc_enc!(DETR.enc["dew"],i,W.bio,PRM_DE.Phi_W[:,i],PRM_DE.a[i])
		fnc_enc!(DETR.enc["dede"],i,DETR.bio,PRM_DE.Phi_DE[:,i],PRM_DE.a[i])
	end
	return PISC,PLAN,DETR
end

###! Total biomass encountered
function sub_ENC(PISC,PLAN,DETR)
    PISC.ENC = vec(sum(PISC.enc["pipi"],1))
              	   + vec(sum(PISC.enc["pipl"],1))
              	   + vec(sum(PISC.enc["pide"],1))
    PLAN.ENC = vec(sum(PLAN.enc["plz"],1))
    DETR.ENC = vec(DETR.enc["dew"]) + vec(sum(DETR.enc["dede"],1))
    return PISC,PLAN,DETR
end

###! Handling times
function sub_tau(PISC,PLAN,DETR,PRM_PI,PRM_PL,PRM_DE)
	@inbounds for i = 1:PRM_PI.N
		fnc_tau!(PISC.tau,i,PISC.met[i])
	end
	@inbounds for i = 1:PRM_PL.N
		fnc_tau!(PLAN.tau,i,PLAN.met[i])
	end
	@inbounds for i = 1:PRM_DE.N
		fnc_tau!(DETR.tau,i,DETR.met[i])
	end
	return PISC,PLAN,DETR
end


###! Consumption
function sub_consume(PISC,PLAN,DETR,W,PRM_PI,PRM_PL,PRM_DE,ZOO,DZc)
	# Calculates TOTAL biomass consumerd g d-1
	PISC.I = zeros(PRM_PI.N)
	PLAN.I = zeros(PRM_PL.N)
	DETR.I = zeros(PRM_DE.N)
	PISC.d = zeros(PRM_PI.N)
	PLAN.d = zeros(PRM_PL.N)
	DETR.d = zeros(PRM_DE.N)
	DZs    = zeros(2) ## zooplankton mortality for offline coupling

	#! piscivore
	@simd for i = 1:PRM_PI.N # pred
		@inbounds for j = 1:PRM_PI.N # prey
		fnc_cons!(PISC.I,PISC.d,i,j,PISC.bio[i],PISC.bio[j],PISC.enc["pipi"][j,i],PISC.ENC[i],PISC.tau[i])
		end
		@inbounds for j = 1:PRM_PL.N
		fnc_cons!(PISC.I,PLAN.d,i,j,PISC.bio[i],PLAN.bio[j],PISC.enc["pipl"][j,i],PISC.ENC[i],PISC.tau[i])
		end
		@inbounds for j = 1:PRM_DE.N
		fnc_cons!(PISC.I,DETR.d,i,j,PISC.bio[i],DETR.bio[j],PISC.enc["pide"][i,j],PISC.ENC[i],PISC.tau[i])
		end
	end

	#! Planktivore
	#! heavy due to COBALT offline coupling
	dzm = zeros(1); dzl = zeros(1);
	con_zm = zeros(PRM_PL.N)
	con_zl = zeros(PRM_PL.N)
	@inbounds for i = 1:PRM_PL.N # pred
		fnc_cons!(con_zm,dzm,i,1,PLAN.bio[i],ZOO[1],PLAN.enc["plz"][1,i],PLAN.ENC[i],PLAN.tau[i])
		fnc_cons!(con_zl,dzl,i,1,PLAN.bio[i],ZOO[2],PLAN.enc["plz"][2,i],PLAN.ENC[i],PLAN.tau[i])
	end
	## Add offline coupling here
	if dzm[1] > DZc[1]
		con_zm = DZc[1] .* (con_zm./sum(con_zm))
	end
	if dzl[1] > DZc[2]
		con_zl = DZc[2] .* (con_zl./sum(con_zl))
	end
	PLAN.I += con_zm + con_zl

	#! Detritivore
	@simd for i = 1:PRM_DE.N # pred
		# detritus feeding
		fnc_cons!(DETR.I,W.d,i,1,DETR.bio[i],W.bio[1],DETR.enc["dew"][i],DETR.ENC[i],DETR.tau[i])
		# canabalism
		@inbounds for j = 1:PRM_DE.N # prey
		fnc_cons!(DETR.I,DETR.d,i,j,DETR.bio[i],DETR.bio[j],DETR.enc["dede"][j,i],DETR.ENC[i],DETR.tau[i])
		end
	end
	return PISC,PLAN,DETR,W
end

###! Nu: total energy avail for somatic growht and reproduction
function sub_nu(PISC,PLAN,DETR,PRM_PI,PRM_PL,PRM_DE)
	@inbounds for i = 1:PRM_PI.N
		fnc_nu!(PISC.nu,i,PISC.bio[i],PISC.I[i],PISC.met[i],PRM_PI.lambda[i])
	end
	@inbounds for i = 1:PRM_PL.N
		fnc_nu!(PLAN.nu,i,PLAN.bio[i],PLAN.I[i],PLAN.met[i],PRM_PL.lambda[i])
	end
	@inbounds for i = 1:PRM_DE.N
		fnc_nu!(DETR.nu,i,DETR.bio[i],DETR.I[i],DETR.met[i],PRM_DE.lambda[i])
	end
	return PISC,PLAN,DETR
end

###! Gamma: energy avail for somatic growth
function sub_gamma(PISC,PLAN,DETR,PRM_PI,PRM_PL,PRM_DE)
	@inbounds for i = 1:PRM_PI.N
		fnc_gamma!(PISC.gamma,i,PISC.nu[i],PISC.bio[i],PISC.d[i],PRM_PI.K[i],PRM_PI.z[i])
	end
	@inbounds for i = 1:PRM_PL.N
		fnc_gamma!(PLAN.gamma,i,PLAN.nu[i],PLAN.bio[i],PLAN.d[i],PRM_PL.K[i],PRM_PL.z[i])
	end
	@inbounds for i = 1:PRM_DE.N
		fnc_gamma!(DETR.gamma,i,DETR.nu[i],DETR.bio[i],DETR.d[i],PRM_DE.K[i],PRM_DE.z[i])
	end
	return PISC,PLAN,DETR
end


###! Biomass reproduced
function sub_rep(PISC,PLAN,DETR,PRM_PI,PRM_PL,PRM_DE)
	#! piscivore
	@inbounds for i = 1:PRM_PI.N
		fnc_rep!(PISC.REP,i,PISC.nu[i],PRM_PI.K[i],PISC.bio[i])
	end
	#! planktivore
	@inbounds for i = 1:PRM_PL.N
		fnc_rep!(PLAN.REP,i,PLAN.nu[i],PRM_PL.K[i],PLAN.bio[i])
	end
	#! detritivore
	@inbounds for i = 1:PRM_DE.N
		fnc_rep!(DETR.REP,i,DETR.nu[i],PRM_DE.K[i],DETR.bio[i])
	end
	return PISC,PLAN,DETR
end

###! Somatic growth
function sub_grw(PISC,PLAN,DETR,PRM_PI,PRM_PL,PRM_DE)
	@inbounds for i = 1:PRM_PI.N
		fnc_grw!(PISC.GRW,i,PISC.nu[i],PISC.bio[i])
	end
	@inbounds for i = 1:PRM_PL.N
		fnc_grw!(PLAN.GRW,i,PLAN.nu[i],PLAN.bio[i])
	end
	@inbounds for i = 1:PRM_DE.N
		fnc_grw!(DETR.GRW,i,DETR.nu[i],DETR.bio[i])
	end
	return PISC,PLAN,DETR
end

###! Maturation
function sub_mat(PISC,PLAN,DETR,PRM_PI,PRM_PL,PRM_DE)
	@inbounds for i = 1:PRM_PI.N
		fnc_mat!(PISC.MAT,i,PISC.gamma[i],PISC.bio[i])
	end
	@inbounds for i = 1:PRM_PL.N
		fnc_mat!(PLAN.MAT,i,PLAN.gamma[i],PLAN.bio[i])
	end
	@inbounds for i = 1:PRM_DE.N
		fnc_mat!(DETR.MAT,i,DETR.gamma[i],DETR.bio[i])
	end
	return PISC,PLAN,DETR
end

###! Update biomass
function sub_update_bio(PISC,PLAN,DETR,W,PRM_PI,PRM_PL,PRM_DE)
	PISC.bio[1] += (sum(PISC.REP) + PISC.GRW[1] - PISC.MAT[1] - PISC.d[1])
	PISC.bio[2:end] += (PISC.MAT[1:end-1] + PISC.GRW[2:end]
				- PISC.REP[2:end] - PISC.MAT[2:end] - PISC.d[2:end])

	PLAN.bio[1] += (sum(PLAN.REP) + PLAN.GRW[1] - PLAN.MAT[1] - PLAN.d[1])
	PLAN.bio[2:end] += (PLAN.MAT[1:end-1] + PLAN.GRW[2:end]
				- PLAN.REP[2:end] - PLAN.MAT[2:end] - PLAN.d[2:end])

	DETR.bio[1] += (sum(DETR.REP) + DETR.GRW[1] - DETR.MAT[1] - DETR.d[1])
	DETR.bio[2:end] += (DETR.MAT[1:end-1] + DETR.GRW[2:end]
				- DETR.REP[2:end] - DETR.MAT[2:end] - DETR.d[2:end])

	W.bio[1] += (W.I[1] - W.d[1])

	return PISC,PLAN,DETR,W
end


###! Forward Euler checks
function sub_euler_checks(PISC,PLAN,DETR,W,PRM_PI,PRM_PL,PRM_DE)
	for i = 1:PRM_PI.N
		if PISC.bio[i] <= 0.
			PISC.bio[i] = eps()
		end
	end
	for i = 1:PRM_PL.N
		if PLAN.bio[i] <= 0.
			PLAN.bio[i] = eps()
		end
	end
	for i = 1:PRM_DE.N
		if DETR.bio[i] <= 0.
			DETR.bio[i] = eps()
		end
	end
	if W.bio[1] <= 0.
		W.bio[1] = eps()
	end
	return PISC,PLAN,DETR,W
end



