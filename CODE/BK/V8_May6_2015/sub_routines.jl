###! Get COBALT data
function get_COBALT!(COBALT,ID,DY,ENV_Tp,ENV_Tb,ENV_Zm,ENV_Zl,ENV_dZm,ENV_dZl,ENV_det)
    ## Get data
    ENV_Tp[:,1]  = COBALT["Tp"][ID,DY]
    ENV_Tb[:,1]  = COBALT["Tb"][ID,DY]
    ENV_Zm[:,1]  = COBALT["Zm"][ID,DY]
    ENV_Zl[:,1]  = COBALT["Zl"][ID,DY]
    ENV_det[:,1] = COBALT["det"][ID,DY]
    ENV_dZm[:,1] = COBALT["dZm"][ID,DY]
    ENV_dZl[:,1] = COBALT["dZl"][ID,DY]
end

###! Temperature multiplier for metabolism
function sub_tmet(tmet,ENV_Tp)
	@devec tmet = exp(0.0548.*ENV_Tp)
end


###! Metabolism
function sub_metabolism_pi!(met,tmet)
    for i = 1:PI_N
        met[i] = PI_bas[i] * tmet * PI_act[i] * 5.258
    end
    nothing
end
function sub_metabolism_pl!(met,tmet)
    for i = 1:PL_N
        met[i] = PL_bas[i] * tmet * PL_act[i] * 5.258
    end
    nothing
end
function sub_metabolism_de!(met,tmet)
    for i = 1:DE_N
        met[i] = DE_bas[i] * tmet * DE_act[i] * 5.258
    end
    nothing
end


###! Fraction of time spent in pelagic (for piscivore)
function sub_tdif(bio1,bio2)
	@devec tdif = sum(bio1) ./ (sum(bio1).+sum(bio2))
end


###! Handling times
function sub_tau_pi!(tau,met)
	for i = 1:PI_N
		tau[i] = 1 / (4*met[i])
	end
	nothing
end
function sub_tau_pl!(tau,met)
	for i = 1:PL_N
		tau[i] = 1 / (4*met[i])
	end
	nothing
end
function sub_tau_de!(tau,met)
	for i = 1:DE_N
		tau[i] = 1 / (4*met[i])
	end
	nothing
end


###! Piscivore encounter rates
function sub_enc_pipi!(enc,prey,tdif)
	for i = 1:PI_N # pred
		for j = 1:PI_N # prey
			enc[j,i] = prey[j]*PI_phi_PI[j,i]*PI_a[i]*tdif
		end
	end
	nothing
end
function sub_enc_pipl!(enc,prey,tdif)
	for i = 1:PI_N # pred
		for j = 1:PL_N # prey
			enc[j,i] = prey[j]*PI_phi_PL[j,i]*PI_a[i]*tdif
		end
	end
	nothing
end
function sub_enc_pide!(enc,prey,tdif)
	for i = 1:PI_N # pred
		for j = 1:DE_N # prey
			enc[j,i] = prey[j]*PI_phi_DE[j,i]*PI_a[i]*(1-tdif)
		end
	end
	nothing
end
function sub_enc_piz!(enc,Zm,Zl,tdif)
	for i = 1:PI_N # pred
		enc[1,i] = Zm*PI_phi_Z[1,i]*PI_a[i]*tdif
		enc[2,i] = Zl*PI_phi_Z[2,i]*PI_a[i]*tdif
	end
	nothing
end

###! Planktivore encounter rates
function sub_enc_plz!(enc,Zm,Zl)
	for i = 1:PL_N # pred
		enc[1,i] = Zm*PL_phi_Z[1,i]*PL_a[i]
		enc[2,i] = Zl*PL_phi_Z[2,i]*PL_a[i]
	end
	nothing
end

###! Planktivore encounter rates
function sub_enc_dede!(enc,prey)
	for i = 1:DE_N # pred
		for j = 1:DE_N # prey
			enc[j,i] = prey[j]*DE_phi_DE[j,i]*DE_a[i]
		end
	end
	nothing
end
function sub_enc_dew!(enc,prey)
	for i = 1:DE_N # pred
		enc[i] = prey[1]*DE_phi_W[i]*DE_a[i]
	end
	nothing
end


###! reset consumption and mortality
function sub_reset_con_fish(I_pi::Array{Any},I_pl::Array{Any},I_de::Array{Any})
	I_pi = I_pi .* 0.0
	I_pl = I_pl .* 0.0
	I_de = I_de .* 0.0
	return I_pi,I_pl,I_de
end
function sub_reset_con_zoo(I_pi::Array{Any},I_pl::Array{Any})
	I_pi = I_pi .* 0.0
	I_pl = I_pl .* 0.0
	return I_pi,I_pl
end
function sub_reset_mort(d_pi::Array{Any},d_pl::Array{Any},d_de::Array{Any},d_w::Array{Any})
	d_pi = d_pi .* 0.0
	d_pl = d_pl .* 0.0
	d_de = d_de .* 0.0
	d_w  = d_w .* 0.0
	return d_pi, d_pl, d_de, d_w
end


###! piscivore consumption
function sub_consume_pipi!(I_pi,d_pi,bio_pi,enc_pipi,ENC_pi,tau_pi)
	for i = 1:PI_N #pred 
		for j = 1:PI_N #prey
			con = (bio_pi[i]*enc_pipi[j,i]) / (1 + (tau_pi[i]*ENC_pi[i]))
			I_pi[i] = d_pi[j] += con
		end
	end
end
function sub_consume_pipl!(I_pi,d_pl,bio_pi,enc_pipl,ENC_pi,tau_pi)
	for i = 1:PI_N #pred 
		for j = 1:PL_N #prey
			con = (bio_pi[i]*enc_pipl[j,i]) / (1 + (tau_pi[i]*ENC_pi[i]))
			I_pi[i] = d_pl[j] += con
		end
	end
end
function sub_consume_pide!(I_pi,d_de,bio_pi,enc_pide,ENC_pi,tau_pi)
	for i = 1:PI_N #pred 
		for j = 1:DE_N #prey
			con = (bio_pi[i]*enc_pide[j,i]) / (1 + (tau_pi[i]*ENC_pi[i]))
			I_pi[i] = d_de[j] += con
		end
	end
end

###! Detrivore consumption
function sub_consume_dede!(I_de,d_de,bio_de,enc_dede,ENC_de,tau_de)
	for i = 1:DE_N #pred 
		for j = 1:DE_N #prey
			con = (bio_de[i]*enc_dede[j,i]) / (1 + (tau_de[i]*ENC_de[i]))
			I_de[i] = d_de[j] += con
		end
	end
end
function sub_consume_dew!(I_de,d_w,bio_de,enc_dew,ENC_de,tau_de)
	for i = 1:DE_N #pred 
		con = (bio_de[i]*enc_dew[i]) / (1 + (tau_de[i]*ENC_de[i]))
		I_de[i] = d_w[1] += con
	end
end

###! Consumption of zooplankton
#! by piscivore
function sub_consume_piz!(I_z,bio_pi,enc_piz,ENC_pi,tau_pi)
	for i = 1:PI_N #pred 
		for j = 1:2 #prey
			con = (bio_pi[i]*enc_piz[j,i]) / (1 + (tau_pi[i]*ENC_pi[i]))
			I_z[j,i] += con
		end
	end
end
function sub_consume_plz!(I_z,bio_pl,enc_plz,ENC_pl,tau_pl)
	for i = 1:PL_N #pred 
		for j = 1:2 #prey
			con = (bio_pl[i]*enc_plz[j,i]) / (1 + (tau_pl[i]*ENC_pl[i]))
			I_z[j,i] += con
		end
	end
end

###! Offline coupling
function sub_offline!(pi_z,pl_z,dZm,dZl)
	ZOO = sum(pi_z,2) + sum(pl_z,2)
	if ZOO[1] > dZm
		pi_z[1,:] = (pi_z[1,:]./ZOO[1]) .* dZm
		pl_z[1,:] = (pl_z[1,:]./ZOO[1]) .* dZm
	end
	if ZOO[2] > dZl
		pi_z[2,:] = (pi_z[2,:]./ZOO[2]) .* dZl
		pl_z[2,:] = (pl_z[2,:]./ZOO[2]) .* dZl
	end
end


###! ENERGY AVAILABLE FOR GROWTH NU
function sub_nu_pi!(nu,bio,I,met)
	for i = 1:PI_N
		nu[i] = (I[i]/bio[i]*PI_lambda[i]) - met[i]
	end
end
function sub_nu_pl!(nu,bio,I,met)
	for i = 1:PL_N
		nu[i] = (I[i]/bio[i]*PL_lambda[i]) - met[i]
	end
end
function sub_nu_de!(nu,bio,I,met)
	for i = 1:DE_N
		nu[i] = (I[i]/bio[i]*DE_lambda[i]) - met[i]
	end
end


###! ENERGY AVAILABLE FOR SOMATIC GROWTH
function sub_gamma_pi!(gamma,nu,bio,d)
    # note: divide by bio to get biomass specific units
    for i = 1:PI_N
		gg = ((PI_K[i]*nu[i]) - (d[i]/bio[i]))/(1-(PI_z[i]^(1-((d[i]/bio[i])/(PI_K[i]*nu[i])))))
		if gg < 0 || isnan(gg)==true
			gamma[i] = 0.0
		else
			gamma[i] = gg
		end
	end
end
function sub_gamma_pl!(gamma,nu,bio,d)
    # note: divide by bio to get biomass specific units
    for i = 1:PL_N
		gg = ((PL_K[i]*nu[i]) - (d[i]/bio[i]))/(1-(PL_z[i]^(1-((d[i]/bio[i])/(PL_K[i]*nu[i])))))
		if gg < 0 || isnan(gg)==true
			gamma[i] = 0.0
		else
			gamma[i] = gg
		end
	end
end
function sub_gamma_de!(gamma,nu,bio,d)
    # note: divide by bio to get biomass specific units
    for i = 1:DE_N
		gg = ((DE_K[i]*nu[i]) - (d[i]/bio[i]))/(1-(DE_z[i]^(1-((d[i]/bio[i])/(DE_K[i]*nu[i])))))
		if gg < 0 || isnan(gg)==true
			gamma[i] = 0.0
		else
			gamma[i] = gg
		end
	end
end


###! TOTAL BIOMASS MADE FROM REPRODUCTION
function sub_rep_pi!(rep,nu,bio)
	for i = 1:PI_N
		if nu[i] > 0.
			rep[i] = (1-PI_K[i]) * nu[i] * bio[i]
		else
			rep[i] = 0.
			#nu[i] = 0.
		end
	end
end
function sub_rep_pl!(rep,nu,bio)
	for i = 1:PL_N
		if nu[i] > 0.
			rep[i] = (1-PL_K[i]) * nu[i] * bio[i]
		else
			rep[i] = 0.
			#nu[i] = 0.
		end
	end
end
function sub_rep_de!(rep,nu,bio)
	for i = 1:DE_N
		if nu[i] > 0.
			rep[i] = (1-DE_K[i]) * nu[i] * bio[i]
		else
			rep[i] = 0.
			#nu[i] = 0.
		end
	end
end


###! TOTAL BIOMASS SOMATIC GROWTH
function sub_grw_pi!(grw,nu,bio)
	for i = 1:PI_N
		grw[i] = nu[i] * bio[i]
	end
end
function sub_grw_pl!(grw,nu,bio)
	for i = 1:PL_N
		grw[i] = nu[i] * bio[i]
	end
end
function sub_grw_de!(grw,nu,bio)
	for i = 1:DE_N
		grw[i] = nu[i] * bio[i]
	end
end


###! TOTAL BIOMASS MATURING
function sub_mat_pi!(mat,gam,bio)
	for i = 1:PI_N
		mat[i] = gam[i] * bio[i]
	end
end
function sub_mat_pl!(mat,gam,bio)
	for i = 1:PL_N
		mat[i] = gam[i] * bio[i]
	end
end
function sub_mat_de!(mat,gam,bio)
	for i = 1:DE_N
		mat[i] = gam[i] * bio[i]
	end
end


###! TOTAL BIOMASS dying from background mortality
function sub_mrt_pi!(mrt,bio)
	for i = 1:PI_N
		mrt[i] = PI_mrt[i] * bio[i]
	end
end
function sub_mrt_pl!(mrt,bio)
	for i = 1:PL_N
		mrt[i] = PL_mrt[i] * bio[i]
	end
end
function sub_mrt_de!(mrt,bio)
	for i = 1:DE_N
		mrt[i] = DE_mrt[i] * bio[i]
	end
end


###! Update biomass
function sub_update_pi!(bio,rep,grw,mat,d,mrt)
	bio[1] += (sum(rep) + grw[1] - mat[1] - d[1] - mrt[1]) * DT
	bio[2:PI_N] += (mat[1:PI_N-1] + grw[2:PI_N] - mrt[2:PI_N] - 
					rep[2:PI_N] - mat[2:PI_N] - d[2:PI_N]) * DT
end
function sub_update_pl!(bio,rep,grw,mat,d,mrt)
	bio[1] += (sum(rep) + grw[1] - mat[1] - d[1] - mrt[1]) * DT
	bio[2:PL_N] += (mat[1:PL_N-1] + grw[2:PL_N] - mrt[2:PL_N] - 
					rep[2:PL_N] - mat[2:PL_N] - d[2:PL_N]) * DT
end
function sub_update_de!(bio,rep,grw,mat,d,mrt)
	bio[1] += (sum(rep) + grw[1] - mat[1] - d[1] - mrt[1]) * DT
	bio[2:DE_N] += (mat[1:DE_N-1] + grw[2:DE_N] - mrt[2:DE_N] - 
					rep[2:DE_N] - mat[2:DE_N] - d[2:DE_N]) * DT
end
function sub_update_w!(bio,d,det)
    bio[1] += (det - d[1] - (bio[1]*0.01)) * DT # with sedimentation
end


###! Forward Euler checks
function sub_check_pi!(bio)
	for i = 1:PI_N
		if bio[i] <= 0.0
			bio[i] = eps()
		end
	end
end
function sub_check_pl!(bio)
	for i = 1:PL_N
		if bio[i] <= 0.0
			bio[i] = eps()
		end
	end
end
function sub_check_de!(bio)
	for i = 1:DE_N
		if bio[i] <= 0.0
			bio[i] = eps()
		end
	end
end
function sub_check_w!(bio)
	if bio[1] <= 0.0
		bio[1] = eps()
	end
end


		
