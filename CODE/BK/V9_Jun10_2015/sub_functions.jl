###! Get COBALT data
function get_COBALT!(COBALT,ID,DY,ENVR)
    ## Get data
    ENVR.Tp[:,1]  = COBALT["Tp"][ID,DY]
    ENVR.Tb[:,1]  = COBALT["Tb"][ID,DY]
    ENVR.Zm[:,1]  = COBALT["Zm"][ID,DY]
    ENVR.Zl[:,1]  = COBALT["Zl"][ID,DY]
    ENVR.det[:,1] = COBALT["det"][ID,DY]
    ENVR.dZm[:,1] = COBALT["dZm"][ID,DY]
    ENVR.dZl[:,1] = COBALT["dZl"][ID,DY]
end

###! Temperature multiplier for metabolism
function sub_tmet(tmet,T)
	tmet = exp(0.0*T)
	#tmet = exp(0.0548.*T)
end


###! Metabolism
function sub_metabolism_pi!(met::Array{Float64},tmet::Float64)
    for i = 1:PI_N
        met[i] = PI_bas[i] * tmet * 5.258
    end
    nothing
end
function sub_metabolism_pl!(met::Array{Float64},tmet::Float64)
    for i = 1:PL_N
        met[i] = PL_bas[i] * tmet * 5.258
    end
    nothing
end
function sub_metabolism_de!(met::Array{Float64},tmet::Float64)
    for i = 1:DE_N
        met[i] = DE_bas[i] * tmet * 5.258
    end
    nothing
end


###! Fraction of time spent in pelagic (for piscivore)
function sub_tdif(Z,bio1::Array{Float64},bio2::Array{Float64})
	if Z < PI_be_cutoff
		tdif = sum(bio1) ./ (sum(bio1).+sum(bio2))
	else
		tdif = 1.0
	end
end

###! Handling times
function sub_tau_pi!(tau::Array{Float64},met::Array{Float64})
	for i = 1:PI_N
		tau[i] = 1 / (4*met[i])
	end
	nothing
end
function sub_tau_pl!(tau::Array{Float64},met::Array{Float64})
	for i = 1:PL_N
		tau[i] = 1 / (4*met[i])
	end
	nothing
end
function sub_tau_de!(tau::Array{Float64},met::Array{Float64})
	for i = 1:DE_N
		tau[i] = 1 / (4*met[i])
	end
	nothing
end


###! Piscivore encounter rates
function sub_enc_pipi!(enc::Array{Float64},prey::Array{Float64},tdif::Float64)
	for i = 1:PI_N # pred
		for j = 1:PI_N # prey
			enc[j,i] = prey[j]*PI_phi_PI[j,i]*PI_a[i]*tdif
		end
	end
	nothing
end

function sub_enc_pipl!(enc::Array{Float64},prey::Array{Float64},tdif::Float64)
	for i = 1:PI_N # pred
		for j = 1:PL_N # prey
			enc[j,i] = prey[j]*PI_phi_PL[j,i]*PI_a[i]*tdif
		end
	end
	nothing
end

function sub_enc_pide!(enc::Array{Float64},prey::Array{Float64},tdif::Float64)
	for i = 1:PI_N # pred
		for j = 1:DE_N # prey
			enc[j,i] = prey[j]*PI_phi_DE[j,i]*PI_a[i]*(1-tdif)
		end
	end
	nothing
end
function sub_enc_piz!(enc::Array{Float64},Zm,Zl,tdif::Float64)
	for i = 1:PI_N # pred
		enc[1,i] = Zm*PI_phi_Z[1,i]*PI_a[i]*tdif
		enc[2,i] = Zl*PI_phi_Z[2,i]*PI_a[i]*tdif
	end
	nothing
end

###! Planktivore encounter rates
function sub_enc_plz!(enc::Array{Float64},Zm,Zl)
	for i = 1:PL_N # pred
		enc[1,i] = Zm*PL_phi_Z[1,i]*PL_a[i]
		enc[2,i] = Zl*PL_phi_Z[2,i]*PL_a[i]
	end
	nothing
end

###! Planktivore encounter rates
function sub_enc_dede!(enc::Array{Float64},prey::Array{Float64})
	for i = 1:DE_N # pred
		for j = 1:DE_N # prey
			enc[j,i] = prey[j]*DE_phi_DE[j,i]*DE_a[i]
		end
	end
	nothing
end
function sub_enc_debe!(enc::Array{Float64},prey::Array{Float64})
	for i = 1:DE_N # pred
		enc[i] = prey[1]*DE_phi_BE[i]*DE_a[i]
	end
	nothing
end

###! Total biomass encountered
function sub_ENC_pi!(ENC,enc_pi,enc_pl,enc_de,enc_z)
	for i = 1:PI_N
		ENC[i] = 0.0
		for j = 1:PI_N
			ENC[i] += enc_pi[j,i]
		end
		for j = 1:PL_N
			ENC[i] += enc_pl[j,i]
		end
		for j = 1:DE_N
			ENC[i] += enc_de[j,i]
		end
		for j = 1:2
			ENC[i] += enc_z[j,i]
		end
	end
end
function sub_ENC_de!(ENC,enc_de,enc_be)
	for i = 1:DE_N
		ENC[i] = 0.0
		for j = 1:DE_N
			ENC[i] += enc_de[j,i]
		end
		ENC[i] += enc_be[i]
	end
end
function sub_ENC_pl!(ENC,enc_z)
	for i = 1:PL_N
		ENC[i] = 0.0
		ENC[i] += enc_z[1,i]
		ENC[i] += enc_z[2,i]
	end
end


###! reset consumption and mortality
function sub_reset_pisc!(I,Iz,d)
	for i = 1:PI_N
		I[i] = 0.0
		Iz[:,i] = [0.,0.];
		d[i] = 0.0
	end
end
function sub_reset_plan!(I,Iz,d)
	for i = 1:PI_N
		I[i] = 0.0
		Iz[:,i] = [0.,0.];
		d[i] = 0.0
	end
end
function sub_reset_detr!(I,d)
	for i = 1:DE_N
		I[i] = 0.0
		d[i] = 0.0
	end
end


###! type II feeding
function sub_typeII(pred,enc,tau,ENC)
	return (pred*enc) / (1 + (tau*ENC))
end

###! piscivore consumption
function sub_consume_pipi!(I_pi,d_pi,bio_pi,enc_pipi,ENC_pi,tau_pi)
	for i = 1:PI_N, j = 1:PI_N#pred 
		con = sub_typeII(bio_pi[i],enc_pipi[j,i],tau_pi[i],ENC_pi[i])
		I_pi[i] += con
		d_pi[j] += con
	end
end
function sub_consume_pipl!(I_pi,d_pl,bio_pi,enc_pipl,ENC_pi,tau_pi)
	for i = 1:PI_N, j = 1:PL_N #pred 
		con = sub_typeII(bio_pi[i],enc_pipl[j,i],tau_pi[i],ENC_pi[i])
		I_pi[i] += con
		d_pl[j] += con
	end
end
function sub_consume_pide!(I_pi,d_de,bio_pi,enc_pide,ENC_pi,tau_pi)
	for i = 1:PI_N, j=1:DE_N #pred 
		con = sub_typeII(bio_pi[i],enc_pide[j,i],tau_pi[i],ENC_pi[i])
		I_pi[i] += con 
		d_de[j] += con
	end
end

###! Detrivore consumption
function sub_consume_dede!(I_de,d_de,bio_de,enc_dede,ENC_de,tau_de)
	for i = 1:DE_N, j = 1:DE_N #pred 
		con = sub_typeII(bio_de[i],enc_dede[j,i],tau_de[i],ENC_de[i])
		I_de[i] += con
		d_de[j] += con
	end
end
function sub_consume_debe!(I_de,d_be,bio_de,enc_debe,ENC_de,tau_de)
	for i = 1:DE_N #pred 
		con = sub_typeII(bio_de[i],enc_debe[1,i],tau_de[i],ENC_de[i])
		I_de[i] += con
		d_be[1] += con
	end
end

###! Consumption of zooplankton
#! by piscivore
function sub_consume_piz!(I_z,bio_pi,enc_piz,ENC_pi,tau_pi)
	for i = 1:PI_N #pred 
		I_z[1,i] += sub_typeII(bio_pi[i],enc_piz[1,i],tau_pi[i],ENC_pi[i])
		I_z[2,i] += sub_typeII(bio_pi[i],enc_piz[2,i],tau_pi[i],ENC_pi[i])
	end
end
function sub_consume_plz!(I_z,bio_pl,enc_plz,ENC_pl,tau_pl)
	for i = 1:PL_N #pred 
		I_z[1,i] += sub_typeII(bio_pl[i],enc_plz[1,i],tau_pl[i],ENC_pl[i])
		I_z[2,i] += sub_typeII(bio_pl[i],enc_plz[2,i],tau_pl[i],ENC_pl[i])
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

###! Add consumed zooplankton to diet
function sub_consume_pizoo!(I,Iz)
	for i = 1:PI_N
		I[i] += (Iz[1,i] + Iz[2,i])
	end
end
function sub_consume_plzoo!(I,Iz)
	for i = 1:PL_N
		I[i] += (Iz[1,i] + Iz[2,i])
	end
end


###! ENERGY AVAILABLE FOR GROWTH NU
function sub_nu_pi!(nu,bio,I,met)
	for i = 1:PI_N
		nu[i] = ((I[i]/bio[i])*PI_lambda[i]) - met[i]
	end
end
function sub_nu_pl!(nu,bio,I,met)
	for i = 1:PL_N
		nu[i] = ((I[i]/bio[i])*PL_lambda[i]) - met[i]
	end
end
function sub_nu_de!(nu,bio,I,met)
	for i = 1:DE_N
		nu[i] = ((I[i]/bio[i])*DE_lambda[i]) - met[i]
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
		end
	end
end
function sub_rep_pl!(rep,nu,bio)
	for i = 1:PL_N
		if nu[i] > 0.
			rep[i] = (1-PL_K[i]) * nu[i] * bio[i]
		else
			rep[i] = 0.
		end
	end
end
function sub_rep_de!(rep,nu,bio)
	for i = 1:DE_N
		if nu[i] > 0.
			rep[i] = (1-DE_K[i]) * nu[i] * bio[i]
		else
			rep[i] = 0.
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
	for i = 2:PI_N
		bio[i] += (mat[i-1] + grw[i] - mrt[i] - rep[i] - mat[i] - d[i]) * DT
	end
end
function sub_update_pl!(bio,rep,grw,mat,d,mrt)
	bio[1] += (sum(rep) + grw[1] - mat[1] - d[1] - mrt[1]) * DT
	for i = 2:PL_N
		bio[i] += (mat[i-1] + grw[i] - mrt[i] - rep[i] - mat[i] - d[i]) * DT
	end
end
function sub_update_de!(bio,rep,grw,mat,d,mrt)
	bio[1] += (sum(rep) + grw[1] - mat[1] - d[1] - mrt[1]) * DT
	for i = 2:DE_N
		bio[i] += (mat[i-1] + grw[i] - mrt[i] - rep[i] - mat[i] - d[i]) * DT
	end
end
function sub_update_be!(bio,d,det)
    bio[1] += (det - d[1] - (bio[1]*0.01)) * DT # with sedimentation
end

###! Fishing
function sub_fishing!(bio_pi,bio_pl,bio_de,AREA)
	ALL_BIO = Array(Float64,NX*(PI_N+PL_N+DE_N))
	ALL_pi = Array(Float64,NX,PI_N)
	ALL_pl = Array(Float64,NX,PI_N)
	ALL_de = Array(Float64,NX,PI_N)

	#! Calc total biomass of fish in the ocean
	bio_pi .*= AREA;
	bio_pl .*= AREA;
	bio_de .*= AREA;

	#for i = 1:NX
	#	ALL_pi[i,:] = bio_pi[i] * AREA[i]
	#	ALL_pl[i,:] = bio_pl[i] * AREA[i]
	#	ALL_de[i,:] = bio_de[i] * AREA[i]
	#end

	#! Total fish biomass
	ALL_BIO[:] = [bio_pi[:]; bio_pl[:]; bio_de[:]]
	#ALL_BIO[:] = [ALL_pi[:]; ALL_pl[:]; ALL_de[:]]

	#! Calc Fishing using Ideal Free Distribution
	#! i.e. fishing in proportion to biomass
	ALL_BIO[:] -= (ALL_BIO/sum(ALL_BIO)) * FISHING # 80MT in kg

	#! Recalc biomasses
	ALL_pi = reshape(ALL_BIO[1:NX*PI_N],(NX,PI_N))
	splice!(ALL_BIO,1:NX*PI_N)
	ALL_pl = reshape(ALL_BIO[1:NX*PL_N],(NX,PL_N))
	splice!(ALL_BIO,1:NX*PI_N)
	ALL_de = reshape(ALL_BIO[1:NX*DE_N],(NX,DE_N))

	#! Calc total biomass of fish in the ocean
	for i = 1:NX
		bio_pi[i] = ALL_pi[i,:] ./ AREA[i]
		bio_pl[i] = ALL_pl[i,:] ./ AREA[i]
		bio_de[i] = ALL_de[i,:] ./ AREA[i]
	end
	
end


###! Forward Euler checks
function sub_check_pi!(bio)
	for i = 1:PI_N
		if bio[i] <= 0.0 || isnan(bio[i])==1
			bio[i] = eps()
		end
	end
end
function sub_check_pl!(bio)
	for i = 1:PL_N
		if bio[i] <= 0.0 || isnan(bio[i])==1
			bio[i] = eps()
		end
	end
end
function sub_check_de!(bio)
	for i = 1:DE_N
		if bio[i] <= 0.0 || isnan(bio[i])==1
			bio[i] = eps()
		end
	end
end
function sub_check_be!(bio)
	if bio[1] <= 0.0 || isnan(bio[1])==1
		bio[1] = eps()
	end
end




		
