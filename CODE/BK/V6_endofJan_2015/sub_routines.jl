###! Get COBALT data
function get_COBALT(COBALT::Dict,ID::Int,X::Int,DY::Int)
	TEMP_p   = COBALT["Tp"][ID,DY]
	TEMP_b   = COBALT["Tb"][ID,DY]
	ZOO      = float64([COBALT["Zm"][ID,DY],COBALT["Zl"][ID,DY]])
	DZc      = float64([COBALT["dZm"][ID,DY],COBALT["dZl"][ID,DY]])
	WI       = float64(COBALT["det"][ID,DY])
	return TEMP_p::Float32, TEMP_b::Float32, 
		   ZOO::Array{Float64,1}, DZc::Array{Float64,1},WI::Float64
end


###! Metabolism
function sub_metabolism(FISH::fish,PRM::param,TEMP_p::Float32)
	@inbounds for i = 1:PRM.N
		fnc_met!(FISH.met,i,PRM.s[i],TEMP_p)
	end
	return FISH::fish
end

function sub_metabolism!(FISH::fish,PRM::param,TEMP_p::Float32)
	@inbounds for i = 1:PRM.N
		fnc_met!(FISH.met,i,PRM.s[i],TEMP_p)
	end
    nothing
end


###! Encounter rate
function sub_enc_pi(PISC::fish,PLAN::fish,DETR::fish,ZOO::Array{Float64,1},PRM_PI::param,X::Int)
	# note: piscivore splits time between necton and benthos
	# in proportion to food availability
	rr = sum(PLAN.bio[:,X]) / (sum(DETR.bio[:,X])+sum(PLAN.bio[:,X]))
	for i = 1:PRM_PI.N
		fnc_enc!(PISC.enc["pipi"],i,PISC.bio[:,X],PRM_PI.Phi_PI[:,i],PRM_PI.a[i]*rr)
		fnc_enc!(PISC.enc["pipl"],i,PLAN.bio[:,X],PRM_PI.Phi_PL[:,i],PRM_PI.a[i]*rr)
		fnc_enc!(PISC.enc["pide"],i,DETR.bio[:,X],PRM_PI.Phi_DE[:,i],PRM_PI.a[i]*(1-rr))
		fnc_enc!(PISC.enc["piz"],i,ZOO,PRM_PI.Phi_Z[:,i],PRM_PI.a[i]*rr)
	end
	return PISC::fish
end

function sub_enc_pl(PLAN::fish,ZOO::Array{Float64,1},PRM_PL::param)
	for i = 1:PRM_PL.N
		fnc_enc!(PLAN.enc["plz"],i,ZOO,PRM_PL.Phi_Z[:,i],PRM_PL.a[i])
	end
	return PLAN::fish
end

function sub_enc_de(DETR::fish,W::Array{Float64,1},PRM_DE::param,X::Int)
	for i = 1:PRM_DE.N
		fnc_enc!(DETR.enc["dew"],i,W,PRM_DE.Phi_W[:,i],PRM_DE.a[i])
		fnc_enc!(DETR.enc["dede"],i,DETR.bio[:,X],PRM_DE.Phi_DE[:,i],PRM_DE.a[i])
	end
	return DETR::fish
end


###! Total biomass encountered
function sub_ENC(PISC::fish,PLAN::fish,DETR::fish)
    PISC.ENC = (sum(PISC.enc["pipi"],1)
              	   + sum(PISC.enc["pipl"],1)
              	   + sum(PISC.enc["pide"],1))
    PLAN.ENC = sum(PLAN.enc["plz"],1)
    DETR.ENC = DETR.enc["dew"] + sum(DETR.enc["dede"],1)
    return PISC::fish,PLAN::fish,DETR::fish
end


###! Handling times
function sub_tau(FISH::fish,PRM::param)
	@inbounds for i = 1:PRM.N
		fnc_tau!(FISH.tau,i,FISH.met[i])
	end
	return FISH::fish
end


###! Consumption
###! Consumption
function sub_consume_pi(PISC::fish,PLAN::fish,DETR::fish,W::detritus,
					 PRM_PI::param,PRM_PL::param,PRM_DE::param,
					 X::Int)
	PISC.I = zeros(PRM_PI.N)
	PISC.d = zeros(PRM_PI.N)
	PLAN.d = zeros(PRM_PL.N)
	DETR.d = zeros(PRM_DE.N)
	for i = 1:PRM_PI.N # pred
		@inbounds for j = 1:PRM_PI.N # prey
		fnc_cons!(PISC.I,PISC.d,i,j,PISC.bio[i,X],PISC.bio[j,X],PISC.enc["pipi"][j,i],PISC.ENC[i],PISC.tau[i])
		end
		@inbounds for j = 1:PRM_PL.N
		fnc_cons!(PISC.I,PLAN.d,i,j,PISC.bio[i,X],PLAN.bio[j,X],PISC.enc["pipl"][j,i],PISC.ENC[i],PISC.tau[i])
		end
		@inbounds for j = 1:PRM_DE.N
		fnc_cons!(PISC.I,DETR.d,i,j,PISC.bio[i,X],DETR.bio[j,X],PISC.enc["pide"][j,i],PISC.ENC[i],PISC.tau[i])
		end
	end
	return PISC::fish,PLAN::fish,DETR::fish
end


function sub_consume_pl(PISC::fish,PLAN::fish,PRM_PI::param,PRM_PL::param,
					 ZOO::Array{Float64,1},DZc::Array{Float64,1},X::Int)
	# Calculates TOTAL biomass consumerd g d-1
	PLAN.I = zeros(PRM_PL.N)
	DZs    = zeros(2) ## zooplankton mortality for offline coupling

	#! Planktivore + Piscivore eating zooplankton
	#! heavy due to COBALT offline coupling
	dzm = zeros(1); dzl = zeros(1);
	con_zm_pl = zeros(PRM_PL.N)
	con_zm_pi = zeros(PRM_PI.N)
	con_zl_pl = zeros(PRM_PL.N)
	con_zl_pi = zeros(PRM_PI.N)
	@inbounds for i = 1:PRM_PL.N
		# planktivore
		fnc_cons!(con_zm_pl,dzm,i,1,PLAN.bio[i,X],ZOO[1],PLAN.enc["plz"][1,i],PLAN.ENC[i],PLAN.tau[i])
		fnc_cons!(con_zl_pl,dzl,i,1,PLAN.bio[i,X],ZOO[2],PLAN.enc["plz"][2,i],PLAN.ENC[i],PLAN.tau[i])
	end
	@inbounds for i = 1:PRM_PI.N
		# piscivore
		fnc_cons!(con_zm_pi,dzm,i,1,PISC.bio[i,X],ZOO[1],PISC.enc["piz"][1,i],PISC.ENC[i],PISC.tau[i])
		fnc_cons!(con_zl_pi,dzl,i,1,PISC.bio[i,X],ZOO[2],PISC.enc["piz"][2,i],PISC.ENC[i],PISC.tau[i])
	end

	## Add offline coupling here
	if dzm[1] > DZc[1]
		tot = sum(con_zm_pl) + sum(con_zm_pi)
		con_zm_pl = DZc[1] .* (con_zm_pl./tot)
		con_zm_pi = DZc[1] .* (con_zm_pi./tot)
	end
	if dzl[1] > DZc[2]
		tot = sum(con_zl_pl) + sum(con_zl_pi)
		con_zl_pl = DZc[2] .* (con_zl_pl./tot)
		con_zl_pi = DZc[2] .* (con_zl_pi./tot)
	end
	PLAN.I += con_zm_pl + con_zl_pl
	PISC.I += con_zm_pi + con_zl_pi
	return PISC::fish, PLAN::fish
end

function sub_consume_de(DETR::fish,W::detritus,PRM_DE::param,X::Int)
	# Calculates TOTAL biomass consumerd g d-1
	DETR.I = zeros(PRM_DE.N)
	DETR.d = zeros(PRM_DE.N)
	W.d    = zeros(Float64,1)
	DZs    = zeros(2) ## zooplankton mortality for offline coupling

    #! Detritivore
    @simd for i = 1:PRM_DE.N # pred
        # detritus feeding
        fnc_cons!(DETR.I,W.d,i,1,DETR.bio[i,X],W.bio[1,X],DETR.enc["dew"][i],DETR.ENC[i],DETR.tau[i])
        # canabalism
        @inbounds for j = 1:PRM_DE.N # prey
        fnc_cons!(DETR.I,DETR.d,i,j,DETR.bio[i,X],DETR.bio[j,X],DETR.enc["dede"][j,i],DETR.ENC[i],DETR.tau[i])
        end
    end
    return DETR::fish,W::detritus
end


###! Nu: total energy avail for somatic growht and reproduction
function sub_nu(FISH::fish,PRM::param,X::Int)
	@inbounds for i = 1:PRM.N
		fnc_nu!(FISH.nu,i,FISH.bio[i,X],FISH.I[i],FISH.met[i],PRM.lambda[i])
	end
	return FISH::fish
end


###! Gamma: energy avail for somatic growth
function sub_gamma(FISH::fish,bio::Vector{Float64},PRM::param)
	@inbounds for i = 1:PRM.N
		fnc_gamma!(FISH.gamma,i,FISH.nu[i],bio[i],FISH.d[i],PRM.K[i],PRM.z[i])
	end
	return FISH::fish
end


###! Biomass reproduced
function sub_rep(FISH::fish,bio::Vector{Float64},PRM::param)
	@inbounds for i = 1:PRM.N
		fnc_rep!(FISH.REP,i,FISH.nu[i],PRM.K[i],bio[i])
	end
	return FISH::fish
end


###! Somatic growth
function sub_grw(FISH::fish,bio::Vector{Float64},PRM::param)
	@inbounds for i = 1:PRM.N
		fnc_grw!(FISH.GRW,i,FISH.nu[i],bio[i])
	end
	return FISH::fish
end


###! Maturation
function sub_mat(FISH::fish,bio::Vector{Float64},PRM::param)
	@inbounds for i = 1:PRM.N
		fnc_mat!(FISH.MAT,i,FISH.gamma[i],bio[i])
	end
	return FISH::fish
end


###! Background mortality
function sub_mrt(FISH::fish,bio::Vector{Float64},PRM::param)
	@inbounds for i = 1:PRM.N
		fnc_mrt!(FISH.MRT,i,bio[i],PRM.mrt[i])
	end
	return FISH::fish
end



###! Update biomass
function sub_update_bio(FISH::fish,PRM::param,X::Int,DT::Float64)
	FISH.bio[1,X] += (sum(FISH.REP)+FISH.GRW[1]-FISH.MAT[1]-FISH.d[1]-FISH.MRT[1]) * DT
	FISH.bio[2:end,X] += (FISH.MAT[1:end-1] + FISH.GRW[2:end] - FISH.MRT[2:end]
				- FISH.REP[2:end] - FISH.MAT[2:end] - FISH.d[2:end]) * DT
	return FISH::fish
end
function sub_update_W(W::detritus,X::Int,DT::Float64)
	W.bio[1,X] += (W.I[1] - W.d[1] - (W.bio[1,X]*0.01)) * DT # with sedimentation
	return W::detritus
end


###! Forward Euler checks
function sub_check_fish(bio::Vector{Float64},PRM)
	for i = 1:PRM.N
		if bio[i] <= 0.
			bio[i] = eps()
		end
	end
	return bio
end

function sub_check_W(bio::Vector{Float64})
	if bio[1] <= 0.
		bio[1] = eps()
	end
	return bio
end


###########!!!!! Iterate Forward In Time !!!!!########
function sub_demog(PISC::fish,PLAN::fish,DETR::fish,W::detritus,
				   PRM_PI::param,PRM_PL::param,PRM_DE::param,
				   TEMP_p::Float32,TEMP_b::Float32,ZOO::Array{Float64,1},
				   Dzc::Array{Float64,1},X::Int64)

	#! metabolism
	PISC = sub_metabolism(PISC,PRM_PI,TEMP_p);
	PLAN = sub_metabolism(PLAN,PRM_PL,TEMP_p);
	DETR = sub_metabolism(DETR,PRM_DE,TEMP_p);

	###! Encounter rates
	PISC = sub_enc_pi(PISC,PLAN,DETR,ZOO,PRM_PI,X);
	PLAN = sub_enc_pl(PLAN,ZOO,PRM_PL);
	DETR = sub_enc_de(DETR,W.bio[:,X],PRM_DE,X);

	###! Total prey biomass encountered
	PISC,PLAN,DETR = sub_ENC(PISC,PLAN,DETR);

	###! Handling times
	PISC = sub_tau(PISC,PRM_PI);
	PLAN = sub_tau(PLAN,PRM_PL);
	DETR = sub_tau(DETR,PRM_DE);

	###! Consumption/Predation
	PISC,PLAN,DETR = sub_consume_pi(PISC,PLAN,DETR,W,PRM_PI,PRM_PL,PRM_DE,X);
	PISC,PLAN = sub_consume_pl(PISC,PLAN,PRM_PI,PRM_PL,ZOO,DZc,X);
	DETR,W = sub_consume_de(DETR,W,PRM_DE,X);

	####! Energy available for growth
	PISC = sub_nu(PISC,PRM_PI,X);
	PLAN = sub_nu(PLAN,PRM_PL,X);
	DETR = sub_nu(DETR,PRM_DE,X);

	#! Energy available for somatic growth
	PISC = sub_gamma(PISC,PISC.bio[:,X],PRM_PI);
	PLAN = sub_gamma(PLAN,PLAN.bio[:,X],PRM_PL);
	DETR = sub_gamma(DETR,DETR.bio[:,X],PRM_DE);

	# Total biomass eggs produced
	PISC = sub_rep(PISC,PISC.bio[:,X],PRM_PI);
	PLAN = sub_rep(PLAN,PLAN.bio[:,X],PRM_PL);
	DETR = sub_rep(DETR,DETR.bio[:,X],PRM_DE);

	#! total biomass to somatic growth
	PISC = sub_grw(PISC,PISC.bio[:,X],PRM_PI);
	PLAN = sub_grw(PLAN,PLAN.bio[:,X],PRM_PL);
	DETR = sub_grw(DETR,DETR.bio[:,X],PRM_DE);

	#! total biomass that matures to next size class
	PISC = sub_mat(PISC,PISC.bio[:,X],PRM_PI);
	PLAN = sub_mat(PLAN,PLAN.bio[:,X],PRM_PL);
	DETR = sub_mat(DETR,DETR.bio[:,X],PRM_DE);

	#! total biomass that dies from background sources
	PISC = sub_mrt(PISC,PISC.bio[:,X],PRM_PI);
	PLAN = sub_mrt(PLAN,PLAN.bio[:,X],PRM_PL);
	DETR = sub_mrt(DETR,DETR.bio[:,X],PRM_DE);

	###! MASS BALANCE
	PISC = sub_update_bio(PISC,PRM_PI,X,DT);
	PLAN = sub_update_bio(PLAN,PRM_PL,X,DT);
	DETR = sub_update_bio(DETR,PRM_DE,X,DT);
	W    = sub_update_W(W,X,DT);

	###! Forward Euler checks
	PISC.bio[:,X] = sub_check_fish(PISC.bio[:,X],PRM_PI);
	PLAN.bio[:,X] = sub_check_fish(PLAN.bio[:,X],PRM_PL);
	DETR.bio[:,X] = sub_check_fish(DETR.bio[:,X],PRM_DE);
	W.bio[:,X]    = sub_check_W(W.bio[:,X]);

	return PISC::fish,PLAN::fish,DETR::fish,W::detritus
end


