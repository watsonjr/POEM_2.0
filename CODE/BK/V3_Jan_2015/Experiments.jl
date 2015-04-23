

#========= LOCAL EXPERIMENT =========#
#! size-based model run at one location
function make_local()
	##! Load COBALT
	COBALT = load("./Data/Data_daily_averages.jld");
	GRD = load("./Data/Data_grid.jld")

	##! Extract a particular location
	XY = zeros(360,200);
	XY[GRD["ID"]] =[1:GRD["N"]]
	id = XY[272,156] # Iberian location
	
	#! get data
	Tp  = squeeze(COBALT["Tp"][id,:],1);
	Tb  = squeeze(COBALT["Tb"][id,:],1);
	Zm  = squeeze(COBALT["Zm"][id,:],1);
	Zl  = squeeze(COBALT["Zl"][id,:],1);
	dZm = squeeze(COBALT["dZm"][id,:],1);
	dZl = squeeze(COBALT["dZl"][id,:],1);
	det = squeeze(COBALT["det"][id,:],1)

	#! Setup data files to save
	Data_PISC = open("./Data/CSV/Data_PISC.csv","w")
	Data_PLAN = open("./Data/CSV/Data_PLAN.csv","w")
	Data_DETR = open("./Data/CSV/Data_DETR.csv","w")
	Data_W    = open("./Data/CSV/Data_W.csv","w")

	##! Setup intermediate variables
	# mass encountered
	enc_pi_pi = Array(Float64,PRM_PI.N,PRM_PI.N) 
	enc_pi_pl = Array(Float64,PRM_PL.N,PRM_PI.N) 
	enc_pi_de = Array(Float64,PRM_DE.N,PRM_PI.N) 
	enc_pi_z  = Array(Float64,2,PRM_PI.N) 
	enc_pl_z  = Array(Float64,2,PRM_PL.N)
	enc_de_de = Array(Float64,PRM_DE.N,PRM_DE.N)
	enc_de_w  = Array(Float64,1,PRM_DE.N)

	# mass metabolized 
	meta_pi = Array(Float64,PRM_PI.N)
	meta_pl = Array(Float64,PRM_PL.N)
	meta_de = Array(Float64,PRM_DE.N)

	# handling times
	tau_pi = Array(Float64,PRM_PI.N)
	tau_pl = Array(Float64,PRM_PL.N)
	tau_de = Array(Float64,PRM_DE.N)

	# mass ingested (I)
	I_pi = Array(Float64,PRM_PI.N)
	I_pl = Array(Float64,PRM_PL.N)
	I_de = Array(Float64,PRM_DE.N)

	# total energy available for growth
	nu_pi = Array(Float64,PRM_PI.N)
	nu_pl = Array(Float64,PRM_PL.N)
	nu_de = Array(Float64,PRM_DE.N)

	# energy available for somatic growth
	gamma_pi = Array(Float64,PRM_PI.N)
	gamma_pl = Array(Float64,PRM_PL.N)
	gamma_de = Array(Float64,PRM_DE.N)

	# mass lost to predation (d)
	# in biomass specific units
	d_pi = Array(Float64,PRM_PI.N)
	d_pl = Array(Float64,PRM_PL.N)
	d_de = Array(Float64,PRM_DE.N)

	# mass lost to predation (d)
	# note; in total biomass units, factoring pred density
	D_PI = Array(Float64,PRM_PI.N)
	D_PL = Array(Float64,PRM_PL.N)
	D_DE = Array(Float64,PRM_DE.N)
	D_W  = Array(Float64,PRM_DE.N)
	
	#! Forward Euler integration of size-based model
	YEARS = 10;
	DAYS = 365;
	for YR = 1:YEARS # years
		for DY = 1:DAYS # days

			#! ticker
			println(YR/YEARS," , ",DY/DAYS)

			#! Write to file
			writecsv(Data_PISC,float32(BIOMASS.PISC'))
			writecsv(Data_PLAN,float32(BIOMASS.PLAN'))
			writecsv(Data_DETR,float32(BIOMASS.DETR'))
			writecsv(Data_W,float32(BIOMASS.W))

			#! COBALT information
			ZM     = Zm[DY]
			ZL     = Zl[DY]
			TEMP_p = Tp[DY]
			TEMP_b = Tb[DY]
			DET    = det[DY]

			#! metabolism (required to calc handling times)
			for i = 1:PRM_PL.N
				meta_pl[i] = fnc_met(PRM_PL.s[i],TEMP_p)
			end
			for i = 1:PRM_PI.N
				meta_pi[i] = fnc_met(PRM_PI.s[i],TEMP_p)
			end
			for i = 1:PRM_DE.N
				meta_de[i] = fnc_met(PRM_DE.s[i],TEMP_p) ## change this when changed bottom temp to account for pressure
			end

			###! Encounter rates
			for i = 1:PRM_PI.N
				enc_pi_pi[:,i] = BIOMASS.PISC.*PRM_PI.Phi_PI[:,i].*PRM_PI.a[i]
				enc_pi_pl[:,i] = BIOMASS.PLAN.*PRM_PI.Phi_PL[:,i].*PRM_PI.a[i]
				enc_pi_z[:,i]  = [ZM;ZL].*PRM_PI.Phi_Z[:,i].*PRM_PI.a[i]
				enc_pi_de[:,i] = BIOMASS.DETR.*PRM_PI.Phi_DE[:,i].*PRM_PI.a[i]
			end
			for i = 1:PRM_PL.N
				enc_pl_z[:,i] = [ZM;ZL].*PRM_PL.Phi_Z[:,i].*PRM_PL.a[i]
			end
			for i = 1:PRM_DE.N
				enc_de_w[i] = BIOMASS.W[1].*PRM_DE.Phi_W[i].*PRM_DE.a[i]
				enc_de_de[:,i] = BIOMASS.DETR.*PRM_DE.Phi_DE[:,i].*PRM_DE.a[i]
			end

			###! Total prey biomass encountered
			ENC_PI = sum(enc_pi_pi,1) + sum(enc_pi_pl,1) + sum(enc_pi_de,1)
			ENC_PL = sum(enc_pl_z,1)
			ENC_DE = enc_de_w + sum(enc_de_de,1)

			###! Handling times
			for i = 1:PRM_PI.N
				tau_pi[i] = fnc_tau(meta_pi[i])
			end
			for i = 1:PRM_PL.N
				tau_pl[i] = fnc_tau(meta_pl[i])
			end
			for i = 1:PRM_DE.N
				tau_de[i] = fnc_tau(meta_de[i])
			end

			###! Consumption/Predation
			##! note, this is total eat/eaten (i.e. factoring in
			##! pred and prey biomasses: NOTE death by predation in total
			##! biomass units, after factoring in pred biomass density
			#! Piscivore
			# piscivore eating piscivore
			for i = 1:PRM_PI.N # prey
				for j = 1:PRM_PI.N # pred
					eat = fnc_cons(BIOMASS.PISC[i],
										   enc_pi_pi[i,j],ENC_PI[j],tau_pi[j])
					I_pi[j] += eat
					d_pi[j] += eat
					D_PI[i] += eat .* BIOMASS.PISC[j]
				end
			end
			# piscivore eating planktivore
			for i = 1:PRM_PL.N # prey
				for j = 1:PRM_PI.N # pred
					eat = fnc_cons(BIOMASS.PLAN[i],
										   enc_pi_pl[i,j],ENC_PI[j],tau_pi[j])
					I_pi[j] += eat
					d_pi[j] += eat
					D_PL[i] += eat .* BIOMASS.PISC[j]
				end
			end
			# eating detritivore
			for i = 1:PRM_DE.N # prey
				for j = 1:PRM_PI.N # pred
					eat = fnc_cons(BIOMASS.DETR[i],
										   enc_pi_de[i,j],ENC_PI[j],tau_pi[j])
					I_pi[j] += eat
					d_pi[j] += eat
					D_DE[i] += eat .* BIOMASS.PISC[j]
				end
			end


			#! Planktivore
			ZOO = [ZM,ZL]
			for i = 1:2
				for j = 1:PRM_PL.N # pred
					eat = fnc_cons(ZOO[i],enc_pl_z[i,j],ENC_PL[j],tau_pl[j])
					I_pl[j] += eat
					# could add COBALT offline check here
					# eaten should = COBALT mortality
				end
			end


			#! Detritivore
			for i = 1:PRM_DE.N # prey
				for j = 1:PRM_DE.N # pred
					eat = fnc_cons(BIOMASS.DETR[i],
										   enc_de_de[i,j],ENC_DE[j],tau_de[j])
					I_de[j] += eat
					d_de[j] += eat
					D_DE[i] += eat .* BIOMASS.DETR[j]
				end
			end
			for i = 1:PRM_DE.N
					eat = fnc_cons(BIOMASS.W[1],enc_de_w[i],ENC_DE[i],tau_de[i])
					I_de[i] += eat
					d_de[i] += eat
					D_W[i] += eat .* BIOMASS.DETR[i]
			end


			####! Energy available for growth
			for i = 1:PRM_PI.N
				nu_pi[i] = fnc_nu(I_pi[i],meta_pi[i],PRM_PI.lambda[i])
			end
			for i = 1:PRM_PL.N
				nu_pl[i] = fnc_nu(I_pl[i],meta_pl[i],PRM_PL.lambda[i])
			end
			for i = 1:PRM_DE.N
				nu_de[i] = fnc_nu(I_de[i],meta_de[i],PRM_DE.lambda[i])
			end


			#! Energy available for somatic growth
			for i = 1:PRM_PI.N
				gamma_pi[i] = fnc_gamma(nu_pi[i],d_pi[i],PRM_PI.K[i],PRM_PI.z[i])
			end
			for i = 1:PRM_PL.N
				gamma_pl[i] = fnc_gamma(nu_pl[i],d_pl[i],PRM_PL.K[i],PRM_PL.z[i])
			end
			for i = 1:PRM_DE.N
				gamma_de[i] = fnc_gamma(nu_de[i],d_de[i],PRM_DE.K[i],PRM_DE.z[i])
			end


			####! TOTAL RATES OF CHANGE (factoring in biomass densities)
			#! total biomass to egg production
			# eggs produced
			REP_PI = Array(Float64,PRM_PI.N)
			REP_PL = Array(Float64,PRM_PL.N)
			REP_DE = Array(Float64,PRM_DE.N)
			#! piscivore
			for i = 1:PRM_PI.N
				REP_PI[i] = fnc_rep(nu_pi[i],PRM_PI.K[i],BIOMASS.PISC[i])
			end
			#! planktivore
			for i = 1:PRM_PL.N
				REP_PL[i] = fnc_rep(nu_pl[i],PRM_PL.K[i],BIOMASS.PLAN[i])
			end
			#! detritivore
			for i = 1:PRM_DE.N
				REP_DE[i] = fnc_rep(nu_de[i],PRM_DE.K[i],BIOMASS.DETR[i])
			end

			#! total biomass to somatic growth
			GRW_PI = Array(Float64,PRM_PI.N)
			GRW_PL = Array(Float64,PRM_PL.N)
			GRW_DE = Array(Float64,PRM_DE.N)
			for i = 1:PRM_PI.N
				GRW_PI[i] = nu_pi[i] * BIOMASS.PISC[i]
			end
			for i = 1:PRM_PL.N
				GRW_PL[i] = nu_pl[i] * BIOMASS.PLAN[i]
			end
			for i = 1:PRM_DE.N
				GRW_DE[i] = nu_de[i] * BIOMASS.DETR[i]
			end

			#! total biomass that matures to next size class
			MAT_PI = Array(Float64,PRM_PI.N)
			MAT_PL = Array(Float64,PRM_PL.N)
			MAT_DE = Array(Float64,PRM_DE.N)
			for i = 1:PRM_PI.N
				MAT_PI[i] = gamma_pi[i] * BIOMASS.PISC[i]
			end
			for i = 1:PRM_PL.N
				MAT_PL[i] = gamma_pl[i] * BIOMASS.PLAN[i]
			end
			for i = 1:PRM_DE.N
				MAT_DE[i] = gamma_de[i] * BIOMASS.DETR[i]
			end

			###! MASS BALANCE
			BIOMASS.PISC[1] += (sum(REP_PI) + GRW_PI[1] - MAT_PI[1] - D_PI[1])
			BIOMASS.PISC[2:end] += (MAT_PI[1:end-1] + GRW_PI[2:end] 
						- REP_PI[2:end] - MAT_PI[2:end] - D_PI[2:end])

			BIOMASS.PLAN[1] += (sum(REP_PL) + GRW_PL[1] - MAT_PL[1] - D_PL[1])
			BIOMASS.PLAN[2:end] += (GRW_PL[1:end-1] - REP_PL[2:end]
						- REP_PL[2:end] - MAT_PL[2:end] - D_PL[2:end])
			
			BIOMASS.DETR[1] += (sum(REP_DE) + GRW_DE[1] - MAT_DE[1] - D_DE[1])
			BIOMASS.DETR[2:end] += (GRW_DE[1:end-1] - REP_DE[2:end]
						- REP_DE[2:end] - MAT_DE[2:end] - D_DE[2:end])

			BIOMASS.W[1] += (DET - sum(D_W))

			###! Forward Euler checks
			for i = 1:PRM_PI.N
				if BIOMASS.PISC[i] <= 0
					BIOMASS.PISC[i] = eps()
				end
			end
			for i = 1:PRM_PL.N
				if BIOMASS.PLAN[i] <= 0
					BIOMASS.PLAN[i] = eps()
				end
			end
			for i = 1:PRM_DE.N
				if BIOMASS.DETR[i] <= 0
					BIOMASS.DETR[i] = eps()
				end
			end
			if BIOMASS.W[1] <= 0
				BIOMASS.W[1] = eps()
			end
		end
	end

	#! close write
	close(Data_PISC)
	close(Data_PLAN)
	close(Data_DETR)
	close(Data_W)
end
