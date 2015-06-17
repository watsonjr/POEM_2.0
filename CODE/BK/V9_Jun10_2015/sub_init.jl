
#============== INITIAL CONDITIONS =============#
function sub_init_env(ID)
	###! Number of spatial cells
	NX = length(ID)
	###! environment
	ENV_Tp = Array(Float64,NX)
	ENV_Tb = Array(Float64,NX)
	ENV_Zm = Array(Float64,NX)
	ENV_Zl = Array(Float64,NX)
	ENV_dZm = Array(Float64,NX)
	ENV_dZl = Array(Float64,NX)
	ENV_det = Array(Float64,NX)
	#return ENV_Tp,ENV_Tb,ENV_Zm,ENV_Zl,ENV_dZm,ENV_dZl,ENV_det
	ENVR = environment(ENV_Tp,ENV_Tb,ENV_Zm,ENV_Zl,ENV_dZm,ENV_dZl,ENV_det)
end

function sub_init_fish(ID)

	#===== VARIABLES =====#
	###! Number of spatial cells
	NX = length(ID)

	###! fish
	#! biomass
	#X = 0.17981663628808964; # mean Zm off Spain
	X = 0.00001; # very small amount
	bio_pi = Array(Any,NX)
	bio_pl = Array(Any,NX)
	bio_de = Array(Any,NX)
	bio_be = Array(Any,NX)
	for i = 1:NX
		bio_pi[i] = ones(Float64,PI_N) .* X
		bio_pl[i] = ones(Float64,PL_N) .* X
		bio_de[i] = ones(Float64,DE_N) .* X
		bio_be[i] = ones(Float64,1) .* X
	end

	# fraction of time spent in the pelagic
	tdif_pi = Array(Float64,NX)
	tdif_pl = Array(Float64,NX)
	tdif_de = Array(Float64,NX)

	# temperature multiplier
	tmet_pi = Array(Float64,NX)
	tmet_pl = Array(Float64,NX)
	tmet_de = Array(Float64,NX)

	# metabolism
	met_pi = Array(Any,NX)
	met_pl = Array(Any,NX)
	met_de = Array(Any,NX)
	for i = 1:NX
		met_pi[i] = Array(Float64,PI_N)
		met_pl[i] = Array(Float64,PL_N)
		met_de[i] = Array(Float64,DE_N)
	end
	
    # handling times
    tau_pi = Array(Any,NX)
    tau_pl = Array(Any,NX)
    tau_de = Array(Any,NX)
	for i = 1:NX
		tau_pi[i] = Array(Float64,PI_N)
		tau_pl[i] = Array(Float64,PL_N)
		tau_de[i] = Array(Float64,DE_N)
	end
	
    # mass ingested (I)
    I_pi = Array(Any,NX)
    I_pl = Array(Any,NX)
    I_de = Array(Any,NX)
    for i = 1:NX
    	I_pi[i] = Array(Float64,PI_N)
    	I_pl[i] = Array(Float64,PL_N)
    	I_de[i] = Array(Float64,DE_N)
    end

	# mass of zooplankton ingested (for COBALT link)
	I_piz = Array(Any,NX)
	I_plz = Array(Any,NX)
	for i = 1:NX
		I_piz[i] = Array(Float64,2,PI_N)
		I_plz[i] = Array(Float64,2,PL_N)
	end

    # mass lost to predation (g d-1)
    d_pi = Array(Any,NX)
    d_pl = Array(Any,NX)
    d_de = Array(Any,NX)
    d_be = Array(Any,NX)
    for i = 1:NX
    	d_pi[i] = Array(Float64,PI_N)
    	d_pl[i] = Array(Float64,PL_N)
    	d_de[i] = Array(Float64,DE_N)
    	d_be[i] = Array(Float64,1,1)
    end

    # total energy available for growtha
    nu_pi = Array(Any,NX)
    nu_pl = Array(Any,NX)
    nu_de = Array(Any,NX)
    for i = 1:NX
    	nu_pi[i] = Array(Float64,PI_N)
    	nu_pl[i] = Array(Float64,PL_N)
    	nu_de[i] = Array(Float64,DE_N)
    end

    # energy available for somatic growth
    gamma_pi = Array(Any,NX)
    gamma_pl = Array(Any,NX)
    gamma_de = Array(Any,NX)
    for i = 1:NX
    	gamma_pi[i] = Array(Float64,PI_N)
    	gamma_pl[i] = Array(Float64,PL_N)
    	gamma_de[i] = Array(Float64,DE_N)
    end

    #! total biomass to reproduction
    REP_PI = Array(Any,NX)
    REP_PL = Array(Any,NX)
    REP_DE = Array(Any,NX)
    for i = 1:NX
		REP_PI[i] = Array(Float64,PI_N)
		REP_PL[i] = Array(Float64,PL_N)
		REP_DE[i] = Array(Float64,DE_N)
    end

    #! total biomass to somatic growth
    GRW_PI = Array(Any,NX)
    GRW_PL = Array(Any,NX)
    GRW_DE = Array(Any,NX)
    for i = 1:NX
    	GRW_PI[i] = Array(Float64,PI_N)
    	GRW_PL[i] = Array(Float64,PL_N)
    	GRW_DE[i] = Array(Float64,DE_N)
    end

    #! total biomass maturing to next sizeclass
    MAT_PI = Array(Any,NX)
    MAT_PL = Array(Any,NX)
    MAT_DE = Array(Any,NX)
    for i = 1:NX
    	MAT_PI[i] = Array(Float64,PI_N)
    	MAT_PL[i] = Array(Float64,PL_N)
    	MAT_DE[i] = Array(Float64,DE_N)
    end

	#! background mortality + fishing
    MRT_PI = Array(Any,NX)
    MRT_PL = Array(Any,NX)
    MRT_DE = Array(Any,NX)
    for i = 1:NX
    	MRT_PI[i] = Array(Float64,PI_N)
    	MRT_PL[i] = Array(Float64,PL_N)
    	MRT_DE[i] = Array(Float64,DE_N)
    end

	#! pisc encounter rates
	enc_pipi = Array(Any,NX)
	enc_pipl = Array(Any,NX)
	enc_pide = Array(Any,NX)
	enc_piz  = Array(Any,NX)
    for i = 1:NX
        enc_pipi[i] = Array(Float64,PI_N,PI_N)
        enc_pipl[i] = Array(Float64,PL_N,PI_N)
        enc_pide[i] = Array(Float64,DE_N,PI_N)
        enc_piz[i] = Array(Float64,2,PI_N)
    end

	#! total pisc ENCounter rates
	ENC_pipi = Array(Any,NX)
	ENC_pipl = Array(Any,NX)
	ENC_pide = Array(Any,NX)
	ENC_piz  = Array(Any,NX)
    for i = 1:NX
        ENC_pipi[i] = Array(Float64,1,PI_N)
        ENC_pipl[i] = Array(Float64,1,PI_N)
        ENC_pide[i] = Array(Float64,1,PI_N)
        ENC_piz[i] = Array(Float64,1,PI_N)
    end

	#! planktivore encounter rates
	enc_plz  = Array(Any,NX)
    for i = 1:NX
        enc_plz[i] = Array(Float64,2,PL_N)
    end
	
	#! detrivore encounter rates
	enc_dede = Array(Any,NX)
	enc_debe = Array(Any,NX)
    for i = 1:NX
        enc_dede[i] = Array(Float64,DE_N,DE_N)
        enc_debe[i] = Array(Float64,1,DE_N)
    end

	#! total detrivore Encounter rates
	ENC_dede = Array(Any,NX)
	ENC_debe = Array(Any,NX)
    for i = 1:NX
        ENC_dede[i] = Array(Float64,1,DE_N)
        ENC_debe[i] = Array(Float64,1,DE_N)
    end

 	# total biomass encoutered
	ENC_pl  = Array(Any,NX)
    ENC_pi = Array(Any,NX)
    ENC_de = Array(Any,NX) 
    for i = 1:NX
    	ENC_pi[i] = Array(Float64,1,PI_N)
		ENC_pl[i] = Array(Float64,1,PL_N)
    	ENC_de[i] = Array(Float64,1,DE_N)
    end

	# assign to types
	PISC = piscivore(bio_pi,tmet_pi,tdif_pi,met_pi,I_pi,I_piz,tau_pi,
                nu_pi,gamma_pi,d_pi,REP_PI,GRW_PI,MAT_PI,MRT_PI,
                enc_pipi,enc_pipl,enc_pide,enc_piz,ENC_pi)
	PLAN = planktivore(bio_pl,tmet_pl,met_pl,I_pl,I_plz,tau_pl,
				nu_pl,gamma_pl,d_pl,REP_PL,GRW_PL,MAT_PL,MRT_PL,
				enc_plz,ENC_pl)
	DETR = detrivore(bio_de,tmet_de,met_de,I_de,tau_de,
                nu_de,gamma_de,d_de,REP_DE,GRW_DE,MAT_DE,MRT_DE,
                enc_dede,enc_debe,ENC_de)

	###! Detritus
    BENT = detritus(bio_be,d_be)

	return PISC,PLAN,DETR,BENT
end


