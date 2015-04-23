
#============== LOCAL EXP INITIAL CONDITIONS =============#
function sub_init_local(PRM_PI,PRM_PL,PRM_DE)

	#===== VARIABLES =====#
	###! fish
	#! biomass
	#X = 0.17981663628808964; # mean Zm off Spain
	X = 0.00001
	#X = eps();
	bio_pi = ones(PRM_PI.N) .* X
	bio_pl = ones(PRM_PL.N) .* X
	bio_de = ones(PRM_DE.N) .* X
	bio_W  = ones(1) .* X;

	# metabolism
	met_pi = Array(Float64,PRM_PI.N)
	met_pl = Array(Float64,PRM_PL.N)
	met_de = Array(Float64,PRM_DE.N)

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

    # mass lost to predation (g d-1)
    d_pi = Array(Float64,PRM_PI.N)
    d_pl = Array(Float64,PRM_PL.N)
    d_de = Array(Float64,PRM_DE.N)

    #! total biomass to reproduction
    REP_PI = Array(Float64,PRM_PI.N)
    REP_PL = Array(Float64,PRM_PL.N)
    REP_DE = Array(Float64,PRM_DE.N)

    #! total biomass to somatic growth
    GRW_PI = Array(Float64,PRM_PI.N)
    GRW_PL = Array(Float64,PRM_PL.N)
    GRW_DE = Array(Float64,PRM_DE.N)

    #! total biomass maturing to next sizeclass
    MAT_PI = Array(Float64,PRM_PI.N)
    MAT_PL = Array(Float64,PRM_PL.N)
    MAT_DE = Array(Float64,PRM_DE.N)

	#! background mortality + fishing
    MRT_PI = Array(Float64,PRM_PI.N)
    MRT_PL = Array(Float64,PRM_PL.N)
    MRT_DE = Array(Float64,PRM_DE.N)

    ###! Preallocate a dictionary for encounter rates
	# encounter rates
    enc_pi = ["pipi"=> Array(Float64,PRM_PI.N,PRM_PI.N),
              "pipl"=> Array(Float64,PRM_PL.N,PRM_PI.N),
              "pide"=> Array(Float64,PRM_DE.N,PRM_PI.N),
              "piz" => Array(Float64,2,PRM_PI.N)]
    enc_pl = ["plz" => Array(Float64,2,PRM_PL.N)]
    enc_de = ["dede"=> Array(Float64,PRM_DE.N,PRM_DE.N),
              "dew" => Array(Float64,1,PRM_DE.N)]
	
	# total biomass encoutered
    ENC_pi = zeros(PRM_PI.N)
    ENC_pl = zeros(PRM_PL.N)
    ENC_de = zeros(PRM_DE.N)
	
	# assign to types
	PISC = fish(bio_pi,met_pi,I_pi,tau_pi,nu_pi,gamma_pi,d_pi,REP_PI,GRW_PI,MAT_PI,MRT_PI,enc_pi,ENC_pi)
	PLAN = fish(bio_pl,met_pl,I_pl,tau_pl,nu_pl,gamma_pl,d_pl,REP_PL,GRW_PL,MAT_PL,MRT_PL,enc_pl,ENC_pl)
	DETR = fish(bio_de,met_de,I_de,tau_de,nu_de,gamma_de,d_de,REP_DE,GRW_DE,MAT_DE,MRT_DE,enc_de,ENC_de)

	###! Detritus
	W = detritus(bio_W,zeros(Float64,1),zeros(Float64,PRM_DE.N))

	return PISC, PLAN, DETR, W
end


