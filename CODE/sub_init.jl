
#============== LOCAL EXP INITIAL CONDITIONS =============#
function sub_init(ID)

	#===== VARIABLES =====#
	###! Number of spatial cells
	NX = length(ID)

	###! fish
	#! biomass
	#X = 0.17981663628808964; # mean Zm off Spain
	X = 0.00001
	#X = eps();
	bio_pi = ones(NX,PI_N) .* X
	bio_pl = ones(NX,PL_N) .* X
	bio_de = ones(NX,DE_N) .* X
	bio_W  = ones(NX,1) .* X;

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
    tau_pi = Array(Float64,PI_N)
    tau_pl = Array(Float64,PL_N)
    tau_de = Array(Float64,DE_N)

    # mass ingested (I)
    I_pi = Array(Float64,PI_N)
    I_pl = Array(Float64,PL_N)
    I_de = Array(Float64,DE_N)

    # total energy available for growth
    nu_pi = Array(Float64,PI_N)
    nu_pl = Array(Float64,PL_N)
    nu_de = Array(Float64,DE_N)

    # energy available for somatic growth
    gamma_pi = Array(Float64,PI_N)
    gamma_pl = Array(Float64,PL_N)
    gamma_de = Array(Float64,DE_N)

    # mass lost to predation (g d-1)
    d_pi = Array(Float64,PI_N)
    d_pl = Array(Float64,PL_N)
    d_de = Array(Float64,DE_N)

    #! total biomass to reproduction
    REP_PI = Array(Float64,PI_N)
    REP_PL = Array(Float64,PL_N)
    REP_DE = Array(Float64,DE_N)

    #! total biomass to somatic growth
    GRW_PI = Array(Float64,PI_N)
    GRW_PL = Array(Float64,PL_N)
    GRW_DE = Array(Float64,DE_N)

    #! total biomass maturing to next sizeclass
    MAT_PI = Array(Float64,PI_N)
    MAT_PL = Array(Float64,PL_N)
	PI_bas = 0.0033*PI_s.^-0.13
    MAT_DE = Array(Float64,DE_N)

	#! background mortality + fishing
    MRT_PI = Array(Float64,PI_N)
    MRT_PL = Array(Float64,PL_N)
    MRT_DE = Array(Float64,DE_N)

    ###! Preallocate a dictionary for encounter rates
	# encounter rates
	enc_pi = zeros(PI_N,PI_N)
	enc_pl = zeros(PL_N,PI_N)
	enc_de = zeros(DE_N,PI_N)
	enc_z  = zeros(2,PI_N)
	enc_w  = zeros(1,DE_N)

 	# total biomass encoutered
    ENC_pi = zeros(PI_N)
    ENC_pl = zeros(PL_N)
    ENC_de = zeros(DE_N)
	
	# assign to types
	PISC = fish(bio_pi,met_pi,I_pi,tau_pi,nu_pi,gamma_pi,d_pi,REP_PI,GRW_PI,MAT_PI,MRT_PI,
				enc_pi,enc_pl,enc_de,enc_z,enc_w,ENC_pi)
	PLAN = fish(bio_pl,met_pl,I_pl,tau_pl,nu_pl,gamma_pl,d_pl,REP_PL,GRW_PL,MAT_PL,MRT_PL,
				enc_pi,enc_pl,enc_de,enc_z,enc_w,ENC_pl)
	DETR = fish(bio_de,met_de,I_de,tau_de,nu_de,gamma_de,d_de,REP_DE,GRW_DE,MAT_DE,MRT_DE,
				enc_pi,enc_pl,enc_de,enc_z,enc_w,ENC_de)

	###! Detritus
	W = detritus(bio_W,zeros(Float64,1),zeros(Float64,1))

	return PISC, PLAN, DETR, W
end


